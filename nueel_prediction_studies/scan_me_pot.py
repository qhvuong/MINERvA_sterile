#!/usr/bin/env python3
"""Inventory data exposure and produced-flux MC normalization POT for ME.

Reads Meta only; applies no event selection and does not infer data-to-MC
playlist mapping. Defaults to this repository's P8 data lists; this does not
claim that they exactly match the published nu-e exposure.
"""
import argparse
import csv
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent
PLAYLISTS=["1A","1B","1C","1D","1E","1F","1G","1L","1M","1N","1O","1P","5A","6A","6B","6C","6D","6E","6F","6G","6H","6I","6J"]
SOURCES=["ME1D","ME1M","ME1N","ME5A","ME6A"]


def read_meta_pot(ROOT,path):
    f=ROOT.TFile.Open(path,"READ")
    if not f or f.IsZombie() or f.TestBit(ROOT.TFile.kRecovered):
        raise RuntimeError("Invalid/recovered input: "+path)
    try:
        meta=f.Get("Meta")
        if not meta:
            raise RuntimeError("Missing Meta tree: "+path)
        leaves={name:meta.GetLeaf(name) for name in ["POT_Used","POT_Total"]}
        if not leaves["POT_Used"]:
            raise RuntimeError("Missing POT_Used: "+path)
        sums={name:0. if leaf else None for name,leaf in leaves.items()}
        total_issues=[]
        if meta.GetEntries()==0:
            raise RuntimeError("Empty Meta tree: "+path)
        for i in range(meta.GetEntries()):
            if meta.GetEntry(i)<=0:
                raise RuntimeError("Unreadable Meta entry: "+path)
            for name,leaf in leaves.items():
                if leaf:
                    if leaf.GetLen()!=1:
                        if name=="POT_Total":
                            total_issues.append(dict(entry=i,value="nonscalar"));continue
                        raise RuntimeError("Expected scalar "+name)
                    value=float(leaf.GetValue())
                    if not math.isfinite(value) or value<0:
                        if name=="POT_Total":
                            total_issues.append(dict(entry=i,value=repr(value)));continue
                        raise RuntimeError("Invalid {} at Meta entry {}: {}".format(name,i,path))
                    sums[name]+=value
        if total_issues:
            sums["POT_Total"]=None
        return dict(meta_entries=int(meta.GetEntries()),**sums,
            POT_Total_status="invalid" if total_issues else ("valid" if leaves["POT_Total"] else "missing"),
            POT_Total_issues=json.dumps(total_issues))
    finally:
        f.Close()


def save(out,name,rows):
    if rows:
        with (out/name).open("w",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--playlists",nargs="+",default=PLAYLISTS,help="Data playlists, e.g. 1D 1M 1N 5A 6A")
    ap.add_argument("--data-list-pattern",default=str(HERE.parent/"file_option/playlist_datame{playlist}p8.txt"))
    ap.add_argument("--flux-dir",type=Path,default=Path("/exp/minerva/data/users/qvuong/flux_studies/producedFluxes_p8_ME_offsetPPFX"))
    ap.add_argument("--max-data-files",type=int,default=0,help="0=all; nonzero is a PARTIAL diagnostic, not a complete exposure")
    ap.add_argument("--out",type=Path,default=Path("me_pot_inventory"))
    ap.add_argument("--resume-dir",type=Path,help="Reuse successful per-file rows from a previous inventory; read only missing files. Assumes input files have not changed.")
    args=ap.parse_args()
    if args.max_data_files<0:
        ap.error("max-data-files must be nonnegative")
    import ROOT
    ROOT.gROOT.SetBatch(True)
    args.out.mkdir(parents=True,exist_ok=True)
    totals=[];files=[];seen=set();failures=[];flux_rows=[];warnings=[];cache={};reused=0
    if args.resume_dir:
        with (args.resume_dir/"data_file_pot.csv").open() as stream:
            for row in csv.DictReader(stream):
                if row["path"] in cache:
                    raise RuntimeError("Duplicate cached data path")
                used=float(row["POT_Used"])
                total=float(row["POT_Total"]) if row.get("POT_Total") else None
                if not math.isfinite(used) or used<0 or int(row["meta_entries"])<=0:
                    raise RuntimeError("Invalid cached POT_Used: "+row["path"])
                if total is not None and (not math.isfinite(total) or total<0):
                    raise RuntimeError("Invalid cached POT_Total: "+row["path"])
                cache[row["path"]]=(row["playlist"],dict(meta_entries=int(row["meta_entries"]),POT_Used=used,POT_Total=total,
                    POT_Total_status=row.get("POT_Total_status") or ("valid" if total is not None else "missing"),
                    POT_Total_issues=row.get("POT_Total_issues") or "[]"))
    for playlist in args.playlists:
        playlist=playlist[2:] if playlist.startswith("ME") else playlist
        filelist=Path(args.data_list_pattern.format(playlist=playlist))
        paths=[s.strip() for s in filelist.read_text().splitlines() if s.strip() and not s.lstrip().startswith("#")]
        requested=paths[:args.max_data_files] if args.max_data_files else paths
        if not requested:
            raise RuntimeError("Empty list: "+str(filelist))
        subtotal={"POT_Used":0.,"POT_Total":0.};missing_total=False;good=0
        print("ME{}: reading Meta from {} of {} data files".format(playlist,len(requested),len(paths)),flush=True)
        for index,path in enumerate(requested):
            # Ignore alternate xrootd doors when detecting duplicate paths.
            identity=path[path.index('/pnfs/'):] if '/pnfs/' in path else path
            identity='/'.join(part for part in identity.split('/') if part)
            if identity in seen:
                raise RuntimeError("Duplicate data input across lists: "+path)
            seen.add(identity)
            try:
                if path in cache:
                    cached_playlist,values=cache[path]
                    if cached_playlist!="ME"+playlist:
                        raise RuntimeError("Cached playlist mismatch: "+path)
                    reused+=1
                else:
                    values=read_meta_pot(ROOT,path)
                if values["POT_Total_status"]!="valid":
                    warnings.append(dict(playlist="ME"+playlist,path=path,field="POT_Total",status=values["POT_Total_status"],details=values["POT_Total_issues"]))
                    print("  WARNING {} POT_Total; valid POT_Used retained: {}".format(values["POT_Total_status"],path),flush=True)
                files.append(dict(playlist="ME"+playlist,path=path,**values));good+=1
                subtotal["POT_Used"]+=values["POT_Used"]
                if values["POT_Total"] is None:
                    missing_total=True
                else:
                    subtotal["POT_Total"]+=values["POT_Total"]
            except RuntimeError as error:
                failures.append(dict(playlist="ME"+playlist,path=path,error=str(error)))
                print("  ERROR "+str(error),flush=True)
            if (index+1)%25==0:
                print("  read {}/{}".format(index+1,len(requested)),flush=True)
        complete=good==len(paths)
        totals.append(dict(playlist="ME"+playlist,beam="FHC" if playlist.startswith("1") else "RHC",filelist=str(filelist),
            files_in_list=len(paths),files_read=good,complete=complete,POT_Used=subtotal["POT_Used"],POT_Total=None if missing_total else subtotal["POT_Total"]))
        print("  POT_Used={:.12e}; {}".format(subtotal["POT_Used"],"COMPLETE" if complete else "INCOMPLETE: do not use as full exposure"),flush=True)
        save(args.out,"data_playlist_pot.csv",totals);save(args.out,"data_file_pot.csv",files)
    for source in SOURCES:
        for flavor in ["nue","nuebar","numu","numubar"]:
            path=args.flux_dir/(source+"_"+flavor+".root")
            f=ROOT.TFile.Open(str(path),"READ")
            if not f or f.IsZombie():
                failures.append(dict(playlist=source,path=str(path),error="Cannot read produced flux metadata"));continue
            try:
                parameter=f.Get("total_POT")
                value=float(parameter.GetVal()) if parameter else None
                flux_rows.append(dict(source=source,flavor=flavor,path=str(path),MC_total_POT=value))
                if value is None or not math.isfinite(value) or value<=0:
                    failures.append(dict(playlist=source,path=str(path),error="Missing/invalid total_POT normalization parameter"))
            finally:
                f.Close()
    save(args.out,"flux_mc_pot.csv",flux_rows)
    report=dict(data_totals={beam:sum(r["POT_Used"] for r in totals if r["beam"]==beam) for beam in ["FHC","RHC"]},
        data_playlists_complete=bool(totals) and all(r["complete"] for r in totals),failures=failures,
        warnings=warnings,reused_data_files=reused,resume_dir=str(args.resume_dir) if args.resume_dir else None,
        note="Data POT weights exposure mixtures. Produced flux total_POT is MC normalization already included in flux per POT; do not apply that normalization again. Data list coverage and mapping to the five MC sources must be confirmed before a paper comparison.")
    (args.out/"inventory.json").write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    return 1 if failures else 0


if __name__=="__main__":
    raise SystemExit(main())
