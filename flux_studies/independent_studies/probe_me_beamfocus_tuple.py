#!/usr/bin/env python3
"""Small-file ME tuple BeamFocus probe; no produced flux files are required.

Defaults to two P9 ME5A files and five previously investigated energy bins.
The existing common-bin CSV optionally supplies P8-produced/legacy RMS values.
This is a weight-scale diagnostic, not a full-sample flux or closure test.
"""
import argparse
import csv
import json
import math
from pathlib import Path


DEFAULT_BINS=["12:1.8:1.9","12:12.5:13","-12:18.5:19","14:9.1:9.2","-14:5.4:5.5"]


def parse_bins(specs):
    bins=[]
    for spec in specs:
        pdg,low,high=spec.split(":")
        pdg,low,high=int(pdg),float(low),float(high)
        if pdg not in [12,-12,14,-14] or not math.isfinite(low+high) or not 0<=low<high:
            raise ValueError("Invalid PDG:LOW:HIGH bin: "+spec)
        item=dict(pdg=pdg,low=low,high=high)
        if any(b["pdg"]==pdg and low<b["high"] and high>b["low"] for b in bins):
            raise ValueError("Overlapping bins for the same species")
        bins.append(item)
    return bins


def decode_scale(leaf_type, override):
    if override is not None:
        if not math.isfinite(override) or override<=0:
            raise ValueError("Weight scale must be finite and positive")
        return override
    if leaf_type in ["Int_t","UInt_t","int","unsigned int","Long64_t","ULong64_t"]:
        return 1e-7
    if leaf_type in ["Float_t","Double_t","float","double"]:
        return 1.
    raise ValueError("Unrecognized BeamFocus leaf type {}; specify --weight-scale after checking the schema".format(leaf_type))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--playlist",default="ME5A")
    ap.add_argument("--production",choices=["p8","p9"],default="p9")
    ap.add_argument("--filelist",type=Path)
    ap.add_argument("--max-files",type=int,default=2,help="0 reads the complete specified list")
    ap.add_argument("--bins",nargs="+",default=DEFAULT_BINS,metavar="PDG:LOW:HIGH")
    ap.add_argument("--selection",choices=["producer","fv-only"],default="producer",
                    help="producer matches the previous P8 carbon-CC FV audit; fv-only removes interaction and target cuts")
    ap.add_argument("--weight-scale",type=float,help="Default: 1e-7 for integer branches, 1 for floating branches; printed and recorded")
    ap.add_argument("--n-universes",type=int,default=1000)
    ap.add_argument("--comparison-csv",type=Path,default=Path("me_flux_common_bin_comparison/weights_by_energy.csv"))
    ap.add_argument("--out",type=Path,default=Path("p9_ME5A_beamfocus_probe"))
    args=ap.parse_args()
    if not args.playlist.startswith(("ME1","ME5","ME6")) or args.max_files<0 or args.n_universes<1:
        ap.error("Require an ME playlist, max-files >=0 and positive universe count")
    try:
        bins=parse_bins(args.bins)
    except ValueError as error:
        ap.error(str(error))
    filelist=args.filelist or Path("/pnfs/minerva/persistent/DataPreservation/{}/FullDetector/{}_Playlist{}_MC.txt".format(args.production,args.production.upper(),args.playlist[2:]))
    paths=[s.strip() for s in filelist.read_text().splitlines() if s.strip() and not s.lstrip().startswith("#")]
    total_files=len(paths)
    if args.max_files:
        paths=paths[:args.max_files]
    if not paths:
        raise RuntimeError("Empty input file list")
    reference=[]
    if args.comparison_csv.is_file():
        with args.comparison_csv.open() as stream:
            reference=list(csv.DictReader(stream))
    import ROOT
    import numpy as np
    ROOT.gROOT.SetBatch(True)
    chain=ROOT.TChain("Truth")
    for path in paths:
        print("INPUT "+path,flush=True)
        if not chain.Add(path):
            raise RuntimeError("Cannot add "+path)
    if not chain.GetEntries():
        raise RuntimeError("No Truth entries")
    chain.LoadTree(0)
    chain.SetBranchStatus("*",0)
    needed=["mc_incomingE","mc_incoming","mc_vtx","mc_cvweight_total"]
    if args.selection=="producer":
        needed += ["mc_current","mc_intType","mc_targetZ"]
    for name in needed+["mc_wgt_Flux_BeamFocus","mc_hornCurrent_cvweight"]:
        if not chain.GetBranch(name):
            raise RuntimeError("Missing branch "+name)
    for name in needed:
        chain.SetBranchStatus(name,1)
    leaf=chain.GetLeaf("mc_wgt_Flux_BeamFocus")
    leaf_type=str(leaf.GetTypeName()) if leaf else "unknown"
    scale=decode_scale(leaf_type,args.weight_scale)
    print("BeamFocus type={}, decode scale={}, expected universes={}".format(leaf_type,scale,args.n_universes),flush=True)
    species=" || ".join("mc_incoming=={}".format(p) for p in sorted({b["pdg"] for b in bins}))
    cut=("({}) && mc_incomingE*0.001>={:.17g} && mc_incomingE*0.001<{:.17g}"
         " && mc_vtx[2]>6117 && mc_vtx[2]<8193 && abs(mc_vtx[0])<850.0"
         " && (abs(mc_vtx[1])<850.0/sqrt(3.0) || abs(mc_vtx[1])<1700.0/sqrt(3.0)-abs(mc_vtx[0]/sqrt(3.0)))").format(
             species,min(b["low"] for b in bins),max(b["high"] for b in bins))
    if args.selection=="producer":
        cut += " && mc_current==1 && mc_intType!=8 && mc_targetZ==6"
    if chain.Draw(">>beamfocus_probe_entries",cut,"goff")<0:
        raise RuntimeError("Truth selection failed")
    entries=ROOT.gDirectory.Get("beamfocus_probe_entries")
    if not entries:
        raise RuntimeError("Missing selected-entry list")
    for name in ["mc_wgt_Flux_BeamFocus","mc_hornCurrent_cvweight"]:
        chain.SetBranchStatus(name,1)
    for b in bins:
        b.update(n=0,skipped=0,sumcv=0.,sumcv2=0.,sums=np.zeros(args.n_universes),first=None,max_event_vector_difference=0.)
    events=[]
    for i in range(entries.GetN()):
        entry=int(entries.GetEntry(i))
        if chain.GetEntry(entry)<=0:
            raise RuntimeError("Could not read selected entry")
        energy=float(chain.mc_incomingE)*.001;pdg=int(chain.mc_incoming)
        b=next((b for b in bins if b["pdg"]==pdg and b["low"]<=energy<b["high"]),None)
        if b is None:
            continue
        cv=float(chain.mc_cvweight_total);horn=float(chain.mc_hornCurrent_cvweight)
        if not math.isfinite(cv) or cv<1e-6 or not math.isfinite(horn) or abs(horn)<1e-12:
            b["skipped"]+=1
            continue
        raw=np.asarray(chain.mc_wgt_Flux_BeamFocus,dtype=float)
        if raw.ndim!=1 or len(raw)!=args.n_universes:
            raise RuntimeError("Unexpected BeamFocus array shape {} at entry {}; expected {}".format(raw.shape,entry,args.n_universes))
        weights=raw*scale/horn
        if not np.all(np.isfinite(weights)):
            raise RuntimeError("Nonfinite weights")
        if b["first"] is None:
            b["first"]=weights.copy()
        b["max_event_vector_difference"]=max(b["max_event_vector_difference"],float(np.max(np.abs(weights-b["first"]))))
        b["sums"]+=cv*weights;b["sumcv"]+=cv;b["sumcv2"]+=cv*cv;b["n"]+=1
        events.append(dict(pdg=pdg,energy_low=b["low"],energy_high=b["high"],entry=entry,energy_GeV=energy,
            total_cv=cv,horn_cv=horn,relative_mean=float(weights.mean()),relative_RMS_about_1=float(np.sqrt(np.mean((weights-1)**2))),
            relative_min=float(weights.min()),relative_max=float(weights.max())))
    args.out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for b in bins:
        key="gen2thin-pdg{}-minervame{}".format(b["pdg"],args.playlist[2:])
        ref=next((r for r in reference if r["case"]==key and r["band"]=="Flux_BeamFocus"
                  and abs(float(r["energy_low"])-b["low"])<1e-8 and abs(float(r["energy_high"])-b["high"])<1e-8),None)
        row=dict(production=args.production,playlist=args.playlist,pdg=b["pdg"],energy_low=b["low"],energy_high=b["high"],
            events=b["n"],skipped=b["skipped"],effective_events=b["sumcv"]**2/b["sumcv2"] if b["n"] else None,
            tuple_RMS_percent=None,tuple_mean=None,tuple_min=None,tuple_max=None,
            max_event_vector_difference=b["max_event_vector_difference"] if b["n"] else None,
            P8_produced_RMS_percent=100*float(ref["produced_weight_rms"]) if ref else None,
            legacy_RMS_percent=100*float(ref["reference_weight_rms"]) if ref else None)
        if b["n"]:
            w=b["sums"]/b["sumcv"]
            row.update(tuple_RMS_percent=100*float(np.sqrt(np.mean((w-1)**2))),tuple_mean=float(w.mean()),tuple_min=float(w.min()),tuple_max=float(w.max()))
            np.savetxt(args.out/"pdg{}_E{:g}_{:g}_weights.csv".format(b["pdg"],b["low"],b["high"]),
                       np.column_stack([np.arange(args.n_universes),w]),delimiter=",",header="universe,tuple_native",comments="")
        rows.append(row)
        print("PDG {} {:.4g}--{:.4g} GeV: N={} | {} RMS={} % | P8-produced={} % | legacy={} %".format(
            b["pdg"],b["low"],b["high"],b["n"],args.production,row["tuple_RMS_percent"],row["P8_produced_RMS_percent"],row["legacy_RMS_percent"]),flush=True)
        if not b["n"]:
            print("  No usable events: increase --max-files for this bin; absence is not evidence of small weights.")
    for name,records in [("summary.csv",rows),("events.csv",events)]:
        if records:
            with (args.out/name).open("w",newline="") as stream:
                writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    metadata=dict(production=args.production,playlist=args.playlist,filelist=str(filelist),files_used=paths,total_files_in_list=total_files,
        selection=cut,selection_mode=args.selection,leaf_type=leaf_type,decode_scale=scale,n_universes=args.n_universes,
        comparison_csv=str(args.comparison_csv),note="Small-file native-weight scale probe; no flux production or full-file closure claim. No universe offset is needed for RMS. Integer conversion is inferred from branch storage type and must be checked if P9 changed its encoding convention.")
    (args.out/"inputs.json").write_text(json.dumps(metadata,indent=2))
    print("Saved diagnostics to "+str(args.out))


if __name__=="__main__":
    main()
