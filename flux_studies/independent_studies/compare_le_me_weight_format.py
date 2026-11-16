#!/usr/bin/env python3
"""Compare LE/ME tuple weight storage and CV conventions on small samples.

Reads both BeamFocus and PPFX. This does not produce flux histograms or assume
that their universe entries are absolute/relative. Integer decoding defaults
to the known P8/P9 1e-7 storage scale; floating weights default to scale 1.
Ensemble means near one or near the component CV are diagnostics, not proof
of the producer's intended convention. No flux files or PlotUtils are needed.
"""
import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

COMPONENTS={"BeamFocus":("mc_wgt_Flux_BeamFocus","mc_hornCurrent_cvweight"),
            "PPFX":("mc_wgt_ppfx1_Total","mc_ppfx1_cvweight")}
INTEGER_TYPES={"Int_t","UInt_t","int","unsigned int","Long64_t","ULong64_t"}
FLOAT_TYPES={"Float_t","Double_t","float","double"}


def scale_for_type(typename,int_scale,float_scale):
    if typename in INTEGER_TYPES:
        return int_scale
    if typename in FLOAT_TYPES:
        return float_scale
    raise ValueError("Unsupported array type: "+typename)


def summarize(rows,np):
    result=dict(events=len(rows),universes=sorted({r["universes"] for r in rows}),
                branch_types=sorted({r["branch_type"] for r in rows}))
    for column in ["component_cv","decoded_mean","decoded_RMS_about_1","divided_mean","divided_RMS_about_1"]:
        values=np.array([r[column] for r in rows],dtype=float)
        values=values[np.isfinite(values)]
        for name,function in [("min",np.min),("median",np.median),("max",np.max)]:
            result[column+"_"+name]=float(function(values)) if len(values) else None
    valid=[r for r in rows if math.isfinite(r["component_cv"]) and abs(r["component_cv"])>1e-12]
    informative=[r for r in valid if abs(r["component_cv"]-1)>1e-6]
    result["invalid_or_zero_component_CV_events"]=len(rows)-len(valid)
    result["nonunit_component_CV_events"]=len(informative)
    if informative:
        result["nonunit_CV_median_abs_decoded_mean_minus_1"]=float(np.median([abs(r["decoded_mean"]-1) for r in informative]))
        result["nonunit_CV_median_abs_divided_mean_minus_1"]=float(np.median([abs(r["divided_mean"]-1) for r in informative]))
    else:
        result["nonunit_CV_median_abs_decoded_mean_minus_1"]=None
        result["nonunit_CV_median_abs_divided_mean_minus_1"]=None
    result["interpretation"]=("CV is unity within tolerance: direct and CV-divided weights cannot distinguish the conventions here."
        if not informative else "Compare the two deviations from one on nonunit-CV events; closeness of ensemble means is diagnostic, not a definitive convention test.")
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample",action="append",required=True,metavar="LABEL=FILELIST",
                    help="Repeat for LE/ME lists, e.g. LE1_p8=../file_option/playlist_mcle1p8.txt")
    ap.add_argument("--max-files",type=int,default=2)
    ap.add_argument("--entries-per-file",type=int,default=10000,help="0 reads every entry in each selected file")
    ap.add_argument("--fv-only",action="store_true",help="Optional producer truth FV; default is no FV or interaction cuts")
    ap.add_argument("--int-scale",type=float,default=1e-7)
    ap.add_argument("--float-scale",type=float,default=1.)
    ap.add_argument("--out",type=Path,default=Path("le_me_weight_format_check"))
    args=ap.parse_args()
    if args.max_files<1 or args.entries_per_file<0 or not all(math.isfinite(x) and x>0 for x in [args.int_scale,args.float_scale]):
        ap.error("Require max-files >=1, entries-per-file >=0 and finite positive decoding scales")
    samples={}
    for item in args.sample:
        if "=" not in item:
            ap.error("Expected LABEL=FILELIST")
        label,name=item.split("=",1)
        if not label or label in samples:
            ap.error("Empty or duplicate sample label")
        samples[label]=Path(name)
    import ROOT
    import numpy as np
    ROOT.gROOT.SetBatch(True)
    args.out.mkdir(parents=True,exist_ok=True)
    groups=defaultdict(list);schemas=[];inputs=[];failures=[]
    for label,filelist in samples.items():
        paths=[s.strip() for s in filelist.read_text().splitlines() if s.strip() and not s.lstrip().startswith("#")][:args.max_files]
        if not paths:
            raise RuntimeError("Empty filelist: "+str(filelist))
        for path in paths:
            print("{} INPUT {}".format(label,path),flush=True)
            f=ROOT.TFile.Open(path,"READ")
            if not f or f.IsZombie():
                raise RuntimeError("Cannot open "+path)
            try:
                tree=f.Get("Truth")
                if not tree:
                    raise RuntimeError("Missing Truth tree in "+path)
                tree.SetBranchStatus("*",0)
                needed=["mc_incoming"]+(["mc_vtx"] if args.fv_only else [])
                for name in needed:
                    if not tree.GetBranch(name):
                        raise RuntimeError("Missing "+name)
                    tree.SetBranchStatus(name,1)
                active={}
                for component,(branch,cvbranch) in COMPONENTS.items():
                    if not tree.GetBranch(branch) or not tree.GetBranch(cvbranch):
                        failures.append(dict(sample=label,file=path,component=component,error="Missing weight array or component CV branch"))
                        continue
                    leaf=tree.GetLeaf(branch)
                    cvleaf=tree.GetLeaf(cvbranch)
                    typename=str(leaf.GetTypeName()) if leaf else "unknown"
                    scale=scale_for_type(typename,args.int_scale,args.float_scale)
                    active[component]=(branch,cvbranch,typename,scale)
                    tree.SetBranchStatus(branch,1);tree.SetBranchStatus(cvbranch,1)
                    schemas.append(dict(sample=label,file=path,component=component,branch=branch,branch_type=typename,
                        leaf_title=str(leaf.GetTitle()),static_length=int(leaf.GetLenStatic()),
                        component_CV_branch=cvbranch,component_CV_type=str(cvleaf.GetTypeName()) if cvleaf else "unknown",decode_scale=scale))
                    print("  {}: {} x {}; CV={}".format(component,typename,scale,cvbranch),flush=True)
                available=int(tree.GetEntries())
                count=min(available,args.entries_per_file) if args.entries_per_file else available
                inputs.append(dict(sample=label,filelist=str(filelist),file=path,available_entries=available,entries_scanned=count))
                for entry in range(count):
                    if tree.GetEntry(entry)<=0:
                        raise RuntimeError("Unreadable entry {} in {}".format(entry,path))
                    pdg=int(tree.mc_incoming)
                    if pdg not in [12,-12,14,-14]:
                        continue
                    if args.fv_only:
                        x,y,z=[float(tree.mc_vtx[i]) for i in range(3)]
                        if not (6117<z<8193 and abs(x)<850 and (abs(y)<850/math.sqrt(3) or abs(y)<1700/math.sqrt(3)-abs(x/math.sqrt(3)))):
                            continue
                    for component,(branch,cvbranch,typename,scale) in active.items():
                        raw=np.asarray(getattr(tree,branch),dtype=float)
                        if raw.ndim!=1 or not len(raw) or not np.all(np.isfinite(raw)):
                            raise RuntimeError("Invalid array in {} entry {}".format(branch,entry))
                        decoded=raw*scale;cv=float(getattr(tree,cvbranch))
                        divided=decoded/cv if math.isfinite(cv) and abs(cv)>1e-12 else None
                        groups[(label,pdg,component)].append(dict(sample=label,file=path,entry=entry,pdg=pdg,component=component,
                            branch_type=typename,universes=len(raw),component_cv=cv,raw_min=float(raw.min()),raw_max=float(raw.max()),
                            decoded_mean=float(decoded.mean()),decoded_RMS_about_1=float(np.sqrt(np.mean((decoded-1)**2))),
                            divided_mean=float(divided.mean()) if divided is not None else float("nan"),
                            divided_RMS_about_1=float(np.sqrt(np.mean((divided-1)**2))) if divided is not None else float("nan")))
            finally:
                f.Close()
    summaries=[]
    for (label,pdg,component),rows in sorted(groups.items()):
        summary=dict(sample=label,pdg=pdg,component=component,**summarize(rows,np));summaries.append(summary)
        print("{} PDG {} {}: N={} universes={} | CV min/median/max={:.6g}/{:.6g}/{:.6g}; nonunit={} | median mean: decoded={:.6g}, divided={:.6g}".format(
            label,pdg,component,len(rows),summary["universes"],summary["component_cv_min"],summary["component_cv_median"],summary["component_cv_max"],
            summary["nonunit_component_CV_events"],summary["decoded_mean_median"],summary["divided_mean_median"] or float("nan")),flush=True)
    def write_csv(name,rows):
        if rows:
            with (args.out/name).open("w",newline="") as stream:
                writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    write_csv("summary.csv",summaries);write_csv("schema.csv",schemas)
    write_csv("events.csv",[r for rows in groups.values() for r in rows])
    (args.out/"inputs.json").write_text(json.dumps(dict(inputs=inputs,failures=failures,FV_only=args.fv_only,
        int_scale=args.int_scale,float_scale=args.float_scale,note="Storage type/length and numerical normalization checks, without a relative/absolute assumption. No total-CV division or physics weighting. Universe means need not equal the CV exactly. Confirm semantics with tuple producer if numerical checks are inconclusive."),indent=2))
    print("Saved diagnostics to "+str(args.out))
    return 1 if failures else 0


if __name__=="__main__":
    raise SystemExit(main())
