#!/usr/bin/env python3
"""Inspect suspect P8 data files and matching P4 files without changing inputs.

Dumps full Meta records, top-level tree counts, and a few event identifiers.
Reads the suspect-file list from the POT inventory; never replaces exposure.
"""
import argparse
import json
import math
import re
from pathlib import Path

HERE=Path(__file__).resolve().parent


def leaf_value(leaf,limit=16):
    typename=str(leaf.GetTypeName())
    size=int(leaf.GetNdata())
    if typename in ["string","std::string"] or (typename in ["Char_t","char"] and size>1):
        return str(leaf.GetValueString())
    integer=typename in ["Int_t","UInt_t","Long64_t","ULong64_t","Short_t","UShort_t","Bool_t","UChar_t"]
    values=[]
    for index in range(min(size,limit)):
        value=int(leaf.GetValueLong64(index)) if integer else float(leaf.GetValue(index))
        values.append(value if integer or math.isfinite(value) else repr(value))
    if size==1:
        return values[0]
    return dict(values=values,length=size,truncated=size>limit)


def inspect_file(ROOT,path,max_meta_entries):
    f=ROOT.TFile.Open(path,"READ")
    if not f or f.IsZombie():
        return dict(path=path,error="Cannot open ROOT file")
    result=dict(path=path,recovered=bool(f.TestBit(ROOT.TFile.kRecovered)),objects=[],trees={})
    try:
        names=[]
        for key in f.GetListOfKeys():
            name=str(key.GetName())
            result["objects"].append(dict(name=name,class_name=str(key.GetClassName()),cycle=int(key.GetCycle())))
            if name not in names:
                names.append(name)
        for name in names:
            tree=f.Get(name)
            if not tree or not tree.InheritsFrom("TTree"):
                continue
            count=int(tree.GetEntries())
            info=dict(entries=count)
            if name=="Meta":
                leaves=list(tree.GetListOfLeaves())
                info["schema"]=[dict(name=str(l.GetName()),type=str(l.GetTypeName()),title=str(l.GetTitle())) for l in leaves]
                limit=count if max_meta_entries==0 else min(count,max_meta_entries)
                info["records_truncated"]=limit<count
                info["records"]=[]
                for entry in range(limit):
                    if tree.GetEntry(entry)<=0:
                        raise RuntimeError("Cannot read Meta entry "+str(entry))
                    row=dict(entry=entry)
                    for leaf in leaves:
                        try:
                            row[str(leaf.GetName())]=leaf_value(leaf)
                        except Exception as error:
                            row[str(leaf.GetName())]=dict(read_error=str(error))
                    info["records"].append(row)
            else:
                # Sample identifiers only; do not read the full event payload.
                tree.SetBranchStatus("*",0)
                leaves=[l for l in tree.GetListOfLeaves() if re.search(r"(^|_)(run|subrun|gate|spill|event|time|gps|is_mc)($|_)",str(l.GetName()),re.I)
                        and int(l.GetLenStatic())==1]
                leaves=leaves[:30]
                for leaf in leaves:
                    tree.SetBranchStatus(str(leaf.GetBranch().GetName()),1)
                info["sampled_identifiers"]=[]
                if leaves:
                    for entry in sorted({0,count//2,count-1}) if count else []:
                        tree.GetEntry(entry)
                        row=dict(entry=entry)
                        for leaf in leaves:
                            try:
                                row[str(leaf.GetName())]=leaf_value(leaf)
                            except Exception as error:
                                row[str(leaf.GetName())]=dict(read_error=str(error))
                        info["sampled_identifiers"].append(row)
            result["trees"][name]=info
        return result
    except Exception as error:
        result["error"]=str(error)
        return result
    finally:
        f.Close()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inventory",type=Path,required=True)
    ap.add_argument("--p4-list-pattern",default=str(HERE.parent/"file_option/playlist_datame{playlist}p4.txt"))
    ap.add_argument("--max-meta-entries",type=int,default=0,help="0 dumps every Meta entry")
    ap.add_argument("--out",type=Path,default=Path("me_zero_pot_inspection"))
    args=ap.parse_args()
    if args.max_meta_entries<0:
        ap.error("Require max-meta-entries >=0")
    manifest=json.loads(args.inventory.read_text())
    targets={}
    for row in manifest.get("warnings",[])+manifest.get("failures",[]):
        if "POT_Total" in row.get("field","")+row.get("error",""):
            targets[row["path"]]=row["playlist"]
    if not targets:
        raise RuntimeError("No suspect POT_Total files found in inventory")
    import ROOT
    ROOT.gROOT.SetBatch(True)
    args.out.mkdir(parents=True,exist_ok=True)
    results=[]
    for path,playlist in targets.items():
        print("Inspecting {} {}".format(playlist,Path(path).name),flush=True)
        p8=inspect_file(ROOT,path,args.max_meta_entries)
        otherlist=Path(args.p4_list_pattern.format(playlist=playlist[2:]))
        matches=[]
        if otherlist.is_file():
            matches=[s.strip() for s in otherlist.read_text().splitlines() if s.strip() and not s.lstrip().startswith("#") and Path(s.strip()).name==Path(path).name]
        p4=[inspect_file(ROOT,p,args.max_meta_entries) for p in matches]
        result=dict(playlist=playlist,p8=p8,p4_matches=p4,p4_list=str(otherlist))
        results.append(result)
        for label,info in [("P8",p8)]+[("P4",x) for x in p4]:
            print("  {} trees: {}".format(label,{k:v["entries"] for k,v in info.get("trees",{}).items()}),flush=True)
            for record in info.get("trees",{}).get("Meta",{}).get("records",[]):
                print("  {} Meta {}".format(label,json.dumps(record)),flush=True)
        if not matches:
            print("  No matching P4 basename in "+str(otherlist),flush=True)
        # Checkpoint each pair in case a later remote file cannot be read.
        (args.out/"file_inspection.json").write_text(json.dumps(dict(inventory=str(args.inventory),results=results,
            note="Read-only inspection. Matching basename is a cross-production diagnostic, not permission to substitute P4 POT. No exposure has been repaired or assigned."),indent=2))
    print("Saved "+str(args.out/"file_inspection.json"))


if __name__=="__main__":
    main()
