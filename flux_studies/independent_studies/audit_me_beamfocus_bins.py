#!/usr/bin/env python3
"""Audit large BeamFocus discrepancies directly in P8 Truth tuples.

Run one playlist at a time. For each species, select the largest discrepant
energy bins plus its reference flux-peak control bin; read the selected events
once and compare the CV-weighted tuple response with the produced flux band.
"""
import argparse
from bisect import bisect_right
import csv
import json
import math
from pathlib import Path


def target_at_energy(targets, lower_edges, energy):
    """Find a requested [lo, hi) bin; gaps in a sparse audit are excluded."""
    if not math.isfinite(energy):
        return None
    index = bisect_right(lower_edges, energy)-1
    if index >= 0 and energy < targets[index]["hi"]:
        return targets[index]
    return None


def energy_envelope_cut(targets):
    """Keep ROOT's formula size constant, regardless of the number of bins."""
    return "mc_incomingE*0.001>={:.17g} && mc_incomingE*0.001<{:.17g}".format(
        targets[0]["lo"], targets[-1]["hi"])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--playlist", required=True)
    ap.add_argument("--pdgs", nargs="+", type=int, default=[12, -12, 14, -14])
    ap.add_argument("--comparison-dir", type=Path, default=Path("me_flux_common_bin_comparison"))
    ap.add_argument("--flux-dir", type=Path, default=Path("/exp/minerva/data/users/qvuong/flux_studies/producedFluxes_p8_ME_offsetPPFX"))
    ap.add_argument("--top-bins", type=int, default=3, help="Largest flagged bins per species; 0 audits all flagged bins")
    ap.add_argument("--all-bins", action="store_true", help="Audit every positive-CV common bin in the energy range, including bins below the peak-fraction threshold")
    ap.add_argument("--selection", choices=["producer", "fv-only"], default="producer",
                    help="producer: carbon CC plus FV; fv-only: truth FV without interaction/target cuts. Both require finite usable CV weights.")
    ap.add_argument("--min-rms-difference-pp", type=float, default=2.)
    ap.add_argument("--min-spread-factor", type=float, default=1.5)
    ap.add_argument("--min-peak-fraction", type=float, default=.01)
    ap.add_argument("--emin", type=float, default=.5)
    ap.add_argument("--emax", type=float, default=20.)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    if not args.playlist.startswith(("ME1", "ME5", "ME6")):
        ap.error("Expected an ME1*, ME5*, or ME6* playlist")
    if args.top_bins < 0 or args.emax <= args.emin:
        ap.error("Require top-bins >= 0 and emax > emin")
    species_names = {12:"nue", -12:"nuebar", 14:"numu", -14:"numubar"}
    if any(pdg not in species_names for pdg in args.pdgs):
        ap.error("PDGs must be one of 12, -12, 14, -14")
    out = args.out or Path(args.playlist+"_beamfocus_bin_audit")
    out.mkdir(parents=True, exist_ok=True)
    with (args.comparison_dir/"weights_by_energy.csv").open() as f:
        comparison = list(csv.DictReader(f))
    filelist = Path("../file_option/playlist_mcme{}p8.txt".format(args.playlist[2:]))
    paths = [x.strip() for x in filelist.read_text().splitlines() if x.strip() and not x.lstrip().startswith("#")]
    if not paths:
        raise RuntimeError("Empty filelist")
    offset = 509 if args.playlist.startswith("ME1") else 811

    import ROOT
    import PlotUtils.LoadPlotUtilsLib
    import numpy as np
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    summaries, event_rows, failed = [], [], []
    manifest = dict(playlist=args.playlist, filelist=str(filelist), files=paths, offset=offset,
        selection=args.selection+"; producer truth FV; total CV >= 1e-6; selected energy bins",
        comparison_dir=str(args.comparison_dir), flux_dir=str(args.flux_dir),
        thresholds=dict(emin=args.emin,emax=args.emax,min_peak_fraction=args.min_peak_fraction,
            min_rms_difference_pp=args.min_rms_difference_pp,min_spread_factor=args.min_spread_factor,
            top_bins=args.top_bins,all_bins=args.all_bins),
        note="Flags prioritize investigation; they do not establish an error. Each species also has a peak-flux control bin. Histogram RMS means sqrt(mean((universe/CV-1)^2)). Effective event count uses CV weights only. No input files are modified.")
    for pdg in args.pdgs:
        case = "gen2thin-pdg{}-minervame{}".format(pdg,args.playlist[2:])
        source_rows = [r for r in comparison if r["case"]==case and r["band"]=="Flux_BeamFocus"]
        if not source_rows:
            raise RuntimeError("No comparison rows for "+case)
        peak_p = max(float(r["produced_cv"]) for r in source_rows)
        peak_r = max(float(r["reference_cv"]) for r in source_rows)
        eligible, flagged = [], []
        for r in source_rows:
            lo, hi = float(r["energy_low"]), float(r["energy_high"])
            pr, rr = float(r["produced_weight_rms"]),float(r["reference_weight_rms"])
            if (lo < args.emin-1e-8 or hi > args.emax+1e-8 or
                float(r["produced_cv"]) <= 0 or float(r["reference_cv"]) <= 0 or
                (not args.all_bins and (float(r["produced_cv"]) < args.min_peak_fraction*peak_p or
                float(r["reference_cv"]) < args.min_peak_fraction*peak_r)) or
                not np.isfinite(pr+rr) or rr<=0):
                continue
            item=dict(lo=lo,hi=hi,produced_rms=pr,reference_rms=rr,
                      delta_pp=100*(pr-rr),spread_factor=pr/rr,reference_cv=float(r["reference_cv"]))
            eligible.append(item)
            factor=max(pr/rr, rr/pr) if pr>0 else float("inf")
            if abs(item["delta_pp"])>=args.min_rms_difference_pp and factor>=args.min_spread_factor:
                flagged.append(item)
        if not eligible:
            raise RuntimeError("No eligible energy bins for "+case)
        ranked=sorted(flagged,key=lambda x:abs(x["delta_pp"]),reverse=True)
        targets=([dict(x,reason="full_energy_scan") for x in eligible] if args.all_bins else
                 [dict(x,reason="large_discrepancy") for x in (ranked[:args.top_bins] if args.top_bins else ranked)])
        control=max(eligible,key=lambda x:x["reference_cv"])
        if any(abs(x["lo"]-control["lo"])<1e-8 for x in targets):
            for t in targets:
                if abs(t["lo"]-control["lo"])<1e-8:
                    t["reason"]+="+reference_peak"
        else:
            targets.append(dict(control,reason="reference_peak_control"))
        targets.sort(key=lambda x:x["lo"])
        print("{} PDG {}: {} flagged bins; auditing {} bins using all {} files".format(
            args.playlist,pdg,len(flagged),len(targets),len(paths)),flush=True)
        fpath=args.flux_dir/(args.playlist+"_"+species_names[pdg]+".root")
        f=ROOT.TFile.Open(str(fpath),"READ")
        if not f or f.IsZombie():
            raise RuntimeError("Cannot open "+str(fpath))
        try:
            h=f.Get("flux_E_cvweighted")
            if not h or not h.HasVertErrorBand("Flux_BeamFocus"):
                raise RuntimeError("Missing CV or BeamFocus band")
            band=h.GetVertErrorBand("Flux_BeamFocus")
            if band.GetNHists()!=1000:
                raise RuntimeError("Expected 1000 universes")
            for t in targets:
                ib=h.FindBin((t["lo"]+t["hi"])/2)
                if not (np.isclose(h.GetXaxis().GetBinLowEdge(ib),t["lo"],rtol=0,atol=1e-8) and
                        np.isclose(h.GetXaxis().GetBinUpEdge(ib),t["hi"],rtol=0,atol=1e-8)):
                    raise RuntimeError("Target is not one native produced bin")
                if h.GetBinContent(ib)<=0:
                    raise RuntimeError("Nonpositive CV in target bin")
                t["saved"]=np.array([band.GetHist(u).GetBinContent(ib)/h.GetBinContent(ib) for u in range(1000)])
                t.update(sums=np.zeros(1000),sumw=0.,sumw2=0.,n=0,skipped=0,fallbacks=0,
                    maxw=0.,maxvec=None,firstvec=None,max_event_vector_difference=0.,parents={})
        finally:
            f.Close()
        chain=ROOT.TChain("Truth")
        for path in paths:
            chain.Add(path)
        if chain.GetEntries()==0:
            raise RuntimeError("Empty Truth chain")
        chain.LoadTree(0)
        chain.SetBranchStatus("*",0)
        for name in ("mc_incomingE","mc_incoming","mc_current","mc_intType","mc_targetZ","mc_vtx","mc_cvweight_total"):
            if not chain.GetBranch(name):
                raise RuntimeError("Missing "+name)
            chain.SetBranchStatus(name,1)
        leaf=chain.GetLeaf("mc_wgt_Flux_BeamFocus")
        if not leaf or str(leaf.GetTypeName()) not in ("Int_t","int"):
            raise RuntimeError("Expected P8 integer BeamFocus weights")
        if not chain.GetBranch("mc_hornCurrent_cvweight"):
            raise RuntimeError("Missing BeamFocus CV")
        lower_edges=[t["lo"] for t in targets]
        energy_cut=energy_envelope_cut(targets)
        cut=("mc_incoming=={}"
            " && mc_vtx[2]>6117 && mc_vtx[2]<8193 && abs(mc_vtx[0])<850.0"
            " && (abs(mc_vtx[1])<850.0/sqrt(3.0) || abs(mc_vtx[1])<1700.0/sqrt(3.0)-abs(mc_vtx[0]/sqrt(3.0)))"
            " && ({})").format(pdg,energy_cut)
        if args.selection=="producer":
            cut += " && mc_current==1 && mc_intType!=8 && mc_targetZ==6"
        manifest.setdefault("selection_expressions",{})[str(pdg)]=cut
        manifest["energy_selection_note"]="ROOT uses the enclosing energy range; exact requested half-open bins and gaps are checked in Python."
        listname="audit_bf_entries_{}".format(str(pdg).replace("-","minus"))
        if chain.Draw(">>"+listname,cut,"goff")<0:
            raise RuntimeError("Selection failed")
        entries=ROOT.gDirectory.Get(listname)
        if not entries:
            raise RuntimeError("Missing selected entry list")
        chain.SetBranchStatus("mc_wgt_Flux_BeamFocus",1)
        chain.SetBranchStatus("mc_hornCurrent_cvweight",1)
        has_parent=bool(chain.GetBranch("mc_fr_nuParentID"))
        if has_parent:
            chain.SetBranchStatus("mc_fr_nuParentID",1)
        for i in range(entries.GetN()):
            entry=int(entries.GetEntry(i))
            if chain.GetEntry(entry)<=0:
                raise RuntimeError("Could not read entry")
            energy=float(chain.mc_incomingE)*.001
            t=target_at_energy(targets,lower_edges,energy)
            if t is None:
                # The enclosing ROOT range can contain gaps between requested
                # bins. Such events do not enter any audit accumulation.
                continue
            cv=float(chain.mc_cvweight_total)
            if math.isnan(cv) or cv<1e-6:
                t["skipped"]+=1
                continue
            denom=float(chain.mc_hornCurrent_cvweight)
            component=denom
            if math.isnan(denom) or abs(denom)<1e-12:
                denom=1.
                t["fallbacks"]+=1
            try:
                w=np.asarray(chain.mc_wgt_Flux_BeamFocus,dtype=float)
                if w.ndim!=1 or len(w)<1000:
                    raise ValueError("Unexpected BeamFocus array shape")
                w=w[:1000]*1e-7/denom
            except (TypeError,ValueError):
                w=np.fromiter((int(chain.mc_wgt_Flux_BeamFocus[u]) for u in range(1000)),dtype=float,count=1000)*1e-7/denom
            if not math.isfinite(cv) or not np.all(np.isfinite(w)):
                raise RuntimeError("Nonfinite weights")
            if t["firstvec"] is None:
                t["firstvec"]=w.copy()
            t["max_event_vector_difference"]=max(t["max_event_vector_difference"],float(np.max(np.abs(w-t["firstvec"]))))
            t["sums"]+=cv*w
            t["sumw"]+=cv
            t["sumw2"]+=cv*cv
            t["n"]+=1
            if cv>t["maxw"]:
                t["maxw"],t["maxvec"]=cv,w.copy()
            parent=int(chain.mc_fr_nuParentID) if has_parent else "unknown"
            t["parents"][str(parent)]=t["parents"].get(str(parent),0)+1
            event_rows.append(dict(playlist=args.playlist,pdg=pdg,energy_low=t["lo"],energy_high=t["hi"],
                chain_entry=entry,energy_GeV=energy,parent_pdg=parent,total_cv=cv,beamfocus_cv=component,
                relative_mean=float(w.mean()),relative_std=float(w.std()),
                relative_RMS_about_1=float(np.sqrt(np.mean((w-1)**2))),relative_min=float(w.min()),relative_max=float(w.max())))
        for t in targets:
            if t["n"]==0:
                failed.append(dict(pdg=pdg,energy_low=t["lo"],energy_high=t["hi"],error="No usable selected events"))
                print("  {}--{} GeV: NO EVENTS".format(t["lo"],t["hi"]),flush=True)
                continue
            native=t["sums"]/t["sumw"]
            shifted=np.roll(native,-offset)
            rms=lambda x:float(np.sqrt(np.mean((x-1)**2)))
            agrees=bool(np.allclose(shifted,t["saved"],rtol=1e-7,atol=1e-9))
            # A changed event selection is a robustness comparison, not closure.
            closure=agrees if args.selection=="producer" else "not_applicable_FV_only"
            leaveout=(t["sums"]-t["maxw"]*t["maxvec"])/(t["sumw"]-t["maxw"]) if t["sumw"]>t["maxw"] else None
            row=dict(playlist=args.playlist,pdg=pdg,energy_low=t["lo"],energy_high=t["hi"],reason=t["reason"],
                events=t["n"],effective_events=t["sumw"]**2/t["sumw2"],
                largest_event_CV_fraction=t["maxw"]/t["sumw"],skipped=t["skipped"],CV_fallbacks=t["fallbacks"],
                tuple_RMS_percent=100*rms(native),produced_RMS_percent=100*rms(t["saved"]),
                reference_RMS_percent=100*t["reference_rms"],
                RMS_without_largest_CV_event_percent=100*rms(leaveout) if leaveout is not None else float("nan"),
                max_absolute_weight_difference=float(np.max(np.abs(shifted-t["saved"]))),
                comparison_CSV_RMS_difference_percent=100*(rms(t["saved"])-t["produced_rms"]),
                closure_match=closure,selection=args.selection,agrees_with_produced=agrees,
                max_event_vector_difference=t["max_event_vector_difference"],
                tuple_q16=float(np.quantile(native,.16)),tuple_median=float(np.median(native)),tuple_q84=float(np.quantile(native,.84)),
                parent_event_counts=json.dumps(t["parents"],sort_keys=True))
            summaries.append(row)
            if args.selection=="producer" and not agrees:
                failed.append(dict(pdg=pdg,energy_low=t["lo"],energy_high=t["hi"],error="Tuple/produced mismatch"))
            np.savetxt(out/"pdg{}_E{:g}_{:g}_weights.csv".format(pdg,t["lo"],t["hi"]),
                np.column_stack([np.arange(1000),native,shifted,t["saved"]]),delimiter=",",
                header="universe,tuple_native,tuple_shifted,produced",comments="")
            print("  {:.4g}--{:.4g} GeV: N={}, Neff={:.2f}, largest event={:.1%}; RMS tuple/produced/reference={:.3f}/{:.3f}/{:.3f}%; {}".format(
                t["lo"],t["hi"],t["n"],row["effective_events"],row["largest_event_CV_fraction"],
                row["tuple_RMS_percent"],row["produced_RMS_percent"],row["reference_RMS_percent"],
                ("MATCH" if agrees else "DIFFERENT")+(" (FV robustness comparison)" if args.selection=="fv-only" else "")),flush=True)
        # Write results after every species so completed checks survive interruption.
        for filename,rows in (("summary.csv",summaries),("events.csv",event_rows)):
            if rows:
                with (out/filename).open("w",newline="") as stream:
                    writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        manifest["failures"]=failed
        (out/"inputs.json").write_text(json.dumps(manifest,indent=2))
        del chain
    print("Completed {} bin audits; {} failures. Outputs: {}".format(len(summaries),len(failed),out))
    return 1 if failed else 0


if __name__=="__main__":
    raise SystemExit(main())
