#!/usr/bin/env python3
"""Check native P8 BeamFocus weights against one produced flux energy bin.

Uses the producer's default carbon CC fiducial selection. The standalone
BeamFocus universe/CV ratio cancels the common bin-width/POT/xsec scaling.
No flux production or input modification is performed.
"""
import argparse
import csv
import json
import math
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--playlist", default="ME5A")
    ap.add_argument("--pdg", type=int, default=12, choices=[12, -12, 14, -14])
    ap.add_argument("--filelist", type=Path)
    ap.add_argument("--energy-low", type=float, default=2.1)
    ap.add_argument("--energy-high", type=float, default=2.2)
    ap.add_argument("--offset", type=int)
    ap.add_argument("--max-files", type=int, default=0, help="0=all files; a subset cannot establish full-file closure")
    ap.add_argument("--flux-dir", type=Path, default=Path("/exp/minerva/data/users/qvuong/flux_studies/producedFluxes_p8_ME_offsetPPFX"))
    ap.add_argument("--out", type=Path, default=Path("me5A_nue_tuple_beamfocus_check"))
    args = ap.parse_args()
    if args.energy_high <= args.energy_low or args.max_files < 0:
        ap.error("Require high > low and max-files >= 0")
    if not args.playlist.startswith(("ME1", "ME5", "ME6")):
        ap.error("Expected an ME1*, ME5*, or ME6* playlist")
    offset = args.offset if args.offset is not None else (509 if args.playlist.startswith("ME1") else 811)
    offset %= 1000
    # Actual naming is lower-case 'me', with the playlist's final letter upper-case.
    filelist = args.filelist or Path("../file_option/playlist_mcme{}p8.txt".format(args.playlist[2:]))
    paths = [s.strip() for s in filelist.read_text().splitlines() if s.strip() and not s.lstrip().startswith("#")]
    all_files = len(paths)
    if args.max_files:
        paths = paths[:args.max_files]
    if not paths:
        raise RuntimeError("Empty filelist")

    import ROOT
    import PlotUtils.LoadPlotUtilsLib
    import numpy as np
    ROOT.gROOT.SetBatch(True)
    chain = ROOT.TChain("Truth")
    for path in paths:
        if chain.Add(path) == 0:
            raise RuntimeError("Could not add "+path)
    if chain.GetEntries() == 0:
        raise RuntimeError("Truth chain is empty")
    chain.LoadTree(0)
    chain.SetBranchStatus("*", 0)
    for name in ("mc_incomingE", "mc_incoming", "mc_current", "mc_intType",
                 "mc_targetZ", "mc_vtx", "mc_cvweight_total"):
        if not chain.GetBranch(name):
            raise RuntimeError("Missing branch "+name)
        chain.SetBranchStatus(name, 1)
    if not chain.GetBranch("mc_wgt_Flux_BeamFocus") or not chain.GetBranch("mc_hornCurrent_cvweight"):
        raise RuntimeError("Missing native BeamFocus weight branches")
    leaf = chain.GetLeaf("mc_wgt_Flux_BeamFocus")
    if not leaf or str(leaf.GetTypeName()) not in ("Int_t", "int"):
        raise RuntimeError("Expected integer-scaled P8 BeamFocus branch; found {}".format(leaf.GetTypeName() if leaf else "no leaf"))
    cut = (
        "mc_incoming=={pdg} && mc_current==1 && mc_intType!=8 && mc_targetZ==6"
        " && mc_vtx[2]>6117 && mc_vtx[2]<8193"
        " && abs(mc_vtx[0])<850.0"
        " && (abs(mc_vtx[1])<850.0/sqrt(3.0)"
        " || abs(mc_vtx[1])<1700.0/sqrt(3.0)-abs(mc_vtx[0]/sqrt(3.0)))"
        " && mc_incomingE*0.001>={low:.17g} && mc_incomingE*0.001<{high:.17g}"
    ).format(pdg=args.pdg, low=args.energy_low, high=args.energy_high)
    print("Selecting {} of {} files, PDG {}, {:.4g}--{:.4g} GeV".format(len(paths), all_files, args.pdg,args.energy_low,args.energy_high), flush=True)
    if chain.Draw(">>bf_check_entries", cut, "goff") < 0:
        raise RuntimeError("Could not evaluate selection")
    entries = ROOT.gDirectory.Get("bf_check_entries")
    if not entries or not entries.GetN():
        raise RuntimeError("No selected events")
    chain.SetBranchStatus("mc_wgt_Flux_BeamFocus", 1)
    chain.SetBranchStatus("mc_hornCurrent_cvweight", 1)
    sums = np.zeros(1000)
    unweighted_sums = np.zeros(1000)
    total_cv = 0.
    events, skipped, fallback = [], 0, 0
    for i in range(entries.GetN()):
        if chain.GetEntry(entries.GetEntry(i)) <= 0:
            raise RuntimeError("Could not read selected entry")
        cv = float(chain.mc_cvweight_total)
        if math.isnan(cv) or cv < 1e-6:
            skipped += 1
            continue
        if not math.isfinite(cv):
            raise RuntimeError("Nonfinite positive total CV weight")
        component_cv = float(chain.mc_hornCurrent_cvweight)
        denominator = component_cv
        if math.isnan(denominator) or abs(denominator)<1e-12:
            denominator = 1.
            fallback += 1
        raw = np.fromiter((int(chain.mc_wgt_Flux_BeamFocus[u]) for u in range(1000)), dtype=float, count=1000)
        decoded = raw * 1e-7
        relative = decoded / denominator
        if not np.all(np.isfinite(relative)):
            raise RuntimeError("Nonfinite relative tuple weights")
        sums += cv * relative
        unweighted_sums += relative
        total_cv += cv
        events.append(dict(chain_entry=int(entries.GetEntry(i)), energy_GeV=float(chain.mc_incomingE)*.001,
            total_cv=cv, beamfocus_cv=component_cv, decoded_mean=float(decoded.mean()),
            relative_mean=float(relative.mean()), relative_std=float(relative.std()),
            relative_RMS_about_1=float(np.sqrt(np.mean((relative-1)**2))),
            relative_min=float(relative.min()), relative_max=float(relative.max())))
    if not events or total_cv<=0:
        raise RuntimeError("No usable events")
    native = sums / total_cv
    shifted = np.roll(native, -offset)
    plain = unweighted_sums / len(events)
    species = {12:"nue", -12:"nuebar", 14:"numu", -14:"numubar"}[args.pdg]
    flux_path = args.flux_dir/(args.playlist+"_"+species+".root")
    f = ROOT.TFile.Open(str(flux_path), "READ")
    if not f or f.IsZombie():
        raise RuntimeError("Cannot open "+str(flux_path))
    try:
        h = f.Get("flux_E_cvweighted")
        if not h or not h.HasVertErrorBand("Flux_BeamFocus"):
            raise RuntimeError("Missing produced flux or BeamFocus band")
        bin_number = h.GetXaxis().FindBin((args.energy_low+args.energy_high)/2)
        if not (np.isclose(h.GetXaxis().GetBinLowEdge(bin_number), args.energy_low,rtol=0,atol=1e-8)
                and np.isclose(h.GetXaxis().GetBinUpEdge(bin_number), args.energy_high,rtol=0,atol=1e-8)):
            raise RuntimeError("Requested range must equal exactly one produced flux bin")
        band = h.GetVertErrorBand("Flux_BeamFocus")
        if band.GetNHists()!=1000 or h.GetBinContent(bin_number)<=0:
            raise RuntimeError("Expected 1000 universes and positive bin CV")
        saved = np.array([band.GetHist(u).GetBinContent(bin_number)/h.GetBinContent(bin_number) for u in range(1000)])
    finally:
        f.Close()
    rms = lambda a: float(np.sqrt(np.mean((a-1)**2)))
    full = len(paths)==all_files
    match = bool(np.allclose(shifted, saved, rtol=1e-7, atol=1e-9))
    report = dict(playlist=args.playlist,pdg=args.pdg,filelist=str(filelist),files_used=paths,
        full_filelist=full,flux_file=str(flux_path),selection=cut,offset=offset,
        events_used=len(events),events_skipped_bad_total_cv=skipped,beamfocus_CV_fallbacks=fallback,
        tuple_aggregate_RMS_about_1=rms(native),shifted_tuple_RMS_about_1=rms(shifted),
        produced_RMS_about_1=rms(saved),tuple_aggregate_std=float(native.std()),
        produced_std=float(saved.std()),unweighted_tuple_aggregate_RMS_about_1=rms(plain),
        max_absolute_tuple_vs_produced_weight_difference=float(np.max(np.abs(shifted-saved))),
        tuple_vs_produced_allclose=match,
        interpretation="Full-file closure" if full else "Subset diagnostic only; not a full-file closure test")
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/"summary.json").write_text(json.dumps(report,indent=2))
    with (args.out/"events.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(events[0]));writer.writeheader();writer.writerows(events)
    np.savetxt(args.out/"universe_weights.csv",np.column_stack([np.arange(1000),native,shifted,saved]),
        delimiter=",",header="universe,tuple_native,tuple_shifted,produced",comments="")
    print("Events used: {}; BeamFocus CV fallbacks: {}".format(len(events),fallback))
    print("Tuple aggregate RMS: {:.6f}% | after offset: {:.6f}% | produced: {:.6f}%".format(100*rms(native),100*rms(shifted),100*rms(saved)))
    print("Maximum absolute weight difference: {:.6g}".format(report["max_absolute_tuple_vs_produced_weight_difference"]))
    print("{}: {}".format(report["interpretation"], "MATCH" if match else "DIFFERENT"))
    print("Saved diagnostics to",args.out)
    return 1 if full and not match else 0


if __name__ == "__main__":
    raise SystemExit(main())
