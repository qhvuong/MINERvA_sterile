#!/usr/bin/env python

import sys
import ROOT
import PlotUtils
import math

ROOT.TH1.AddDirectory(False)

def zero_hist_bins(h):
    for i in range(h.GetSize()):
        h.SetBinContent(i, 0.0)
        h.SetBinError(i, 0.0)

def add_hist_bins_into(target, source):
    """
    Add source into target bin-by-bin, with errors in quadrature.
    Works for TH1D/TH2D-style global bin indexing.
    """
    for i in range(target.GetSize()):
        c = target.GetBinContent(i) + source.GetBinContent(i)
        e2 = target.GetBinError(i)**2 + source.GetBinError(i)**2
        target.SetBinContent(i, c)
        target.SetBinError(i, math.sqrt(e2))

def sync_band_cvs(hist):
    """
    Copy parent/main CV into every vertical-band CV.
    Does not touch universe histograms.
    """
    if hist.InheritsFrom("PlotUtils::MnvH1D"):
        cv = ROOT.TH1D(hist)
    elif hist.InheritsFrom("PlotUtils::MnvH2D"):
        cv = ROOT.TH2D(hist)
    else:
        return

    cv.SetDirectory(0)

    for bandname in hist.GetVertErrorBandNames():
        band = hist.GetVertErrorBand(str(bandname))
        if band:
            cv.Copy(band)

def repair_one_hist(out_hist, input_hists):
    """
    Replace out_hist contents with explicit sum of input_hists:
      parent CV
      each vertical-band CV
      each universe histogram
    """

    # Parent/main CV
    zero_hist_bins(out_hist)
    for h in input_hists:
        add_hist_bins_into(out_hist, h)

    # Vertical bands
    for bandname in out_hist.GetVertErrorBandNames():
        bandname = str(bandname)
        out_band = out_hist.GetVertErrorBand(bandname)
        if not out_band:
            continue

        input_bands = []
        for h in input_hists:
            if h.HasVertErrorBand(bandname):
                input_bands.append(h.GetVertErrorBand(bandname))

        if not input_bands:
            continue

        # Band CV
        zero_hist_bins(out_band)
        for b in input_bands:
            add_hist_bins_into(out_band, b)

        # Universe histograms
        n_univ = out_band.GetNHists()
        for iu in range(n_univ):
            out_u = out_band.GetHist(iu)
            zero_hist_bins(out_u)

            for b in input_bands:
                if iu < b.GetNHists():
                    add_hist_bins_into(out_u, b.GetHist(iu))

    # Final bookkeeping sync: band CV should equal parent CV
    sync_band_cvs(out_hist)

def repair_file(out_path, input_paths):
    print("[OPEN OUT]", out_path)
    fout = ROOT.TFile.Open(out_path, "UPDATE")
    if not fout or fout.IsZombie():
        raise RuntimeError("Could not open output file: {}".format(out_path))

    input_files = []
    for p in input_paths:
        print("[OPEN IN ]", p)
        f = ROOT.TFile.Open(p, "READ")
        if not f or f.IsZombie():
            raise RuntimeError("Could not open input file: {}".format(p))
        input_files.append(f)

    keys = fout.GetListOfKeys().Clone()
    n_repaired = 0
    n_skipped = 0

    for key in keys:
        name = key.GetName()
        out_obj = fout.Get(name)

        if not out_obj:
            n_skipped += 1
            continue

        is_mnv = (
            out_obj.InheritsFrom("PlotUtils::MnvH1D")
            or out_obj.InheritsFrom("PlotUtils::MnvH2D")
        )

        if not is_mnv:
            n_skipped += 1
            continue

        input_objs = []
        missing = False
        for f in input_files:
            obj = f.Get(name)
            if not obj:
                missing = True
                break
            input_objs.append(obj)

        if missing:
            print("[SKIP missing in at least one input]", name)
            n_skipped += 1
            continue

        print("[REPAIR]", name)
        repair_one_hist(out_obj, input_objs)

        fout.cd()
        out_obj.Write(name, ROOT.TObject.kOverwrite)
        n_repaired += 1

    for f in input_files:
        f.Close()

    fout.Close()

    print("[DONE] repaired:", n_repaired)
    print("[DONE] skipped :", n_skipped)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage:")
        print("  python repair_fhc_mnv_universes.py OUT_FHC.root IN1.root IN2.root IN3.root IN4.root")
        sys.exit(1)

    repair_file(sys.argv[1], sys.argv[2:])