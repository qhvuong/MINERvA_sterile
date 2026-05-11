#!/usr/bin/env python

import sys
import math
import ROOT
import PlotUtils

ROOT.TH1.AddDirectory(False)

DEFAULT_FILE = "/exp/minerva/data/users/qvuong/nu_e/kin_dist_mcleFHC_CCnue_allSystematics_newFlux_MAD.root"

def is_mnvh1d(obj):
    return obj and obj.InheritsFrom("PlotUtils::MnvH1D")

def is_mnvh2d(obj):
    return obj and obj.InheritsFrom("PlotUtils::MnvH2D")

def get_hist_names(f, include_2d=False, name_filter=None):
    names = []

    for key in f.GetListOfKeys():
        name = key.GetName()
        obj = key.ReadObj()

        if name_filter and name_filter not in name:
            continue

        if is_mnvh1d(obj):
            names.append(name)
        elif include_2d and is_mnvh2d(obj):
            names.append(name)

    return sorted(names)

def summarize_band_inventory(h):
    rows = []

    for bandname in h.GetVertErrorBandNames():
        bandname = str(bandname)
        band = h.GetVertErrorBand(bandname)
        if not band:
            continue

        n = band.GetNHists()
        bcv = band.Integral()
        cv = h.Integral()
        ratio = bcv / cv if cv else 0.0

        rows.append((bandname, n, ratio))

    return rows

FOCUS_BANDS = [
    "Target_Mass_CH",
    "elE_ECAL",
    "elE_HCAL",
    "elE_Tracker",
    "response_em",
    "response_p",
    "response_meson",
    "response_other",
    "response_xtalk",

    # A few GENIE / model dials that looked reasonably centered
    "GENIE_D2_MaRES",
    "GENIE_D2_NormCCRES",
    "GENIE_MFP_pi",
    "RPA_LowQ2",
    "RPA_HighQ2",
]

KEY_HISTS = [
    "EN4",
    "EN4_CCNuEQE",
    "EN4_CCNuEDelta",
    "EN4_NCPi0",
    "EN4_CCNuEDIS",
    "EN4_CCPi",
    "EN4_CCPi0",
    "EN4_NCPi",
    "EN4_NCDiff",
    "EN4_NCCohPi0",
    "EN4_NCOther",
    "EN4_NuEElastic",
]


def check_two_universe_band(h, bandname, frac_offset_warn=0.05):
    band = h.GetVertErrorBand(bandname)
    if not band or band.GetNHists() != 2:
        return None

    u0 = band.GetHist(0)
    u1 = band.GetHist(1)

    nbins = h.GetNbinsX()

    n_nonzero = 0
    n_same_side = 0
    n_large_offset = 0
    max_abs_frac_offset = 0.0
    max_offset_over_halfspread = 0.0

    for b in range(1, nbins + 1):
        cv = h.GetBinContent(b)
        v0 = u0.GetBinContent(b)
        v1 = u1.GetBinContent(b)

        if cv == 0 and v0 == 0 and v1 == 0:
            continue

        n_nonzero += 1

        midpoint = 0.5 * (v0 + v1)
        halfspread = 0.5 * abs(v0 - v1)
        offset = midpoint - cv

        frac_offset = offset / cv if cv else 0.0
        offset_over_halfspread = abs(offset) / halfspread if halfspread else float("inf")
        same_side = ((v0 - cv) * (v1 - cv)) > 0

        if same_side:
            n_same_side += 1
        if abs(frac_offset) > frac_offset_warn:
            n_large_offset += 1

        max_abs_frac_offset = max(max_abs_frac_offset, abs(frac_offset))
        if math.isfinite(offset_over_halfspread):
            max_offset_over_halfspread = max(max_offset_over_halfspread, offset_over_halfspread)

    same_side_frac = n_same_side / float(n_nonzero) if n_nonzero else 0.0

    return {
        "band": bandname,
        "n_nonzero": n_nonzero,
        "n_same_side": n_same_side,
        "same_side_frac": same_side_frac,
        "n_large_offset": n_large_offset,
        "max_abs_frac_offset": max_abs_frac_offset,
        "max_offset_over_halfspread": max_offset_over_halfspread,
    }

def analyze_focused_bands(h, hname):
    print("\n" + "=" * 120)
    print("Focused ±1σ dial check:", hname)
    print("Integral:", h.Integral())

    for bandname in FOCUS_BANDS:
        if not h.HasVertErrorBand(bandname):
            print(f"  [missing] {bandname}")
            continue

        band = h.GetVertErrorBand(bandname)
        if band.GetNHists() != 2:
            print(f"  [skip] {bandname}: nUniverses={band.GetNHists()}, not 2")
            continue

        result = check_two_universe_band(h, bandname, frac_offset_warn=0.05)

        status = "OK"
        if result["max_abs_frac_offset"] > 0.05:
            status = "CHECK"
        elif result["max_abs_frac_offset"] > 0.02:
            status = "WARN"

        print(
            "  {:25s} {:5s}  "
            "same_side={:2d}/{:2d} ({:.3f})  "
            "large_offset_bins={:2d}  "
            "max|offset/CV|={:.5g}  "
            "max|offset|/halfspread={:.5g}".format(
                bandname,
                status,
                result["n_same_side"],
                result["n_nonzero"],
                result["same_side_frac"],
                result["n_large_offset"],
                result["max_abs_frac_offset"],
                result["max_offset_over_halfspread"],
            )
        )

def print_detailed_two_universe_check(h, bandname):
    band = h.GetVertErrorBand(bandname)
    if not band or band.GetNHists() != 2:
        return

    u0 = band.GetHist(0)
    u1 = band.GetHist(1)

    print("\n" + "-" * 100)
    print("Detailed two-universe check:", h.GetName(), bandname)

    for b in range(1, h.GetNbinsX() + 1):
        cv = h.GetBinContent(b)
        v0 = u0.GetBinContent(b)
        v1 = u1.GetBinContent(b)

        if cv == 0 and v0 == 0 and v1 == 0:
            continue

        midpoint = 0.5 * (v0 + v1)
        halfspread = 0.5 * abs(v0 - v1)
        offset = midpoint - cv

        frac_offset = offset / cv if cv else 0.0
        frac_halfspread = halfspread / cv if cv else 0.0
        offset_over_halfspread = abs(offset) / halfspread if halfspread else float("inf")
        same_side = ((v0 - cv) * (v1 - cv)) > 0

        print(
            "bin {:2d} x=[{:.3f},{:.3f}] "
            "CV={:.6g} u0={:.6g} u1={:.6g} "
            "mid/CV={:.6g} halfspread/CV={:.6g} "
            "offset/CV={:.6g} |offset|/halfspread={:.6g} same_side={}"
            .format(
                b,
                h.GetXaxis().GetBinLowEdge(b),
                h.GetXaxis().GetBinUpEdge(b),
                cv,
                v0,
                v1,
                midpoint / cv if cv else 0.0,
                frac_halfspread,
                frac_offset,
                offset_over_halfspread,
                same_side,
            )
        )

def analyze_hist(f, hname, detailed_bands=None):
    h = f.Get(hname)
    if not h or not is_mnvh1d(h):
        return []

    print("\n" + "=" * 120)
    print("Histogram:", hname)
    print("Title:", h.GetTitle())
    print("Integral:", h.Integral())

    print("\nBand inventory:")
    for bandname, n, ratio in summarize_band_inventory(h):
        kind = "candidate ±1σ" if n == 2 else ("many-universe" if n > 2 else "other")
        print("  {:30s} n={:3d} bandCV/main={:.8g} {}".format(bandname, n, ratio, kind))

    flagged = []

    print("\nTwo-universe summary:")
    for bandname, n, _ in summarize_band_inventory(h):
        if n != 2:
            continue

        result = check_two_universe_band(h, bandname)
        if not result:
            continue

        status = "OK"
        if result["same_side_frac"] > 0.5 or result["n_large_offset"] > 0:
            status = "CHECK"
            flagged.append((hname, result))

        print(
            "  {:30s} status={:5s} "
            "same_side={:2d}/{:2d} ({:.3f}) "
            "large_offset_bins={:2d} "
            "max|offset/CV|={:.5g} "
            "max|offset|/halfspread={:.5g}"
            .format(
                bandname,
                status,
                result["n_same_side"],
                result["n_nonzero"],
                result["same_side_frac"],
                result["n_large_offset"],
                result["max_abs_frac_offset"],
                result["max_offset_over_halfspread"],
            )
        )

    if detailed_bands:
        for bandname in detailed_bands:
            if h.HasVertErrorBand(bandname):
                print_detailed_two_universe_check(h, bandname)

    return flagged

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_FILE

    # Optional args:
    #   --filter EN4
    #   --detail electron_scale,LowQ2Pi
    #   --include-2d
    name_filter = None
    detailed_bands = []
    include_2d = False

    args = sys.argv[2:]
    i = 0
    while i < len(args):
        if args[i] == "--filter":
            name_filter = args[i + 1]
            i += 2
        elif args[i] == "--detail":
            detailed_bands = [x.strip() for x in args[i + 1].split(",") if x.strip()]
            i += 2
        elif args[i] == "--include-2d":
            include_2d = True
            i += 1
        else:
            i += 1

    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        raise RuntimeError("Could not open file: {}".format(path))

    hist_names = get_hist_names(f, include_2d=include_2d, name_filter=name_filter)

    print("Opened:", path)
    print("Found {} histograms matching filter={!r}".format(len(hist_names), name_filter))

    all_flagged = []

    # for hname in hist_names:
    #     # This version only analyzes MnvH1D in detail.
    #     h = f.Get(hname)
    #     if not is_mnvh1d(h):
    #         continue

    #     flagged = analyze_hist(f, hname, detailed_bands=detailed_bands)
    #     all_flagged.extend(flagged)

    for hname in KEY_HISTS:
        h = f.Get(hname)
        if not is_mnvh1d(h):
            print("[missing or not MnvH1D]", hname)
            continue

        analyze_focused_bands(h, hname)

        if detailed_bands:
            for bandname in detailed_bands:
                if h.HasVertErrorBand(bandname):
                    print_detailed_two_universe_check(h, bandname)

    print("\n" + "#" * 120)
    print("SUMMARY OF FLAGGED TWO-UNIVERSE BANDS")
    print("#" * 120)

    if not all_flagged:
        print("No flagged two-universe bands.")
    else:
        for hname, r in all_flagged:
            print(
                "{:35s} {:30s} "
                "same_side={:2d}/{:2d} ({:.3f}) "
                "large_offset_bins={:2d} "
                "max|offset/CV|={:.5g} "
                "max|offset|/halfspread={:.5g}"
                .format(
                    hname,
                    r["band"],
                    r["n_same_side"],
                    r["n_nonzero"],
                    r["same_side_frac"],
                    r["n_large_offset"],
                    r["max_abs_frac_offset"],
                    r["max_offset_over_halfspread"],
                )
            )

    f.Close()

if __name__ == "__main__":
    main()