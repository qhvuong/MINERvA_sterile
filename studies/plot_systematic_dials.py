#!/usr/bin/env python3

import os
import sys
import math
import ROOT
import PlotUtils

ROOT.TH1.AddDirectory(False)

# DEFAULT_FILE = "/exp/minerva/data/users/qvuong/nu_e/kin_dist_mcleFHC_CCnue_allSystematics_newFlux_MAD.root"
DEFAULT_FILE = "/exp/minerva/data/users/qvuong/nu_e/kin_dist_mcle1_p6_test_MAD_0.root"
DEFAULT_HIST = "EN4"
DEFAULT_OUTDIR = "plots/systematic_dials"

# Start focused. You can replace this with all bands later.
DEFAULT_BANDS = [
    "RPA_LowQ2",
    "RPA_HighQ2",
    "Low_Recoil_2p2h_Tune",
    "LowQ2Pi",
    "fsi_weight",
    "SuSA_Valencia_Weight",
    "MK_model",
]

def safe_name(s):
    return (
        s.replace("/", "_")
         .replace(" ", "_")
         .replace("(", "")
         .replace(")", "")
         .replace(",", "_")
    )

def clone_th1(obj, name):
    h = ROOT.TH1D(obj)
    h.SetDirectory(0)
    h.SetName(name)
    h.SetStats(0)
    return h

def get_band_names(h):
    return [str(x) for x in h.GetVertErrorBandNames()]

def make_universe_stats(h, bandname):
    band = h.GetVertErrorBand(bandname)
    n = band.GetNHists()

    h_mean = ROOT.TH1D(h)
    h_min  = ROOT.TH1D(h)
    h_max  = ROOT.TH1D(h)

    h_mean.SetName(f"{h.GetName()}_{bandname}_mean")
    h_min.SetName(f"{h.GetName()}_{bandname}_min")
    h_max.SetName(f"{h.GetName()}_{bandname}_max")

    h_mean.Reset()
    h_min.Reset()
    h_max.Reset()

    h_mean.SetDirectory(0)
    h_min.SetDirectory(0)
    h_max.SetDirectory(0)

    for b in range(0, h.GetNbinsX() + 2):
        vals = [band.GetHist(i).GetBinContent(b) for i in range(n)]
        if not vals:
            continue

        h_mean.SetBinContent(b, sum(vals) / len(vals))
        h_min.SetBinContent(b, min(vals))
        h_max.SetBinContent(b, max(vals))

    return h_mean, h_min, h_max

def compute_asymmetry_summary(h, bandname):
    band = h.GetVertErrorBand(bandname)
    n = band.GetNHists()

    rows = []

    max_abs_mean_shift = 0.0
    max_abs_mid_shift = 0.0
    max_abs_pull = 0.0
    max_envelope_width = 0.0
    n_same_side = 0
    n_nonzero = 0

    for b in range(1, h.GetNbinsX() + 1):
        cv = h.GetBinContent(b)
        bcv = band.GetBinContent(b)
        vals = [band.GetHist(i).GetBinContent(b) for i in range(n)]

        if cv == 0 and all(v == 0 for v in vals):
            continue

        n_nonzero += 1

        ratios = [(v / cv if cv else 0.0) for v in vals]
        mean_ratio = sum(ratios) / len(ratios) if ratios else 0.0
        min_ratio = min(ratios) if ratios else 0.0
        max_ratio = max(ratios) if ratios else 0.0

        mean_shift = mean_ratio - 1.0
        envelope_width = max_ratio - min_ratio

        max_abs_mean_shift = max(max_abs_mean_shift, abs(mean_shift))
        max_envelope_width = max(max_envelope_width, envelope_width)
        max_abs_pull = max(max_abs_pull, max(abs(min_ratio - 1.0), abs(max_ratio - 1.0)))

        same_side = False
        mid_shift = None
        halfspread = None
        offset_over_halfspread = None

        if n == 2:
            u0 = ratios[0]
            u1 = ratios[1]
            midpoint = 0.5 * (u0 + u1)
            halfspread = 0.5 * abs(u0 - u1)
            mid_shift = midpoint - 1.0
            same_side = ((u0 - 1.0) * (u1 - 1.0)) > 0

            if same_side:
                n_same_side += 1

            max_abs_mid_shift = max(max_abs_mid_shift, abs(mid_shift))
            offset_over_halfspread = abs(mid_shift) / halfspread if halfspread else float("inf")

        rows.append({
            "bin": b,
            "xlo": h.GetXaxis().GetBinLowEdge(b),
            "xhi": h.GetXaxis().GetBinUpEdge(b),
            "cv": cv,
            "bandcv_ratio": bcv / cv if cv else 0.0,
            "mean_ratio": mean_ratio,
            "min_ratio": min_ratio,
            "max_ratio": max_ratio,
            "mean_shift": mean_shift,
            "envelope_width": envelope_width,
            "mid_shift": mid_shift,
            "halfspread": halfspread,
            "offset_over_halfspread": offset_over_halfspread,
            "same_side": same_side,
        })

    same_side_frac = n_same_side / float(n_nonzero) if n_nonzero else 0.0

    return {
        "band": bandname,
        "n_universes": n,
        "n_nonzero": n_nonzero,
        "n_same_side": n_same_side,
        "same_side_frac": same_side_frac,
        "max_abs_mean_shift": max_abs_mean_shift,
        "max_abs_mid_shift": max_abs_mid_shift,
        "max_abs_pull": max_abs_pull,
        "max_envelope_width": max_envelope_width,
        "rows": rows,
    }

def print_asymmetry_summary(summary):
    band = summary["band"]
    n = summary["n_universes"]

    print("\n" + "=" * 110)
    print(f"Asymmetry summary: {band}  nUniverses={n}")
    print("=" * 110)

    if n == 2:
        print(
            "nBins={}  sameSide={}/{} ({:.3f})  "
            "max|mid/CV - 1|={:.5g}  max|mean/CV - 1|={:.5g}  "
            "max pull from CV={:.5g}  max envelope width={:.5g}".format(
                summary["n_nonzero"],
                summary["n_same_side"],
                summary["n_nonzero"],
                summary["same_side_frac"],
                summary["max_abs_mid_shift"],
                summary["max_abs_mean_shift"],
                summary["max_abs_pull"],
                summary["max_envelope_width"],
            )
        )
    else:
        print(
            "nBins={}  max|mean/CV - 1|={:.5g}  "
            "max pull from CV={:.5g}  max envelope width={:.5g}".format(
                summary["n_nonzero"],
                summary["max_abs_mean_shift"],
                summary["max_abs_pull"],
                summary["max_envelope_width"],
            )
        )

    print(
        "{:>4s} {:>15s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s} {:>10s}".format(
            "bin", "x-range", "bandCV/CV", "mean/CV", "min/CV", "max/CV", "midShift", "sameSide"
        )
    )

    for r in summary["rows"]:
        mid = r["mid_shift"]
        mid_str = "{:.5g}".format(mid) if mid is not None else "NA"
        print(
            "{:4d} [{:5.2f},{:5.2f}] {:12.5g} {:12.5g} {:12.5g} {:12.5g} {:12s} {:>10s}".format(
                r["bin"],
                r["xlo"],
                r["xhi"],
                r["bandcv_ratio"],
                r["mean_ratio"],
                r["min_ratio"],
                r["max_ratio"],
                mid_str,
                str(r["same_side"]),
            )
        )

def make_ratio_hist(num, den, name):
    """
    Build a plain TH1D ratio histogram num/den bin-by-bin.
    This avoids PyROOT overload issues between MnvH1D and TH1D.
    """
    r = ROOT.TH1D(num)
    r.SetDirectory(0)
    r.SetName(name)
    r.SetStats(0)
    r.Reset()

    for b in range(0, num.GetNbinsX() + 2):
        n = num.GetBinContent(b)
        d = den.GetBinContent(b)

        if d != 0:
            r.SetBinContent(b, n / d)
            # no need to propagate errors for this diagnostic plot
            r.SetBinError(b, 0.0)
        else:
            r.SetBinContent(b, 0.0)
            r.SetBinError(b, 0.0)

    return r

def plot_band(h, bandname, outdir):
    if not h.HasVertErrorBand(bandname):
        print(f"[missing] {h.GetName()} has no band {bandname}")
        return None

    band = h.GetVertErrorBand(bandname)
    n = band.GetNHists()

    os.makedirs(outdir, exist_ok=True)

    h_cv = clone_th1(h, f"{h.GetName()}_{bandname}_CV")
    h_bandcv = clone_th1(band, f"{h.GetName()}_{bandname}_bandCV")
    h_mean, h_min, h_max = make_universe_stats(h, bandname)

    h_cv.SetLineColor(ROOT.kBlack)
    h_cv.SetLineWidth(4)

    h_bandcv.SetLineColor(ROOT.kBlue)
    h_bandcv.SetLineWidth(2)
    h_bandcv.SetLineStyle(2)

    h_mean.SetLineColor(ROOT.kRed + 1)
    h_mean.SetLineWidth(3)

    h_min.SetLineColor(ROOT.kGreen + 2)
    h_min.SetLineWidth(2)
    h_min.SetLineStyle(3)

    h_max.SetLineColor(ROOT.kGreen + 2)
    h_max.SetLineWidth(2)
    h_max.SetLineStyle(3)

    univs = []
    ymax = max(h_cv.GetMaximum(), h_bandcv.GetMaximum(), h_mean.GetMaximum(), h_max.GetMaximum())

    for i in range(n):
        hu = clone_th1(band.GetHist(i), f"{h.GetName()}_{bandname}_u{i}")
        hu.SetLineColor(ROOT.kGray + 1)
        hu.SetLineWidth(1)
        univs.append(hu)
        ymax = max(ymax, hu.GetMaximum())

    c = ROOT.TCanvas(f"c_{safe_name(h.GetName())}_{safe_name(bandname)}", "", 1200, 1000)
    c.Divide(1, 2)

    # Top pad
    top = c.cd(1)
    top.SetPad(0.0, 0.30, 1.0, 1.0)
    top.SetBottomMargin(0.02)

    h_cv.SetMaximum(1.25 * ymax if ymax > 0 else 1.0)
    h_cv.GetYaxis().SetTitle("Events")
    h_cv.GetXaxis().SetLabelSize(0)
    h_cv.SetTitle(f"{h.GetName()} : {bandname} ({n} universes)")
    h_cv.Draw("HIST")

    for hu in univs:
        hu.Draw("HIST SAME")

    h_min.Draw("HIST SAME")
    h_max.Draw("HIST SAME")
    h_mean.Draw("HIST SAME")
    h_bandcv.Draw("HIST SAME")
    h_cv.Draw("HIST SAME")

    leg = ROOT.TLegend(0.55, 0.62, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(h_cv, "CV", "l")
    leg.AddEntry(h_bandcv, "Band CV", "l")
    leg.AddEntry(h_mean, "Universe mean", "l")
    leg.AddEntry(h_min, "Universe min/max", "l")
    if univs:
        leg.AddEntry(univs[0], "Universes", "l")
    leg.Draw()

    # Bottom pad
    bot = c.cd(2)
    bot.SetPad(0.0, 0.0, 1.0, 0.30)
    bot.SetTopMargin(0.03)
    bot.SetBottomMargin(0.30)

    ratios = []
    for i, hu in enumerate(univs):
        rr = make_ratio_hist(hu, h_cv, f"{hu.GetName()}_ratio")
        rr.SetLineColor(ROOT.kGray + 1)
        rr.SetLineWidth(1)
        ratios.append(rr)

    r_mean = make_ratio_hist(h_mean, h_cv, f"{h_mean.GetName()}_ratio")
    r_min = make_ratio_hist(h_min, h_cv, f"{h_min.GetName()}_ratio")
    r_max = make_ratio_hist(h_max, h_cv, f"{h_max.GetName()}_ratio")
    r_bandcv = make_ratio_hist(h_bandcv, h_cv, f"{h_bandcv.GetName()}_ratio")

    r_mean.SetLineColor(ROOT.kRed + 1)
    r_mean.SetLineWidth(3)
    r_min.SetLineColor(ROOT.kGreen + 2)
    r_min.SetLineStyle(3)
    r_min.SetLineWidth(2)
    r_max.SetLineColor(ROOT.kGreen + 2)
    r_max.SetLineStyle(3)
    r_max.SetLineWidth(2)
    r_bandcv.SetLineColor(ROOT.kBlue)
    r_bandcv.SetLineStyle(2)
    r_bandcv.SetLineWidth(2)

    # Auto y-range from ratios
    vals = []
    for rr in ratios + [r_mean, r_min, r_max]:
        for b in range(1, rr.GetNbinsX() + 1):
            v = rr.GetBinContent(b)
            if math.isfinite(v) and v != 0:
                vals.append(v)

    ymin = min(vals) if vals else 0.8
    ymaxr = max(vals) if vals else 1.2
    pad = 0.15 * max(abs(ymaxr - 1.0), abs(1.0 - ymin), 0.05)
    ymin = min(ymin - pad, 0.95)
    ymaxr = max(ymaxr + pad, 1.05)

    frame = h_cv.Clone(f"{h_cv.GetName()}_ratio_frame")
    frame.Reset()
    frame.SetDirectory(0)
    frame.SetStats(0)
    frame.GetYaxis().SetTitle("Universe / CV")
    frame.GetXaxis().SetTitle(h_cv.GetXaxis().GetTitle())
    frame.GetYaxis().SetTitleSize(0.09)
    frame.GetYaxis().SetLabelSize(0.08)
    frame.GetXaxis().SetTitleSize(0.10)
    frame.GetXaxis().SetLabelSize(0.09)
    frame.GetYaxis().SetNdivisions(505)
    frame.SetMinimum(ymin)
    frame.SetMaximum(ymaxr)
    frame.Draw("AXIS")

    for rr in ratios:
        rr.Draw("HIST SAME")

    r_min.Draw("HIST SAME")
    r_max.Draw("HIST SAME")
    r_mean.Draw("HIST SAME")
    r_bandcv.Draw("HIST SAME")

    line = ROOT.TLine(h_cv.GetXaxis().GetXmin(), 1.0, h_cv.GetXaxis().GetXmax(), 1.0)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw("SAME")

    c.cd()
    c.Update()

    # Keep references alive
    c._objs = [h_cv, h_bandcv, h_mean, h_min, h_max, leg, frame, line, r_mean, r_min, r_max, r_bandcv] + univs + ratios

    outname = os.path.join(outdir, f"{safe_name(h.GetName())}_{safe_name(bandname)}.png")
    c.SaveAs(outname)
    print("Saved", outname)

    return compute_asymmetry_summary(h, bandname)

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_FILE
    hname = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_HIST
    outdir = sys.argv[3] if len(sys.argv) > 3 else DEFAULT_OUTDIR

    # Optional:
    #   --all-bands
    #   --bands LowQ2Pi,fsi_weight,MK_model
    all_bands = False
    selected_bands = None

    args = sys.argv[4:]
    i = 0
    while i < len(args):
        if args[i] == "--all-bands":
            all_bands = True
            i += 1
        elif args[i] == "--bands":
            selected_bands = [x.strip() for x in args[i + 1].split(",") if x.strip()]
            i += 2
        else:
            i += 1

    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open {path}")

    h = f.Get(hname)
    if not h:
        raise RuntimeError(f"Could not find histogram {hname}")

    print("Opened:", path)
    print("Histogram:", hname)
    print("Integral:", h.Integral())

    if all_bands:
        bands = get_band_names(h)
    elif selected_bands is not None:
        bands = selected_bands
    else:
        bands = DEFAULT_BANDS

    print("Bands to plot:")
    for b in bands:
        print(" ", b)

    summaries = []

    for bandname in bands:
        if not h.HasVertErrorBand(bandname):
            print("[missing]", bandname)
            continue

        summary = plot_band(h, bandname, outdir)
        if summary:
            summaries.append(summary)
            print_asymmetry_summary(summary)

    print("\n" + "#" * 120)
    print("GLOBAL ASYMMETRY RANKING")
    print("#" * 120)

    summaries.sort(key=lambda s: s["max_abs_pull"], reverse=True)

    print(
        "{:30s} {:>5s} {:>12s} {:>14s} {:>14s} {:>14s} {:>12s}".format(
            "band", "nU", "sameSide", "max|mean-1|", "max|mid-1|", "max pull", "max width"
        )
    )

    for s in summaries:
        print(
            "{:30s} {:5d} {:>12s} {:14.5g} {:14.5g} {:14.5g} {:14.5g}".format(
                s["band"],
                s["n_universes"],
                f"{s['n_same_side']}/{s['n_nonzero']}" if s["n_universes"] == 2 else "NA",
                s["max_abs_mean_shift"],
                s["max_abs_mid_shift"] if s["n_universes"] == 2 else 0.0,
                s["max_abs_pull"],
                s["max_envelope_width"],
            )
        )

    f.Close()

if __name__ == "__main__":
    main()
    # f = ROOT.TFile.Open("/exp/minerva/data/users/qvuong/nu_e/kin_dist_mcle1_p6_test_MAD_0.root")
    # h = f.Get("EN4")
    # b = h.GetVertErrorBand("electron_scale")

    # print("CV integral:", h.Integral(0, h.GetNbinsX()+1))
    # print("band CV integral:", b.Integral(0, h.GetNbinsX()+1))

    # for i in range(b.GetNHists()):
    #     hu = b.GetHist(i)
    #     print(i, hu.GetName(), hu.Integral(0, hu.GetNbinsX()+1))
    #     print(
    #         "univ", i,
    #         "integral:", hu.Integral(0, hu.GetNbinsX()+1),
    #         "ratio:", hu.Integral(0, hu.GetNbinsX()+1) / h.Integral(0, h.GetNbinsX()+1)
    #     )