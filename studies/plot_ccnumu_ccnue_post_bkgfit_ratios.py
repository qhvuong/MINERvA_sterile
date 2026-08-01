import ROOT
import os
import argparse

ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetOptStat(0)


BASE_DIR = "/exp/minerva/data/users/qvuong"


CONFIG = {
    "fhc": {
        "num_data": BASE_DIR + "/nu_mu/kin_dist_dataleFHC_p8Tuples_CCnumu_allSystematics_MAD.root",
        "num_mc":   BASE_DIR + "/nu_mu/kin_dist_mcleFHC_p8Tuples_CCnumu_allSystematics_MAD.root",

        "den_tuned": BASE_DIR + "/nu_e/bkgfit_leFHC_N4_tune_p8Tuples_CCnue_allSystematics_MAD.root",

        "title": "FHC CC#nu_{#mu}/CC#nu_{e} post-bkg-tune",
        "data_label": "Data CC#nu_{#mu}/CC#nu_{e}",
        "mc_label": "MC CC#nu_{#mu}/post-tune CC#nu_{e}",
        "ytitle": "CC#nu_{#mu}/CC#nu_{e} ratio",
    },

    "rhc": {
        "num_data": BASE_DIR + "/antinu_mu/kin_dist_datale5_p8Tuples_CCnumubar_allSystematics_MAD.root",
        "num_mc":   BASE_DIR + "/antinu_mu/kin_dist_mcle5_p8Tuples_CCnumubar_allSystematics_MAD.root",

        "den_tuned": BASE_DIR + "/antinu_e/bkgfit_le5_N4_tune_p8Tuples_CCnuebar_allSystematics_MAD.root",

        "title": "RHC CC#bar{#nu}_{#mu}/CC#bar{#nu}_{e} post-bkg-tune",
        "data_label": "Data CC#bar{#nu}_{#mu}/CC#bar{#nu}_{e}",
        "mc_label": "MC CC#bar{#nu}_{#mu}/post-tune CC#bar{#nu}_{e}",
        "ytitle": "CC#bar{#nu}_{#mu}/CC#bar{#nu}_{e} ratio",
    },
}


def get_pot_from_meta(path):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError("Could not open " + path)

    meta = f.Get("Meta")
    if not meta:
        f.Close()
        raise RuntimeError("Could not find Meta tree in " + path)

    total_pot = 0.0
    for i in range(meta.GetEntries()):
        meta.GetEntry(i)
        total_pot += float(meta.POT_Used)

    f.Close()
    return total_pot


def print_matching_keys(f, hname):
    print("Possible keys containing '{}':".format(hname))
    found = False
    for key in f.GetListOfKeys():
        name = key.GetName()
        if hname in name:
            print("  " + name)
            found = True
    if not found:
        print("  none")


def get_hist(path, hname):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError("Could not open " + path)

    h = f.Get(hname)

    if not h:
        print("\nCould not find exact hist:", hname)
        print("In file:", path)
        print_matching_keys(f, hname)
        raise RuntimeError("Missing histogram " + hname)

    h = h.Clone(os.path.basename(path).replace(".root", "_" + hname))
    h.SetDirectory(0)
    f.Close()
    return h


def divide(num, den, name):
    r = num.Clone(name)
    r.SetDirectory(0)
    r.Divide(num, den)
    return r


def make_bin_number_hist(h_energy, name):
    nbins = h_energy.GetNbinsX()
    h_bin = ROOT.TH1D(name, name, nbins, 0.5, nbins + 0.5)
    h_bin.SetDirectory(0)

    for i in range(1, nbins + 1):
        h_bin.SetBinContent(i, h_energy.GetBinContent(i))
        h_bin.SetBinError(i, h_energy.GetBinError(i))

    return h_bin


def draw_pair(h_data, h_mc, title, xtitle, ytitle, data_label, mc_label, outname, ymin=None, ymax=None):
    c = ROOT.TCanvas("c_" + os.path.basename(outname).replace(".png", ""), "c", 900, 700)

    h_data.SetLineColor(ROOT.kBlack)
    h_data.SetMarkerColor(ROOT.kBlack)
    h_data.SetMarkerStyle(20)
    h_data.SetMarkerSize(1.1)
    h_data.SetLineWidth(2)

    h_mc.SetLineColor(ROOT.kRed + 1)
    h_mc.SetMarkerColor(ROOT.kRed + 1)
    h_mc.SetMarkerStyle(24)
    h_mc.SetMarkerSize(1.0)
    h_mc.SetLineWidth(2)

    h_data.SetTitle(title + ";" + xtitle + ";" + ytitle)

    if ymin is None or ymax is None:
        maxval = max(h_data.GetMaximum(), h_mc.GetMaximum())
        if maxval <= 0:
            maxval = 1.0
        ymin_plot = 0.0
        ymax_plot = 1.25 * maxval
    else:
        ymin_plot = ymin
        ymax_plot = ymax

    h_data.GetYaxis().SetRangeUser(ymin_plot, ymax_plot)

    h_data.Draw("E1")
    h_mc.Draw("E1 SAME")

    leg = ROOT.TLegend(0.45, 0.72, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(h_data, data_label, "lep")
    leg.AddEntry(h_mc, mc_label, "lep")
    leg.Draw()

    c.Print(outname)
    c.Print(outname.replace(".png", ".pdf"))

    print("Wrote", outname)
    print("Wrote", outname.replace(".png", ".pdf"))


def check_binning(h1, h2, name1, name2):
    print("\nChecking binning:")
    print("  {} nbins = {}".format(name1, h1.GetNbinsX()))
    print("  {} nbins = {}".format(name2, h2.GetNbinsX()))

    if h1.GetNbinsX() != h2.GetNbinsX():
        print("  ERROR: different number of bins")
        return False

    ok = True
    for i in range(1, h1.GetNbinsX() + 1):
        lo1 = h1.GetBinLowEdge(i)
        hi1 = lo1 + h1.GetBinWidth(i)
        lo2 = h2.GetBinLowEdge(i)
        hi2 = lo2 + h2.GetBinWidth(i)

        if abs(lo1 - lo2) > 1e-6 or abs(hi1 - hi2) > 1e-6:
            print(
                "  ERROR bin {} mismatch: {} [{:.6g}, {:.6g}] vs {} [{:.6g}, {:.6g}]".format(
                    i, name1, lo1, hi1, name2, lo2, hi2
                )
            )
            ok = False

    if ok:
        print("  binning OK")

    return ok


def print_bin_table(label, h):
    print("\n--- {} ---".format(label))
    print("Integral =", h.Integral())
    for i in range(1, h.GetNbinsX() + 1):
        lo = h.GetBinLowEdge(i)
        hi = lo + h.GetBinWidth(i)
        print(
            "bin {:2d}: [{:.3f}, {:.3f}]  {:.6g}".format(
                i, lo, hi, h.GetBinContent(i)
            )
        )


def process_mode(mode, cfg, data_hname, mc_hname, outdir, ymin, ymax):
    print("\n===== {} =====".format(mode))
    for k in ["num_data", "num_mc", "den_tuned"]:
        print(k + ":", cfg[k])

    # Numerator: CCnumu/CCnumubar has no background tuning
    # Use normal selected data and MC histograms.
    h_num_data = get_hist(cfg["num_data"], "EN4")
    h_num_mc = get_hist(cfg["num_mc"], "EN4")

    # Scale CCnumu/CCnumubar MC to data POT using Meta tree
    data_pot = get_pot_from_meta(cfg["num_data"])
    mc_pot = get_pot_from_meta(cfg["num_mc"])

    if mc_pot <= 0:
        raise RuntimeError("MC POT is zero or negative for " + cfg["num_mc"])

    num_mc_pot_scale = data_pot / mc_pot

    print("num data POT =", data_pot)
    print("num MC POT   =", mc_pot)
    print("num MC POT scale =", num_mc_pot_scale)

    h_num_mc.Scale(num_mc_pot_scale)


    # Denominator: CCnue/CCnuebar has background tuning
    # Use the same post-tune objects used by the oscillation fit.
    h_den_data = get_hist(cfg["den_tuned"], data_hname)
    h_den_mc = get_hist(cfg["den_tuned"], mc_hname)

    check_binning(h_num_data, h_den_data, "num data", "den data")
    check_binning(h_num_mc, h_den_mc, "num MC", "den MC")

    print_bin_table("CCnumu data EN4", h_num_data)
    print_bin_table("CCnumu MC EN4", h_num_mc)
    print_bin_table("CCnue data_bkgSubbed", h_den_data)
    print_bin_table("CCnue predicted_Signal", h_den_mc)

    r_data_energy = divide(
        h_num_data,
        h_den_data,
        mode + "_data_numu_raw_over_nue_bkgSubbed_ratio_energy",
    )

    r_mc_energy = divide(
        h_num_mc,
        h_den_mc,
        mode + "_mc_numu_raw_over_nue_predicted_signal_ratio_energy",
    )

    r_data_bin = make_bin_number_hist(
        r_data_energy,
        mode + "_data_numu_raw_over_nue_bkgSubbed_ratio_bin",
    )

    r_mc_bin = make_bin_number_hist(
        r_mc_energy,
        mode + "_mc_numu_raw_over_nue_predicted_signal_ratio_bin",
    )

    draw_pair(
        r_data_energy,
        r_mc_energy,
        cfg["title"],
        "E_{#nu}^{reco} [GeV]",
        cfg["ytitle"],
        cfg["data_label"],
        cfg["mc_label"],
        os.path.join(outdir, mode + "_ccnumu_over_ccnue_mixed_postbkgfit_energy.png"),
        ymin,
        ymax,
    )

    draw_pair(
        r_data_bin,
        r_mc_bin,
        cfg["title"],
        "Analysis bin number",
        cfg["ytitle"],
        cfg["data_label"],
        cfg["mc_label"],
        os.path.join(outdir, mode + "_ccnumu_over_ccnue_mixed_postbkgfit_bin_number.png"),
        ymin,
        ymax,
    )

    return {
        "num_data_raw": h_num_data,
        "num_mc_raw": h_num_mc,
        "den_data_bkgSubbed": h_den_data,
        "den_mc_predicted_signal": h_den_mc,
        "r_data_energy": r_data_energy,
        "r_mc_energy": r_mc_energy,
        "r_data_bin": r_data_bin,
        "r_mc_bin": r_mc_bin,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-hist", default="EN4_data_bkgSubbed")
    parser.add_argument("--mc-hist", default="EN4_predicted_Signal")
    parser.add_argument("--outdir", default="ccnumu_ccnue_postbkgfit_ratio_plots")
    parser.add_argument("--modes", default="fhc,rhc", help="fhc,rhc or one of them")
    parser.add_argument("--ymin", type=float, default=None)
    parser.add_argument("--ymax", type=float, default=None)
    args = parser.parse_args()

    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    modes = [m.strip().lower() for m in args.modes.split(",") if m.strip()]

    fout = ROOT.TFile(os.path.join(args.outdir, "ccnumu_ccnue_postbkgfit_selection_ratios.root"), "RECREATE")

    for mode in modes:
        if mode not in CONFIG:
            raise RuntimeError("Unknown mode: " + mode)

        out = process_mode(
            mode,
            CONFIG[mode],
            args.data_hist,
            args.mc_hist,
            args.outdir,
            args.ymin,
            args.ymax,
        )

        fout.cd()
        for name, hist in out.items():
            hist.Write(mode + "_" + name)

    fout.Close()
    print("\nWrote", os.path.join(args.outdir, "ccnumu_ccnue_postbkgfit_selection_ratios.root"))


if __name__ == "__main__":
    main()