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
        "den_data": BASE_DIR + "/nu_e/kin_dist_dataleFHC_p8Tuples_CCnue_allSystematics_MAD.root",
        "den_mc":   BASE_DIR + "/nu_e/kin_dist_mcleFHC_p8Tuples_CCnue_allSystematics_MAD.root",
        "title": "FHC CC#nu_{#mu}/CC#nu_{e}",
        "data_label": "Data CC#nu_{#mu}/CC#nu_{e}",
        "mc_label": "MC CC#nu_{#mu}/CC#nu_{e}",
    },
    "rhc": {
        "num_data": BASE_DIR + "/antinu_mu/kin_dist_datale5_p8Tuples_CCnumubar_allSystematics_MAD.root",
        "num_mc":   BASE_DIR + "/antinu_mu/kin_dist_mcle5_p8Tuples_CCnumubar_allSystematics_MAD.root",
        "den_data": BASE_DIR + "/antinu_e/kin_dist_datale5_p8Tuples_CCnuebar_allSystematics_MAD.root",
        "den_mc":   BASE_DIR + "/antinu_e/kin_dist_mcle5_p8Tuples_CCnuebar_allSystematics_MAD.root",
        "title": "RHC CC#bar{#nu}_{#mu}/CC#bar{#nu}_{e}",
        "data_label": "Data CC#bar{#nu}_{#mu}/CC#bar{#nu}_{e}",
        "mc_label": "MC CC#bar{#nu}_{#mu}/CC#bar{#nu}_{e}",
    },
}


def get_hist(path, hname):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError("Could not open " + path)

    h = f.Get(hname)
    if not h:
        print("\nCould not find exact hist:", hname)
        print("In file:", path)
        print("Possible matches:")
        for key in f.GetListOfKeys():
            name = key.GetName()
            if hname in name:
                print("  ", name)
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


def draw_pair(h_data, h_mc, title, xtitle, data_label, mc_label, outname, ymin=None, ymax=None):
    c = ROOT.TCanvas("c_" + outname.replace(".png", ""), "c", 900, 700)

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

    h_data.SetTitle(title + ";" + xtitle + ";CC#nu_{#mu}/CC#nu_{e} ratio")

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

    leg = ROOT.TLegend(0.50, 0.72, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(h_data, data_label, "lep")
    leg.AddEntry(h_mc, mc_label, "lep")
    leg.Draw()

    c.Print(outname)
    c.Print(outname.replace(".png", ".pdf"))

    print("Wrote", outname)
    print("Wrote", outname.replace(".png", ".pdf"))


def process_mode(mode, cfg, hname, outdir, ymin, ymax):
    print("\n===== {} =====".format(mode))
    for k in ["num_data", "den_data", "num_mc", "den_mc"]:
        print(k + ":", cfg[k])

    h_num_data = get_hist(cfg["num_data"], hname)
    h_den_data = get_hist(cfg["den_data"], hname)
    h_num_mc = get_hist(cfg["num_mc"], hname)
    h_den_mc = get_hist(cfg["den_mc"], hname)

    r_data_energy = divide(h_num_data, h_den_data, mode + "_data_ccnumu_over_ccnue_energy")
    r_mc_energy = divide(h_num_mc, h_den_mc, mode + "_mc_ccnumu_over_ccnue_energy")

    r_data_bin = make_bin_number_hist(r_data_energy, mode + "_data_ccnumu_over_ccnue_bin")
    r_mc_bin = make_bin_number_hist(r_mc_energy, mode + "_mc_ccnumu_over_ccnue_bin")

    # Titles avoid .format() to keep ROOT #nu_{...} braces safe.
    draw_pair(
        r_data_energy,
        r_mc_energy,
        cfg["title"],
        "E_{#nu}^{reco} [GeV]",
        cfg["data_label"],
        cfg["mc_label"],
        os.path.join(outdir, mode + "_ccnumu_over_ccnue_energy.png"),
        ymin,
        ymax,
    )

    draw_pair(
        r_data_bin,
        r_mc_bin,
        cfg["title"],
        "Analysis bin number",
        cfg["data_label"],
        cfg["mc_label"],
        os.path.join(outdir, mode + "_ccnumu_over_ccnue_bin_number.png"),
        ymin,
        ymax,
    )

    return {
        "num_data": h_num_data,
        "den_data": h_den_data,
        "num_mc": h_num_mc,
        "den_mc": h_den_mc,
        "r_data_energy": r_data_energy,
        "r_mc_energy": r_mc_energy,
        "r_data_bin": r_data_bin,
        "r_mc_bin": r_mc_bin,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hist", default="EN4")
    parser.add_argument("--outdir", default="ccnumu_ccnue_ratio_plots")
    parser.add_argument("--modes", default="fhc,rhc", help="fhc,rhc or one of them")
    parser.add_argument("--ymin", type=float, default=None)
    parser.add_argument("--ymax", type=float, default=None)
    args = parser.parse_args()

    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    modes = [m.strip().lower() for m in args.modes.split(",") if m.strip()]

    fout = ROOT.TFile(os.path.join(args.outdir, "ccnumu_ccnue_selection_ratios.root"), "RECREATE")

    for mode in modes:
        if mode not in CONFIG:
            raise RuntimeError("Unknown mode: " + mode)

        out = process_mode(mode, CONFIG[mode], args.hist, args.outdir, args.ymin, args.ymax)

        fout.cd()
        for name, hist in out.items():
            hist.Write(mode + "_" + name)

    fout.Close()
    print("\nWrote", os.path.join(args.outdir, "ccnumu_ccnue_selection_ratios.root"))


if __name__ == "__main__":
    main()