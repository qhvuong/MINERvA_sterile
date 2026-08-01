import ROOT
import os
import argparse

ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetOptStat(0)


BASE_DIR = "/exp/minerva/data/users/qvuong"


SAMPLES = {
    "fhc_ccnue": {
        "data": BASE_DIR + "/nu_e/kin_dist_dataleFHC_p8Tuples_CCnue_allSystematics_MAD.root",
        "mc":   BASE_DIR + "/nu_e/kin_dist_mcleFHC_p8Tuples_CCnue_allSystematics_MAD.root",
        "label": "FHC CC#nu_{e}",
    },
    "fhc_ccnumu": {
        "data": BASE_DIR + "/nu_mu/kin_dist_dataleFHC_p8Tuples_CCnumu_allSystematics_MAD.root",
        "mc":   BASE_DIR + "/nu_mu/kin_dist_mcleFHC_p8Tuples_CCnumu_allSystematics_MAD.root",
        "label": "FHC CC#nu_{#mu}",
    },
    "rhc_ccnuebar": {
        "data": BASE_DIR + "/antinu_e/kin_dist_datale5_p8Tuples_CCnuebar_allSystematics_MAD.root",
        "mc":   BASE_DIR + "/antinu_e/kin_dist_mcle5_p8Tuples_CCnuebar_allSystematics_MAD.root",
        "label": "RHC CC#bar{#nu}_{e}",
    },
    "rhc_ccnumubar": {
        "data": BASE_DIR + "/antinu_mu/kin_dist_datale5_p8Tuples_CCnumubar_allSystematics_MAD.root",
        "mc":   BASE_DIR + "/antinu_mu/kin_dist_mcle5_p8Tuples_CCnumubar_allSystematics_MAD.root",
        "label": "RHC CC#bar{#nu}_{#mu}",
    },
}


def get_hist(path, hname):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError("Could not open " + path)

    h = f.Get(hname)

    if not h:
        print("\nCould not find exact hist:", hname)
        print("Searching keys containing:", hname)
        matches = []
        for key in f.GetListOfKeys():
            name = key.GetName()
            if hname in name:
                matches.append(name)

        if matches:
            print("Possible matches:")
            for m in matches:
                print("  ", m)
        else:
            print("No matching keys found.")
        raise RuntimeError("Missing histogram " + hname + " in " + path)

    h = h.Clone(os.path.basename(path).replace(".root", "_" + hname))
    h.SetDirectory(0)
    f.Close()
    return h


def make_ratio(data, mc, name):
    r = data.Clone(name)
    r.SetDirectory(0)
    r.Divide(data, mc)
    return r


def make_bin_number_ratio(r_energy, name):
    nbins = r_energy.GetNbinsX()
    r_bin = ROOT.TH1D(name, name, nbins, 0.5, nbins + 0.5)
    r_bin.SetDirectory(0)

    for i in range(1, nbins + 1):
        r_bin.SetBinContent(i, r_energy.GetBinContent(i))
        r_bin.SetBinError(i, r_energy.GetBinError(i))

    return r_bin


def draw_ratio(ratio, title, xtitle, outname, ymin=0.0, ymax=2.0):
    c = ROOT.TCanvas("c_" + outname.replace(".png", ""), "c", 900, 700)

    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetMarkerColor(ROOT.kBlack)
    ratio.SetMarkerStyle(20)
    ratio.SetMarkerSize(1.1)
    ratio.SetLineWidth(2)

    ratio.SetTitle(title + ";" + xtitle + ";Data / MC")
    ratio.GetYaxis().SetRangeUser(ymin, ymax)
    ratio.Draw("E1")

    line = ROOT.TLine(
        ratio.GetXaxis().GetXmin(),
        1.0,
        ratio.GetXaxis().GetXmax(),
        1.0,
    )
    line.SetLineColor(ROOT.kRed + 1)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw("SAME")

    c.Print(outname)
    c.Print(outname.replace(".png", ".pdf"))

    print("Wrote", outname)
    print("Wrote", outname.replace(".png", ".pdf"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hist", default="EN4")
    parser.add_argument("--outdir", default="selection_ratio_plots")
    parser.add_argument("--ymin", type=float, default=0.0)
    parser.add_argument("--ymax", type=float, default=2.0)
    parser.add_argument(
        "--samples",
        default="all",
        help="Comma-separated sample keys or all. Options: " + ",".join(SAMPLES.keys()),
    )
    args = parser.parse_args()

    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    if args.samples == "all":
        sample_keys = list(SAMPLES.keys())
    else:
        sample_keys = [s.strip() for s in args.samples.split(",") if s.strip()]

    fout = ROOT.TFile(os.path.join(args.outdir, "selection_datamc_ratios.root"), "RECREATE")

    for key in sample_keys:
        if key not in SAMPLES:
            raise RuntimeError("Unknown sample key: " + key)

        cfg = SAMPLES[key]
        print("\n===== " + key + " =====")
        print("data:", cfg["data"])
        print("mc:  ", cfg["mc"])

        h_data = get_hist(cfg["data"], args.hist)
        h_mc = get_hist(cfg["mc"], args.hist)

        r_energy = make_ratio(h_data, h_mc, key + "_datamc_ratio_energy")
        r_bin = make_bin_number_ratio(r_energy, key + "_datamc_ratio_bin_number")

        title = cfg["label"] + " " + args.hist + " Data/MC"

        draw_ratio(
            r_energy,
            title,
            "E_{#nu}^{reco} [GeV]",
            os.path.join(args.outdir, key + "_datamc_ratio_energy.png"),
            ymin=args.ymin,
            ymax=args.ymax,
        )

        draw_ratio(
            r_bin,
            title,
            "Analysis bin number",
            os.path.join(args.outdir, key + "_datamc_ratio_bin_number.png"),
            ymin=args.ymin,
            ymax=args.ymax,
        )

        fout.cd()
        h_data.Write(key + "_data_" + args.hist)
        h_mc.Write(key + "_mc_" + args.hist)
        r_energy.Write()
        r_bin.Write()

    fout.Close()
    print("\nWrote", os.path.join(args.outdir, "selection_datamc_ratios.root"))


if __name__ == "__main__":
    main()