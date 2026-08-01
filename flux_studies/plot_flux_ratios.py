import ROOT
import argparse
import os

ROOT.TH1.AddDirectory(False)
ROOT.gStyle.SetOptStat(0)

def get_hist(path, hname="flux_E_cvweighted"):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError("Could not open {}".format(path))

    h = f.Get(hname)
    if not h:
        raise RuntimeError("Could not find {} in {}".format(hname, path))

    h = h.Clone(os.path.basename(path).replace(".root", "_" + hname))
    h.SetDirectory(0)
    f.Close()
    return h

def make_ratio(num, den, name):
    r = num.Clone(name)
    r.SetDirectory(0)
    r.Divide(num, den, 1.0, 1.0, "B")
    return r

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--flux-dir",
        default="/exp/minerva/app/users/qvuong/MAT_AL9/CC-NuE-XSec/custom_plotutils/data/flux",
    )
    parser.add_argument(
        "--tag",
        default="minerva1",
        help="Example: minerva1, minerva5, minervame1A, etc.",
    )
    parser.add_argument(
        "--hname",
        default="flux_E_cvweighted",
    )
    parser.add_argument(
        "--out",
        default=None,
    )
    args = parser.parse_args()

    if args.out is None:
        args.out = "flux_ratios_{}.png".format(args.tag)

    f_nue     = "{}/flux-gen2thin-pdg12-{}.root".format(args.flux_dir, args.tag)
    f_numu    = "{}/flux-gen2thin-pdg14-{}.root".format(args.flux_dir, args.tag)
    f_nuebar  = "{}/flux-gen2thin-pdg-12-{}.root".format(args.flux_dir, args.tag)
    f_numubar = "{}/flux-gen2thin-pdg-14-{}.root".format(args.flux_dir, args.tag)

    h_nue     = get_hist(f_nue, args.hname)
    h_numu    = get_hist(f_numu, args.hname)
    h_nuebar  = get_hist(f_nuebar, args.hname)
    h_numubar = get_hist(f_numubar, args.hname)

    r_numu_nue = make_ratio(h_numu, h_nue, "ratio_numu_over_nue")
    r_numubar_nuebar = make_ratio(h_numubar, h_nuebar, "ratio_numubar_over_nuebar")

    r_numu_nue.SetLineColor(ROOT.kBlue + 1)
    r_numu_nue.SetLineWidth(3)
    r_numubar_nuebar.SetLineColor(ROOT.kRed + 1)
    r_numubar_nuebar.SetLineWidth(3)

    r_numu_nue.SetTitle("Flux ratios: {};E_{{#nu}} [GeV];Flux ratio".format(args.tag))
    r_numu_nue.GetXaxis().SetRangeUser(0, 20)

    ymax = max(r_numu_nue.GetMaximum(), r_numubar_nuebar.GetMaximum())
    r_numu_nue.GetYaxis().SetRangeUser(0, 1.2 * ymax)

    c = ROOT.TCanvas("c", "c", 900, 700)
    r_numu_nue.Draw("HIST")
    r_numubar_nuebar.Draw("HIST SAME")

    leg = ROOT.TLegend(0.55, 0.72, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(r_numu_nue, "#nu_{#mu} / #nu_{e}", "l")
    leg.AddEntry(r_numubar_nuebar, "#bar{#nu}_{#mu} / #bar{#nu}_{e}", "l")
    leg.Draw()

    c.Print(args.out)

    fout = ROOT.TFile("flux_ratios_{}.root".format(args.tag), "RECREATE")
    r_numu_nue.Write()
    r_numubar_nuebar.Write()
    fout.Close()

    print("Wrote", args.out)
    print("Wrote flux_ratios_{}.root".format(args.tag))

if __name__ == "__main__":
    main()