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


# def make_ratio(num, den, name):
#     r = num.Clone(name)
#     r.SetDirectory(0)
#     r.Divide(num, den, 1.0, 1.0, "B")
#     return r
def make_ratio(num, den, name):
    r = num.Clone(name)
    r.SetDirectory(0)
    r.Divide(num, den, 1.0, 1.0)
    return r

def flux_path(flux_dir, pdg, tag):
    return "{}/flux-gen2thin-pdg{}-{}.root".format(flux_dir, pdg, tag)


def make_right_sign_ratio(flux_dir, tag, mode, hname):
    """
    mode = "fhc": numu / nue
    mode = "rhc": numubar / nuebar
    """
    if mode == "fhc":
        num_pdg = 14
        den_pdg = 12
        ratio_name = "ratio_numu_over_nue_{}".format(tag)
    elif mode == "rhc":
        num_pdg = -14
        den_pdg = -12
        ratio_name = "ratio_numubar_over_nuebar_{}".format(tag)
    else:
        raise RuntimeError("Unknown mode: {}".format(mode))

    h_num = get_hist(flux_path(flux_dir, num_pdg, tag), hname)
    h_den = get_hist(flux_path(flux_dir, den_pdg, tag), hname)

    return make_ratio(h_num, h_den, ratio_name)


def draw_compare(r_le, r_me, title, le_label, me_label, out_png, out_pdf=None, xmax=20.0):
    r_le.SetLineColor(ROOT.kBlue + 1)
    r_le.SetLineWidth(3)

    r_me.SetLineColor(ROOT.kRed + 1)
    r_me.SetLineWidth(3)

    r_le.SetTitle("{};E_{{#nu}} [GeV];Flux ratio".format(title))
    r_le.GetXaxis().SetRangeUser(0, xmax)

    ymax = max(r_le.GetMaximum(), r_me.GetMaximum())
    if ymax <= 0:
        ymax = 1.0
    r_le.GetYaxis().SetRangeUser(0, 1.2 * ymax)

    c = ROOT.TCanvas("c_" + out_png.replace(".png", ""), "c", 900, 700)
    r_le.Draw("HIST")
    r_me.Draw("HIST SAME")

    leg = ROOT.TLegend(0.50, 0.72, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(r_le, le_label, "l")
    leg.AddEntry(r_me, me_label, "l")
    leg.Draw()

    c.Print(out_png)
    if out_pdf:
        c.Print(out_pdf)

    print("Wrote", out_png)
    if out_pdf:
        print("Wrote", out_pdf)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--flux-dir",
        default="/exp/minerva/app/users/qvuong/MAT_AL9/CC-NuE-XSec/custom_plotutils/data/flux",
    )
    parser.add_argument("--hname", default="flux_E_cvweighted")
    parser.add_argument("--xmax", default=20.0, type=float)

    parser.add_argument("--le-fhc-tag", default="minerva1")
    parser.add_argument("--me-fhc-tag", default="minervame1N")

    parser.add_argument("--le-rhc-tag", default="minerva5")
    parser.add_argument("--me-rhc-tag", default="minervame5A")

    parser.add_argument("--out-prefix", default="le_me_flux_ratio_compare")
    args = parser.parse_args()

    # FHC: numu / nue
    r_fhc_le = make_right_sign_ratio(args.flux_dir, args.le_fhc_tag, "fhc", args.hname)
    r_fhc_me = make_right_sign_ratio(args.flux_dir, args.me_fhc_tag, "fhc", args.hname)

    draw_compare(
        r_fhc_le,
        r_fhc_me,
        "FHC right-sign flux ratio: #nu_{{#mu}}/#nu_{{e}}",
        "LE FHC #nu_{{#mu}}/#nu_{{e}} ({})".format(args.le_fhc_tag),
        "ME FHC #nu_{{#mu}}/#nu_{{e}} ({})".format(args.me_fhc_tag),
        "{}_fhc.png".format(args.out_prefix),
        "{}_fhc.pdf".format(args.out_prefix),
        xmax=args.xmax,
    )

    # RHC: numubar / nuebar
    r_rhc_le = make_right_sign_ratio(args.flux_dir, args.le_rhc_tag, "rhc", args.hname)
    r_rhc_me = make_right_sign_ratio(args.flux_dir, args.me_rhc_tag, "rhc", args.hname)

    draw_compare(
        r_rhc_le,
        r_rhc_me,
        "RHC right-sign flux ratio: #bar{{#nu}}_{{#mu}}/#bar{{#nu}}_{{e}}",
        "LE RHC #bar{{#nu}}_{{#mu}}/#bar{{#nu}}_{{e}} ({})".format(args.le_rhc_tag),
        "ME RHC #bar{{#nu}}_{{#mu}}/#bar{{#nu}}_{{e}} ({})".format(args.me_rhc_tag),
        "{}_rhc.png".format(args.out_prefix),
        "{}_rhc.pdf".format(args.out_prefix),
        xmax=args.xmax,
    )

    fout = ROOT.TFile("{}.root".format(args.out_prefix), "RECREATE")
    r_fhc_le.Write("fhc_le_numu_over_nue")
    r_fhc_me.Write("fhc_me_numu_over_nue")
    r_rhc_le.Write("rhc_le_numubar_over_nuebar")
    r_rhc_me.Write("rhc_me_numubar_over_nuebar")
    fout.Close()

    print("Wrote {}.root".format(args.out_prefix))


if __name__ == "__main__":
    main()