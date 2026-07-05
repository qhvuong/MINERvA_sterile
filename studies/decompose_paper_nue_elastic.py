#!/usr/bin/env python

import argparse
import array
import math
import ROOT
import PlotUtils

ROOT.TH1.AddDirectory(False)

PAPER_NUE_ELASTIC_POT = 3.43e20
FHC_ANALYSIS_POT = 3.323050142731963e20

# -------------------------
# Published MINERvA nu+e table
# -------------------------

PAPER_EDGES = [0.8, 2.0, 3.0, 5.0, 7.0, 9.0, 20.0]

PAPER_COUNTS = [48.7, 14.4, 20.5, 18.1, 11.9, 21.6]

PAPER_COV = [
    [98.7,   1.22,  1.72,  1.38,  0.420, -0.269],
    [1.22,  27.3,   1.63,  1.14,  0.340,  0.755],
    [1.72,   1.63, 40.1,   1.88,  0.596,  1.35 ],
    [1.38,   1.14,  1.88, 34.7,   0.448,  0.968],
    [0.420,  0.340, 0.596, 0.448,18.9,    0.778],
    [-0.269, 0.755, 1.35,  0.968, 0.778, 59.5  ],
]


SPECIES_1D = {
    "nue":     "Eel_NuEElasticE",
    "nuebar":  "Eel_NuEElasticEBar",
    "numu":    "Eel_NuEElasticMu",
    "numubar": "Eel_NuEElasticMuBar",
}

SPECIES_2D = {
    "nue":     "ElepReco_LE_NuEElasticE",
    "nuebar":  "ElepReco_LE_NuEElasticEBar",
    "numu":    "ElepReco_LE_NuEElasticMu",
    "numubar": "ElepReco_LE_NuEElasticMuBar",
}

def get_pot_from_meta(fin):
    """
    Read POT from the Meta tree in the input ROOT file.
    Tries common branch names and falls back to printing available branches.
    """
    meta = fin.Get("Meta")
    if not meta:
        raise RuntimeError("Input file has no Meta tree, so source POT must be provided manually.")

    branch_names = [b.GetName() for b in meta.GetListOfBranches()]
    print("\n=== Meta tree branches ===")
    for b in branch_names:
        print(" ", b)

    # Try common POT branch names.
    candidates = [
        "POT_Used",
        "POTUsed",
        "used_POT",
        "UsedPOT",
        "POT",
        "pot",
        "POT_Total",
        "total_pot",
    ]

    for name in candidates:
        if name in branch_names:
            total = 0.0
            for i in range(meta.GetEntries()):
                meta.GetEntry(i)
                total += float(getattr(meta, name))
            print("\nRead source POT from Meta branch:", name)
            print("  source POT =", total)
            return total

    raise RuntimeError(
        "Could not find a POT branch in Meta. Available branches: {}".format(branch_names)
    )

def make_h1(name, title, edges):
    arr = array.array("d", edges)
    h = PlotUtils.MnvH1D(name, title, len(edges) - 1, arr)
    h.SetDirectory(0)
    return h


def get_required(fin, name):
    obj = fin.Get(name)
    if not obj:
        raise RuntimeError("Missing object in input file: {}".format(name))
    obj = obj.Clone(name + "_clone")
    obj.SetDirectory(0)
    return obj


def rebin_by_bin_center(h_in, out_name, edges):
    """
    Rebin h_in into arbitrary paper bins by assigning each input-bin center
    to an output bin.

    This is fine if the source bins are finer than the paper bins.
    If the input bins are coarse and do not align with the paper bins,
    replace this with an overlap-based rebin.
    """
    h_out = make_h1(
        out_name,
        "{};Reco E_{{e}} [GeV];Events".format(out_name),
        edges,
    )

    nb = h_out.GetNbinsX()

    for i in range(1, h_in.GetNbinsX() + 1):
        x = h_in.GetBinCenter(i)
        y = h_in.GetBinContent(i)
        e = h_in.GetBinError(i)

        j = h_out.FindBin(x)
        if 1 <= j <= nb:
            old = h_out.GetBinContent(j)
            olde = h_out.GetBinError(j)

            h_out.SetBinContent(j, old + y)
            h_out.SetBinError(j, math.sqrt(olde * olde + e * e))

    return h_out


def make_paper_total_hist(edges, pot_scale=1.0):
    h = make_h1(
        "paper_total_nue_elastic",
        "Published MINERvA #nu+e elastic;Reco E_{e} [GeV];Events",
        edges,
    )

    for i, val in enumerate(PAPER_COUNTS, start=1):
        h.SetBinContent(i, val * pot_scale)
        h.SetBinError(i, math.sqrt(PAPER_COV[i - 1][i - 1]) * pot_scale)

    return h


def make_cov_matrix(pot_scale=1.0):
    n = len(PAPER_COUNTS)
    m = ROOT.TMatrixD(n, n)
    for i in range(n):
        for j in range(n):
            m[i][j] = PAPER_COV[i][j] * pot_scale * pot_scale
    return m


def decompose_1d(fin, fout, edges, paper_pot_scale=1.0, source_pot_scale=1.0):
    """
    Use the ME->LE weighted 1D species histograms to decompose the
    published paper total into nue/nuebar/numu/numubar components.
    """

    h_paper = make_paper_total_hist(edges, pot_scale=paper_pot_scale)

    # Read and rebin species MC/reference hists
    h_species_rebinned = {}
    for sp, hname in SPECIES_1D.items():
        h_raw = get_required(fin, hname)
        h_raw.Scale(source_pot_scale)
        h_species_rebinned[sp] = rebin_by_bin_center(
            h_raw,
            "source_{}_paperbins".format(sp),
            edges,
        )

    # Total source prediction in paper bins
    h_source_total = make_h1(
        "source_total_paperbins",
        "ME#rightarrowLE weighted source total;Reco E_{e} [GeV];Events",
        edges,
    )

    for sp, h in h_species_rebinned.items():
        h_source_total.Add(h)

    # Build fraction hists and paper-normalized species hists
    h_frac = {}
    h_decomp = {}

    for sp in SPECIES_1D:
        h_frac[sp] = make_h1(
            "frac_{}".format(sp),
            "{} fraction;Reco E_{{e}} [GeV];Fraction".format(sp),
            edges,
        )
        h_decomp[sp] = make_h1(
            "paper_decomp_{}".format(sp),
            "Paper-normalized {};Reco E_{{e}} [GeV];Events".format(sp),
            edges,
        )

    print("\n=== Paper decomposition ===")
    print(
        "{:<8s} {:>12s} {:>12s} {:>12s} {:>12s} {:>12s}".format(
            "bin", "paper", "nue", "nuebar", "numu", "numubar"
        )
    )

    for i in range(1, h_paper.GetNbinsX() + 1):
        n_paper = h_paper.GetBinContent(i)
        n_source = h_source_total.GetBinContent(i)

        if n_source <= 0:
            print("WARNING: source total is <= 0 in paper bin", i)
            continue

        row = [i, n_paper]

        for sp in SPECIES_1D:
            n_sp_source = h_species_rebinned[sp].GetBinContent(i)
            frac = n_sp_source / n_source
            n_sp_paper = n_paper * frac

            h_frac[sp].SetBinContent(i, frac)
            h_frac[sp].SetBinError(i, 0.0)

            h_decomp[sp].SetBinContent(i, n_sp_paper)

            # Approximate error assignment: same fractional error as paper total.
            # The full covariance belongs to the total table spectrum.
            paper_err = h_paper.GetBinError(i)
            h_decomp[sp].SetBinError(i, paper_err * frac)

            row.append(n_sp_paper)

        print(
            "{:<8d} {:12.4f} {:12.4f} {:12.4f} {:12.4f} {:12.4f}".format(
                row[0], row[1], row[2], row[3], row[4], row[5]
            )
        )

    # Check closure
    h_decomp_sum = make_h1(
        "paper_decomp_sum",
        "Sum of decomposed species;Reco E_{e} [GeV];Events",
        edges,
    )

    for sp in SPECIES_1D:
        h_decomp_sum.Add(h_decomp[sp])

    print("\n=== Closure check: sum(species) - paper ===")
    for i in range(1, h_paper.GetNbinsX() + 1):
        diff = h_decomp_sum.GetBinContent(i) - h_paper.GetBinContent(i)
        print("bin {:d}: diff = {:.6e}".format(i, diff))

    # Write output
    fout.cd()

    h_paper.Write()
    h_source_total.Write()
    h_decomp_sum.Write()

    for sp in SPECIES_1D:
        h_species_rebinned[sp].Write()
        h_frac[sp].Write()
        h_decomp[sp].Write()

    cov = make_cov_matrix(pot_scale=paper_pot_scale)
    cov.Write("paper_covariance_matrix")

    return h_paper, h_source_total, h_frac, h_decomp


def write_2d_templates_scaled_to_fhc(fin, fout, source_pot_scale=1.0):
    """
    Write ME->LE 2D templates scaled only by LElike source POT -> target FHC POT.

    These are NOT paper-normalized. They keep the source 2D L/E shape and
    normalization at the target FHC exposure.
    """
    for sp, hname in SPECIES_2D.items():
        h2_raw = get_required(fin, hname)

        h2_scaled = h2_raw.Clone("source_fhc_2d_{}".format(sp))
        h2_scaled.SetDirectory(0)
        h2_scaled.Scale(source_pot_scale)

        fout.cd()
        h2_scaled.Write()

        print(
            "Wrote source_fhc_2d_{}: raw integral = {}, scaled integral = {}".format(
                sp, h2_raw.Integral(), h2_scaled.Integral()
            )
        )


def main():
    parser = argparse.ArgumentParser(
        description="Decompose published MINERvA nu+e table into species using ME->LE weighted histograms."
    )
    parser.add_argument(
        "-i",
        "--input",
        default="/exp/minerva/data/users/qvuong/elastic_nue/kin_dist_mcmeFHC_p6_scattering_LElike_MAD.root",
        help="Input ROOT file with ME->LE weighted species histograms.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="paper_nue_elastic_decomposed.root",
        help="Output ROOT file.",
    )
    parser.add_argument(
        "--make-2d",
        action="store_true",
        help="Also write source 2D templates scaled only to target FHC POT.",
    )
    parser.add_argument(
        "--source-pot",
        type=float,
        default=None,
        help="POT/exposure represented by the input LElike ROOT file. If omitted, read from Meta tree.",
    )
    parser.add_argument(
        "--target-pot",
        type=float,
        default=FHC_ANALYSIS_POT,
        help="Target FHC POT to scale both paper and LElike source to.",
    )
    args = parser.parse_args()

    fin = ROOT.TFile.Open(args.input)
    if not fin or fin.IsZombie():
        raise RuntimeError("Could not open input file: {}".format(args.input))

    fout = ROOT.TFile.Open(args.output, "RECREATE")
    if not fout or fout.IsZombie():
        raise RuntimeError("Could not create output file: {}".format(args.output))

    source_pot = args.source_pot
    if source_pot is None:
        source_pot = get_pot_from_meta(fin)

    paper_pot_scale = args.target_pot / PAPER_NUE_ELASTIC_POT
    source_pot_scale = args.target_pot / source_pot

    print("\n=== POT scaling ===")
    print("paper POT       =", PAPER_NUE_ELASTIC_POT)
    print("source POT      =", source_pot)
    print("target FHC POT  =", args.target_pot)
    print("paper scale     =", paper_pot_scale)
    print("source scale    =", source_pot_scale)

    h_paper, h_source_total, h_frac, h_decomp = decompose_1d(
        fin,
        fout,
        PAPER_EDGES,
        paper_pot_scale=paper_pot_scale,
        source_pot_scale=source_pot_scale,
    )

    if args.make_2d:
        write_2d_templates_scaled_to_fhc(
            fin,
            fout,
            source_pot_scale=source_pot_scale,
        )

    fout.Close()
    fin.Close()

    print("\nWrote:", args.output)


if __name__ == "__main__":
    main()