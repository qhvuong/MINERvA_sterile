#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path
from collections import OrderedDict

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.AddDirectory(False)


# ============================================================
# INPUT CONFIGURATION
# ============================================================

SAMPLES = {
    "FHC": {
        "directory": "/exp/minerva/data/users/qvuong/nu_mu",
        "playlist": "leFHC",
        "sample": "CCnumu",
    },
    "RHC": {
        "directory": "/exp/minerva/data/users/qvuong/antinu_mu",
        "playlist": "le5",
        "sample": "CCnumubar",
    },
}


# ============================================================
# EXACT TUNE TAGS
# ============================================================

TAGS = [
    "v1_D2_JOINT_LECOH_DiffCOH",
    "v1_D2_LECOH_DiffCOH",
    "v4_D2_MENU1PI_LECOH_DiffCOH",
    "v4_D2_MENU1PI_MECOH_DiffCOH",
    "v4p3p0_D2_MENU1PI_LECOH_DiffCOH",
    "v4p3p0_D2_MENU1PI_MECOH_DiffCOH",
]

REFERENCE = "v4p3p0_D2_MENU1PI_LECOH_DiffCOH"

PAIRS = [
    (
        "v1_D2_JOINT_LECOH_DiffCOH",
        "v1_D2_LECOH_DiffCOH",
    ),
    (
        "v4_D2_MENU1PI_LECOH_DiffCOH",
        "v4_D2_MENU1PI_MECOH_DiffCOH",
    ),
    (
        "v4p3p0_D2_MENU1PI_LECOH_DiffCOH",
        "v4p3p0_D2_MENU1PI_MECOH_DiffCOH",
    ),
    (
        "v4_D2_MENU1PI_LECOH_DiffCOH",
        "v4p3p0_D2_MENU1PI_LECOH_DiffCOH",
    ),
    (
        "v4_D2_MENU1PI_MECOH_DiffCOH",
        "v4p3p0_D2_MENU1PI_MECOH_DiffCOH",
    ),
]


# ============================================================
# HISTOGRAMS
# ============================================================

HISTS_1D = OrderedDict([
    ("true_n_piplus", "True #pi^{+} multiplicity"),
    ("true_n_piminus", "True #pi^{-} multiplicity"),
    ("true_n_pi0", "True #pi^{0} multiplicity"),
    ("true_n_pions", "Total true pion multiplicity"),
    ("true_n_beamsign_pion", "Beam-sign charged pion multiplicity"),

    ("Q2", "True Q^{2} (GeV^{2})"),

    ("mc_w_raw", "Raw true W (GeV)"),
    ("trueW", "True W (GeV)"),
    ("trueWExp", "True W exp (GeV)"),

    ("mc_intType", "GENIE interaction type"),

    ("MENU1PI_eligible", "MENU1PI eligibility"),
    ("JOINT_eligible", "JOINT eligibility"),

    ("MENU1PI_cv_weight", "MENU1PI CV weight"),
    ("JOINT_cv_weight", "JOINT CV weight"),

    ("Flux", "True neutrino energy (GeV)"),
    ("E_avail", "Reconstructed available energy (GeV)"),
    ("EN4", "Reconstructed neutrino energy estimator (GeV)"),
])

HISTS_2D = OrderedDict([
    ("trueEnu_trueQ2", "True E_{#nu} vs Q^{2}"),
    ("trueW_trueQ2", "True W vs Q^{2}"),

    ("MENU1PIWeight_trueQ2", "MENU1PI weight vs Q^{2}"),
    ("JOINTWeight_trueQ2", "JOINT weight vs Q^{2}"),

    ("MENU1PIWeight_trueEnu", "MENU1PI weight vs E_{#nu}"),
    ("JOINTWeight_trueEnu", "JOINT weight vs E_{#nu}"),
])


# ============================================================
# STYLE
# ============================================================

COLORS = [
    ROOT.kMagenta + 2,
    ROOT.kViolet + 1,
    ROOT.kOrange + 7,
    ROOT.kRed + 1,
    ROOT.kBlack,
    ROOT.kAzure + 2,
]

TUNE_COLORS = dict(zip(TAGS, COLORS))


# ============================================================
# INPUT FILES
# ============================================================

def get_input_path(mode, tag):

    config = SAMPLES[mode]

    filename = (
        f"kin_dist_mc{config['playlist']}_"
        f"Validation_{config['sample']}_{tag}_MAD.root"
    )

    return Path(config["directory"]) / filename


def open_files(mode):

    files = OrderedDict()

    print(f"\n{'=' * 80}")
    print(f"Opening {mode} files")
    print(f"{'=' * 80}")

    for tag in TAGS:

        path = get_input_path(mode, tag)

        if not path.exists():
            print(f"[MISSING] {path}")
            continue

        f = ROOT.TFile.Open(str(path), "READ")

        if not f or f.IsZombie():
            print(f"[INVALID] {path}")
            continue

        files[tag] = f
        print(f"[FOUND] {path.name}")

    return files


# ============================================================
# POT NORMALIZATION
# ============================================================

def get_pot(f):

    meta = f.Get("Meta")

    if not meta or not meta.GetBranch("POT_Used"):
        return None

    pot = 0.0

    for event in meta:
        pot += float(event.POT_Used)

    return pot if pot > 0 else None


def get_scale_factors(files, pot_scale):

    pots = {
        tag: get_pot(f)
        for tag, f in files.items()
    }

    factors = {tag: 1.0 for tag in files}

    if pot_scale:

        reference_pot = pots.get(REFERENCE)

        if reference_pot is None:
            print("[WARNING] Reference POT unavailable; no POT scaling")
        else:
            for tag in files:

                if pots[tag] is not None:
                    factors[tag] = reference_pot / pots[tag]
                else:
                    print(f"[WARNING] Missing POT for {tag}")

    return pots, factors


# ============================================================
# HISTOGRAM ACCESS
# ============================================================

def get_histogram(f, name, dimension=1, scale=1.0):

    obj = f.Get(name)

    if not obj:
        return None

    if hasattr(obj, "GetCVHistoWithStatError"):
        hist = obj.GetCVHistoWithStatError()
    else:
        hist = obj.Clone()

    if not hist:
        return None

    actual_dimension = 2 if hist.InheritsFrom("TH2") else 1

    if actual_dimension != dimension:
        return None

    hist.SetDirectory(0)

    if scale != 1.0:
        hist.Scale(scale)

    return hist


def compatible(h1, h2):

    if h1.GetNbinsX() != h2.GetNbinsX():
        return False

    for i in range(1, h1.GetNbinsX() + 2):
        if abs(
            h1.GetXaxis().GetBinLowEdge(i)
            - h2.GetXaxis().GetBinLowEdge(i)
        ) > 1e-6:
            return False

    if h1.InheritsFrom("TH2") or h2.InheritsFrom("TH2"):

        if not (h1.InheritsFrom("TH2") and h2.InheritsFrom("TH2")):
            return False

        if h1.GetNbinsY() != h2.GetNbinsY():
            return False

        for i in range(1, h1.GetNbinsY() + 2):
            if abs(
                h1.GetYaxis().GetBinLowEdge(i)
                - h2.GetYaxis().GetBinLowEdge(i)
            ) > 1e-6:
                return False

    return True


# ============================================================
# DRAW 1D OVERLAY + RATIO
# ============================================================

def draw_overlay(name, title, hists, reference_tag, output, shape=False):

    if reference_tag not in hists:
        return

    reference_original = hists[reference_tag]

    drawing_hists = OrderedDict()

    for tag, h in hists.items():

        if not compatible(h, reference_original):
            print(f"[BINNING MISMATCH] {name}: {tag}")
            continue

        clone = h.Clone(f"{name}_{tag}_drawing")
        clone.SetDirectory(0)

        if shape and clone.Integral() > 0:
            clone.Scale(1.0 / clone.Integral())

        clone.SetLineColor(TUNE_COLORS[tag])
        clone.SetLineWidth(3)
        clone.SetLineStyle(1 if tag == reference_tag else 2)

        drawing_hists[tag] = clone

    if len(drawing_hists) < 2:
        return

    reference = drawing_hists[reference_tag]

    canvas = ROOT.TCanvas("canvas", "canvas", 1400, 900)

    upper = ROOT.TPad("upper", "upper", 0, 0.30, 1, 1)
    lower = ROOT.TPad("lower", "lower", 0, 0, 1, 0.30)

    upper.SetLeftMargin(0.12)
    upper.SetRightMargin(0.05)
    upper.SetBottomMargin(0.03)

    lower.SetLeftMargin(0.12)
    lower.SetRightMargin(0.05)
    lower.SetTopMargin(0.04)
    lower.SetBottomMargin(0.28)

    upper.Draw()
    lower.Draw()

    # --------------------------------------------------------
    # Upper panel
    # --------------------------------------------------------

    upper.cd()

    reference.SetTitle(title)

    ymax = max(h.GetMaximum() for h in drawing_hists.values())

    reference.SetMinimum(0)
    reference.SetMaximum(ymax * 1.40 if ymax > 0 else 1)

    reference.GetYaxis().SetTitle(
        "Normalized fraction" if shape else "Weighted events / bin"
    )

    reference.GetXaxis().SetLabelSize(0)

    reference.Draw("HIST")

    for tag, h in drawing_hists.items():
        if tag != reference_tag:
            h.Draw("HIST SAME")

    legend = ROOT.TLegend(0.45, 0.59, 0.95, 0.91)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.024)

    for tag, h in drawing_hists.items():
        legend.AddEntry(h, tag, "l")

    legend.Draw()

    # MENU1PI W applicability boundary
    w_line = None

    if name in ("mc_w_raw", "trueW", "trueWExp"):

        w_line = ROOT.TLine(
            1.4,
            0,
            1.4,
            reference.GetMaximum()
        )

        w_line.SetLineStyle(3)
        w_line.SetLineWidth(2)
        w_line.Draw()

    # --------------------------------------------------------
    # Lower ratio panel
    # --------------------------------------------------------

    lower.cd()

    frame = reference.Clone("ratio_frame")
    frame.Reset("ICES")

    frame.SetTitle("")
    frame.SetMinimum(0)
    frame.SetMaximum(2)

    frame.GetYaxis().SetTitle("Tune / ref")
    frame.GetYaxis().SetTitleSize(0.11)
    frame.GetYaxis().SetLabelSize(0.09)
    frame.GetYaxis().SetTitleOffset(0.47)

    frame.GetXaxis().SetTitle(title)
    frame.GetXaxis().SetTitleSize(0.11)
    frame.GetXaxis().SetLabelSize(0.10)

    frame.Draw("AXIS")

    unity = ROOT.TLine(
        frame.GetXaxis().GetXmin(),
        1,
        frame.GetXaxis().GetXmax(),
        1,
    )

    unity.SetLineStyle(2)
    unity.Draw()

    graphs = []

    for tag, h in drawing_hists.items():

        if tag == reference_tag:
            continue

        points = []

        for i in range(1, h.GetNbinsX() + 1):

            denominator = reference.GetBinContent(i)

            if denominator <= 0:
                continue

            numerator = h.GetBinContent(i)

            points.append((
                h.GetBinCenter(i),
                numerator / denominator,
            ))

        if not points:
            continue

        graph = ROOT.TGraph(len(points))

        for i, (x, y) in enumerate(points):
            graph.SetPoint(i, x, y)

        graph.SetLineColor(TUNE_COLORS[tag])
        graph.SetLineWidth(3)
        graph.Draw("L SAME")

        graphs.append(graph)

    output.parent.mkdir(parents=True, exist_ok=True)

    canvas.Print(str(output))
    canvas.Close()


# ============================================================
# DRAW 2D HISTOGRAMS
# ============================================================

def draw_2d(hist, title, output, reference=None):

    h = hist.Clone("hist2d_drawing")
    h.SetDirectory(0)

    if reference is not None:

        if not compatible(h, reference):
            return

        for ix in range(1, h.GetNbinsX() + 1):
            for iy in range(1, h.GetNbinsY() + 1):

                numerator = h.GetBinContent(ix, iy)
                denominator = reference.GetBinContent(ix, iy)

                if denominator > 0:
                    value = numerator / denominator - 1.0
                else:
                    value = 0.0

                h.SetBinContent(ix, iy, value)
                h.SetBinError(ix, iy, 0)

        h.SetMinimum(-0.5)
        h.SetMaximum(0.5)

    canvas = ROOT.TCanvas("canvas2d", "canvas2d", 1100, 850)

    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.17)
    canvas.SetBottomMargin(0.12)

    h.SetTitle(title)
    h.Draw("COLZ")

    output.parent.mkdir(parents=True, exist_ok=True)

    canvas.Print(str(output))
    canvas.Close()


# ============================================================
# COMPARE ONE BEAM MODE
# ============================================================

def process_mode(mode, args):

    output = Path(args.output) / mode
    output.mkdir(parents=True, exist_ok=True)

    files = open_files(mode)

    if REFERENCE not in files:
        print(f"[ERROR] Reference tune missing for {mode}")
        for f in files.values():
            f.Close()
        return

    try:

        pots, factors = get_scale_factors(files, args.pot_scale)

        # ----------------------------------------------------
        # POT summary
        # ----------------------------------------------------

        with open(output / "input_pot.csv", "w", newline="") as fp:

            writer = csv.writer(fp)

            writer.writerow([
                "tag",
                "POT_Used",
                "scale_factor",
            ])

            for tag in files:
                writer.writerow([
                    tag,
                    pots[tag],
                    factors[tag],
                ])

        # ----------------------------------------------------
        # Histogram inventory
        # ----------------------------------------------------

        with open(
            output / "histogram_inventory.csv",
            "w",
            newline=""
        ) as fp:

            writer = csv.writer(fp)

            writer.writerow([
                "tag",
                "histogram",
                "dimension",
                "found",
            ])

            for tag, f in files.items():

                for dimension, histograms in (
                    (1, HISTS_1D),
                    (2, HISTS_2D),
                ):

                    for name in histograms:

                        writer.writerow([
                            tag,
                            name,
                            dimension,
                            bool(f.Get(name)),
                        ])

        # ----------------------------------------------------
        # 1D histograms
        # ----------------------------------------------------

        with open(
            output / "hist_yields.csv",
            "w",
            newline=""
        ) as fp:

            writer = csv.writer(fp)

            writer.writerow([
                "histogram",
                "tag",
                "visible_yield",
                "yield_including_flow",
            ])

            for name, title in HISTS_1D.items():

                hists = OrderedDict()

                for tag, f in files.items():

                    h = get_histogram(
                        f,
                        name,
                        dimension=1,
                        scale=factors[tag],
                    )

                    if h is not None:
                        hists[tag] = h

                if REFERENCE not in hists:
                    print(f"[MISSING 1D HIST] {mode}: {name}")
                    continue

                for tag, h in hists.items():

                    writer.writerow([
                        name,
                        tag,
                        h.Integral(),
                        h.Integral(0, h.GetNbinsX() + 1),
                    ])

                # All six tune comparisons
                draw_overlay(
                    name,
                    title,
                    hists,
                    REFERENCE,
                    output / "overlays" / f"{name}.png",
                    shape=args.shape,
                )

                # Pairwise comparisons
                if args.pairs:

                    for first, second in PAIRS:

                        if first not in hists or second not in hists:
                            continue

                        pair_hists = OrderedDict([
                            (first, hists[first]),
                            (second, hists[second]),
                        ])

                        draw_overlay(
                            name,
                            title,
                            pair_hists,
                            first,
                            output
                            / "pairs"
                            / f"{second}_over_{first}"
                            / f"{name}.png",
                            shape=args.shape,
                        )

                # Eligibility fraction
                if name in ("MENU1PI_eligible", "JOINT_eligible"):

                    for tag, h in hists.items():

                        n0 = h.GetBinContent(h.FindBin(0))
                        n1 = h.GetBinContent(h.FindBin(1))

                        total = n0 + n1

                        if total > 0:
                            fraction = n1 / total

                            print(
                                f"[{mode}] {name} | {tag} | "
                                f"eligible fraction = {fraction:.6f}"
                            )

        # ----------------------------------------------------
        # 2D histograms
        # ----------------------------------------------------

        for name, title in HISTS_2D.items():

            hists = OrderedDict()

            for tag, f in files.items():

                h = get_histogram(
                    f,
                    name,
                    dimension=2,
                    scale=factors[tag],
                )

                if h is not None:
                    hists[tag] = h

            if REFERENCE not in hists:
                print(f"[MISSING 2D HIST] {mode}: {name}")
                continue

            for tag, h in hists.items():

                draw_2d(
                    h,
                    f"{title} | {tag}",
                    output / "maps" / name / f"{tag}.png",
                )

                if args.map_ratios and tag != REFERENCE:

                    draw_2d(
                        h,
                        f"{title} | {tag} / {REFERENCE} - 1",
                        output / "map_ratios" / name / f"{tag}.png",
                        reference=hists[REFERENCE],
                    )

        print(f"\n[DONE] {mode}: {output}")

    finally:

        for f in files.values():
            f.Close()


# ============================================================
# MAIN
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--mode",
        choices=["FHC", "RHC", "all"],
        default="all",
    )

    parser.add_argument(
        "--output",
        default="/exp/minerva/data/users/qvuong/studies/mnv_tune_comparisons",
    )

    parser.add_argument(
        "--pot-scale",
        action="store_true",
    )

    parser.add_argument(
        "--shape",
        action="store_true",
        help="Normalize each 1D histogram to unit area",
    )

    parser.add_argument(
        "--pairs",
        action="store_true",
    )

    parser.add_argument(
        "--map-ratios",
        action="store_true",
    )

    args = parser.parse_args()

    modes = (
        ["FHC", "RHC"]
        if args.mode == "all"
        else [args.mode]
    )

    for mode in modes:
        process_mode(mode, args)


if __name__ == "__main__":
    main()