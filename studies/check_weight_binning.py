#!/usr/bin/env python3

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

SAMPLES = {
    "FHC": (
        "/exp/minerva/data/users/qvuong/nu_mu",
        "leFHC",
        "CCnumu",
    ),
    "RHC": (
        "/exp/minerva/data/users/qvuong/antinu_mu",
        "le5",
        "CCnumubar",
    ),
}

BASELINE = "v1_D2_LECOH_DiffCOH"

TESTS = {
    "JOINT": (
        "v1_D2_JOINT_LECOH_DiffCOH",
        "JOINTWeight_trueEnu",
    ),
    "MENU1PI": (
        "v4p3p0_D2_MENU1PI_LECOH_DiffCOH",
        "MENU1PIWeight_trueEnu",
    ),
}


def open_file(mode, tag):
    directory, playlist, sample = SAMPLES[mode]

    path = (
        f"{directory}/kin_dist_mc{playlist}_"
        f"Validation_{sample}_{tag}_MAD.root"
    )

    f = ROOT.TFile.Open(path)

    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open {path}")

    return f


def get_hist(f, name):
    obj = f.Get(name)

    if not obj:
        raise RuntimeError(
            f"Missing histogram {name} in {f.GetName()}"
        )

    if hasattr(obj, "GetCVHistoWithStatError"):
        h = obj.GetCVHistoWithStatError()
    else:
        h = obj.Clone()

    h.SetDirectory(0)
    return h


def integrate_energy(h, low=None, high=None):
    total = 0.0

    for ix in range(1, h.GetNbinsX() + 1):
        xlow = h.GetXaxis().GetBinLowEdge(ix)
        xhigh = h.GetXaxis().GetBinUpEdge(ix)

        if low is not None and xlow < low - 1e-7:
            continue

        if high is not None and xhigh > high + 1e-7:
            continue

        total += h.GetBinContent(ix)

    return total


def audit(mode, test, tag, map_name):
    fbase = open_file(mode, BASELINE)
    ftune = open_file(mode, tag)

    try:
        hbase = get_hist(fbase, "Flux")
        htune = get_hist(ftune, "Flux")
        hmap = get_hist(fbase, map_name)

        yaxis = hmap.GetYaxis()
        unity_bin = yaxis.FindBin(1.0)

        if not 1 <= unity_bin <= hmap.GetNbinsY():
            raise RuntimeError("Unity falls outside weight axis")

        unity_low = yaxis.GetBinLowEdge(unity_bin)
        unity_high = yaxis.GetBinUpEdge(unity_bin)
        unity_center = yaxis.GetBinCenter(unity_bin)

        print("\n" + "=" * 90)
        print(mode, test)
        print("=" * 90)

        print("Weight map:", map_name)
        print("Number of Y bins:", hmap.GetNbinsY())
        print("Weight axis range:", yaxis.GetXmin(), yaxis.GetXmax())

        print("\nWeight bin containing 1.0:")
        print("  Bin:", unity_bin)
        print("  Low edge:", unity_low)
        print("  High edge:", unity_high)
        print("  Center:", unity_center)
        print("  Center - 1:", unity_center - 1.0)

        # Count baseline-weighted yield in the unity bin.
        unity_yield = 0.0
        total_map_yield = 0.0

        for ix in range(1, hmap.GetNbinsX() + 1):
            for iy in range(1, hmap.GetNbinsY() + 1):
                value = hmap.GetBinContent(ix, iy)
                total_map_yield += value

                if iy == unity_bin:
                    unity_yield += value

        print("\nPopulation in bin containing unity:")
        print(f"  Weighted yield: {unity_yield:.3f}")
        print(
            f"  Fraction of map: "
            f"{100 * unity_yield / total_map_yield:.3f}%"
        )

        naive_unity_bias = unity_yield * (unity_center - 1.0)

        print(
            f"  Potential bin-center contribution: "
            f"{naive_unity_bias:.3f}"
        )

        regions = [
            ("Full visible spectrum", None, None),
            ("Flux peak 2.5-3.8 GeV", 2.5, 3.8),
            ("2.0-5.0 GeV", 2.0, 5.0),
            ("Above 5.0 GeV", 5.0, None),
        ]

        for region_name, emin, emax in regions:

            baseline = integrate_energy(hbase, emin, emax)
            observed = integrate_energy(htune, emin, emax)

            naive_delta = 0.0
            snapped_delta = 0.0
            unity_region_yield = 0.0

            map_region_yield = 0.0

            for ix in range(1, hmap.GetNbinsX() + 1):

                xlow = hmap.GetXaxis().GetBinLowEdge(ix)
                xhigh = hmap.GetXaxis().GetBinUpEdge(ix)

                if emin is not None and xlow < emin - 1e-7:
                    continue

                if emax is not None and xhigh > emax + 1e-7:
                    continue

                for iy in range(1, hmap.GetNbinsY() + 1):

                    content = hmap.GetBinContent(ix, iy)
                    weight = yaxis.GetBinCenter(iy)

                    map_region_yield += content

                    naive_delta += content * (weight - 1.0)

                    if iy == unity_bin:
                        unity_region_yield += content
                        snapped_weight = 1.0
                    else:
                        snapped_weight = weight

                    snapped_delta += content * (
                        snapped_weight - 1.0
                    )

            observed_delta = observed - baseline

            naive_residual = observed_delta - naive_delta
            snapped_residual = observed_delta - snapped_delta

            correction = (
                unity_region_yield * (unity_center - 1.0)
            )

            print("\n" + "-" * 75)
            print(region_name)
            print("-" * 75)

            print(f"Baseline yield:       {baseline:15.3f}")
            print(f"Observed yield:       {observed:15.3f}")

            print(f"Observed delta:       {observed_delta:15.3f}")
            print(f"Naive predicted delta:{naive_delta:15.3f}")
            print(f"Snap predicted delta: {snapped_delta:15.3f}")

            print(f"Naive residual:       {naive_residual:15.3f}")
            print(f"Snap residual:        {snapped_residual:15.3f}")

            print(f"Unity-bin correction: {correction:15.3f}")
            print(f"Unity-bin yield:      {unity_region_yield:15.3f}")

            if baseline:
                print(
                    "Naive residual / baseline: "
                    f"{100 * naive_residual / baseline:9.4f}%"
                )
                print(
                    "Snap residual / baseline:  "
                    f"{100 * snapped_residual / baseline:9.4f}%"
                )

            if observed_delta:
                print(
                    "Naive residual / observed delta: "
                    f"{100 * naive_residual / observed_delta:9.3f}%"
                )
                print(
                    "Snap residual / observed delta:  "
                    f"{100 * snapped_residual / observed_delta:9.3f}%"
                )

            print(
                "Map coverage: "
                f"{100 * map_region_yield / baseline:.5f}%"
                if baseline else "Map coverage: N/A"
            )

    finally:
        fbase.Close()
        ftune.Close()


for mode in ("FHC", "RHC"):
    for test, (tag, map_name) in TESTS.items():
        audit(mode, test, tag, map_name)