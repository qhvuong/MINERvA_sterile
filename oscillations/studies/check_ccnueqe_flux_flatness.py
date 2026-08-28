#!/usr/bin/env python3

import ROOT
import argparse
import math


def check_hist(filename, histname="EN4_CCNuEQE",
               bandname="Flux", universes=(0, 1, 2, 3, 4)):

    f = ROOT.TFile.Open(filename)
    if not f or f.IsZombie():
        raise RuntimeError("Could not open {}".format(filename))

    h = f.Get(histname)
    if not h:
        raise RuntimeError(
            "Could not find histogram '{}' in {}".format(
                histname, filename
            )
        )

    if not h.HasVertErrorBand(bandname):
        raise RuntimeError(
            "Histogram '{}' has no vertical error band '{}'".format(
                histname, bandname
            )
        )

    band = h.GetVertErrorBand(bandname)

    print("\n============================================================")
    print("FILE      :", filename)
    print("HIST      :", histname)
    print("BAND      :", bandname)
    print("N UNIV    :", band.GetNHists())
    print("CV INT    :", h.Integral())
    print("============================================================")

    for u in universes:

        if u >= band.GetNHists():
            print("\nUniverse {} does not exist, skipping.".format(u))
            continue

        hu = band.GetHist(u)

        ratios = []

        print("\n------------------------------------------------------------")
        print("UNIVERSE {}".format(u))
        print("------------------------------------------------------------")
        print(
            "{:>4s} {:>10s} {:>10s} {:>14s} {:>14s} {:>14s}".format(
                "bin",
                "low",
                "high",
                "CV",
                "universe",
                "U/CV",
            )
        )

        for b in range(1, h.GetNbinsX() + 1):

            cv = h.GetBinContent(b)
            uv = hu.GetBinContent(b)

            low = h.GetXaxis().GetBinLowEdge(b)
            high = h.GetXaxis().GetBinUpEdge(b)

            if cv != 0.0:
                ratio = uv / cv
                ratios.append(ratio)

                print(
                    "{:4d} {:10.4f} {:10.4f} "
                    "{:14.8g} {:14.8g} {:14.8f}".format(
                        b,
                        low,
                        high,
                        cv,
                        uv,
                        ratio,
                    )
                )

            elif uv != 0.0:
                print(
                    "{:4d} {:10.4f} {:10.4f} "
                    "{:14.8g} {:14.8g} {:>14s}".format(
                        b,
                        low,
                        high,
                        cv,
                        uv,
                        "CV=0",
                    )
                )

        if len(ratios) > 1:

            mean = sum(ratios) / len(ratios)

            variance = sum(
                (x - mean) ** 2
                for x in ratios
            ) / len(ratios)

            std = math.sqrt(variance)

            rmin = min(ratios)
            rmax = max(ratios)

            print("\nSummary for universe {}:".format(u))
            print("  mean(U/CV) = {:.10f}".format(mean))
            print("  std bins   = {:.10e}".format(std))
            print("  min        = {:.10f}".format(rmin))
            print("  max        = {:.10f}".format(rmax))
            print("  max-min    = {:.10e}".format(rmax - rmin))

            if std < 1.0e-8:
                print("  RESULT     = FLAT")
            else:
                print("  RESULT     = NOT FLAT")

    f.Close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=(
            "Check universe/CV flatness for a raw CCNuEQE histogram."
        )
    )

    parser.add_argument(
        "file",
        help="ROOT selection output file"
    )

    parser.add_argument(
        "--hist",
        default="EN4_CCNuEQE",
        help="Histogram name (default: EN4_CCNuEQE)"
    )

    parser.add_argument(
        "--band",
        default="Flux",
        help="Error-band name in selection output (default: Flux)"
    )

    parser.add_argument(
        "--universes",
        nargs="+",
        type=int,
        default=[0, 1, 2, 3, 4],
        help="Universe indices to inspect"
    )

    args = parser.parse_args()

    check_hist(
        args.file,
        histname=args.hist,
        bandname=args.band,
        universes=args.universes,
    )