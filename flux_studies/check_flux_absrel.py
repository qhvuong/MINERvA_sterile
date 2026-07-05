#!/usr/bin/env python3

import argparse
import math
import ROOT

ROOT.TH1.AddDirectory(False)


def is_bad(x):
    try:
        return not math.isfinite(float(x))
    except Exception:
        return True


def open_file(path):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError("Could not open file: {}".format(path))
    return f


def list_keys(f, prefix=""):
    out = []
    for key in f.GetListOfKeys():
        name = key.GetName()
        obj = key.ReadObj()
        full = prefix + name

        out.append((full, obj.ClassName()))

        if obj.InheritsFrom("TDirectory"):
            out.extend(list_keys(obj, full + "/"))

    return out


def get_obj(f, name):
    obj = f.Get(name)
    if not obj:
        raise RuntimeError("Could not find object '{}' in {}".format(name, f.GetName()))
    return obj


def get_bin_table(cv_hist, univ_hist, emin, emax, label=""):
    rows = []

    nb = cv_hist.GetNbinsX()
    for i in range(1, nb + 1):
        lo = cv_hist.GetBinLowEdge(i)
        hi = lo + cv_hist.GetBinWidth(i)
        width = cv_hist.GetBinWidth(i)

        if hi <= emin or lo >= emax:
            continue

        cv = cv_hist.GetBinContent(i)
        univ = univ_hist.GetBinContent(i)

        ratio = univ / cv if cv != 0 else float("nan")

        rows.append({
            "bin": i,
            "lo": lo,
            "hi": hi,
            "width": width,
            "cv": cv,
            "univ": univ,
            "ratio": ratio,
        })

    return rows


def median(vals):
    vals = [abs(v) for v in vals if not is_bad(v) and v != 0]
    if not vals:
        return None
    vals.sort()
    return vals[len(vals)//2]


def summarize_absrel(rows, label):
    cvs = [r["cv"] for r in rows]
    univs = [r["univ"] for r in rows]
    ratios = [r["ratio"] for r in rows if not is_bad(r["ratio"])]

    cv_med = median(cvs)
    univ_med = median(univs)
    ratio_med = median(ratios)

    print("\n[SUMMARY] {}".format(label))
    print("  median |CV|       =", cv_med)
    print("  median |universe| =", univ_med)
    print("  median |univ/CV|  =", ratio_med)

    if cv_med is None or univ_med is None:
        print("  classification    = UNKNOWN, missing useful nonzero values")
        return

    # Heuristic only. The bin table is the real evidence.
    if 0.2 <= univ_med <= 5.0 and not (0.2 <= cv_med <= 5.0):
        print("  classification    = universe looks RELATIVE-ish, CV not O(1)")
        print("  warning           = if you divide this universe by CV, you may create nonsense")
    elif cv_med > 0 and 0.05 <= (univ_med / cv_med) <= 20.0:
        print("  classification    = universe looks ABSOLUTE-ish, same scale as CV")
        print("  expected handling = relative universe should be universe / matching CV")
    else:
        print("  classification    = AMBIGUOUS")
        print("  note              = inspect bin table and compare to units / expected flux scale")


def print_rows(rows, label):
    print("\n[BINS] {}".format(label))
    print(
        "{:>5s} {:>11s} {:>11s} {:>11s} {:>16s} {:>16s} {:>16s}".format(
            "bin", "lo", "hi", "width", "CV", "universe", "univ/CV"
        )
    )

    for r in rows:
        print(
            "{bin:5d} {lo:11.5g} {hi:11.5g} {width:11.5g} "
            "{cv:16.8e} {univ:16.8e} {ratio:16.8e}".format(**r)
        )


def check_separate_hists(args):
    f = open_file(args.file)

    cv = get_obj(f, args.cv_hist)
    univs = [get_obj(f, u) for u in args.univ_hist]

    print("[FILE]", args.file)
    print("[CV]", args.cv_hist)

    for idx, uh in enumerate(univs):
        label = "separate_hist_univ{}_{}".format(idx, args.univ_hist[idx])
        rows = get_bin_table(cv, uh, args.emin, args.emax, label)
        print_rows(rows, label)
        summarize_absrel(rows, label)


def check_mnv_band(args):
    f = open_file(args.file)

    h = get_obj(f, args.hist)

    if not hasattr(h, "HasVertErrorBand") or not h.HasVertErrorBand(args.band):
        raise RuntimeError(
            "Histogram '{}' does not have vertical error band '{}'".format(
                args.hist, args.band
            )
        )

    band = h.GetVertErrorBand(args.band)
    nh = band.GetNHists()

    print("[FILE]", args.file)
    print("[MNV HIST]", args.hist)
    print("[BAND]", args.band)
    print("[N UNIVERSES]", nh)

    max_univ = min(args.nuniv, nh)

    # Check main CV vs band CV too. This catches CV mismatch / unsynced band CV.
    rows_cv = get_bin_table(h, band, args.emin, args.emax, "mainCV_vs_bandCV")
    print_rows(rows_cv, "mainCV_vs_bandCV")
    summarize_absrel(rows_cv, "mainCV_vs_bandCV")

    for iu in range(max_univ):
        uh = band.GetHist(iu)
        label = "{}_univ{}".format(args.band, iu)
        rows = get_bin_table(h, uh, args.emin, args.emax, label)
        print_rows(rows, label)
        summarize_absrel(rows, label)


def compare_two_mnv_files(args):
    f_old = open_file(args.old_file)
    f_new = open_file(args.new_file)

    h_old = get_obj(f_old, args.hist)
    h_new = get_obj(f_new, args.hist)

    print("[COMPARE]")
    print("  old =", args.old_file)
    print("  new =", args.new_file)
    print("  hist =", args.hist)

    print(
        "{:>5s} {:>11s} {:>11s} {:>16s} {:>16s} {:>16s}".format(
            "bin", "lo", "hi", "old_CV", "new_CV", "new/old"
        )
    )

    for i in range(1, h_old.GetNbinsX() + 1):
        lo = h_old.GetBinLowEdge(i)
        hi = lo + h_old.GetBinWidth(i)

        if hi <= args.emin or lo >= args.emax:
            continue

        old = h_old.GetBinContent(i)
        new = h_new.GetBinContent(i)
        ratio = new / old if old != 0 else float("nan")

        print(
            "{:5d} {:11.5g} {:11.5g} {:16.8e} {:16.8e} {:16.8e}".format(
                i, lo, hi, old, new, ratio
            )
        )

    if args.band:
        if not h_old.HasVertErrorBand(args.band) or not h_new.HasVertErrorBand(args.band):
            print("[WARN] one file does not have band {}; skipping universe compare".format(args.band))
            return

        b_old = h_old.GetVertErrorBand(args.band)
        b_new = h_new.GetVertErrorBand(args.band)

        nh = min(args.nuniv, b_old.GetNHists(), b_new.GetNHists())

        for iu in range(nh):
            u_old = b_old.GetHist(iu)
            u_new = b_new.GetHist(iu)

            print("\n[COMPARE_UNIV] universe", iu)
            print(
                "{:>5s} {:>11s} {:>11s} {:>16s} {:>16s} {:>16s}".format(
                    "bin", "lo", "hi", "old_univ", "new_univ", "new/old"
                )
            )

            for i in range(1, h_old.GetNbinsX() + 1):
                lo = h_old.GetBinLowEdge(i)
                hi = lo + h_old.GetBinWidth(i)

                if hi <= args.emin or lo >= args.emax:
                    continue

                old = u_old.GetBinContent(i)
                new = u_new.GetBinContent(i)
                ratio = new / old if old != 0 else float("nan")

                print(
                    "{:5d} {:11.5g} {:11.5g} {:16.8e} {:16.8e} {:16.8e}".format(
                        i, lo, hi, old, new, ratio
                    )
                )


def main():
    parser = argparse.ArgumentParser(
        description="Check whether flux universe hists look absolute or relative, and inspect 5-7 GeV region."
    )

    parser.add_argument("--file", help="ROOT file to inspect")
    parser.add_argument("--list", action="store_true", help="List ROOT keys and exit")

    parser.add_argument("--cv-hist", help="CV histogram name for separate-hist mode")
    parser.add_argument("--univ-hist", nargs="*", default=[], help="Universe histogram names for separate-hist mode")

    parser.add_argument("--hist", help="MnvH1D histogram name for band mode")
    parser.add_argument("--band", default="Flux", help="Vertical error band name, default Flux")
    parser.add_argument("--nuniv", type=int, default=3, help="Number of universes to print")

    parser.add_argument("--old-file", help="Old/reference produced ROOT file")
    parser.add_argument("--new-file", help="New/P8 produced ROOT file")

    parser.add_argument("--emin", type=float, default=4.5)
    parser.add_argument("--emax", type=float, default=7.5)

    args = parser.parse_args()

    if args.list:
        if not args.file:
            raise RuntimeError("--list requires --file")
        f = open_file(args.file)
        for name, cls in list_keys(f):
            print("{:<70s} {}".format(name, cls))
        return

    if args.old_file and args.new_file:
        if not args.hist:
            raise RuntimeError("--old-file/--new-file compare requires --hist")
        compare_two_mnv_files(args)
        return

    if not args.file:
        raise RuntimeError("Need --file")

    if args.cv_hist and args.univ_hist:
        check_separate_hists(args)
        return

    if args.hist:
        check_mnv_band(args)
        return

    raise RuntimeError(
        "Choose one mode:\n"
        "  --list --file FILE\n"
        "  --file FILE --cv-hist CV --univ-hist U0 U1 ...\n"
        "  --file FILE --hist HIST --band Flux\n"
        "  --old-file OLD.root --new-file NEW.root --hist HIST --band Flux"
    )


if __name__ == "__main__":
    main()