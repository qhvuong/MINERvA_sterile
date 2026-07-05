#!/usr/bin/env python3

import argparse
import math
import ROOT

ROOT.TH1.AddDirectory(False)


def finite(x):
    try:
        return math.isfinite(float(x))
    except Exception:
        return False


def get_tree(f, tree_name=None):
    if tree_name:
        t = f.Get(tree_name)
        if not t:
            raise RuntimeError("Could not find tree '{}'".format(tree_name))
        return t

    # Auto-find first TTree
    for key in f.GetListOfKeys():
        obj = key.ReadObj()
        if obj.InheritsFrom("TTree"):
            print("[AUTO_TREE] using tree:", key.GetName())
            return obj

    raise RuntimeError("No TTree found. Use --list to inspect file.")


def list_file(path):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError("Could not open {}".format(path))

    print("[KEYS]")
    for key in f.GetListOfKeys():
        obj = key.ReadObj()
        print("  {:50s} {}".format(key.GetName(), obj.ClassName()))

    print("\n[BRANCHES]")
    for key in f.GetListOfKeys():
        obj = key.ReadObj()
        if obj.InheritsFrom("TTree"):
            print("\nTREE:", key.GetName())
            for br in obj.GetListOfBranches():
                print("  ", br.GetName())


def get_branch_value(tree, branch_name):
    """
    Returns branch value from current tree entry.
    Works for scalar branches and vector/array-like branches.
    """
    try:
        return getattr(tree, branch_name)
    except Exception:
        return None


def as_float(x):
    try:
        return float(x)
    except Exception:
        return None


def get_universe_value(uobj, universe_index):
    """
    Extract one universe value from a scalar/vector/array branch object.
    """
    # Scalar branch
    val = as_float(uobj)
    if val is not None:
        return val

    # std::vector / C array / Python proxy
    try:
        return float(uobj[universe_index])
    except Exception:
        return None


def classify_values(cv_vals, univ_vals, ratio_vals):
    abs_cv = [abs(x) for x in cv_vals if finite(x) and x != 0]
    abs_univ = [abs(x) for x in univ_vals if finite(x) and x != 0]
    abs_ratio = [abs(x) for x in ratio_vals if finite(x) and x != 0]

    def median(v):
        if not v:
            return None
        v = sorted(v)
        return v[len(v)//2]

    cv_med = median(abs_cv)
    univ_med = median(abs_univ)
    ratio_med = median(abs_ratio)

    print("\n[SUMMARY]")
    print("  median |CV|       =", cv_med)
    print("  median |universe| =", univ_med)
    print("  median |univ/CV|  =", ratio_med)

    print("\n[INTERPRETATION]")
    if cv_med is None or univ_med is None:
        print("  Not enough valid values.")
        return

    if 0.2 <= univ_med <= 5.0 and (cv_med < 0.05 or cv_med > 20.0):
        print("  Universe branch looks RELATIVE-ish.")
        print("  If you divide universe by CV, that is probably wrong.")
    elif cv_med > 0 and 0.05 <= univ_med / cv_med <= 20.0:
        print("  Universe branch looks ABSOLUTE-ish, same scale as CV.")
        print("  Making relative as universe/CV is probably correct, assuming CV is the matching denominator.")
    else:
        print("  Ambiguous. Inspect event-by-event values below.")


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("files", nargs="+", help="Input ntuple ROOT files")
    ap.add_argument("--tree", default=None, help="Tree name. If omitted, auto-uses first TTree.")
    ap.add_argument("--list", action="store_true", help="List keys/branches and exit")

    ap.add_argument("--energy-branch", default="mc_incomingE")
    ap.add_argument("--pdg-branch", default="mc_incoming")

    ap.add_argument("--cv-branch", required=False, help="CV flux weight branch")
    ap.add_argument("--univ-branch", required=False, help="Universe flux weight branch or vector branch")

    ap.add_argument("--univ-index", type=int, default=0, help="Universe index if --univ-branch is vector-like")
    ap.add_argument("--univ-scale", type=float, default=1.0, help="Scale factor applied to universe branch value")
    ap.add_argument("--max-events", type=int, default=50000)
    ap.add_argument("--print-max", type=int, default=200)

    ap.add_argument("--emin", type=float, default=4.5)
    ap.add_argument("--emax", type=float, default=7.5)
    ap.add_argument("--pdg", type=int, default=None, help="Only inspect this incoming neutrino PDG, e.g. 14 or 12")

    args = ap.parse_args()

    if args.list:
        for path in args.files:
            print("\n==========", path, "==========")
            list_file(path)
        return

    if not args.cv_branch or not args.univ_branch:
        raise RuntimeError("Need --cv-branch and --univ-branch, or use --list first.")

    chain = ROOT.TChain(args.tree if args.tree else "")

    if args.tree:
        for path in args.files:
            chain.Add(path)
        tree = chain
    else:
        # Auto tree mode only supports one file cleanly
        if len(args.files) != 1:
            raise RuntimeError("Auto tree mode supports one file. Use --tree for multiple files.")
        f = ROOT.TFile.Open(args.files[0])
        if not f or f.IsZombie():
            raise RuntimeError("Could not open {}".format(args.files[0]))
        tree = get_tree(f, None)

    n = min(int(tree.GetEntries()), args.max_events)

    print("[CONFIG]")
    print("  entries      =", tree.GetEntries())
    print("  scanning     =", n)
    print("  E range GeV  =", args.emin, args.emax)
    print("  CV branch    =", args.cv_branch)
    print("  univ branch  =", args.univ_branch)
    print("  univ index   =", args.univ_index)
    print("  PDG filter   =", args.pdg)

    cv_vals = []
    univ_vals = []
    ratio_vals = []

    n_pass = 0
    n_print = 0

    for i in range(n):
        tree.GetEntry(i)

        e_raw = get_branch_value(tree, args.energy_branch)
        e = as_float(e_raw)
        if e is None:
            continue

        # Most MINERvA ntuples store mc_incomingE in MeV
        e_gev = e * 1e-3 if e > 100.0 else e

        if not (args.emin <= e_gev < args.emax):
            continue

        if args.pdg is not None:
            pdg_val = get_branch_value(tree, args.pdg_branch)
            pdg_val = as_float(pdg_val)
            if pdg_val is None or int(pdg_val) != args.pdg:
                continue

        cv_obj = get_branch_value(tree, args.cv_branch)
        univ_obj = get_branch_value(tree, args.univ_branch)

        cv = as_float(cv_obj)
        univ = get_universe_value(univ_obj, args.univ_index)
        if univ is not None:
            univ *= args.univ_scale

        if cv is None or univ is None:
            continue

        ratio = univ / cv if cv != 0 else float("nan")

        cv_vals.append(cv)
        univ_vals.append(univ)
        ratio_vals.append(ratio)

        n_pass += 1

        if n_print < args.print_max:
            print(
                "[EVT] entry={:8d} E={:8.4f} cv={:16.8e} "
                "univ={:16.8e} univ/cv={:16.8e}".format(
                    i, e_gev, cv, univ, ratio
                )
            )
            n_print += 1

    print("\n[COUNTS]")
    print("  passed events =", n_pass)

    classify_values(cv_vals, univ_vals, ratio_vals)


if __name__ == "__main__":
    main()