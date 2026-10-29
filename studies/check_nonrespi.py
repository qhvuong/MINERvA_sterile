
#!/usr/bin/env python3
"""Audit the exact v1 NonResPi condition on a MINERvA MC ntuple."""

import argparse
from collections import defaultdict, Counter
import ROOT

ROOT.gROOT.SetBatch(True)

LABELS = {
    1: "QE",
    2: "RES",
    3: "DIS",
    4: "COH",
    8: "2p2h",
}


def find_tree(root_file, requested=None):
    if requested:
        tree = root_file.Get(requested)
        if not tree or not tree.InheritsFrom("TTree"):
            raise RuntimeError(f"Cannot find tree: {requested}")
        return tree

    # Find the tree containing the required GENIE branches.
    for key in root_file.GetListOfKeys():
        obj = key.ReadObj()
        if not obj.InheritsFrom("TTree"):
            continue
        names = {b.GetName() for b in obj.GetListOfBranches()}
        required = {
            "mc_intType",
            "truth_genie_wgt_Rvn1pi",
            "truth_genie_wgt_Rvp1pi",
        }
        if required.issubset(names):
            return obj

    raise RuntimeError("Cannot find a tree with the required branches")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", required=True)
    parser.add_argument("--tree", default=None)
    parser.add_argument("--max-events", type=int, default=100000)
    parser.add_argument("--examples", type=int, default=5)
    args = parser.parse_args()

    f = ROOT.TFile.Open(args.file)
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open {args.file}")

    try:
        tree = find_tree(f, args.tree)
        names = {b.GetName() for b in tree.GetListOfBranches()}

        has_pdg = "mc_incoming" in names
        has_current = "mc_current" in names
        has_w = "mc_w" in names

        nentries = tree.GetEntries()
        n = min(nentries, args.max_events)

        print(f"\nTree: {tree.GetName()}")
        print(f"Total entries: {nentries}")
        print(f"Entries to scan: {n}")

        stats = defaultdict(lambda: {
            "n": 0,
            "triggered": 0,
            "sum_factor": 0.0,
            "rvn_below": 0,
            "rvp_below": 0,
            "both_below": 0,
            "short_vectors": 0,
            "examples": [],
            "rvn_values": Counter(),
            "rvp_values": Counter(),
        })

        for i in range(n):
            tree.GetEntry(i)

            itype = int(tree.mc_intType)

            # Focus on neutrino charged-current interactions
            # when those identifying branches are available.
            if has_current and int(tree.mc_current) != 1:
                continue
            if has_pdg and int(tree.mc_incoming) != 14:
                continue

            s = stats[itype]
            s["n"] += 1

            rvn = getattr(tree, "truth_genie_wgt_Rvn1pi")
            rvp = getattr(tree, "truth_genie_wgt_Rvp1pi")

            if len(rvn) <= 2 or len(rvp) <= 2:
                s["short_vectors"] += 1
                continue

            a = float(rvn[2])
            b = float(rvp[2])

            # EXACT condition from GetStandardWeightV1.
            hit_a = a < 1.0
            hit_b = b < 1.0
            triggered = hit_a or hit_b

            factor = 0.43 if triggered else 1.0

            s["triggered"] += int(triggered)
            s["rvn_below"] += int(hit_a)
            s["rvp_below"] += int(hit_b)
            s["both_below"] += int(hit_a and hit_b)
            s["sum_factor"] += factor

            s["rvn_values"][round(a, 4)] += 1
            s["rvp_values"][round(b, 4)] += 1

            if len(s["examples"]) < args.examples:
                w = float(tree.mc_w) if has_w else None
                s["examples"].append(
                    (i, a, b, triggered, factor, w)
                )

        print("\n========== NONRESPI APPLICABILITY ==========")
        print(
            f"{'Type':<10} {'Events':>10} "
            f"{'Flagged':>10} {'Fraction':>11} "
            f"{'Mean factor':>13}"
        )

        for itype in sorted(stats):
            s = stats[itype]
            count = s["n"]
            if count == 0:
                continue

            print(
                f"{LABELS.get(itype, str(itype)):<10} "
                f"{count:>10d} "
                f"{s['triggered']:>10d} "
                f"{s['triggered']/count:>11.4f} "
                f"{s['sum_factor']/count:>13.4f}"
            )

        for itype in sorted(stats):
            s = stats[itype]
            if not s["n"]:
                continue

            print(f"\n========== {LABELS.get(itype, itype)} ==========")
            print("Rvn1pi[2] < 1:", s["rvn_below"])
            print("Rvp1pi[2] < 1:", s["rvp_below"])
            print("Both < 1:    ", s["both_below"])
            print("Short vectors:", s["short_vectors"])

            print("Most common Rvn1pi[2] values:")
            print(s["rvn_values"].most_common(10))

            print("Most common Rvp1pi[2] values:")
            print(s["rvp_values"].most_common(10))

            print("Example events:")
            for example in s["examples"]:
                idx, a, b, triggered, factor, w = example
                print(
                    f"  entry={idx} "
                    f"rvn={a:.6g} rvp={b:.6g} "
                    f"flag={triggered} "
                    f"factor={factor} mc_w={w}"
                )

    finally:
        f.Close()


if __name__ == "__main__":
    main()
