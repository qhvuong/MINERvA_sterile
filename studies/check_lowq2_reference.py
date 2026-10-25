#!/usr/bin/env python3

import os
import math
import ROOT
import PlotUtils

ROOT.gROOT.SetBatch(True)

root_dir = os.environ.get("MPARAMFILESROOT")
if not root_dir:
    raise RuntimeError("MPARAMFILESROOT is not set")

path = os.path.join(
    root_dir, "data", "Reweight", "lowQ2pi_weights.root"
)

f = ROOT.TFile.Open(path)
if not f or f.IsZombie():
    raise RuntimeError(f"Cannot open {path}")

print("Reference ROOT file:", path)
weighter = PlotUtils.weightLowQ2Pi()

print("Loaded low-Q2 pion weighter:", type(weighter))

q2_values = [
    -0.01,
    0.00,
    0.01,
    0.05,
    0.10,
    0.15,
    0.20,
    0.30,
    0.40,
    0.50,
    0.60,
    0.69,
    0.699999,
    0.70,
    0.700001,
    0.80,
    1.00,
    1.50,
    2.00,
]

failures = 0

for channel in ("JOINT", "MENU1PI"):
    print("\n" + "=" * 85)
    print(channel)
    print("=" * 85)

    graph = f.Get(f"{channel}_weights")

    if not graph:
        raise RuntimeError(
            f"Missing graph {channel}_weights"
        )

    print("Graph points:", graph.GetN())

    print(
        f"{'Q2':>10s}"
        f"{'Graph Eval':>16s}"
        f"{'C++ Weight':>16s}"
        f"{'Expected':>16s}"
        f"{'Difference':>16s}"
        f"{'Result':>10s}"
    )

    for q2 in q2_values:
        graph_value = float(graph.Eval(q2))

        actual = float(
            weighter.getWeight(q2, channel, 0, 1000060120)
        )

        if channel == "JOINT" and (q2 < 0 or q2 >= 0.7):
            expected = 1.0
        else:
            expected = graph_value

        difference = actual - expected

        passed = (
            math.isfinite(actual)
            and math.isclose(
                actual, expected,
                rel_tol=1e-10,
                abs_tol=1e-10,
            )
        )

        if not passed:
            failures += 1

        print(
            f"{q2:10.6f}"
            f"{graph_value:16.8f}"
            f"{actual:16.8f}"
            f"{expected:16.8f}"
            f"{difference:16.3e}"
            f"{'PASS' if passed else 'FAIL':>10s}"
        )

    print("\nHydrogen target check:")

    for q2 in (0.05, 0.2, 0.5, 0.8, 1.0):
        actual = float(
            weighter.getWeight(q2, channel, 0, 2212)
        )

        if channel == "MENU1PI":
            expected = 1.0
        elif q2 < 0 or q2 >= 0.7:
            expected = 1.0
        else:
            expected = float(graph.Eval(q2))

        passed = math.isclose(
            actual, expected,
            rel_tol=1e-10,
            abs_tol=1e-10,
        )

        if not passed:
            failures += 1

        print(
            f"  Q2={q2:.3f}: "
            f"actual={actual:.8f}, "
            f"expected={expected:.8f}, "
            f"{'PASS' if passed else 'FAIL'}"
        )

print("\n" + "=" * 85)
print("FINAL RESULT")
print("=" * 85)
print("Failures:", failures)

if failures:
    raise SystemExit(1)

print("All reference-graph comparisons passed.")

f.Close()