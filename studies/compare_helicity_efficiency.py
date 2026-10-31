
#!/usr/bin/env python3

import os
import math
import ROOT
import numpy as np
import matplotlib.pyplot as plt

ROOT.gROOT.SetBatch(True)

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

DIRECTORY = "/exp/minerva/data/users/qvuong/antinu_mu"

TAG_BEFORE = (
    "ValidateCuts_CCnumubar_FV_"
    "InverseHasNoBackExitingTracks_DeadTime_LeptonAngle_v1"
)

TAG_AFTER = "ValidateCuts_CCnumubar_Inclusive_v1"

HIST_NAME = "Eel"
OUTPUT = "studies/helicity_efficiency"

os.makedirs(OUTPUT, exist_ok=True)


def load_hist(sample, tag):
    filename = os.path.join(
        DIRECTORY,
        f"kin_dist_{sample}le5_{tag}_MAD.root"
    )

    f = ROOT.TFile.Open(filename)

    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open {filename}")

    h = f.Get(HIST_NAME)

    if not h:
        raise RuntimeError(f"Missing {HIST_NAME}: {filename}")

    nbins = h.GetNbinsX()

    edges = np.array([
        h.GetXaxis().GetBinLowEdge(i)
        for i in range(1, nbins + 2)
    ])

    values = np.array([
        h.GetBinContent(i)
        for i in range(1, nbins + 1)
    ])

    variances = np.array([
        h.GetBinError(i)**2
        for i in range(1, nbins + 1)
    ])

    f.Close()

    return edges, values, variances


# ------------------------------------------------------------
# Efficiency and uncertainty
# ------------------------------------------------------------

def efficiency(before, after, var_before, var_after):
    """
    after is a subset of before.

    For weighted MC:
        Var(before - after) = Var(before) - Var(after)

    This assumes matching event weights between the two runs.
    """

    eff = np.full(len(before), np.nan)
    err = np.full(len(before), np.nan)

    for i in range(len(before)):
        B = before[i]
        A = after[i]

        if B <= 0:
            continue

        e = A / B

        # Approximate sumw2 variance of rejected events.
        v_pass = var_after[i]
        v_fail = var_before[i] - var_after[i]

        # Small numerical negatives can result from rounding.
        v_fail = max(0.0, v_fail)

        variance = (
            (1.0 - e)**2 * v_pass
            + e**2 * v_fail
        ) / B**2

        eff[i] = e
        err[i] = math.sqrt(max(0.0, variance))

    return eff, err


# ------------------------------------------------------------
# Read histograms
# ------------------------------------------------------------

edges_d0, d0, vd0 = load_hist("data", TAG_BEFORE)
edges_d1, d1, vd1 = load_hist("data", TAG_AFTER)

edges_m0, m0, vm0 = load_hist("mc", TAG_BEFORE)
edges_m1, m1, vm1 = load_hist("mc", TAG_AFTER)

for edges in [edges_d1, edges_m0, edges_m1]:
    if not np.allclose(edges_d0, edges):
        raise RuntimeError("Histogram binning mismatch")

if not all(np.all(np.isfinite(x)) for x in
           [d0, d1, m0, m1, vd0, vd1, vm0, vm1]):
    raise RuntimeError("Nonfinite histogram contents/errors")

centers = (edges_d0[:-1] + edges_d0[1:]) / 2
xerr = (edges_d0[1:] - edges_d0[:-1]) / 2

# Verify subset behavior in data.
if np.any(d1 > d0 + 1e-8):
    raise RuntimeError("After-cut data exceeds before-cut data")

# Weighted MC may have negative individual weights;
# this check is a warning, not a universal requirement.
if np.any(m1 > m0 + 1e-8):
    print("WARNING: weighted MC after > before in some bins")

# ------------------------------------------------------------
# Calculate efficiencies
# ------------------------------------------------------------

eff_d, err_d = efficiency(d0, d1, vd0, vd1)
eff_m, err_m = efficiency(m0, m1, vm0, vm1)

ratio = np.full(len(centers), np.nan)
ratio_err = np.full(len(centers), np.nan)

for i in range(len(centers)):
    if (
        not np.isfinite(eff_d[i])
        or not np.isfinite(eff_m[i])
        or eff_d[i] <= 0
        or eff_m[i] <= 0
    ):
        continue

    ratio[i] = eff_d[i] / eff_m[i]

    ratio_err[i] = ratio[i] * math.sqrt(
        (err_d[i] / eff_d[i])**2
        + (err_m[i] / eff_m[i])**2
    )

# ------------------------------------------------------------
# Integrated results: include underflow/overflow separately
# if needed; these totals use visible bins only.
# ------------------------------------------------------------

total_d0 = d0.sum()
total_d1 = d1.sum()
total_m0 = m0.sum()
total_m1 = m1.sum()

int_eff_d = total_d1 / total_d0
int_eff_m = total_m1 / total_m0
int_ratio = int_eff_d / int_eff_m

print("\n========== INTEGRATED ==========")
print(f"Data before: {total_d0:.3f}")
print(f"Data after:  {total_d1:.3f}")
print(f"MC before:   {total_m0:.3f}")
print(f"MC after:    {total_m1:.3f}")
print(f"Data efficiency: {int_eff_d:.5f}")
print(f"MC efficiency:   {int_eff_m:.5f}")
print(f"Double ratio:    {int_ratio:.5f}")

print("\n========== PER-BIN RESULTS ==========")
print("E_low E_high DataEff MCeff DoubleRatio Error")

for i in range(len(centers)):
    print(
        f"{edges_d0[i]:6.2f} "
        f"{edges_d0[i+1]:6.2f} "
        f"{eff_d[i]:8.4f} "
        f"{eff_m[i]:8.4f} "
        f"{ratio[i]:11.4f} "
        f"{ratio_err[i]:8.4f}"
    )

# ------------------------------------------------------------
# Plot data and MC helicity survival fractions
# ------------------------------------------------------------

plt.figure(figsize=(9, 6))

plt.errorbar(
    centers, eff_d, xerr=xerr, yerr=err_d,
    fmt="o", capsize=2, label="Data"
)

plt.errorbar(
    centers, eff_m, xerr=xerr, yerr=err_m,
    fmt="s", capsize=2, label="MC (weighted)"
)

plt.xlabel("Reconstructed muon energy (GeV)")
plt.ylabel("Helicity-cut survival fraction")
plt.title("RHC: Helicity selection efficiency")
plt.legend()
plt.grid(alpha=0.25)
plt.tight_layout()

plt.savefig(
    os.path.join(OUTPUT, "helicity_efficiency.png"),
    dpi=180
)
plt.close()

# ------------------------------------------------------------
# Plot data/MC double survival ratio
# ------------------------------------------------------------

plt.figure(figsize=(9, 6))

plt.errorbar(
    centers, ratio, xerr=xerr, yerr=ratio_err,
    fmt="o", capsize=2
)

plt.axhline(
    1.0, linestyle="--", color="black",
    label="Equal data/MC survival"
)

plt.axhline(
    int_ratio, linestyle=":", color="gray",
    label=f"Integrated ratio = {int_ratio:.3f}"
)

plt.xlabel("Reconstructed muon energy (GeV)")
plt.ylabel("Data / MC helicity survival ratio")
plt.title("RHC: Relative helicity-cut acceptance")
plt.legend()
plt.grid(alpha=0.25)
plt.tight_layout()

plt.savefig(
    os.path.join(OUTPUT, "helicity_double_ratio.png"),
    dpi=180
)
plt.close()

print("\nSaved plots to:", OUTPUT)
