#!/usr/bin/env python3
"""Read-only closure for the CC electron inputs written to stitch."""
import argparse
from contextlib import ExitStack
from pathlib import Path

import numpy as np


SIGNAL = ("CCNuE", "CCNuEQE", "CCNuEDelta", "CCNuEDIS", "CCNuE2p2h")
BACKGROUND = (
    "NCDiff", "CCPi0", "NCCohPi0", "NCPi0", "CCPi", "NCPi",
    "CCOther", "NCOther", "Other",
)
FIXED_OR_WRONG_SIGN = ("CCNuEWrongSign", "NuEElastic")
PRODUCTION_SIGNAL = SIGNAL + ("CCNuEWrongSign",)
INPUTS = {
    "fhc": ("nu_e", "leFHC", "CCnue", 0.15031290821716903),
    "rhc": ("antinu_e", "le5", "CCnuebar", 0.1095830686235881),
}


def read_hist(file, key, edges=None):
    hist = file.Get(key)
    if not hist or hist.GetDimension() != 1:
        raise RuntimeError(f"{file.GetName()}: missing 1D histogram {key}")
    axis = hist.GetXaxis()
    n = hist.GetNbinsX()
    bins = np.array([axis.GetBinLowEdge(i) for i in range(1, n + 2)])
    if edges is not None and not np.array_equal(bins, edges):
        raise RuntimeError(f"{file.GetName()}: incompatible bin edges in {key}")
    return bins, np.array([hist.GetBinContent(i) for i in range(1, n + 1)])


def read_region(data, mc, key, pot_scale, edges):
    _, observed = read_hist(data, key, edges)
    _, total = read_hist(mc, key, edges)
    # Match HistHolder's truth-category lookup. A shared prefix also matches
    # diagnostic subsets such as EN4_CCNuEQE_PDG12.
    category_names = SIGNAL + BACKGROUND + FIXED_OR_WRONG_SIGN
    categories = {
        name: read_hist(mc, f"{key}_{name}", edges)[1] * pot_scale
        for name in category_names
    }
    signal = sum(categories[name] for name in SIGNAL)
    wrong_sign = categories["CCNuEWrongSign"]
    elastic = categories["NuEElastic"]
    background = sum(categories[name] for name in BACKGROUND)
    total *= pot_scale
    category_sum = signal + wrong_sign + elastic + background
    residual = total - category_sum
    ignored = sorted({k.GetName() for k in mc.GetListOfKeys()}
                     - {f"{key}_{name}" for name in category_names})
    ignored = [name for name in ignored if name.startswith(key + "_")
               and not (key == "EN4" and name.startswith("EN4_dEdX"))]
    print(f"{key}: max |raw MC total - category sum| = {np.max(np.abs(residual)):.8g}")
    if ignored:
        print(f"  Extra keys excluded from category sum: {', '.join(ignored)}")
    if not np.allclose(total, category_sum, rtol=1e-8, atol=1e-8):
        print("  WARNING: category sum does not close. See MC_gap by bin.")
    return observed, total, signal, background, wrong_sign, elastic, residual


def solve(sr, sb, wrong_sign_treatment):
    d_sr, _, s_sr, b_sr, w_sr, e_sr, _ = sr
    d_sb, _, s_sb, b_sb, w_sb, e_sb, _ = sb
    if wrong_sign_treatment == "floating":
        b_sr, b_sb = b_sr + w_sr, b_sb + w_sb
        w_sr, w_sb = 0.0, 0.0
    elif wrong_sign_treatment == "signal":
        s_sr, s_sb = s_sr + w_sr, s_sb + w_sb
        w_sr, w_sb = 0.0, 0.0
    x_sr, x_sb = d_sr - w_sr - e_sr, d_sb - w_sb - e_sb
    det = s_sr * b_sb - s_sb * b_sr
    with np.errstate(divide="ignore", invalid="ignore"):
        alpha = (x_sr * b_sb - x_sb * b_sr) / det
        beta = (s_sr * x_sb - s_sb * x_sr) / det
        separation = det / ((s_sr + b_sr) * (s_sb + b_sb))
    fitted_bkg = beta * b_sr + w_sr + e_sr
    return alpha, beta, d_sr - fitted_bkg, separation, s_sr


def data_only_sigma(sr, sb):
    d_sr, _, s_sr, b_sr, w_sr, _, _ = sr
    d_sb, _, s_sb, b_sb, w_sb, _, _ = sb
    s_sr, s_sb = s_sr + w_sr, s_sb + w_sb
    det = s_sr * b_sb - s_sb * b_sr
    return np.sqrt((s_sr * b_sb / det)**2 * d_sr + (s_sr * b_sr / det)**2 * d_sb)


def stat_errors(hist):
    return np.array([hist.GetBinError(i) for i in range(1, hist.GetNbinsX() + 1)])


def template_stat_errors(mc, key, pot_scale):
    groups = (SIGNAL, BACKGROUND, ("CCNuEWrongSign",), ("NuEElastic",))
    return tuple(pot_scale * np.sqrt(sum(stat_errors(mc.Get(f"{key}_{cat}"))**2 for cat in cats))
                 for cats in groups)


def fluctuation_check(sr, sb, sr_errors, sb_errors, n_toys, seed):
    rng = np.random.default_rng(seed)
    toy_sr, toy_sb = list(sr), list(sb)
    for toy in (toy_sr, toy_sb):
        toy[0] = rng.poisson(toy[0], size=(n_toys, len(toy[0])))
    data_only = solve(toy_sr, toy_sb, "signal")[2]
    # One draw per region, bin and component. The same draw enters the
    # scale solve and the subtraction. Histograms provide marginal sumw2.
    for toy, errors in ((toy_sr, sr_errors), (toy_sb, sb_errors)):
        for index, error in enumerate(errors, start=2):
            toy[index] = rng.normal(toy[index], error, size=toy[0].shape)
    _, _, data_mc, sep, _ = solve(toy_sr, toy_sb, "signal")
    nominal_sep = solve(sr, sb, "signal")[3]
    valid = np.isfinite(data_mc) & np.isfinite(sep)
    crossed = np.mean(valid & (np.sign(sep) != np.sign(nominal_sep)), axis=0)
    data_mc = np.where(valid, data_mc, np.nan)
    q_low, q_high = np.nanquantile(data_mc, [0.158655, 0.841345], axis=0)
    return (data_only.std(axis=0, ddof=1), np.nanstd(data_mc, axis=0, ddof=1),
            (q_high - q_low) / 2, crossed, np.sum(~valid, axis=0))


def flux_prediction_check(mc, prediction, pot_scale):
    sources = [mc.Get(f"EN4_{cat}") for cat in PRODUCTION_SIGNAL]
    if any("Flux" not in h.GetVertErrorBandNames() for h in [prediction] + sources):
        print("Flux prediction closure: FAIL, missing Flux band")
        return False
    target = prediction.GetVertErrorBand("Flux")
    bands = [h.GetVertErrorBand("Flux") for h in sources]
    count = target.GetNHists()
    if count == 0 or any(b.GetNHists() != count for b in bands):
        print("Flux prediction closure: FAIL, inconsistent universe counts")
        return False

    def contents(hist):
        return np.array([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)])

    def compare(hist, components):
        observed = contents(hist)
        expected = pot_scale * sum(map(contents, components))
        return np.isclose(observed, expected, rtol=1e-8, atol=1e-8), np.abs(observed - expected)

    cv_match, _ = compare(target, bands)
    failures = 0
    maximum = 0.0
    for i in range(count):
        match, residual = compare(target.GetHist(i), [b.GetHist(i) for b in bands])
        failures += np.count_nonzero(~match)
        maximum = max(maximum, float(np.max(residual))) if np.all(np.isfinite(residual)) else float("inf")
    ok = cv_match.all() and failures == 0
    print(f"Flux prediction closure: {'PASS' if ok else 'FAIL'}, {count} universes, "
          f"{failures} mismatched universe bins, max absolute residual {maximum:.8g}, "
          f"band CV {'matches' if cv_match.all() else 'differs'}")
    return ok


def audit(beam, base, root, uncertainty_toys=0, seed=431):
    folder, playlist, selection, pot_scale = INPUTS[beam]
    stem = f"{selection}_FullCuts_v4p3p1_MAD.root"
    paths = (
        base / folder / f"kin_dist_data{playlist}_{stem}",
        base / folder / f"kin_dist_mc{playlist}_{stem}",
        base / folder / f"bkgfit_{playlist}_N4_tune_{stem}",
    )
    with ExitStack() as stack:
        files = []
        print(f"\n{beam.upper()}  MC-to-data POT scale: {pot_scale:.15g}")
        for path in paths:
            print(path)
            file = root.TFile.Open(str(path), "READ")
            if not file or file.IsZombie():
                raise RuntimeError(f"Cannot open {path}")
            stack.callback(file.Close)
            files.append(file)
        data, mc, fit = files
        edges, _ = read_hist(data, "EN4")
        sr = read_region(data, mc, "EN4", pot_scale, edges)
        sb = read_region(data, mc, "EN4_dEdX", pot_scale, edges)
        saved = tuple(read_hist(fit, key, edges)[1] for key in (
            "signal_Scale_Factor", "background_Scale_Factor",
            "EN4_data_bkgSubbed", "EN4_predicted_Signal",
        ))
        uncertainty_ok = True
        if uncertainty_toys:
            uncertainty_ok = flux_prediction_check(mc, fit.Get("EN4_predicted_Signal"), pot_scale)
            production = solve(sr, sb, "signal")
            if not all(np.allclose(a, b, rtol=1e-8, atol=1e-8)
                       for a, b in zip((*production[:3], production[4]), saved)):
                print("Statistical closure: skipped, production central values differ")
                uncertainty_ok = False
            else:
                saved_stat = stat_errors(fit.Get("EN4_data_bkgSubbed").GetCVHistoWithStatError())
                exact = data_only_sigma(sr, sb)
                sr_errors = template_stat_errors(mc, "EN4", pot_scale)
                sb_errors = template_stat_errors(mc, "EN4_dEdX", pot_scale)
                toy = fluctuation_check(sr, sb, sr_errors, sb_errors, uncertainty_toys, seed)
                print(f"\nStatistical closure ({uncertainty_toys} toys, seed {seed})")
                print("bin   saved_stat  data_exact   data_toy  data+MC_SD  data+MC_68half  sep_flip_%  nonfinite")
                for i, row in enumerate(zip(saved_stat, exact, *toy), start=1):
                    saved_err, expected, data_sd, mc_sd, width, crossed, bad = row
                    print(f"{i:3d} {saved_err:12.5g} {expected:11.5g} {data_sd:10.5g}"
                          f" {mc_sd:11.5g} {width:15.5g} {100 * crossed:11.5g} {bad:10d}")

    matches = []
    for label, treatment in (
        ("fixed wrong-sign", "fixed"),
        ("floating wrong-sign", "floating"),
        ("wrong-sign included in signal", "signal"),
    ):
        alpha, beta, subtracted, separation, nominal_signal = solve(sr, sb, treatment)
        calculated = (alpha, beta, subtracted, nominal_signal)
        matching = all(np.allclose(a, b, rtol=1e-8, atol=1e-8) for a, b in zip(calculated, saved))
        print(f"\n{label}: {'MATCH' if matching else 'MISMATCH'}")
        for key, a, b in zip(("alpha", "beta", "subtracted", "nominal signal"), calculated, saved):
            print(f"  {key:15s} max |reconstructed - saved| = {np.max(np.abs(a - b)):.8g}")
        if matching:
            matches.append(label)
        print("bin       E range      D_SR   MC_SR   S_SR   B_SR   WS_SR   el_SR   MC_gap_SR"
              "    D_SB    S_SB    B_SB   WS_SB   el_SB   MC_gap_SB   alpha    beta"
              "     sep  raw_D-MC  bkg_shift   sub_D-S   fit_S_SR")
        for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
            raw_delta = sr[0][i] - sr[1][i]
            sub_delta = subtracted[i] - nominal_signal[i]
            values = [v[i] for v in sr] + [v[i] for j, v in enumerate(sb) if j != 1]
            bkg_shift = sub_delta - raw_delta - sr[-1][i]
            values += [alpha[i], beta[i], separation[i], raw_delta, bkg_shift, sub_delta,
                       nominal_signal[i]]
            print(f"{i + 1:3d} [{lo:5.2f},{hi:5.2f}) " + " ".join(f"{v:8.3f}" for v in values))
    return bool(matches) and uncertainty_ok


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--beam", choices=("fhc", "rhc", "both"), default="both")
    parser.add_argument("--base", type=Path, default=Path("/exp/minerva/data/users/qvuong"))
    parser.add_argument("--uncertainty-toys", type=int, default=0,
                        help="Also audit Flux and run this many statistical toys per beam")
    parser.add_argument("--seed", type=int, default=431)
    args = parser.parse_args()
    if args.uncertainty_toys < 0 or args.uncertainty_toys == 1:
        parser.error("--uncertainty-toys must be zero or at least two")
    import ROOT
    import PlotUtils  # Load the MnvH1D dictionaries before opening files.

    ROOT.gROOT.SetBatch(True)
    beams = INPUTS if args.beam == "both" else (args.beam,)
    results = [audit(beam, args.base, ROOT, args.uncertainty_toys, args.seed + i)
               for i, beam in enumerate(beams)]
    return 0 if all(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
