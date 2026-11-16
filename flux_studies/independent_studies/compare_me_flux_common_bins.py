#!/usr/bin/env python3
"""Read-only ME flux comparison on exact shared bin boundaries.

Inputs are flux densities per GeV by default. Only whole native bins are
combined; no interpolation, fractional bin splitting, or playlist averaging.
Universe ensembles are compared without pairing their universe indices.
"""
import argparse
import csv
import json
from pathlib import Path

np = None


def common_edges(a, b, emin, emax):
    edges = [float(x) for x in a if emin - 1e-8 <= x <= emax + 1e-8
             and np.any(np.isclose(b, x, rtol=0, atol=1e-8))]
    if len(edges) < 2:
        raise ValueError("Fewer than two common boundaries in requested range")
    return np.array(edges)


def rebin_density(values, old_edges, new_edges, convention="density"):
    indices = []
    for edge in new_edges:
        hits = np.flatnonzero(np.isclose(old_edges, edge, rtol=0, atol=1e-8))
        if len(hits) != 1:
            raise ValueError("Ambiguous/nonmatching common boundary: {}".format(edge))
        indices.append(int(hits[0]))
    integrated = values * np.diff(old_edges) if convention == "density" else values
    cumulative = np.concatenate([
        np.zeros(values.shape[:-1] + (1,)), np.cumsum(integrated, axis=-1)
    ], axis=-1)
    return np.diff(cumulative[..., indices], axis=-1) / np.diff(new_edges)


def read_hist(ROOT, path, hname):
    f = ROOT.TFile.Open(str(path), "READ")
    if not f or f.IsZombie():
        raise RuntimeError("Cannot open {}".format(path))
    try:
        h = f.Get(hname)
        if not h or not hasattr(h, "GetVertErrorBandNames"):
            raise RuntimeError("Missing MnvH1D {} in {}".format(hname, path))
        n = h.GetNbinsX()
        arr = lambda obj: np.array([obj.GetBinContent(i) for i in range(1, n+1)])
        bands = {}
        band_cvs = {}
        for name in h.GetVertErrorBandNames():
            name = str(name)
            band = h.GetVertErrorBand(name)
            bands[name] = np.array([arr(band.GetHist(u)) for u in range(band.GetNHists())])
            band_cvs[name] = arr(band)
        result = dict(
            edges=np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1, n+2)]),
            cv=arr(h), bands=bands, band_cvs=band_cvs,
            x_title=str(h.GetXaxis().GetTitle()), y_title=str(h.GetYaxis().GetTitle()),
            recovered=bool(f.TestBit(ROOT.TFile.kRecovered)),
            underflow=float(h.GetBinContent(0)), overflow=float(h.GetBinContent(n+1)))
        if not np.all(np.isfinite(result["cv"])):
            raise ValueError("Nonfinite CV: {}".format(path))
        return result
    finally:
        f.Close()


def weight_stats(universes, cv, widths):
    if not np.all(np.isfinite(universes)) or len(universes) == 0:
        raise ValueError("Empty or nonfinite universe ensemble")
    weights = np.full_like(universes, np.nan, dtype=float)
    np.divide(universes, cv, out=weights, where=cv[None, :] > 0)
    good = cv > 0
    mean = np.full(cv.shape, np.nan)
    std = mean.copy()
    rms = mean.copy()
    quantiles = np.full((3, len(cv)), np.nan)
    mean[good] = weights[:, good].mean(axis=0)
    std[good] = weights[:, good].std(axis=0, ddof=0)
    rms[good] = np.sqrt(np.mean((weights[:, good]-1)**2, axis=0))
    quantiles[:, good] = np.quantile(weights[:, good], [.16, .5, .84], axis=0)
    total = np.sum(cv * widths)
    if total <= 0:
        raise ValueError("Nonpositive integrated CV flux")
    integrated_weights = np.sum(universes * widths, axis=1) / total
    return dict(weights=weights, mean=mean, std=std, rms=rms,
                q16=quantiles[0], median=quantiles[1], q84=quantiles[2],
                integrated_weights=integrated_weights)


def write_csv(path, rows):
    if rows:
        with path.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def plot_case(path, label, edges, p, r, stats, plot_max):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    colors = {"produced": "#007f8b", "reference": "#cc4c33"}
    for side, cv in (("produced", p), ("reference", r)):
        axs[0, 0].stairs(cv, edges, label=side, color=colors[side])
    axs[0, 0].set(yscale="log", ylabel="Flux / GeV (native area and POT units)")
    ratio = np.divide(p, r, out=np.full_like(p, np.nan), where=r > 0)
    axs[0, 1].stairs(ratio, edges, color=colors["produced"])
    axs[0, 1].axhline(1, color="gray", ls=":")
    axs[0, 1].set_ylabel("Produced CV / reference CV")
    styles = ["-", "--", ":", "-."]
    for i, (band, sides) in enumerate(stats.items()):
        for side, values in sides.items():
            axs[1, 0].stairs(100*values["rms"], edges, color=colors[side],
                             linestyle=styles[i % 4], label=side+": "+band)
    axs[1, 0].set_ylabel("Universe/CV RMS about 1 (%)")
    band = "Flux" if "Flux" in stats else next(iter(stats), None)
    if band:
        arrays = [v["integrated_weights"] for v in stats[band].values()]
        bins = np.histogram_bin_edges(np.concatenate(arrays), bins=45)
        for side, values in stats[band].items():
            axs[1, 1].hist(values["integrated_weights"], bins=bins, density=True,
                            histtype="step", label=side, color=colors[side])
        axs[1, 1].set(xlabel="Integrated flux universe / CV", ylabel="Ensemble density",
                      title=band+": {:.4g}--{:.4g} GeV".format(edges[0], edges[-1]))
    for ax in (axs[0, 0], axs[0, 1], axs[1, 0]):
        ax.set_xlim(edges[0], min(plot_max, edges[-1]))
        ax.set_xlabel("Neutrino energy (GeV)")
    for ax in axs.flat:
        ax.grid(alpha=.15)
    axs[0, 0].legend()
    axs[1, 0].legend(fontsize=7)
    if band:
        axs[1, 1].legend()
    fig.suptitle(label+" | whole-bin aggregation, independent ensemble comparison")
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--produced-dir", default="/exp/minerva/app/users/qvuong/MAT_AL9/CC-NuE-XSec/custom_plotutils/data/flux")
    parser.add_argument("--reference-dir", default="/exp/minerva/app/users/qvuong/MAT_AL9/opt/lib/data/flux")
    parser.add_argument("--playlists", nargs="+", default=["1D", "1M", "1N", "5A", "6A"])
    parser.add_argument("--pdgs", nargs="+", type=int, default=[12, -12, 14, -14])
    parser.add_argument("--prefixes", nargs="+", choices=["gen2thin", "g4numiv6"], default=["gen2thin", "g4numiv6"])
    parser.add_argument("--produced-content", choices=["density", "integral"], default="density")
    parser.add_argument("--reference-content", choices=["density", "integral"], default="density")
    parser.add_argument("--emin", type=float, default=0.)
    parser.add_argument("--emax", type=float, default=100.)
    parser.add_argument("--plot-max", type=float, default=20.)
    parser.add_argument("--plots", action="store_true")
    parser.add_argument("--save-arrays", action="store_true")
    parser.add_argument("--out", type=Path, default=Path("me_flux_common_bin_comparison"))
    args = parser.parse_args()
    # ROOT must initialize before NumPy/Matplotlib in the GPVM environment.
    import ROOT
    import PlotUtils.LoadPlotUtilsLib
    global np
    import numpy as np
    ROOT.gROOT.SetBatch(True)
    args.out.mkdir(parents=True, exist_ok=True)
    summaries, bands_summary, bins_summary, metadata, failures = [], [], [], [], []
    for playlist in args.playlists:
        for pdg in args.pdgs:
            for prefix in args.prefixes:
                label = "{}-pdg{}-minervame{}".format(prefix, pdg, playlist)
                name = "flux-"+label+".root"
                hname = "flux_E_cvweighted" if prefix == "gen2thin" else "flux_E_unweighted"
                try:
                    a = read_hist(ROOT, Path(args.produced_dir)/name, hname)
                    b = read_hist(ROOT, Path(args.reference_dir)/name, hname)
                    edges = common_edges(a["edges"], b["edges"], args.emin, args.emax)
                    widths = np.diff(edges)
                    p = rebin_density(a["cv"], a["edges"], edges, args.produced_content)
                    r = rebin_density(b["cv"], b["edges"], edges, args.reference_content)
                    ip, ir = float(p@widths), float(r@widths)
                    if ip <= 0 or ir <= 0:
                        raise ValueError("Nonpositive CV integral")
                    core = r > 1e-4*np.max(r)
                    summary = dict(case=label, produced_bins=len(a["cv"]), reference_bins=len(b["cv"]),
                        common_bins=len(p), emin=float(edges[0]), emax=float(edges[-1]),
                        produced_integral=ip, reference_integral=ir, integral_ratio=ip/ir,
                        fractional_L1_difference=float(np.sum(np.abs(p-r)*widths)/ir),
                        shape_L1_difference=float(np.sum(np.abs(p/ip-r/ir)*widths)),
                        max_abs_CV_ratio_minus_1_core=float(np.max(np.abs(p[core]/r[core]-1))))
                    stats, arrays = {}, dict(edges=edges, produced_cv=p, reference_cv=r)
                    for band in sorted(set(a["bands"]) & set(b["bands"])):
                        stats[band] = {}
                        for side, source, cv, convention in [
                            ("produced", a, p, args.produced_content),
                            ("reference", b, r, args.reference_content)]:
                            u = rebin_density(source["bands"][band], source["edges"], edges, convention)
                            s = weight_stats(u, cv, widths)
                            stats[band][side] = s
                            bc = rebin_density(source["band_cvs"][band], source["edges"], edges, convention)
                            iw = s["integrated_weights"]
                            bands_summary.append(dict(case=label, band=band, source=side, universes=len(u),
                                integrated_weight_mean=float(iw.mean()), integrated_weight_std=float(iw.std()),
                                integrated_weight_RMS_about_1=float(np.sqrt(np.mean((iw-1)**2))),
                                integrated_weight_q16=float(np.quantile(iw,.16)),
                                integrated_weight_median=float(np.median(iw)),
                                integrated_weight_q84=float(np.quantile(iw,.84)),
                                band_CV_max_absolute_difference=float(np.max(np.abs(bc-cv)))))
                            arrays[side+"_"+band+"_universes"] = u
                        for i in range(len(p)):
                            row = dict(case=label, band=band, energy_low=edges[i], energy_high=edges[i+1],
                                produced_cv=p[i], reference_cv=r[i],
                                cv_ratio=p[i]/r[i] if r[i]>0 else float("nan"))
                            for side in ("produced", "reference"):
                                for key in ("mean", "std", "rms", "q16", "median", "q84"):
                                    row[side+"_weight_"+key] = stats[band][side][key][i]
                            bins_summary.append(row)
                    metadata.append(dict(case=label, common_edges=edges.tolist(),
                        produced=dict(edges=a["edges"].tolist(), x_title=a["x_title"], y_title=a["y_title"],
                            recovered=a["recovered"], underflow=a["underflow"], overflow=a["overflow"]),
                        reference=dict(edges=b["edges"].tolist(), x_title=b["x_title"], y_title=b["y_title"],
                            recovered=b["recovered"], underflow=b["underflow"], overflow=b["overflow"]),
                        produced_only_bands=sorted(set(a["bands"])-set(b["bands"])),
                        reference_only_bands=sorted(set(b["bands"])-set(a["bands"]))))
                    if args.plots:
                        plot_case(args.out/(label+".png"), label, edges, p, r, stats, args.plot_max)
                    if args.save_arrays:
                        np.savez_compressed(args.out/(label+".npz"), **arrays)
                    summaries.append(summary)
                    print("{}: {} / {} -> {} bins; {:.4g}--{:.4g} GeV; integral ratio={:.6g}; shape L1={:.4g}".format(
                        label,len(a["cv"]),len(b["cv"]),len(p),edges[0],edges[-1],ip/ir,summary["shape_L1_difference"]), flush=True)
                except Exception as exc:
                    failures.append(dict(case=label,error=str(exc)))
                    print("FAILED {}: {}".format(label,exc), flush=True)
    write_csv(args.out/"cv_summary.csv", summaries)
    write_csv(args.out/"weight_summary.csv", bands_summary)
    write_csv(args.out/"weights_by_energy.csv", bins_summary)
    (args.out/"binning_and_inputs.json").write_text(json.dumps(dict(
        produced_dir=args.produced_dir, reference_dir=args.reference_dir,
        produced_content=args.produced_content, reference_content=args.reference_content,
        note="Regular bins only. Common energy boundaries; no interpolation or playlist merging. Weight statistics use each input's own main CV. STD is mean-centered; RMS is about unity. Inputs must have matching area/POT units. These are flux-histogram weights, not individual event-weight distributions.",
        cases=metadata, failures=failures), indent=2))
    print("Completed {} comparisons; {} failures. Outputs: {}".format(len(summaries),len(failures),args.out))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
