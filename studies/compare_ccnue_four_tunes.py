
#!/usr/bin/env python3
"""Read-only CCnue / CCnuebar tune comparison."""

import argparse
import csv
import glob
from pathlib import Path

# ROOT must load before matplotlib in the FNAL environment.
import ROOT
ROOT.gROOT.SetBatch(True)

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


BASE = "/exp/minerva/data/users/qvuong"

SERIES = {
    "old": {
        "v4.3.0": "ValidateCuts_{sample}_FullCuts_noEavail_v4p3p0",
        "v4.3.0 LECOH": "ValidateCuts_{sample}_FullCuts_noEavail_v4p3p0_LECOH",
        "v4.4.1": "ValidateCuts_{sample}_FullCuts_noEavail_v4p4p1",
    },
    "new": {
        "v1": "ValidateTune_{sample}_noEavail_v1",
        "v4.3.1": "ValidateTune_{sample}_noEavail_v4p3p1",
        "v4.3.1 no MENU1PI": "ValidateTune_{sample}_noEavail_v4p3p1_noMENU1PI",
        "v4.4.1": "ValidateTune_{sample}_noEavail_v4p4p1",
    },
}

SAMPLES = {
    "CCnue": "nu_e",
    "CCnuebar": "antinu_e",
}

ALIASES = {
    "Eavail": "E_avail",
    "Visible Energy": "E_avail",
    "Lepton Energy": "Eel",
    "Lepton Pt": "Lepton_Pt",
    "Lepton Ppar": "Lepton_Ppar",
}

TWO_D = {
    "Elepton_Eavail": ("Electron energy", "Available energy"),
    "Lepton_Pt_vs_Eavail": ("Electron pT", "Available energy"),
    "Lepton_Ppar_vs_Eavail": ("Electron p-parallel", "Available energy"),
}

COMPONENTS = [
    "CCNuE", "CCNuE2p2h", "CCNuEDIS", "CCNuEDelta",
    "CCNuEQE", "CCNuEWrongSign", "CCOther", "CCPi",
    "CCPi0", "NCCohPi0", "NCDiff", "NCOther",
    "NCPi", "NCPi0", "NuEElastic", "Other",
]


def open_file(path):
    f = ROOT.TFile.Open(str(path), "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open {path}")
    return f


def find_file(base, folder, kind, playlist, tag):
    path = (
        Path(base) / folder /
        f"kin_dist_{kind}{playlist}_{tag}_MAD.root"
    )
    if path.is_file():
        return str(path)

    # Fallback only if there is exactly one match.
    pattern = str(
        Path(base) / folder /
        f"kin_dist_{kind}*_{tag}_MAD.root"
    )
    matches = sorted(glob.glob(pattern))

    if len(matches) == 1:
        print("[FALLBACK]", matches[0])
        return matches[0]

    if matches:
        print("[AMBIGUOUS]", kind, tag)
        for m in matches:
            print("  ", m)
    else:
        print("[MISSING]", path)

    return None


def get_pot(f):
    meta = f.Get("Meta")
    if not meta or not meta.GetBranch("POT_Used"):
        raise RuntimeError(
            f"Missing Meta.POT_Used in {f.GetName()}"
        )
    return sum(float(e.POT_Used) for e in meta)


def read_hist(f, name, dim=1):
    h = f.Get(name)
    if not h or not h.InheritsFrom("TH1"):
        return None
    if h.GetDimension() != dim:
        return None

    x = np.array([
        h.GetXaxis().GetBinLowEdge(i)
        for i in range(1, h.GetNbinsX() + 2)
    ])

    if dim == 1:
        v = np.array([
            h.GetBinContent(i)
            for i in range(1, h.GetNbinsX() + 1)
        ], dtype=float)
        e = np.array([
            h.GetBinError(i)
            for i in range(1, h.GetNbinsX() + 1)
        ], dtype=float)
        return x, v, e

    y = np.array([
        h.GetYaxis().GetBinLowEdge(j)
        for j in range(1, h.GetNbinsY() + 2)
    ])
    v = np.array([
        [
            h.GetBinContent(i, j)
            for j in range(1, h.GetNbinsY() + 1)
        ]
        for i in range(1, h.GetNbinsX() + 1)
    ], dtype=float)

    return x, y, v


def ratio(a, b):
    return np.divide(
        a, b,
        out=np.full_like(a, np.nan, dtype=float),
        where=b != 0
    )


def save_fig(fig, path):
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(str(path) + "." + ext, dpi=180)
    plt.close(fig)


def plot_1d(
    name, data, predictions, out, title, xlabel
):
    edges, d, de = data
    centers = 0.5 * (edges[:-1] + edges[1:])

    fig, (ax, ra) = plt.subplots(
        2, 1, figsize=(10, 7),
        sharex=True,
        gridspec_kw={"height_ratios": [2.5, 1.2]}
    )

    ax.errorbar(
        centers, d, yerr=de,
        fmt="ko", markersize=3.5,
        label="Data", zorder=10
    )

    for label, m in predictions.items():
        ax.stairs(m, edges, label=label, linewidth=1.7)
        ra.plot(
            centers, ratio(d, m),
            marker=".", markersize=5, label=label
        )

    ra.axhline(1, color="gray", linestyle="--")
    ax.set_ylabel("Events / bin")
    ra.set_ylabel("Data / MC")
    ra.set_xlabel(xlabel)
    ax.set_title(title)
    ax.legend(fontsize=9, ncol=2)
    ax.grid(alpha=0.2)
    ra.grid(alpha=0.2)

    save_fig(fig, out / f"{name}_datamc")

    if len(predictions) > 1:
        ref = next(iter(predictions))
        fig, ax = plt.subplots(figsize=(10, 4.5))
        for label, m in predictions.items():
            if label == ref:
                continue
            ax.plot(
                centers, ratio(m, predictions[ref]),
                marker=".", label=f"{label} / {ref}"
            )
        ax.axhline(1, color="gray", linestyle="--")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("MC tune ratio")
        ax.set_title(title)
        ax.grid(alpha=0.2)
        ax.legend()
        save_fig(fig, out / f"{name}_tune_ratios")


def analyze_1d(sample, series, records, hist, out, components):
    hist = ALIASES.get(hist, hist)

    ref = next(
        (r for r in records.values() if r["data"]),
        None
    )
    if ref is None:
        return

    data = read_hist(ref["data"], hist)
    if data is None:
        print("[SKIP] missing data histogram:", sample, hist)
        return

    edges, d, de = data
    predictions = {}
    valid = {}

    for label, r in records.items():
        h = read_hist(r["mc"], hist)
        if h is None:
            print("[SKIP] missing MC histogram", label, hist)
            continue

        e, m, me = h
        if len(e) != len(edges) or not np.allclose(e, edges):
            print("[BINNING MISMATCH]", sample, label, hist)
            print("  data:", edges.tolist())
            print("  MC:  ", e.tolist())
            continue

        predictions[label] = m * r["scale"]
        valid[label] = r

    if not predictions:
        return

    out.mkdir(parents=True, exist_ok=True)

    print(f"\n{sample} {series} {hist}")
    print(f"  DATA: {d.sum():.3f}")

    for label, m in predictions.items():
        print(
            f"  {label:<24} MC={m.sum():.3f}"
            f"  D/MC={d.sum()/m.sum():.5f}"
        )

    columns = [
        "bin", "low", "high", "data", "data_error"
    ]
    for label in predictions:
        columns += [f"mc_{label}", f"data_mc_{label}"]

    with (out / f"{hist}_bins.csv").open(
        "w", newline=""
    ) as f:
        w = csv.DictWriter(f, fieldnames=columns)
        w.writeheader()
        for i in range(len(d)):
            row = {
                "bin": i + 1,
                "low": edges[i],
                "high": edges[i+1],
                "data": d[i],
                "data_error": de[i],
            }
            for label, m in predictions.items():
                row[f"mc_{label}"] = m[i]
                row[f"data_mc_{label}"] = (
                    d[i]/m[i] if m[i] else float("nan")
                )
            w.writerow(row)

    plot_1d(
        hist, data, predictions, out,
        f"{sample} {series}: {hist} (POT normalized)",
        hist + " (GeV)"
    )

    if not components:
        return

    rows = []
    for component in COMPONENTS:
        key = f"{hist}_{component}"
        yields = {}
        for label, r in valid.items():
            h = read_hist(r["mc"], key)
            if h is not None:
                yields[label] = float(h[1].sum() * r["scale"])
        if yields:
            rows.append({"component": component, **yields})

    with (out / f"{hist}_components.csv").open(
        "w", newline=""
    ) as f:
        w = csv.DictWriter(
            f,
            fieldnames=["component", *valid.keys()],
            extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(rows)

    print("  Component yields:")
    for row in rows:
        print("   ", row)


def analyze_2d(sample, records, out):
    ref = next(
        (r for r in records.values() if r["data"]),
        None
    )
    if ref is None:
        return

    for name, (xlabel, ylabel) in TWO_D.items():
        data = read_hist(ref["data"], name, dim=2)

        if data is None:
            print("[2D SKIP]", sample, name)
            continue

        x, y, d = data
        preds = {}

        for label, r in records.items():
            h = read_hist(r["mc"], name, dim=2)
            if h is None:
                continue
            xx, yy, m = h

            if not (
                np.array_equal(x, xx)
                and np.array_equal(y, yy)
            ):
                print("[2D BINNING MISMATCH]", label, name)
                continue

            preds[label] = m * r["scale"]

        if not preds:
            continue

        rows = []
        for ix in range(len(x) - 1):
            ds = d[ix]
            curves = {
                label: m[ix]
                for label, m in preds.items()
            }

            if ds.sum() == 0 and all(
                m.sum() == 0 for m in curves.values()
            ):
                continue

            for label, m in curves.items():
                rows.append({
                    "variable": name,
                    "x_low": x[ix],
                    "x_high": x[ix+1],
                    "tune": label,
                    "data": ds.sum(),
                    "mc": m.sum(),
                    "data_mc": (
                        ds.sum()/m.sum()
                        if m.sum() else float("nan")
                    ),
                })

            de = np.sqrt(np.maximum(ds, 0))
            plot_1d(
                f"{name}_slice_{ix+1:02d}",
                (y, ds, de),
                curves,
                out,
                f"{sample}: {xlabel} [{x[ix]:g}, {x[ix+1]:g}]",
                ylabel + " (GeV)"
            )

        out.mkdir(parents=True, exist_ok=True)
        with (out / f"{name}_slice_summary.csv").open(
            "w", newline=""
        ) as f:
            fields = [
                "variable", "x_low", "x_high",
                "tune", "data", "mc", "data_mc"
            ]
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)

        print("[2D]", sample, name, len(rows), "rows")


def analyze_sample(args, sample, series):
    folder = SAMPLES[sample]
    playlist = (
        args.fhc_playlist if sample == "CCnue" else "le5"
    )

    records = {}
    handles = []

    for label, template in SERIES[series].items():
        tag = template.format(sample=sample)

        mcpath = find_file(
            args.base, folder, "mc", playlist, tag
        )
        datapath = find_file(
            args.base, folder, "data", playlist, tag
        )

        if not mcpath:
            continue

        mc = open_file(mcpath)
        data = open_file(datapath) if datapath else None
        handles.append(mc)
        if data:
            handles.append(data)

        mcpot = get_pot(mc)
        dpot = get_pot(data) if data else None

        records[label] = {
            "mc": mc,
            "data": data,
            "mcpot": mcpot,
            "datapot": dpot,
        }

        print(
            f"{sample} {series} {label}:"
            f" MC POT={mcpot}, DATA POT={dpot}"
        )

    if not records:
        return

    dataref = next(
        (r for r in records.values() if r["data"]),
        None
    )
    if dataref is None:
        print("[ERROR] No data file:", sample, series)
        return

    dpot = dataref["datapot"]
    if not dpot or dpot <= 0:
        raise RuntimeError("Invalid data POT")

    for label, r in records.items():
        if r["datapot"] is not None and not np.isclose(
            r["datapot"], dpot, rtol=1e-6
        ):
            print("[WARNING] Data POT differs:", label)
        if not r["mcpot"] or r["mcpot"] <= 0:
            raise RuntimeError(f"Invalid MC POT for {label}")
        r["scale"] = dpot / r["mcpot"]

    out = Path(args.outdir) / series / sample

    for hist in args.hists:
        analyze_1d(
            sample, series, records, hist,
            out, args.components
        )

    if series == "new" and args.two_d:
        analyze_2d(sample, records, out)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--base", default=BASE)
    p.add_argument(
        "--outdir",
        default=BASE + "/studies/ccnue_tune_comparison"
    )
    p.add_argument(
        "--series",
        choices=["old", "new", "both"],
        default="both"
    )
    p.add_argument(
        "--samples", nargs="+",
        choices=list(SAMPLES),
        default=list(SAMPLES)
    )
    p.add_argument(
        "--hists", nargs="+",
        default=["E_avail", "Eel", "Lepton_Pt"]
    )
    p.add_argument("--components", action="store_true")
    p.add_argument("--two-d", action="store_true")
    p.add_argument(
        "--fhc-playlist",
        choices=["le1", "leFHC"],
        default="le1"
    )
    args = p.parse_args()

    series_list = (
        ["old", "new"]
        if args.series == "both"
        else [args.series]
    )

    for series in series_list:
        for sample in args.samples:
            analyze_sample(args, sample, series)


if __name__ == "__main__":
    main()
