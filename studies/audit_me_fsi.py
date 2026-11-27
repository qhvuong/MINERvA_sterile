#!/usr/bin/env python3
import argparse
import csv
import json
import math
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--base", type=Path,
        default=Path("/exp/minerva/data/users/qvuong"),
    )
    parser.add_argument("--out", type=Path, default=Path("ME_FSI_audit"))
    parser.add_argument(
        "--hist", action="append",
        help="Histogram name; repeat as needed. Default: discover histograms with DIS categories.",
    )
    parser.add_argument("--dis-category", default="CCNuMuDIS")
    parser.add_argument("--list", action="store_true", help="List histogram keys and exit.")
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument(
        "--pot-json", type=Path,
        help="Optional JSON mapping ROOT filenames to explicit POT overrides.",
    )
    args = parser.parse_args()

    # Load ROOT before NumPy/Matplotlib in the usual analysis environment.
    import ROOT
    import PlotUtils
    import numpy as np

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)

    tunes = ("v4p3p0", "v4p3p1")
    categories = (
        "NCDiff", "NuEElastic", "NC", "CCNuMuWrongSign",
        "CCNuMuQE", "CCNuMuDelta", "CCNuMuDIS",
        "CCNuMu2p2h", "CCNuMu", "CC",
    )
    overrides = json.loads(args.pot_json.read_text()) if args.pot_json else {}
    args.out.mkdir(parents=True, exist_ok=True)

    files = {}
    paths = {}
    for beam, directory, selection in (
        ("FHC", "nu_mu", "ME_CCnumu"),
        ("RHC", "antinu_mu", "ME_CCnumubar"),
    ):
        for sample in ("data", "mc"):
            for tune in tunes:
                name = (
                    f"kin_dist_{sample}me{beam}_"
                    f"{selection}_RyansCuts_sameOffset_{tune}_MAD.root"
                )
                path = args.base / directory / name
                fin = ROOT.TFile.Open(str(path), "READ")
                if not fin or fin.IsZombie():
                    raise RuntimeError(f"Cannot open {path}")
                if fin.TestBit(ROOT.TFile.kRecovered):
                    raise RuntimeError(f"Recovered ROOT file: {path}")
                key = (beam, sample, tune)
                files[key], paths[key] = fin, path

    def histograms(fin):
        result = {}
        for name in sorted({k.GetName() for k in fin.GetListOfKeys()}):
            obj = fin.Get(name)
            if obj and obj.InheritsFrom("TH1"):
                result[name] = obj
        return result

    inventory = {key: histograms(fin) for key, fin in files.items()}

    if args.list:
        for key, hists in inventory.items():
            print("\n", key, paths[key])
            for name, hist in hists.items():
                print(f"  {name}: dimension={hist.GetDimension()}")
        return

    def cv(hist):
        h = (
            hist.GetCVHistoWithStatError()
            if hasattr(hist, "GetCVHistoWithStatError")
            else hist.Clone()
        )
        h.SetDirectory(0)
        return h

    def cells(hist, errors=False):
        getter = hist.GetBinError if errors else hist.GetBinContent
        return np.array([getter(i) for i in range(hist.GetNcells())])

    def axes(hist):
        getters = (hist.GetXaxis, hist.GetYaxis, hist.GetZaxis)
        return [
            [
                float(axis.GetBinLowEdge(i))
                for i in range(1, axis.GetNbins() + 2)
            ]
            for axis in (g() for g in getters[:hist.GetDimension()])
        ]

    def compare_hist(a, b):
        """Compare every cell, including underflow/overflow."""
        if a.GetDimension() != b.GetDimension() or axes(a) != axes(b):
            return {"equal": False, "reason": "dimension/binning differs"}
        result = {"equal": True}
        for label, errors in (("contents", False), ("stat_errors", True)):
            x, y = cells(a, errors), cells(b, errors)
            bad = np.flatnonzero(~np.isclose(x, y, rtol=0, atol=0))
            finite = bool(np.all(np.isfinite(x)) and np.all(np.isfinite(y)))
            result[label] = {
                "different_cells": int(len(bad)),
                "finite": finite,
                "max_abs_difference": float(np.max(np.abs(y - x))) if finite else None,
                "examples": [
                    {"global_bin": int(i), "old": float(x[i]), "new": float(y[i])}
                    for i in bad[:10]
                ],
            }
            result["equal"] &= not len(bad) and finite
        result["entries_equal"] = bool(a.GetEntries() == b.GetEntries())
        result["equal"] &= result["entries_equal"]
        return result

    def audit_data(a, b):
        result = {"nominal": compare_hist(cv(a), cv(b))}
        for kind in ("Vert", "Lat"):
            get_names = f"Get{kind}ErrorBandNames"
            get_band = f"Get{kind}ErrorBand"
            na = set(map(str, getattr(a, get_names)())) if hasattr(a, get_names) else set()
            nb = set(map(str, getattr(b, get_names)())) if hasattr(b, get_names) else set()
            result[kind + "_band_names_equal"] = na == nb
            for name in sorted(na & nb):
                ba, bb = getattr(a, get_band)(name), getattr(b, get_band)(name)
                ok = ba.GetNHists() == bb.GetNHists()
                failures = []
                if ok:
                    for u in range(ba.GetNHists()):
                        check = compare_hist(ba.GetHist(u), bb.GetHist(u))
                        if not check["equal"]:
                            failures.append({"universe": u, "comparison": check})
                result[f"{kind}:{name}"] = {
                    "equal": ok and not failures,
                    "failures": failures[:10],
                }
        result["equal"] = (
            result["nominal"]["equal"]
            and result["Vert_band_names_equal"]
            and result["Lat_band_names_equal"]
            and all(
                value["equal"]
                for name, value in result.items()
                if ":" in name
            )
        )
        return result

    def pot(key):
        path = paths[key]
        override = overrides.get(str(path), overrides.get(path.name))
        if override is not None:
            value, method = float(override), "explicit override"
        else:
            fin = files[key]
            meta = fin.Get("Meta")
            used = None
            if meta and meta.InheritsFrom("TTree") and meta.GetLeaf("POT_Used"):
                used = 0.0
                for entry in range(meta.GetEntries()):
                    meta.GetEntry(entry)
                    used += float(meta.GetLeaf("POT_Used").GetValue())

            value, method = None, None
            if meta:
                try:
                    chain = ROOT.TChain("Meta")
                    chain.Add(str(path))
                    candidate = float(
                        ROOT.PlotUtils.POTCounter().getPOTfromTChain(chain)
                    )
                    if math.isfinite(candidate) and candidate > 0:
                        value, method = candidate, "PlotUtils.POTCounter"
                except Exception:
                    pass

            if used is not None:
                if value is not None and not math.isclose(
                    value, used, rel_tol=1e-9, abs_tol=0
                ):
                    raise RuntimeError(
                        f"POTCounter and sum(Meta.POT_Used) disagree: {path}"
                    )
                value, method = used, "sum(Meta.POT_Used)"

            if value is None:
                for name in ("POT_Used", "used_POT"):
                    obj = fin.Get(name)
                    if obj and hasattr(obj, "GetVal"):
                        value, method = float(obj.GetVal()), name
                        break

        if value is None or not math.isfinite(value) or value <= 0:
            raise RuntimeError(
                f"Cannot establish used POT for {path}. "
                "Supply --pot-json; do not assume unit normalization."
            )
        return value, method

    exposure = {key: pot(key) for key in files}
    report = {
        "files": {
            ":".join(key): {
                "path": str(paths[key]),
                "POT": exposure[key][0],
                "POT_source": exposure[key][1],
            }
            for key in files
        },
        "data_checks": {},
        "yields": [],
    }
    per_bin = []
    data_ok = True

    def values(hist):
        return np.array([
            hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)
        ])

    def edges(hist):
        return np.array([
            hist.GetBinLowEdge(i) for i in range(1, hist.GetNbinsX() + 2)
        ])

    def divide(a, b):
        return np.divide(
            a, b, out=np.full_like(np.asarray(a, dtype=float), np.nan),
            where=np.asarray(b) != 0,
        )

    for beam in ("FHC", "RHC"):
        dkeys = [(beam, "data", t) for t in tunes]
        mkeys = [(beam, "mc", t) for t in tunes]
        d0, d1 = (inventory[k] for k in dkeys)
        m0, m1 = (inventory[k] for k in mkeys)
        p0, p1 = (exposure[k][0] for k in dkeys)

        check = {
            "same_histogram_keys": set(d0) == set(d1),
            "same_POT": math.isclose(p0, p1, rel_tol=1e-10, abs_tol=0),
            "histograms": {},
        }
        for name in sorted(set(d0) & set(d1)):
            check["histograms"][name] = audit_data(d0[name], d1[name])
        check["equal"] = (
            check["same_histogram_keys"] and check["same_POT"]
            and all(x["equal"] for x in check["histograms"].values())
        )
        data_ok &= check["equal"]
        report["data_checks"][beam] = check
        print(f"\n{beam}: DATA IDENTITY {'PASS' if check['equal'] else 'FAIL'}")
        print(f"  data POT: {p0:.12g}, {p1:.12g}")

        names = args.hist or sorted(
            name for name in set(d0) & set(d1) & set(m0) & set(m1)
            if d0[name].GetDimension() == 1
            and name + "_" + args.dis_category in m0
            and name + "_" + args.dis_category in m1
        )
        if not names:
            raise RuntimeError(f"No DIS-split observables found for {beam}; run --list.")

        # Both MC tunes are compared at the same reference data exposure.
        scales = [p0 / exposure[k][0] for k in mkeys]
        print(f"  MC POT: {[exposure[k][0] for k in mkeys]}")
        print(f"  MC scale to common data POT: {scales}")

        for name in names:
            for collection in (d0, d1, m0, m1):
                if name not in collection or collection[name].GetDimension() != 1:
                    raise RuntimeError(f"Missing 1D histogram {beam}/{name}")

            data = cv(d0[name])
            data_new = cv(d1[name])
            h0, h1 = cv(m0[name]), cv(m1[name])
            if not (axes(data) == axes(data_new) == axes(h0) == axes(h1)):
                raise RuntimeError(f"Binning mismatch: {beam}/{name}")

            xedges = edges(data)
            data_values = values(data)
            component_values = {}
            available = [
                c for c in categories
                if name + "_" + c in m0 and name + "_" + c in m1
            ]
            for category in ("TOTAL",) + tuple(available):
                key = name if category == "TOTAL" else name + "_" + category
                ha, hb = cv(m0[key]), cv(m1[key])
                if axes(ha) != axes(data) or axes(hb) != axes(data):
                    raise RuntimeError(f"Category binning mismatch: {key}")
                raw0, raw1 = values(ha), values(hb)
                if not np.all(np.isfinite(raw0)) or not np.all(np.isfinite(raw1)):
                    raise RuntimeError(f"Nonfinite MC content: {beam}/{key}")
                norm0, norm1 = raw0 * scales[0], raw1 * scales[1]
                component_values[category] = (norm0, norm1)
                total0, total1 = float(norm0.sum()), float(norm1.sum())
                ratio = total1 / total0 if total0 else None
                row = {
                    "beam": beam, "histogram": name, "category": category,
                    "raw_MC_430": float(raw0.sum()),
                    "raw_MC_431": float(raw1.sum()),
                    "normalized_MC_430": total0,
                    "normalized_MC_431": total1,
                    "MC_431_over_430": ratio,
                    "change_percent": 100 * (ratio - 1) if ratio is not None else None,
                    "underflow_430": ha.GetBinContent(0) * scales[0],
                    "underflow_431": hb.GetBinContent(0) * scales[1],
                    "overflow_430": ha.GetBinContent(ha.GetNbinsX() + 1) * scales[0],
                    "overflow_431": hb.GetBinContent(hb.GetNbinsX() + 1) * scales[1],
                }
                report["yields"].append(row)
                print(f"  {name:20s} {category:18s}: "
                      f"{total0:.7g} -> {total1:.7g}; ratio={ratio}")
                for i in range(len(norm0)):
                    per_bin.append({
                        "beam": beam, "histogram": name, "category": category,
                        "bin": i + 1, "low": xedges[i], "high": xedges[i + 1],
                        "MC_430": norm0[i], "MC_431": norm1[i],
                        "MC_431_over_430": norm1[i] / norm0[i] if norm0[i] else None,
                        "data_430": data_values[i] if category == "TOTAL" else None,
                    })

            if args.no_plots:
                continue
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots(2, 2, figsize=(12, 8))
            a, b = component_values["TOTAL"]
            centers = (xedges[:-1] + xedges[1:]) / 2
            errors = np.array([
                data.GetBinError(i) for i in range(1, data.GetNbinsX() + 1)
            ])
            for vals, tune, color in ((a, "4.3.0", "C0"), (b, "4.3.1", "C1")):
                ax[0, 0].stairs(vals, xedges, label=f"MC {tune}", color=color)
                ax[1, 0].stairs(
                    divide(vals, data_values), xedges,
                    label=f"MC {tune} / data 4.3.0", color=color,
                )
            ax[0, 0].errorbar(
                centers, data_values, yerr=errors, fmt="ko", label="Data 4.3.0"
            )
            if not check["equal"]:
                ax[0, 0].errorbar(
                    centers, values(data_new) * p0 / p1, fmt="rx",
                    label="Data 4.3.1, common POT",
                )
            ax[1, 1].stairs(divide(b, a), xedges, label="Total MC")
            dis = component_values.get(args.dis_category)
            if dis is not None:
                da, db = dis
                ax[0, 1].stairs(da, xedges, label="DIS 4.3.0", color="C0")
                ax[0, 1].stairs(db, xedges, label="DIS 4.3.1", color="C1")
                ax[1, 1].stairs(
                    divide(db, da), xedges, label="Right-sign CC DIS"
                )
            ax[0, 0].set_title(f"{beam}: selected total")
            ax[0, 1].set_title(f"{beam}: {args.dis_category}")
            ax[0, 0].set_ylabel("Events / bin")
            ax[0, 1].set_ylabel("Events / bin")
            ax[1, 0].set_ylabel("Total MC / data")
            ax[1, 1].set_ylabel("MC 4.3.1 / MC 4.3.0")
            for panel in ax.flat:
                panel.set_xlabel(data.GetXaxis().GetTitle() or name)
                panel.grid(alpha=.2)
                panel.legend(fontsize=8)
            for panel in ax[1]:
                panel.axhline(1, color="black", linestyle="--", linewidth=1)
            fig.suptitle(
                "Same reference POT; no shape normalization. "
                "Tune-ratio curves are CV diagnostics without paired-MC error bars.",
                fontsize=10,
            )
            fig.tight_layout()
            safe = "".join(c if c.isalnum() or c in "_-" else "_" for c in name)
            for extension in ("png", "pdf"):
                fig.savefig(args.out / f"{beam}_{safe}.{extension}")
            plt.close(fig)

    (args.out / "audit.json").write_text(json.dumps(report, indent=2) + "\n")
    for filename, rows in (
        ("integrated_yields.csv", report["yields"]),
        ("bin_yields.csv", per_bin),
    ):
        if rows:
            with (args.out / filename).open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
    for fin in files.values():
        fin.Close()
    print("\nOverall data identity:", "PASS" if data_ok else "FAIL")
    print("Outputs:", args.out.resolve())
    if not data_ok:
        raise SystemExit(2)


if __name__ == "__main__":
    main()