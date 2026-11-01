#!/usr/bin/env python3
"""Compare MINERvA V1 vs v4.3.0 MC ROOT histograms (read-only).

Designed for the CCnumu FHC / CCnumubar LE5 selection outputs.
Requires PyROOT. Avoids TFile.Close()/Clone() object ownership issues.
"""
import argparse
import csv
import math
import os

DEFAULTS = {
    "FHC": ("/exp/minerva/data/users/qvuong/nu_mu", "leFHC", "ValidateCuts_CCnumu_MyCuts"),
    "RHC": ("/exp/minerva/data/users/qvuong/antinu_mu", "le5", "ValidateCuts_CCnumubar_MyCuts"),
}


def open_file(ROOT, path, opened):
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        raise RuntimeError("Cannot open " + path)
    opened.append(f)
    return f


def get_hist(f, key):
    h = f.Get(key)
    if not h or not h.InheritsFrom("TH1") or h.GetDimension() != 1:
        raise RuntimeError("Missing 1D histogram {} in {}".format(key, f.GetName()))
    return h


def edges(h):
    a = h.GetXaxis()
    return [a.GetBinLowEdge(i) for i in range(1, h.GetNbinsX() + 2)]


def values(h):
    return [float(h.GetBinContent(i)) for i in range(1, h.GetNbinsX() + 1)]


def ratio(n, d):
    return n / d if d != 0 else float("nan")


def compare(ROOT, beam, folder, playlist, stem, args, writer, opened):
    paths = {}
    for tune in ("v1", "v4p3p0"):
        for typ in ("data", "mc"):
            paths[(tune, typ)] = os.path.join(
                folder,
                "kin_dist_{}{}_{}_{}_MAD.root".format(typ, playlist, stem, tune),
            )
    files = {k: open_file(ROOT, p, opened) for k, p in paths.items()}
    print("\n" + "=" * 72 + "\n" + beam)
    for k, p in paths.items():
        print("{} {}: {}".format(*k, p))

    hdata1 = get_hist(files[("v1", "data")], args.hist)
    hdata4 = get_hist(files[("v4p3p0", "data")], args.hist)
    hmc1 = get_hist(files[("v1", "mc")], args.hist)
    hmc4 = get_hist(files[("v4p3p0", "mc")], args.hist)
    hs = (hdata1, hdata4, hmc1, hmc4)
    ref_edges = edges(hdata1)
    for h in hs[1:]:
        if len(edges(h)) != len(ref_edges) or any(abs(a-b) > 1e-8 * max(1, abs(a), abs(b)) for a, b in zip(ref_edges, edges(h))):
            raise RuntimeError("Histogram bin edges differ for {}".format(beam))

    d1, d4, m1, m4 = [values(h) for h in hs]
    data_diff = max((abs(a-b) for a,b in zip(d1,d4)), default=0.)
    print("Data V1-v4 max bin difference:", data_diff,
          "(MATCH)" if data_diff < 1e-8 else "(WARNING: data changed)")
    # Optionally pass actual data / MC POT ratio. Never assume shared normalization.
    scale1 = args.scale_v1
    scale4 = args.scale_v4
    print("MC scale V1 = {}, v4.3.0 = {}".format(scale1, scale4))
    print("{:>13} {:>12} {:>12} {:>11} {:>10} {:>10} {:>10}".format(
        "E range", "MC V1", "MC v4", "v4/V1", "D/MC V1", "D/MC v4", "data"))
    for i, (x, y, a, b, c, d) in enumerate(zip(ref_edges[:-1], ref_edges[1:], d1, d4, m1, m4)):
        print("{:>5.1f}-{:>5.1f} {:12.1f} {:12.1f} {:11.4f} {:10.4f} {:10.4f} {:10.1f}".format(
            x, y, c, d, ratio(d,c), ratio(a, scale1*c), ratio(b,scale4*d),a))
        writer.writerow({"beam":beam,"hist":args.hist,"low_GeV":x,"high_GeV":y,
                         "data_v1":a,"data_v4":b,"mc_v1":c,"mc_v4":d,
                         "mc_v4_over_v1":ratio(d,c),"data_over_mc_v1":ratio(a,scale1*c),
                         "data_over_mc_v4":ratio(b,scale4*d)})
    print("Visible-bin sums: D V1={:.3f}, D v4={:.3f}, MC V1={:.3f}, MC v4={:.3f}, MC v4/V1={:.5f}".format(
        sum(d1),sum(d4),sum(m1),sum(m4),ratio(sum(m4),sum(m1))))
    if args.categories:
        keys1 = {k.GetName() for k in files[("v1", "mc")].GetListOfKeys()}
        keys4 = {k.GetName() for k in files[("v4p3p0", "mc")].GetListOfKeys()}
        cats = sorted(k for k in (keys1 & keys4) if k.startswith(args.hist + "_") and
                      "_CV_WithErr" not in k)
        print("\nMatching truth-category histograms:")
        if not cats: print("(none found)")
        for cat in cats:
            try:
                c1=get_hist(files[("v1","mc")],cat)
                c4=get_hist(files[("v4p3p0","mc")],cat)
                print("{:45s} V1={:12.2f} v4={:12.2f} ratio={:8.4f}".format(
                    cat, sum(values(c1)),sum(values(c4)),ratio(sum(values(c4)),sum(values(c1)))))
            except RuntimeError:
                continue


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hist",default="Eel",help="1D histogram name")
    ap.add_argument("--scale-v1",type=float,default=1.0,
                    help="MC-to-data POT scale for V1; default 1 (raw ratio; D/MC not physical until correct)")
    ap.add_argument("--scale-v4",type=float,default=1.0,
                    help="MC-to-data POT scale for v4.3.0")
    ap.add_argument("--categories",action="store_true")
    ap.add_argument("--out",default="v1_v4p3p0_comparison.csv")
    ap.add_argument("--beam",choices=["both","FHC","RHC"],default="both")
    args=ap.parse_args()
    import ROOT
    ROOT.gROOT.SetBatch(True)
    opened=[]
    fields=["beam","hist","low_GeV","high_GeV","data_v1","data_v4","mc_v1","mc_v4",
            "mc_v4_over_v1","data_over_mc_v1","data_over_mc_v4"]
    with open(args.out,"w",newline="") as output:
        writer=csv.DictWriter(output,fieldnames=fields)
        writer.writeheader()
        for beam in ([args.beam] if args.beam!="both" else ["FHC","RHC"]):
            compare(ROOT,beam,*DEFAULTS[beam],args,writer,opened)
    print("\nWrote",os.path.abspath(args.out))
    # Do not manually close files: GPVM ROOT/cppyy versions can crash on close.

if __name__=="__main__":
    main()
