#!/usr/bin/env python3
"""Check shared numerator statistics using the nominal ratio and direct ROOT files."""
import argparse
from pathlib import Path

import numpy as np


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ratio-file", type=Path, required=True)
    p.add_argument("--direct-file", type=Path, required=True)
    args = p.parse_args()

    import ROOT
    import PlotUtils  # registers MnvH1D with PyROOT

    def open_file(path):
        f = ROOT.TFile.Open(str(path), "READ")
        if not f or f.IsZombie():
            raise OSError(f"Cannot open {path}")
        return f

    def get(f, name):
        h = f.Get(name)
        if not h or not hasattr(h, "GetStatErrorMatrix"):
            raise ValueError(f"{f.GetName()}: missing MnvH1D {name}")
        return h

    def bin_matches(h1, h2, i):
        return np.allclose([h1.GetBinLowEdge(i), h1.GetBinLowEdge(i+1)],
                           [h2.GetBinLowEdge(i), h2.GetBinLowEdge(i+1)], rtol=0, atol=1e-9)

    ratio = open_file(args.ratio_file)
    direct = open_file(args.direct_file)
    failures = 0
    for beam, ratio_start, mu_start in (("fhc", 7, 31), ("rhc", 19, 43)):
        print(f"\n{beam.upper()} (v4.3.1 stitched layout; ROOT bin numbers)")
        for kind in ("data", "mc"):
            hr = get(ratio, f"{kind}_{beam}_ratio")
            hn = get(ratio, f"{kind}_{beam}_numu_selection")
            hd = get(direct, f"{kind}_{beam}_nue_selection")
            hg = get(ratio, f"sample_{kind}")
            stat = hg.GetStatErrorMatrix()
            print(f"{kind}: bin   muon       electron   ratio     expected Cov_stat(R,N)   stored   implied extra ratio stat var")
            for i in range(1, hr.GetNbinsX()+1):
                if not bin_matches(hr, hn, i) or not bin_matches(hr, hd, i):
                    raise ValueError(f"{beam} {kind} bin {i}: mismatched edges")
                n,d,r = (h.GetBinContent(i) for h in (hn,hd,hr))
                a,b = ratio_start+i-1,mu_start+i-1
                if d <= 0 or not np.allclose([r, r, n], [n/d, hg.GetBinContent(a),
                                                         hg.GetBinContent(b)], rtol=1e-8, atol=1e-10):
                    raise ValueError(f"{beam} {kind} bin {i}: component or stitched CV does not close")
                # Shared numerator term, assuming independent numerator/denominator
                # counting statistics; systematic correlations are separate.
                expected = hn.GetBinError(i)**2/d
                observed = float(stat(a, b))
                extra = float(stat(a, a)) - hn.GetBinError(i)**2/d**2
                failures += not np.isclose(expected, observed, rtol=1e-6, atol=1e-10)
                print(f"  {i:2d}  {n:10.4f} {d:10.4f} {r:9.4f}"
                      f" {expected:18.6g} {observed:12.6g} {extra:17.6g}")
    print(f"\nEntries differing from the shared-numerator statistical term: {failures}")
    print("This checks the stored statistical matrix only, not correlations in systematic bands.")


if __name__ == "__main__":
    main()
