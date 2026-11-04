#!/usr/bin/env python3
"""Read-only check of the generated mu/e flux pairs used by SwappedFlux."""
import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools.SwappedFlux import _DEFAULT_FLUX_DIR, _get_flux_pair, _flux_value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--flux-dir', default=str(_DEFAULT_FLUX_DIR))
    parser.add_argument('--playlists', nargs='+', default=['minerva1', 'minerva13', 'minerva5'],
                        choices=['minerva1', 'minerva13', 'minerva5'])
    parser.add_argument('--energies', nargs='+', type=float, default=[1., 2., 3., 5., 10., 20.])
    args = parser.parse_args()
    directory = str(Path(args.flux_dir).expanduser().resolve())
    failures = 0
    for tag in args.playlists:
        for source_pdg, target_pdg in [(12, 14), (-12, -14)]:
            print('\n=== {}: generated PDG {} / {} ==='.format(tag, target_pdg, source_pdg))
            for pdg in [source_pdg, target_pdg]:
                print('file:', Path(directory) / 'flux-g4numiv5-pdg{}-{}.root'.format(pdg, tag))
            try:
                source, target = _get_flux_pair(directory, tag, source_pdg)
            except (RuntimeError, ValueError) as error:
                print('FAIL:', error)
                failures += 1
                continue
            n = source.GetNbinsX()
            axis = source.GetXaxis()
            widths = [axis.GetBinWidth(i) for i in range(1, n+1)]
            print('PASS: all {} edges match (including final upper edge).'.format(n+1))
            print('bins={} range=[{}, {}) width_min={} width_max={}'.format(
                n, axis.GetXmin(), axis.GetXmax(), min(widths), max(widths)))
            print('source title={!r}; x-axis={!r}'.format(source.GetTitle(), axis.GetTitle()))
            print('target title={!r}; x-axis={!r}'.format(target.GetTitle(), target.GetXaxis().GetTitle()))
            print('E_GeV  e_bin mu_bin   source_bin_content target_bin_content   ratio_bin   ratio_helper')
            for energy in args.energies:
                if not math.isfinite(energy) or not axis.GetXmin() <= energy < axis.GetXmax():
                    print('{}: outside histogram range; skipped'.format(energy))
                    continue
                se = source.GetXaxis().FindFixBin(energy)
                mu = target.GetXaxis().FindFixBin(energy)
                den = float(source.GetBinContent(se))
                num = float(target.GetBinContent(mu))
                raw = num/den if den > 0 else float('nan')
                d = _flux_value(source, energy)
                v = _flux_value(target, energy)
                interpolated = v/d if d > 0 else float('nan')
                print('{:6.2f} {:6d} {:6d} {:18.10g} {:18.10g} {:11.6g} {:14.6g}'.format(
                    energy, se, mu, den, num, raw, interpolated))
    print('\nSUMMARY: {} failed pairs.'.format(failures))
    print('Matching edges checks binning, not absolute flux units/POT normalization.')
    return int(failures > 0)


if __name__ == '__main__':
    sys.exit(main())
