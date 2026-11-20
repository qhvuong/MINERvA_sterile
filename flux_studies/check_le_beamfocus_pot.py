#!/usr/bin/env python3
"""Read-only check of recovered LE BeamFocus weights and embedded POT.

Reconstruct every bin/universe from responses and the supplied pull CSV.
Writes an optional JSON report, never changes the ROOT input or pull table.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

DIALS = ['Horn1Offset', 'BaffleScraping', 'POT', 'HornIMiscal',
         'Horn1InnerConductor', 'TargetZ', 'WaterLayer']


def reconstruct(np, responses, pulls, historical):
    weights = np.ones((len(pulls), responses.shape[1]))
    for k, dial in enumerate(DIALS):
        twist = pulls[:, k] - (1. if historical and dial == 'Horn1InnerConductor' else 0.)
        weights *= 1. + twist[:, None] * responses[k][None, :]
    return np.maximum(weights, 0.)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--weights-file', type=Path, default=Path(
        '/exp/minerva/app/users/qvuong/FluxAna/LE_BeamFocus_universes_nue.root'))
    ap.add_argument('--pulls', type=Path, default=Path(
        '/exp/minerva/app/users/qvuong/FluxAna/LE_BeamFocus_dial_pulls_1000.csv'))
    ap.add_argument('--output', type=Path)
    args = ap.parse_args()
    import ROOT
    import numpy as np
    ROOT.gROOT.SetBatch(True)
    for p in (args.weights_file, args.pulls):
        if not p.is_file():
            raise FileNotFoundError(p)
    with args.pulls.open() as stream:
        reader = csv.DictReader(stream)
        missing = set(['universe'] + DIALS) - set(reader.fieldnames or [])
        if missing:
            raise ValueError('Missing CSV columns: ' + ', '.join(sorted(missing)))
        rows = sorted(reader, key=lambda row: int(row['universe']))
    if [int(row['universe']) for row in rows] != list(range(1000)):
        raise ValueError('CSV must contain each universe 0..999 exactly once')
    pulls = np.array([[float(row[dial]) for dial in DIALS] for row in rows])
    if not np.isfinite(pulls).all():
        raise ValueError('Nonfinite pulls')
    f = ROOT.TFile.Open(str(args.weights_file), 'READ')
    if not f or f.IsZombie() or f.TestBit(ROOT.TFile.kRecovered):
        raise RuntimeError('Invalid/recovered ROOT file: ' + str(args.weights_file))
    try:
        h = f.Get('BeamFocusWeights')
        if not h or not h.InheritsFrom('TH2'):
            raise ValueError('Missing TH2 BeamFocusWeights')
        nx, nu = h.GetNbinsX(), h.GetNbinsY()
        if nu != 1000 or nx < 1:
            raise ValueError(f'Unexpected dimensions: {nx} energy bins x {nu} universes')
        np.testing.assert_allclose([h.GetYaxis().GetBinCenter(j+1) for j in range(nu)],
                                   np.arange(nu), rtol=0, atol=1e-10)
        edges = np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1,nx+2)])
        observed = np.array([[h.GetBinContent(i+1,j+1) for i in range(nx)] for j in range(nu)])
        responses = []
        for dial in DIALS:
            r = f.Get('responses/' + dial)
            if not r or r.GetNbinsX() != nx:
                raise ValueError('Missing/incompatible saved response: responses/' + dial)
            np.testing.assert_allclose([r.GetXaxis().GetBinLowEdge(i) for i in range(1,nx+2)],
                                       edges, rtol=0, atol=1e-10)
            responses.append([r.GetBinContent(i) for i in range(1,nx+1)])
        responses = np.asarray(responses)
        if not np.isfinite(responses).all() or not np.isfinite(observed).all():
            raise ValueError('Nonfinite ROOT weights/responses')
        meta = {}
        for name in ['ResponseFile', 'PullTable', 'Mode', 'NonPOTResponsePrefix',
                     'POTResponsePrefix', 'BeamTag']:
            obj = f.Get(name)
            meta[name] = str(obj.GetString().Data()) if obj else None
        geometry = f.Get('HistoricalGeometryOffset')
        geometry_value = int(geometry.GetVal()) if geometry else None
        candidates = {}
        for historical in (False, True):
            rebuilt = reconstruct(np,responses,pulls,historical)
            candidates[historical] = dict(
                max_absolute_difference=float(np.max(np.abs(rebuilt-observed))),
                matches=bool(np.allclose(rebuilt,observed,rtol=1e-9,atol=1e-10)))
        if geometry_value is not None:
            if geometry_value not in (0,1):
                raise ValueError('Invalid HistoricalGeometryOffset metadata')
            matched = candidates[bool(geometry_value)]['matches']
            selected_geometry = bool(geometry_value)
        else:
            selected_geometry = min(candidates,key=lambda k:candidates[k]['max_absolute_difference'])
            matched = candidates[selected_geometry]['matches']
        pot_response = responses[DIALS.index('POT')]
        pot_is_two_percent = bool(np.allclose(pot_response,.02,rtol=0,atol=1e-12))
        factors = 1 + pulls[:,DIALS.index('POT')] * .02
        report = dict(weights_file=str(args.weights_file),pulls_file=str(args.pulls),
            pulls_sha256=hashlib.sha256(args.pulls.read_bytes()).hexdigest(),
            energy_bins=nx,universes=nu,metadata=meta,
            historical_geometry_metadata=geometry_value,
            matched_geometry='historical (-1+pull)' if selected_geometry else 'nominal-relative (pull)',
            reconstruction_checks={str(k):v for k,v in candidates.items()},
            all_stored_weights_match=matched,
            pot_response_min=float(pot_response.min()),pot_response_max=float(pot_response.max()),
            pot_response_is_flat_two_percent=pot_is_two_percent,
            pot_factor_min=float(factors.min()),pot_factor_max=float(factors.max()),
            pot_factor_population_std=float(factors.std()),
            stored_zero_weights=int(np.count_nonzero(observed==0)),
            status='PASS' if matched and pot_is_two_percent and np.all(factors>0) else 'FAIL')
        print(json.dumps(report,indent=2))
        if args.output:
            args.output.parent.mkdir(parents=True,exist_ok=True)
            args.output.write_text(json.dumps(report,indent=2)+'\n')
        if report['status'] != 'PASS':
            raise SystemExit(1)
    finally:
        f.Close()


if __name__ == '__main__':
    main()
