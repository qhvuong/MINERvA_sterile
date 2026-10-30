#!/usr/bin/env python3
"""Inspect raw MINERvA MasterAnaDev_nuHelicity in data and MC.

Run from the CC-NuE-XSec repository root with your ROOT/PlotUtils environment loaded.
This audits RAW ntuples before reconstruction selections. No selection files are modified.
"""
import argparse
from collections import Counter, defaultdict
import ROOT

ROOT.gROOT.SetBatch(True)


def scan(label, playlist, tag, nfiles, limit, beam):
    from tools import Utilities
    wrapper = Utilities.fileChain(playlist, label, tag, None, 0, nfiles)
    tree = wrapper.GetChain()
    names = {b.GetName() for b in tree.GetListOfBranches()}
    branch = 'MasterAnaDev_nuHelicity'
    if branch not in names:
        raise RuntimeError(f'{label}: missing branch {branch}; sample branches: {sorted(names)[:25]}')
    true_branch = next((b for b in ('mc_incoming', 'MasterAnaDev_mc_incoming') if b in names), None)
    total = tree.GetEntries()
    n = total if limit <= 0 else min(total, limit)
    counts = Counter()
    accepted = Counter()
    pdg_counts = defaultdict(Counter)
    expected_sign = -1 if beam == 'fhc' else 1
    for i in range(n):
        tree.GetEntry(i)
        v = int(getattr(tree, branch))
        counts[v] += 1
        mapped = -1 if v == 1 else 1
        if mapped == expected_sign:
            accepted[v] += 1
        if label == 'mc' and true_branch:
            pdg_counts[v][int(getattr(tree, true_branch))] += 1
    print(f'\n=== {label.upper()} | {beam.upper()} | {playlist} ===')
    print(f'Entries scanned: {n:,} of {total:,}; helicity cut accepted: {sum(accepted.values()):,} ({sum(accepted.values()) / n:.3%} of raw entries)' if n else 'No entries')
    print('raw_value | count | fraction | mapped_hel | current_cut | true PDG distribution (MC)')
    for v, cnt in sorted(counts.items()):
        m = -1 if v == 1 else 1
        outcome = 'PASS' if v in accepted else 'FAIL'
        pdgs = dict(sorted(pdg_counts[v].items())) if label == 'mc' else ''
        print(f'{v:>9} | {cnt:>8} | {cnt/n:>8.3%} | {m:>10} | {outcome:>11} | {pdgs}')
    print('NOTE: these are RAW event fractions, not weighted selection efficiencies.')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--playlist', default='le5_p8')
    p.add_argument('--ntuple-tag', default='MAD')
    p.add_argument('--files', type=int, default=1)
    p.add_argument('--max-events', type=int, default=200000, help='<=0 scans all events in chosen files')
    p.add_argument('--beam', choices=('fhc', 'rhc'), default='rhc')
    a = p.parse_args()
    for typ in ('data', 'mc'):
        scan(typ, a.playlist, a.ntuple_tag, a.files, a.max_events, a.beam)


if __name__ == '__main__':
    main()
