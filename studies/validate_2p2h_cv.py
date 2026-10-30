#!/usr/bin/env python3
"""Run normal eventSelection with the low-recoil 2p2h CV factor ON or unity.
Run from the usual CC-NuE-XSec directory with CCAntiNuMu/CCNuMu config active.
Does not edit repository or existing ROOT files.
"""
import argparse
import math
import runpy
import sys
from collections import Counter

ap = argparse.ArgumentParser()
ap.add_argument('--mode', choices=['on', 'off'], required=True)
ap.add_argument('--max-print', type=int, default=15)
args, selection_args = ap.parse_known_args()
if selection_args and selection_args[0] == '--':
    selection_args = selection_args[1:]
if not selection_args:
    ap.error('Pass eventSelection.py arguments after --')

# AnalysisConfig parses sys.argv at import time. Strip wrapper-only args first.
sys.argv = ['eventSelection.py'] + selection_args

from tools import SystematicsUniverse as su

original = su.ROOT.PythonMinervaUniverse.GetLowRecoil2p2hWeight
counts = Counter()
examples = []

def checked_weight(self, *a, **kw):
    weight = float(original(self, *a, **kw))
    if getattr(self, 'mc_intType', None) == 8 and getattr(self, 'ShortName', lambda: '')() == 'cv':
        counts['calls'] += 1
        if not math.isfinite(weight):
            counts['nonfinite'] += 1
        elif weight == 1:
            counts['unity'] += 1
        elif weight > 1:
            counts['above_one'] += 1
        else:
            counts['below_one'] += 1
        if len(examples) < args.max_print:
            examples.append((getattr(self, 'mc_incoming', None), getattr(self, 'mc_incomingE', None), weight))
    return weight if args.mode == 'on' else 1.0

su.CVSystematicUniverse.GetLowRecoil2p2hWeight = checked_weight
print('[2P2H CV AUDIT] patched CVSystematicUniverse; original ROOT method:', original, flush=True)
print('[2P2H CV AUDIT] mode:', args.mode, 'selection args:', selection_args, flush=True)
try:
    sys.argv = ['eventSelection.py'] + selection_args
    runpy.run_path('selection/eventSelection.py', run_name='__main__')
finally:
    print('[2P2H CV AUDIT] MEC cv call counts:', dict(counts), flush=True)
    print('[2P2H CV AUDIT] first examples (PDG, incoming E [MeV], weight):', examples, flush=True)
