#!/usr/bin/env python3
"""Run existing eventSelection.py with original or unity MINOS CV factor.

Execute from CC-NuE-XSec directory with its standard experiment environment.
Does not edit production files. Use two DIFFERENT output directories.
"""
import argparse
import os
import runpy
import sys

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--mode', choices=('on', 'off'), required=True)
parser.add_argument('--selection-script', default='selection/eventSelection.py')
args, remaining = parser.parse_known_args()
if remaining and remaining[0] == '--':
    remaining = remaining[1:]
script = os.path.abspath(args.selection_script)
if not os.path.isfile(script):
    parser.error('Not found: ' + script)
if not remaining:
    parser.error('Provide eventSelection arguments after --')
if '--mc_only' not in remaining and '--mc-only' not in remaining:
    parser.error('For this test use --mc_only; do not run on data')
if '--output' not in remaining and '-o' not in remaining:
    parser.error('Provide --output pointing to a separate temporary directory')

# AnalysisConfig parses argv during its import. Set it BEFORE any analysis import.
sys.argv = [script] + remaining

# Ensure the same imports are resolved as normal eventSelection.
sys.path.insert(0, os.getcwd())
from tools import SystematicsUniverse as su

orig = su.CVPythonUniverse.GetMyMinosEfficiencyWeight
if args.mode == 'off':
    def unity_minos(self):
        return 1.0
    su.CVPythonUniverse.GetMyMinosEfficiencyWeight = unity_minos

print('[MINOS CV VALIDATION] mode=%s; MINOS weight %s; all other weight logic unchanged' %
      (args.mode, 'ORIGINAL' if args.mode == 'on' else 'FORCED TO 1.0'), flush=True)
print('[MINOS CV VALIDATION] eventSelection args:', ' '.join(remaining), flush=True)
runpy.run_path(script, run_name='__main__')
