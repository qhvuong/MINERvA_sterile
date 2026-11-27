#!/usr/bin/env python3
"""Compare two completed component dumps. Requires only standard Python."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from tools.FSIEventPair import compare_runs

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--v430',required=True)
    parser.add_argument('--v431',required=True)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    summary=compare_runs(args.v430,args.v431,args.output)
    print('Comparison: {} ({} matched events)'.format(summary['status'],summary['matched_events']))
    for row in summary['totals']:
        print(row)
    sys.exit(0 if summary['status']=='passed' else 1)
