#!/usr/bin/env python3
"""Run an unchanged MINERvA selection with either LE or ME1A MHRW renorm.

Usage (from CC-NuE-XSec):
 python studies/run_selection_mhrw_variant.py --variant le -- \
    --playlist le1_p8 --ntuple_tag MAD --use-sideband --truth --cal_POT \
    --selection_tag MHRW_LE --exclude-universes all --mc_only --test

Each variant MUST run in a fresh Python process. No MAT library changes.
"""
import argparse
import os
import runpy
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
FILENAME = 'Renorm_Kine_Truth_Minerva1_v10_HD_withElastics_RwgtNeutronCV_z5991-8408_a850.root'
LE = REPO / 'custom_mhrwKineRenorm' / FILENAME
ME = (REPO.parent / 'opt/lib/data/mhrwKineRenorm' /
      FILENAME.replace('Minerva1_', 'MinervaME1A_'))
SELECTION = REPO / 'selection/eventSelection.py'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--variant', required=True, choices=['le', 'me'])
    parser.add_argument('--selection-script', type=Path, default=SELECTION)
    parser.add_argument('selection_args', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    forwarded = args.selection_args
    if forwarded and forwarded[0] == '--':
        forwarded = forwarded[1:]

    source = LE if args.variant == 'le' else ME
    if not source.is_file():
        parser.error(f'Missing renormalization ROOT file: {source}')
    selection = args.selection_script.resolve()
    if not selection.is_file():
        parser.error(f'Missing selection script: {selection}')

    # The loader always constructs a filename based on the *analysis* playlist.
    # Symlink ME1A contents under that expected name, without editing the ROOT file.
    variant_dir = REPO / 'studies/mhrw_variant_files' / args.variant
    variant_dir.mkdir(parents=True, exist_ok=True)
    link = variant_dir / FILENAME
    if link.is_symlink():
        if link.resolve() != source.resolve():
            link.unlink()
    elif link.exists():
        if link.resolve() != source.resolve():
            parser.error(f'Refusing to overwrite an existing file: {link}')
    if not link.exists():
        link.symlink_to(source.resolve())

    # Some MINERvA configuration modules parse sys.argv at IMPORT time.
    # Forward selection arguments before importing any analysis module.
    sys.argv = [str(selection), *forwarded]
    sys.path.insert(0, str(REPO))
    sys.path.insert(0, str(selection.parent))
    import tools.SystematicsUniverse as SU

    # Must happen BEFORE the first GetGeantHadronWeight(). The ordinary
    # GetAllSystematicsUniverses() custom override is skipped via its guard.
    SU.CVSystematicUniverse.SetMHRWeightFilename(
        'Renorm_Kine_Truth', str(variant_dir)
    )
    SU.CVPythonUniverse._mhrw_filename_set = True

    print(f'[MHRW VARIANT] {args.variant.upper()}')
    print(f'[MHRW VARIANT] Requested renorm: {source}', flush=True)
    print(f'[MHRW VARIANT] Loader path:      {link}', flush=True)
    print(f'[MHRW VARIANT] Selection:        {selection}', flush=True)
    runpy.run_path(str(selection), run_name='__main__')


if __name__ == '__main__':
    main()
