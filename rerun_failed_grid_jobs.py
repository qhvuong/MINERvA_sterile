#!/usr/bin/env python3
"""Generate rerun commands for missing or invalid ROOT outputs in MINERvA combine logs.

Run from CC-NuE-XSec/selection on a GPVM; no jobs are submitted.
Reads Log_<stamp>.txt plus the matching <stamp>.txt submission log.
Copies original gridSelection arguments (truth, POT, sidebands, exclusions,
selection, tarball, etc.) and changes only the job slice to one failed job.
"""
import argparse
import re
import shlex
import sys
from collections import defaultdict
from pathlib import Path

SKIP = re.compile(r'^\s*\[SKIP\]\s+(?:Invalid ROOT file:|ROOT open/check timed out[^:]*:|.*?Timeout[^:]*:)\s*(\S+?\.root)(?:\s|$)')
FILE = re.compile(r'^kin_dist_(mc|data)(.+)_([^_]+)_(\d+)\.root$')
STAMP = re.compile(r'^Diagnostic log generated from runningNotes stamp:\s*(\S+)')
TAG = re.compile(r'^Selection tag \(from file\):\s*(\S+)')
RUN = re.compile(r'^\s*Running command:\s*$')
OUTPUT = re.compile(r'^\[GRID_OUTPUT_HISTS\]\s+(\S+)')
TARBALL = re.compile(r'^Using shared tarball:\s+(\S+)')
JOBSUB = re.compile(r'^\[GRID_JOBSUB_CMD\].*?\s-N\s+(\d+)(?:\s|$)')

# Options replaced with job-specific settings, rather than copied from the
# original submission. Preserve all other original analysis flags.
REPLACE = {'--playlist', '--playlists', '-p', '--start', '--end', '--njobs',
           '--count', '--tarball', '--mc_only', '--mc-only', '--data_only',
           '--data-only'}
WITH_VALUE = {'--playlist', '--playlists', '-p', '--start', '--end', '--njobs',
              '--count', '--tarball'}


def read_submission(path):
    lines = path.read_text(errors='replace').splitlines()
    tar = next((m.group(1) for line in lines if (m := TARBALL.match(line))), None)
    records = []
    command = None
    current_record = None
    for index, line in enumerate(lines):
        if RUN.match(line):
            command = None
            for candidate in lines[index+1:index+5]:
                candidate = candidate.strip()
                if candidate.startswith('python ') or candidate.startswith('python3 '):
                    command = candidate
                    break
        elif (m := OUTPUT.match(line)):
            if command:
                argv = shlex.split(command)
                # The same log can contain MC, data, and special MC submissions.
                current_record = {'dir': m.group(1).rstrip('/'), 'argv': argv, 'njobs': None}
                records.append(current_record)
                command = None
        elif (m := JOBSUB.match(line)) and current_record is not None:
            current_record['njobs'] = int(m.group(1))
            current_record = None
    return tar, records


def get_opt(argv, names, default=None):
    for i, arg in enumerate(argv):
        if arg in names and i + 1 < len(argv):
            return argv[i+1]
    return default


def clean_original(argv):
    # Starting from the script arguments, keep all original analysis switches.
    try:
        i = next(i for i, x in enumerate(argv) if x.endswith('gridSelection.py')) + 1
    except StopIteration:
        raise ValueError('Submission command does not invoke gridSelection.py')
    original = argv[i:]
    result = []
    j = 0
    while j < len(original):
        flag = original[j]
        if flag in REPLACE:
            j += 2 if flag in WITH_VALUE else 1
        else:
            result.append(flag)
            j += 1
    return result


def find_option_value(argv, name, default):
    return int(get_opt(argv, [name], str(default)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('combine_log', type=Path, help='Log_<stamp>.txt (combining diagnostic log)')
    p.add_argument('--submission-log', type=Path, help='Original submission log (auto-located otherwise)')
    p.add_argument('--exclude-universes', nargs='+', metavar='NAME', help='Override original excluded universe names, e.g. all')
    p.add_argument('--sideband', nargs='+', metavar='NAME', help='Override original sidebands; use None to force none')
    p.add_argument('--tarball', help='Override the exact tarball in the submission log')
    p.add_argument('--extra', default='', help='Extra arguments appended to each command')
    args = p.parse_args()

    combine_lines = args.combine_log.read_text(errors='replace').splitlines()
    stamp = next((m.group(1) for line in combine_lines if (m := STAMP.match(line))), None)
    tag = next((m.group(1) for line in combine_lines if (m := TAG.match(line))), None)
    if not stamp or not tag:
        p.error('Combine log missing runningNotes stamp or selection tag header')
    submit_path = args.submission_log or args.combine_log.with_name(stamp + '.txt')
    if not submit_path.is_file():
        p.error(f'Submission log not found: {submit_path}. Supply --submission-log PATH')

    tarball, records = read_submission(submit_path)
    tarball = args.tarball or tarball
    if not tarball:
        p.error('Submission log lacks "Using shared tarball:"; specify --tarball')
    if not records:
        p.error('No submission commands paired with [GRID_OUTPUT_HISTS] found')
    by_dir = {r['dir']: r for r in records}
    failures = defaultdict(set)
    broken = defaultdict(set)
    for line in combine_lines:
        m = SKIP.match(line)
        if not m:
            continue
        path = Path(m.group(1))
        fm = FILE.match(path.name)
        if not fm:
            print(f'# WARNING: unrecognized skipped file: {path}', file=sys.stderr)
            continue
        kind, body, ntuple, start = fm.groups()
        suffix = '_' + tag
        if not body.endswith(suffix):
            print(f'# WARNING: selection mismatch: {path}', file=sys.stderr)
            continue
        playlist = body[:-len(suffix)]
        record = by_dir.get(str(path.parent))
        if not record:
            print(f'# WARNING: cannot find submission for directory {path.parent}; skipping', file=sys.stderr)
            continue
        original = record['argv']
        original_playlist = get_opt(original, ['--playlists', '--playlist', '-p'])
        original_kind = 'mc' if any(x in original for x in ('--mc_only', '--mc-only')) else 'data' if any(x in original for x in ('--data_only', '--data-only')) else None
        original_ntuple = get_opt(original, ['--ntuple_tag', '--ntuple-tag'])
        if original_playlist != playlist or original_kind != kind or original_ntuple != ntuple:
            print(f'# WARNING: submission/output mismatch for {path}; skipping', file=sys.stderr)
            continue
        key = (str(path.parent), playlist, kind, ntuple)
        failures[key].add(int(start))
        broken[key].add(int(start))

    # Actual job count is taken from the -N argument in [GRID_JOBSUB_CMD].
    missing = defaultdict(set)
    for record in records:
        hist_dir = record['dir']
        original = record['argv']
        playlist = get_opt(original, ['--playlists', '--playlist', '-p'])
        ntuple = get_opt(original, ['--ntuple_tag', '--ntuple-tag'])
        kind = ('mc' if any(x in original for x in ('--mc_only', '--mc-only'))
                else 'data' if any(x in original for x in ('--data_only', '--data-only'))
                else None)
        if not (playlist and ntuple and kind):
            print(f'# WARNING: cannot identify playlist/type for {hist_dir}', file=sys.stderr)
            continue
        n_jobs = record['njobs']
        if n_jobs is None:
            print(f'# WARNING: missing [GRID_JOBSUB_CMD] -N for {hist_dir}; cannot check missing files', file=sys.stderr)
            continue
        count = find_option_value(original, '--count', 1 if kind == 'mc' else 1000)
        first = find_option_value(original, '--start', 0)
        key = (hist_dir, playlist, kind, ntuple)
        directory = Path(hist_dir)
        if not directory.is_dir():
            p.error(f'Histogram directory inaccessible: {directory}. Run this on GPVM with PNFS mounted.')
        observed = set()
        prefix = f'kin_dist_{kind}{playlist}_{tag}_{ntuple}_'
        for entry in directory.iterdir():
            if entry.name.startswith(prefix) and entry.name.endswith('.root'):
                offset = entry.name[len(prefix):-5]
                if offset.isdigit():
                    observed.add(int(offset))
        expected = {first + j * count for j in range(n_jobs)}
        missing[key] = expected - observed
        failures[key].update(missing[key])

    print('#!/bin/bash')
    print('set -e')
    print(f'# Combining log: {args.combine_log}')
    print(f'# Submission log: {submit_path}')
    print(f'# Exact original tarball: {tarball}')
    print('# Commands must be run in CC-NuE-XSec/selection on GPVM.')
    print('# New grid outputs must be validated before replacing failed originals.')
    print(f'ls -lh {shlex.quote(tarball)}')
    total = 0
    for (hist_dir, playlist, kind, ntuple), starts in sorted(failures.items()):
        original = by_dir[hist_dir]['argv']
        count = find_option_value(original, '--count', 1 if kind == 'mc' else 1000)
        original_start = find_option_value(original, '--start', 0)
        original_selection = get_opt(original, ['--selection_tag', '--selection-tag'])
        if original_selection != tag:
            print(f'# WARNING: selection mismatch for {playlist}; skipping', file=sys.stderr)
            continue
        if count <= 0:
            print(f'# WARNING: invalid count for {playlist}', file=sys.stderr)
            continue
        print(f'\n# {playlist} {kind}: missing {sorted(missing[(hist_dir, playlist, kind, ntuple)])}; ' +
              f'broken {sorted(broken[(hist_dir, playlist, kind, ntuple)])}; ' +
              f'original count={count}; original start={original_start}')
        for start in sorted(starts):
            if start < original_start or (start - original_start) % count:
                print(f'# WARNING: unexpected job offset {start} for {playlist}; skipping', file=sys.stderr)
                continue
            cmd = ['python', 'selection/gridSelection.py', '--playlists', playlist]
            retained = clean_original(original)
            # Only override analysis options explicitly requested by user.
            if args.exclude_universes is not None:
                retained = remove_multivalue(retained, {'--exclude-universes', '--exclude_universes'})
                retained += ['--exclude-universes'] + args.exclude_universes
            if args.sideband is not None:
                retained = remove_multivalue(retained, {'--use-sideband', '--use_sideband'})
                retained += ['--use-sideband'] + ([] if args.sideband == ['None'] else args.sideband)
            cmd += retained
            cmd += ['--mc_only' if kind == 'mc' else '--data_only', '--start', str(start),
                    '--count', str(count), '--njobs', '1', '--tarball', tarball]
            if args.extra:
                cmd += shlex.split(args.extra)
            print(shlex.join(cmd))
            total += 1
    print(f'\n# Total: {total} rerun commands')
    if not total:
        print('# No missing outputs or recognized broken files found')


def remove_multivalue(argv, flags):
    """Remove nargs=* flags and their non-option values, leaving all other args."""
    out = []
    i = 0
    while i < len(argv):
        if argv[i] in flags:
            i += 1
            while i < len(argv) and not argv[i].startswith('-'):
                i += 1
        else:
            out.append(argv[i])
            i += 1
    return out


if __name__ == '__main__':
    main()
