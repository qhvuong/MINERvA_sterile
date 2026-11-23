#!/usr/bin/env python3
"""Read-only ME POT audit. Run with PyROOT on GPVM; no analysis config imports.

Reads official PNFS text lists directly, checks every source Meta tree, compares
POT with per-playlist/final selection files, and compares EN4 combination methods.
Never edits ROOT inputs, configs/POT.json, or existing selection outputs.
"""
import argparse
import csv
import hashlib
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
from collections import Counter

PLAYLISTS = {
    'FHC': ['1A','1B','1C','1D','1E','1F','1G','1L','1M','1N','1O','1P'],
    'RHC': ['5A','6A','6B','6C','6D','6E','6F','6G','6H','6I','6J'],
}


def read_root(task):
    path, include_hist = task
    result = {'path': path, 'ok': False}
    f = None
    try:
        import ROOT
        ROOT.PyConfig.IgnoreCommandLineOptions = True
        ROOT.gROOT.SetBatch(True)
        # Load dictionaries only for selection histograms.
        if include_hist:
            import PlotUtils  # noqa: F401
        f = ROOT.TFile.Open(path, 'READ')
        if not f or f.IsZombie():
            raise RuntimeError('Cannot open ROOT file')
        if f.TestBit(ROOT.TFile.kRecovered):
            raise RuntimeError('Recovered/incomplete ROOT file')
        tree = f.Get('Meta')
        if not tree or not tree.GetBranch('POT_Used'):
            raise RuntimeError('Missing Meta tree or POT_Used branch')
        tree.SetBranchStatus('*', 0)
        tree.SetBranchStatus('POT_Used', 1)
        values = []
        for i in range(tree.GetEntries()):
            if tree.GetEntry(i) <= 0:
                raise RuntimeError(f'Cannot read Meta entry {i}')
            value = float(tree.GetLeaf('POT_Used').GetValue())
            if not math.isfinite(value) or value < 0:
                raise RuntimeError(f'Invalid POT_Used at entry {i}: {value}')
            values.append(value)
        if not values:
            raise RuntimeError('Empty Meta tree')
        result.update(ok=True, pot=math.fsum(values), meta_entries=len(values))
        if include_hist:
            h = f.Get('EN4')
            if h:
                result['edges'] = [h.GetBinLowEdge(i) for i in range(1, h.GetNbinsX()+2)]
                result['counts'] = [h.GetBinContent(i) for i in range(1, h.GetNbinsX()+1)]
            else:
                result['hist_warning'] = 'EN4 histogram missing; POT check still available'
    except Exception as error:
        result.update(ok=False, error=str(error))
    finally:
        if f:
            f.Close()
    return result


def listing_path(directory, playlist, kind):
    wanted = f'P8_Playlist{playlist}_{kind}.txt'.lower()
    matches = [p for p in directory.iterdir() if p.name.lower() == wanted]
    if len(matches) != 1:
        raise RuntimeError(f'Expected one {wanted} in {directory}, found {len(matches)}')
    return matches[0]


def load_list(path):
    entries = []
    for line in path.read_text().splitlines():
        line = line.split('#', 1)[0].strip()
        if not line:
            continue
        if len(line.split()) != 1 or '.root' not in line:
            raise RuntimeError(f'Unrecognized list entry in {path}: {line!r}')
        entries.append(line)
    if not entries:
        raise RuntimeError(f'Empty input list: {path}')
    return entries


def close(a, b):
    return math.isclose(a, b, rel_tol=1e-8, abs_tol=1e-9)


def predict(rows):
    """Compare (sum data POT / sum MC POT)*sum MC with sum (D_j/M_j)*MC_j."""
    pairs = [(r['data']['selection'], r['mc']['selection']) for r in rows]
    if not all(d.get('ok') and m.get('ok') and 'counts' in d and 'counts' in m
               and d['pot'] > 0 and m['pot'] > 0 for d, m in pairs):
        return {'status': 'UNAVAILABLE', 'reason': 'Missing/invalid selection POT or EN4'}
    edges = pairs[0][0]['edges']
    if any(d['edges'] != edges or m['edges'] != edges for d, m in pairs):
        return {'status': 'UNAVAILABLE', 'reason': 'Inconsistent EN4 bin edges'}
    scale = math.fsum(d['pot'] for d, _ in pairs) / math.fsum(m['pot'] for _, m in pairs)
    bins = []
    for i in range(len(edges)-1):
        data = math.fsum(d['counts'][i] for d, _ in pairs)
        global_mc = scale * math.fsum(m['counts'][i] for _, m in pairs)
        per_playlist_mc = math.fsum(d['pot']/m['pot']*m['counts'][i] for d, m in pairs)
        bins.append(dict(low=edges[i], high=edges[i+1], data=data,
                         global_mc=global_mc, per_playlist_mc=per_playlist_mc,
                         per_playlist_over_global=per_playlist_mc/global_mc if global_mc else None,
                         data_over_global=data/global_mc if global_mc else None,
                         data_over_per_playlist=data/per_playlist_mc if per_playlist_mc else None))
    return dict(status='COMPUTED', global_pot_scale=scale, bins=bins)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--beam', choices=['FHC','RHC','both'], default='both')
    p.add_argument('--pnfs-base', type=Path, default=Path('/pnfs/minerva/persistent/DataPreservation/p8'))
    p.add_argument('--analysis-root', type=Path, default=Path.cwd())
    p.add_argument('--selection-base', type=Path, default=Path('/exp/minerva/data/users/qvuong'))
    p.add_argument('--fhc-tag', default='ME_CCnumu_FullCuts_sameOffset_v4p3p1')
    p.add_argument('--rhc-tag', default='ME_CCnumubar_FullCuts_sameOffset_v4p3p1')
    p.add_argument('--workers', type=int, default=4)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--list-only', action='store_true', help='Validate list discovery and membership without opening ROOT files')
    a = p.parse_args()
    if a.workers < 1:
        p.error('--workers must be positive')
    a.output_dir.mkdir(parents=True, exist_ok=True)
    beams = ['FHC','RHC'] if a.beam == 'both' else [a.beam]
    pot_config_path = a.analysis_root/'configs/POT.json'
    pot_config = json.loads(pot_config_path.read_text()) if pot_config_path.is_file() else {}
    manifest, unique_sources = {}, set()
    problems = []
    for beam in beams:
        manifest[beam] = []
        for pl in PLAYLISTS[beam]:
            row = {'playlist': 'me'+pl}
            for kind, folder in [('data','Data'), ('mc','FullDetector')]:
                txt = listing_path(a.pnfs_base/folder, pl, 'Data' if kind == 'data' else 'MC')
                entries = load_list(txt)
                duplicate_count = len(entries)-len(set(entries))
                if duplicate_count:
                    problems.append(f'{beam} me{pl} {kind}: {duplicate_count} duplicate list entries')
                item = dict(list=str(txt), sha256=hashlib.sha256(txt.read_bytes()).hexdigest(),
                            files=entries, duplicate_entries=duplicate_count)
                # Read existing mapping without importing AnalysisConfig or changing POT.json.
                mapping = pot_config.get('me'+pl+'_p8', {}).get(kind, {}).get('MAD', {})
                local = mapping.get('playlist_location')
                if local:
                    local = Path(local)
                    if not local.is_absolute():
                        local = a.analysis_root/local
                    item['configured_list'] = str(local)
                    if local.is_file():
                        actual = load_list(local)
                        item['configured_files_absent_from_pnfs_list'] = sorted(set(actual)-set(entries))
                        item['pnfs_files_absent_from_configured_list'] = sorted(set(entries)-set(actual))
                        item['configured_duplicate_entries'] = len(actual)-len(set(actual))
                        item['configured_same_order'] = actual == entries
                        if Counter(actual) != Counter(entries):
                            problems.append(f'{beam} me{pl} {kind}: configured list differs from PNFS list')
                    else:
                        item['configured_list_warning'] = 'Configured list not accessible'
                tag = a.fhc_tag if beam == 'FHC' else a.rhc_tag
                subdir = 'nu_mu' if beam == 'FHC' else 'antinu_mu'
                item['selection_path'] = str(a.selection_base/subdir/f'kin_dist_{kind}me{pl}_{tag}_MAD.root')
                row[kind] = item
                unique_sources.update(entries)
                print(f'{beam} me{pl:<3} {kind:<4}: {len(entries):5d} PNFS files; duplicates={duplicate_count}', flush=True)
            manifest[beam].append(row)
    # A shared file between playlists would count the same source exposure twice.
    owners = {}
    for beam, rows in manifest.items():
        for row in rows:
            for kind in ['data','mc']:
                for path in set(row[kind]['files']):
                    owners.setdefault(path, []).append(f'{beam}/{row["playlist"]}/{kind}')
    overlaps = {path: names for path, names in owners.items() if len(names) > 1}
    if overlaps:
        problems.append(f'{len(overlaps)} source files shared across playlist/type lists')
    report = dict(manifest=manifest, source_overlaps=overlaps, problems=problems,
                  note='POT agreement checks exposure bookkeeping, not full event-loop completion.')
    destination = a.output_dir/f'me_pot_audit_{a.beam}.json'
    if a.list_only:
        report['status'] = 'LIST_CHECK_ONLY'
        destination.write_text(json.dumps(report, indent=2, allow_nan=False))
        print(f'List-only report: {destination}', flush=True)
        for problem in problems:
            print('[WARN] '+problem, flush=True)
        return
    # Save the manifest before the expensive reads so failures leave a useful record.
    destination.write_text(json.dumps(report, indent=2, allow_nan=False))
    source_results = {}
    with (a.output_dir/f'source_reads_{a.beam}.jsonl').open('w') as log:
        with mp.get_context('spawn').Pool(a.workers) as pool:
            tasks = [(path, False) for path in sorted(unique_sources)]
            for i, result in enumerate(pool.imap_unordered(read_root, tasks, chunksize=1), 1):
                source_results[result['path']] = result
                log.write(json.dumps(result, allow_nan=False)+'\n')
                log.flush()
                if not result['ok']:
                    print('[ERROR] '+result['path']+': '+result['error'], flush=True)
                if i % 100 == 0 or i == len(tasks):
                    print(f'Source metadata: {i}/{len(tasks)} files checked', flush=True)
    csv_rows = []
    for beam, rows in manifest.items():
        for row in rows:
            for kind in ['data','mc']:
                item = row[kind]
                bad = [path for path in item['files'] if not source_results[path]['ok']]
                item['bad_source_files'] = bad
                item['source_pot_readable_files_only'] = math.fsum(source_results[path]['pot'] for path in item['files'] if source_results[path]['ok'])
                item['source_meta_entries_readable_files_only'] = sum(source_results[path]['meta_entries'] for path in item['files'] if source_results[path]['ok'])
                item['selection'] = selected = read_root((item['selection_path'], True))
                ratio = selected['pot']/item['source_pot_readable_files_only'] if selected['ok'] and item['source_pot_readable_files_only'] > 0 else None
                matched = (not bad and not item['duplicate_entries'] and selected['ok']
                           and close(selected['pot'], item['source_pot_readable_files_only'])
                           and selected['meta_entries'] == item['source_meta_entries_readable_files_only'])
                item['pot_status'] = 'MATCH' if matched else 'CHECK'
                if not matched:
                    problems.append(f'{beam} {row["playlist"]} {kind}: source/selection POT or Meta entries do not match, or source audit incomplete')
                csv_rows.append(dict(beam=beam, playlist=row['playlist'], kind=kind,
                                     source_files=len(item['files']), bad_source_files=len(bad),
                                     source_pot=item['source_pot_readable_files_only'],
                                     selection_pot=selected.get('pot'), selection_over_source=ratio,
                                     status=item['pot_status']))
                print(f'{beam} {row["playlist"]:<5} {kind:<4}: selection/source={ratio if ratio is not None else "unavailable"} [{item["pot_status"]}]', flush=True)
        tag = a.fhc_tag if beam == 'FHC' else a.rhc_tag
        subdir = 'nu_mu' if beam == 'FHC' else 'antinu_mu'
        final = {}
        for kind in ['data','mc']:
            path = a.selection_base/subdir/f'kin_dist_{kind}me{beam}_{tag}_MAD.root'
            result = read_root((str(path), True))
            selections = [row[kind]['selection'] for row in rows]
            ready = all(x['ok'] for x in selections)
            expected = math.fsum(x['pot'] for x in selections) if ready else None
            result['playlist_pot_sum'] = expected
            result['matches_playlist_sum'] = bool(result['ok'] and ready and close(result['pot'], expected)
                and result['meta_entries'] == sum(x['meta_entries'] for x in selections))
            if not result['matches_playlist_sum']:
                problems.append(f'{beam} {kind}: final merged POT/Meta does not match per-playlist sum')
            if result['ok'] and 'counts' in result and all('counts' in x and x['edges'] == result['edges'] for x in selections):
                summed = [math.fsum(x['counts'][i] for x in selections) for i in range(len(result['counts']))]
                result['en4_matches_playlist_sum'] = all(close(x,y) for x,y in zip(result['counts'], summed))
                if not result['en4_matches_playlist_sum']:
                    problems.append(f'{beam} {kind}: final EN4 does not match per-playlist sum')
            final[kind] = result
            print(f'{beam} final {kind}: POT={result.get("pot")} matches_playlist_sum={result["matches_playlist_sum"]}', flush=True)
        report[beam] = dict(final=final, combination_comparison=predict(rows))
        comparison = report[beam]['combination_comparison']
        if comparison['status'] == 'COMPUTED':
            with (a.output_dir/f'en4_combination_{beam}.csv').open('w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=list(comparison['bins'][0]))
                writer.writeheader(); writer.writerows(comparison['bins'])
            print(f'{beam}: per-playlist/global MC in 2-8 GeV bins:', flush=True)
            for b in comparison['bins']:
                if b['low'] < 8 and b['high'] > 2:
                    print(f'  {b["low"]:g}-{b["high"]:g}: {b["per_playlist_over_global"]}', flush=True)
    with (a.output_dir/f'pot_summary_{a.beam}.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(csv_rows[0]))
        writer.writeheader(); writer.writerows(csv_rows)
    report['status'] = 'CHECK' if problems else 'POT_BOOKKEEPING_MATCH'
    destination.write_text(json.dumps(report, indent=2, allow_nan=False))
    print(f'\nStatus: {report["status"]}; report: {destination}', flush=True)
    for problem in problems:
        print('[WARN] '+problem, flush=True)


if __name__ == '__main__':
    main()
