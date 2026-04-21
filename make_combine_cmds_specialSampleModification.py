#!/usr/bin/env python3
import os
import re
import sys
import argparse
import subprocess

argv_backup = sys.argv[:]
sys.argv = [sys.argv[0]]
from config.AnalysisConfig import AnalysisConfig
sys.argv = argv_backup

from datetime import datetime
from typing import Dict, Optional, List, Tuple


RE_RUN_CMD = re.compile(r'^\s*Running command:\s*$')
RE_SEL_CMD = re.compile(r'^\s*python\s+selection/gridSelection\.py\b(.*)$')
RE_PLAYLIST = re.compile(r'--playlist\s+(\S+)')
RE_TAG = re.compile(r'--selection_tag\s+(\S+)')
RE_TAR = re.compile(r'qvuong-CCNUE_selection_(\d{4}-\d{2}-\d{2}-\d{6})\.tar\.gz')
RE_JOBS = re.compile(r'^\s*(\d+)\s+job\(s\)\s+submitted\b')
RE_NOTES_STAMP = re.compile(r'(.+_\d{4}-\d{2}-\d{2}_\d{6})\.txt$')

# Detect which submission type we're in by wrapper name
RE_WRAPPER_TYPE = re.compile(r'CCNuE-[^-]+\-(mc|data)_wrapper\.sh\b')

# Heuristic token match for ROOT output type
RE_MC_TOKEN = re.compile(r'(^|[._-])mc([._-]|$)', re.IGNORECASE)
RE_DATA_TOKEN = re.compile(r'(^|[._-])data([._-]|$)', re.IGNORECASE)

def strip_p6(s: str) -> str:
    return s[:-3] if s.endswith("_p6") else s

def hist_dir_from_stamp(stamp: str) -> str:
    return f"/pnfs/minerva/persistent/users/qvuong/CCNUE_selection_{stamp}_hists/"

def count_root_files_for_playlist(path: str, playlist_raw: str) -> Dict[str, int]:
    counts = {"mc": 0, "data": 0, "other": 0, "ambiguous": 0, "total": 0}

    mc_prefix = f"kin_dist_mc{playlist_raw}_"
    data_prefix = f"kin_dist_data{playlist_raw}_"

    try:
        for fn in os.listdir(path):
            if not fn.endswith(".root"):
                continue

            is_mc = fn.startswith(mc_prefix)
            is_data = fn.startswith(data_prefix)

            if is_mc or is_data:
                counts["total"] += 1

            if is_mc and not is_data:
                counts["mc"] += 1
            elif is_data and not is_mc:
                counts["data"] += 1
            elif is_mc and is_data:
                counts["ambiguous"] += 1
            else:
                counts["other"] += 1
    except FileNotFoundError:
        pass

    return counts

def format_block_log(block: dict, file_counts: Dict[str, int]) -> List[str]:
    playlist_raw = block["playlist"]
    tag = block["selection_tag"]
    stamp = block["stamp"]
    exp_mc = block["expected_mc_jobs"]
    exp_data = block["expected_data_jobs"]
    exp_total = exp_mc + exp_data

    playlist = strip_p6(playlist_raw) if playlist_raw else None
    in_dir = f"/pnfs/minerva/persistent/users/qvuong/CCNUE_selection_{stamp}_hists/" if stamp else "(unknown)"

    seen_mc = file_counts["mc"]
    seen_data = file_counts["data"]
    seen_total = file_counts["total"]

    lines = []
    lines.append("==================================================")
    lines.append(f"Playlist           : {playlist} (raw: {playlist_raw})")
    lines.append(f"Selection tag      : {tag}")
    lines.append(f"Stamp              : {stamp}")
    lines.append(f"Input dir          : {in_dir}")
    lines.append(f"Expected MC jobs   : {exp_mc}")
    lines.append(f"Expected DATA jobs : {exp_data}")
    lines.append(f"Expected TOTAL     : {exp_total}")
    lines.append(f"Seen MC ROOT files : {seen_mc}")
    lines.append(f"Seen DATA ROOT files: {seen_data}")
    lines.append(f"Seen TOTAL ROOT files: {seen_total}")
    if file_counts["ambiguous"] or file_counts["other"]:
        lines.append(f"Seen ambiguous ROOT files: {file_counts['ambiguous']}")
        lines.append(f"Seen other ROOT files    : {file_counts['other']}")

    # Warnings (per-type and total)
    if seen_mc != exp_mc:
        lines.append(f"[WARN] MC mismatch (expected {exp_mc}, saw {seen_mc})")
    else:
        lines.append("[OK] MC count matches")

    if seen_data != exp_data:
        lines.append(f"[WARN] DATA mismatch (expected {exp_data}, saw {seen_data})")
    else:
        lines.append("[OK] DATA count matches")

    if seen_total != exp_total:
        lines.append(f"[WARN] TOTAL mismatch (expected {exp_total}, saw {seen_total})")
    else:
        lines.append("[OK] TOTAL count matches")

    lines.append("==================================================")
    lines.append("")
    return lines

def block_ok(block: dict, file_counts: Dict[str, int]) -> Tuple[bool, bool]:
    ok_mc = (file_counts["mc"] == block["expected_mc_jobs"])
    ok_data = (file_counts["data"] == block["expected_data_jobs"])
    return ok_mc, ok_data

def handle_block(block: dict, dry_run: bool, log_sink: Optional[List[str]] = None, special_sample_dir: Optional[str] = None):
    if not block:
        return

    playlist_raw = block["playlist"]
    tag = block["selection_tag"]
    stamp = block["stamp"]

    if not playlist_raw or not tag or not stamp:
        msg = f"[WARN] Incomplete block, skipping: {block}"
        print(msg, file=sys.stderr)
        if log_sink is not None:
            log_sink.append(msg)
            log_sink.append("")
        return

    playlist = strip_p6(playlist_raw)
    in_dir = f"/pnfs/minerva/persistent/users/qvuong/CCNUE_selection_{stamp}_hists/"
    file_counts = count_root_files_for_playlist(in_dir, playlist_raw)
    ok_mc, ok_data = block_ok(block, file_counts)
    # should_run = ok_mc and ok_data
    should_run = True

    # Console output (similar to your old behavior, but richer)
    print("\n".join(format_block_log(block, file_counts)))

    # Log output
    if log_sink is not None:
        log_sink.extend(format_block_log(block, file_counts))

    cmd = [
        "python", "combine_file_specialSampleModification.py",
        "--playlist", playlist,
        "--i", in_dir,
        "--cal_POT",
        "--ntuple_tag", "MAD",
        "--selection_tag", tag,
    ]

    print("Command:")
    print(" ", " ".join(cmd))
    if log_sink is not None:
        log_sink.append("Command:")
        log_sink.append("  " + " ".join(cmd))
        log_sink.append("")

    if dry_run:
        if should_run:
            print("[DRY-RUN] would execute combine_file_specialSampleModification.py (counts OK for MC and DATA)")
            if log_sink is not None:
                log_sink.append("[DRY-RUN] would execute combine_file_specialSampleModification.py (counts OK for MC and DATA)")
                log_sink.append("")
        else:
            print("[DRY-RUN] would SKIP combine_file_specialSampleModification.py (count mismatch)")
            if log_sink is not None:
                log_sink.append("[DRY-RUN] would SKIP combine_file_specialSampleModification.py (count mismatch)")
                log_sink.append("")
    else:
        if should_run:
            print("[RUNNING] counts OK for MC and DATA")
            if log_sink is not None:
                log_sink.append("[RUNNING] counts OK for MC and DATA")
                log_sink.append("")

            if special_sample_dir:
                special_input = special_sample_dir + "\n\n"
                print(f"[INFO] Supplying special MC sample dir: {special_sample_dir}")
                if log_sink is not None:
                    log_sink.append(f"[INFO] Supplying special MC sample dir: {special_sample_dir}")
                    log_sink.append("")
            else:
                special_input = "\n"

            subprocess.run(cmd, input=special_input, text=True)
        else:
            print("[SKIP] Not running combine_file_specialSampleModification.py because counts did not match for MC and/or DATA")
            if log_sink is not None:
                log_sink.append("[SKIP] Not running combine_file_specialSampleModification.py because counts did not match for MC and/or DATA")
                log_sink.append("")

EXPECTED_MC_PLAYLISTS = ["le1", "le7", "le9", "le13C"]
EXPECTED_DATA_PLAYLISTS = ["le1", "le7", "le9", "le13A", "le13B", "le13C", "le13D", "le13E"]

def combined_output_path(playlist: str, selection_tag: str, is_data: bool, ntuple_tag: str = "MAD") -> str:
    old_playlist = getattr(AnalysisConfig, "playlist", None)
    old_selection_tag = getattr(AnalysisConfig, "selection_tag", None)
    old_ntuple_tag = getattr(AnalysisConfig, "ntuple_tag", None)

    try:
        AnalysisConfig.playlist = playlist
        AnalysisConfig.selection_tag = selection_tag
        AnalysisConfig.ntuple_tag = ntuple_tag
        return AnalysisConfig.SelectionHistoPath(playlist, is_data)
    finally:
        AnalysisConfig.playlist = old_playlist
        AnalysisConfig.selection_tag = old_selection_tag
        AnalysisConfig.ntuple_tag = old_ntuple_tag


def maybe_hadd_fhc(selection_tag: str, dry_run: bool, log_sink=None):
    mc_files = [
        combined_output_path(pl, selection_tag, False)
        for pl in EXPECTED_MC_PLAYLISTS
    ]
    data_files = [
        combined_output_path(pl, selection_tag, True)
        for pl in EXPECTED_DATA_PLAYLISTS
    ]

    # keep only the standard combined MC outputs
    mc_files = [
        f for f in mc_files
        if "_no2p2h" not in os.path.basename(f) and "_only2p2h" not in os.path.basename(f)
    ]

    mc_dir = os.path.dirname(mc_files[0]) if mc_files else "."
    data_dir = os.path.dirname(data_files[0]) if data_files else "."

    out_mc = os.path.join(mc_dir, f"kin_dist_mcleFHC_{selection_tag}_MAD.root")
    out_data = os.path.join(data_dir, f"kin_dist_dataleFHC_{selection_tag}_MAD.root")

    missing_mc = [f for f in mc_files if not os.path.exists(f)]
    missing_data = [f for f in data_files if not os.path.exists(f)]

    lines = []
    lines.append("################ FHC MERGE CHECK ################")
    lines.append(f"Selection tag: {selection_tag}")
    lines.append(f"MC output dir : {mc_dir}")
    lines.append(f"DATA output dir: {data_dir}")
    lines.append(f"Found MC      : {len(mc_files) - len(missing_mc)}/{len(mc_files)}")
    lines.append(f"Found DATA    : {len(data_files) - len(missing_data)}/{len(data_files)}")

    if missing_mc:
        lines.append("[WARN] Missing MC files:")
        lines.extend([f"  {x}" for x in missing_mc])
    else:
        lines.append("[OK] All MC playlist outputs are present")

    if missing_data:
        lines.append("[WARN] Missing DATA files:")
        lines.extend([f"  {x}" for x in missing_data])
    else:
        lines.append("[OK] All DATA playlist outputs are present")

    print("\n".join(lines))
    if log_sink is not None:
        log_sink.extend(lines)
        log_sink.append("")

    cmd_mc = ["hadd", "-f", out_mc] + mc_files
    cmd_data = ["hadd", "-f", out_data] + data_files

    if not missing_mc:
        print("FHC MC hadd command:")
        print(" ", " ".join(cmd_mc))
        if log_sink is not None:
            log_sink.append("FHC MC hadd command:")
            log_sink.append("  " + " ".join(cmd_mc))
            log_sink.append("")

        if dry_run:
            print("[DRY-RUN] would run FHC MC hadd")
            if log_sink is not None:
                log_sink.append("[DRY-RUN] would run FHC MC hadd")
                log_sink.append("")
        else:
            print("[RUNNING] hadd for FHC MC")
            subprocess.run(cmd_mc, check=False)
    else:
        print("[SKIP] Not running MC hadd because some expected MC playlist outputs are missing")
        if log_sink is not None:
            log_sink.append("[SKIP] Not running MC hadd because some expected MC playlist outputs are missing")
            log_sink.append("")

    if not missing_data:
        print("FHC DATA hadd command:")
        print(" ", " ".join(cmd_data))
        if log_sink is not None:
            log_sink.append("FHC DATA hadd command:")
            log_sink.append("  " + " ".join(cmd_data))
            log_sink.append("")

        if dry_run:
            print("[DRY-RUN] would run FHC DATA hadd")
            if log_sink is not None:
                log_sink.append("[DRY-RUN] would run FHC DATA hadd")
                log_sink.append("")
        else:
            print("[RUNNING] hadd for FHC DATA")
            subprocess.run(cmd_data, check=False)
    else:
        print("[SKIP] Not running DATA hadd because some expected DATA playlist outputs are missing")
        if log_sink is not None:
            log_sink.append("[SKIP] Not running DATA hadd because some expected DATA playlist outputs are missing")
            log_sink.append("")

# def maybe_hadd_fhc(selection_tag: str, dry_run: bool, log_sink=None):
#     mc_files = [
#         combined_output_path(pl, selection_tag, False)
#         for pl in EXPECTED_MC_PLAYLISTS
#     ]
#     data_files = [
#         combined_output_path(pl, selection_tag, True)
#         for pl in EXPECTED_DATA_PLAYLISTS
#     ]

#     # keep only the standard combined MC outputs
#     mc_files = [
#         f for f in mc_files
#         if "_no2p2h" not in os.path.basename(f) and "_only2p2h" not in os.path.basename(f)
#     ]

#     mc_dir = os.path.dirname(mc_files[0]) if mc_files else "."
#     data_dir = os.path.dirname(data_files[0]) if data_files else "."

#     out_mc = os.path.join(mc_dir, f"kin_dist_mcleFHC_{selection_tag}_MAD.root")
#     out_data = os.path.join(data_dir, f"kin_dist_dataleFHC_{selection_tag}_MAD.root")

#     missing_mc = [f for f in mc_files if not os.path.exists(f)]
#     missing_data = [f for f in data_files if not os.path.exists(f)]

#     lines = []
#     lines.append("################ FHC MERGE CHECK ################")
#     lines.append(f"Selection tag: {selection_tag}")
#     lines.append(f"MC output dir : {mc_dir}")
#     lines.append(f"DATA output dir: {data_dir}")
#     lines.append(f"Found MC      : {len(mc_files) - len(missing_mc)}/{len(mc_files)}")
#     lines.append(f"Found DATA    : {len(data_files) - len(missing_data)}/{len(data_files)}")

#     if missing_mc:
#         lines.append("[WARN] Missing MC files:")
#         lines.extend([f"  {x}" for x in missing_mc])
#     else:
#         lines.append("[OK] All MC playlist outputs are present")

#     if missing_data:
#         lines.append("[WARN] Missing DATA files:")
#         lines.extend([f"  {x}" for x in missing_data])
#     else:
#         lines.append("[OK] All DATA playlist outputs are present")

#     print("\n".join(lines))
#     if log_sink is not None:
#         log_sink.extend(lines)
#         log_sink.append("")

#     if missing_mc or missing_data:
#         msg = "[SKIP] Not running hadd because not all expected merged playlist outputs exist"
#         print(msg)
#         if log_sink is not None:
#             log_sink.append(msg)
#             log_sink.append("")
#         return

#     cmd_mc = ["hadd", "-f", out_mc] + mc_files
#     cmd_data = ["hadd", "-f", out_data] + data_files

#     print("FHC MC hadd command:")
#     print(" ", " ".join(cmd_mc))
#     print("FHC DATA hadd command:")
#     print(" ", " ".join(cmd_data))

#     if log_sink is not None:
#         log_sink.append("FHC MC hadd command:")
#         log_sink.append("  " + " ".join(cmd_mc))
#         log_sink.append("FHC DATA hadd command:")
#         log_sink.append("  " + " ".join(cmd_data))
#         log_sink.append("")

#     if dry_run:
#         print("[DRY-RUN] would run FHC MC/DATA hadd")
#         if log_sink is not None:
#             log_sink.append("[DRY-RUN] would run FHC MC/DATA hadd")
#             log_sink.append("")
#     else:
#         print("[RUNNING] hadd for FHC MC")
#         subprocess.run(cmd_mc, check=False)
#         print("[RUNNING] hadd for FHC DATA")
#         subprocess.run(cmd_data, check=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("runningNotes", help="runningNotes.txt file")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands only, do not execute combine_file_specialSampleModification.py",
    )
    parser.add_argument(
        "--max-blocks",
        type=int,
        default=None,
        help="Maximum number of submission blocks to process (default: all)",
    )
    args = parser.parse_args()

    all_blocks: List[dict] = []

    notes_basename = os.path.basename(args.runningNotes)

    m_notes = RE_NOTES_STAMP.match(notes_basename)
    notes_stamp = m_notes.group(1) if m_notes else "UNKNOWN"

    current = None
    waiting_for_cmd = False
    blocks_processed = 0

    # For per-block job-type parsing
    last_submit_type: Optional[str] = None  # 'mc' or 'data'

    # Build log lines in-memory, then write once at end
    log_lines: List[str] = []
    first_tag_seen: Optional[str] = None
    # log_timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    tags_seen = set()

    def start_new_block():
        nonlocal current, last_submit_type
        last_submit_type = None
        current = {
            "playlist": None,
            "selection_tag": None,
            "stamp": None,
            "expected_mc_jobs": 0,
            "expected_data_jobs": 0,
        }

    def finalize_block():
        nonlocal blocks_processed, first_tag_seen
        if current is None:
            return
        all_blocks.append(current.copy())

    with open(args.runningNotes, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            # Unified start-of-block: either marker or direct command
            m_cmd = RE_SEL_CMD.match(line)

            if RE_RUN_CMD.match(line) or m_cmd:
                # finalize previous block
                if current is not None:
                    finalize_block()
                    if args.max_blocks is not None and blocks_processed >= args.max_blocks:
                        print(f"[INFO] Reached --max-blocks={args.max_blocks}, stopping.")
                        break

                # start new block
                start_new_block()

                # If the command is on the same line, parse it immediately
                if m_cmd:
                    args_line = m_cmd.group(1)
                    mp = RE_PLAYLIST.search(args_line)
                    mt = RE_TAG.search(args_line)
                    current["playlist"] = mp.group(1) if mp else None
                    current["selection_tag"] = mt.group(1) if mt else None
                    waiting_for_cmd = False
                else:
                    # marker line: next line will contain the command
                    waiting_for_cmd = True
                continue

            if waiting_for_cmd:
                m = RE_SEL_CMD.match(line)
                if m:
                    waiting_for_cmd = False
                    args_line = m.group(1)
                    mp = RE_PLAYLIST.search(args_line)
                    mt = RE_TAG.search(args_line)
                    current["playlist"] = mp.group(1) if mp else None
                    current["selection_tag"] = mt.group(1) if mt else None
                continue

            # Accumulate stamp + job counts inside the current block
            if current is not None:
                mtar = RE_TAR.search(line)
                if mtar and current["stamp"] is None:
                    current["stamp"] = mtar.group(1)

                # Detect submission type from wrapper copy lines
                mw = RE_WRAPPER_TYPE.search(line)
                if mw:
                    last_submit_type = mw.group(1)  # 'mc' or 'data'

                mj = RE_JOBS.search(line)
                if mj:
                    n = int(mj.group(1))
                    # Attribute this job count to the most recently seen wrapper type
                    if last_submit_type == "mc":
                        current["expected_mc_jobs"] += n
                    elif last_submit_type == "data":
                        current["expected_data_jobs"] += n
                    else:
                        # Fallback: if we can't tell, treat as DATA (safer for data_only)
                        current["expected_data_jobs"] += n

    # Handle last block (if we didn't break early due to max-blocks)
    if current is not None and (args.max_blocks is None or blocks_processed < args.max_blocks):
        finalize_block()
                
    special_block_map = {}
    for block in all_blocks:
        pl = block.get("playlist")
        tag = block.get("selection_tag")
        stamp = block.get("stamp")
        if pl == "le13C_2p2hp6" and tag and stamp:
            special_block_map[("le13C_p6", tag)] = hist_dir_from_stamp(stamp)

    for block in all_blocks:
        playlist = block.get("playlist")
        tag = block.get("selection_tag")

        if playlist == "le13C_2p2hp6":
            if args.dry_run:
                print("[DRY-RUN] Showing diagnostics for special-sample block")
                log_lines.append("[DRY-RUN] Showing diagnostics for special-sample block")
                log_lines.append("")
                handle_block(block, True, log_lines, special_sample_dir=None)
            else:
                print("[SKIP] Special-sample submission block is used only as input to le13C_p6 combine step")
                log_lines.append("[SKIP] Special-sample submission block is used only as input to le13C_p6 combine step")
                log_lines.append("")
            continue

        special_dir = None
        if playlist == "le13C_p6" and tag:
            special_dir = special_block_map.get((playlist, tag))

        handle_block(block, args.dry_run, log_lines, special_sample_dir=special_dir)

        blocks_processed += 1
        if first_tag_seen is None and tag:
            first_tag_seen = tag
        if tag:
            tags_seen.add(tag)

        if args.max_blocks is not None and blocks_processed >= args.max_blocks:
            print(f"[INFO] Reached --max-blocks={args.max_blocks}, stopping.")
            break

    for tag in sorted(tags_seen):
        maybe_hadd_fhc(tag, args.dry_run, log_lines)
        
    # Write diagnostic log next to runningNotes (ALSO in dry-run)
    out_dir = os.path.dirname(os.path.abspath(args.runningNotes)) or "."
    tag_for_name = first_tag_seen or "UNKNOWN"
    # log_name = f"Log_{tag_for_name}_{notes_stamp}.txt"
    log_name = f"Log_{notes_stamp}.txt"
    log_path = os.path.join(out_dir, log_name)

    header = [
        f"Diagnostic log generated from runningNotes stamp: {notes_stamp}",
        f"Mode: {'DRY-RUN' if args.dry_run else 'RUN'}",
        f"Input runningNotes: {os.path.abspath(args.runningNotes)}",
        f"Selection tag (from file): {tag_for_name}",
        "",
    ]
    try:
        with open(log_path, "w", encoding="utf-8") as lf:
            lf.write("\n".join(header + log_lines))
        print(f"[INFO] Wrote diagnostic log: {log_path}")
    except OSError as e:
        print(f"[ERROR] Failed to write log file {log_path}: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()






# import os
# import re
# import sys
# import shlex
# import argparse
# import subprocess

# argv_backup = sys.argv[:]
# sys.argv = [sys.argv[0]]
# from config.AnalysisConfig import AnalysisConfig
# sys.argv = argv_backup

# from typing import Dict, Optional, List, Tuple


# RE_SUBMITTING = re.compile(r'^\s*Submitting:\s*(\S+)\s*$')
# RE_PLAYLIST_LINE = re.compile(r'^\s*playlist\s*=\s*(\S+)\s*$')
# RE_SAMPLE_LINE = re.compile(r'^\s*sample\s*=\s*(mc|data)\s*$', re.IGNORECASE)
# RE_NJOBS_LINE = re.compile(r'^\s*njobs\s*=\s*(\d+)\s*$')
# RE_OUTTAG_LINE = re.compile(r'^\s*out_tag\s*=\s*(\S+)\s*$')
# RE_RUNNING_COMMAND = re.compile(r'^\s*jobsub_submit\b')
# RE_NOTES_STAMP = re.compile(r'(.+_\d{4}-\d{2}-\d{2}_\d{6})\.txt$')

# EXPECTED_MC_PLAYLISTS = ["le1", "le7", "le9", "le13C"]
# EXPECTED_DATA_PLAYLISTS = ["le1", "le7", "le9", "le13A", "le13B", "le13C", "le13D", "le13E"]


# def strip_p6(s: str) -> str:
#     return s[:-3] if s.endswith("_p6") else s


# def count_root_files_for_playlist(path: str, playlist_raw: str) -> Dict[str, int]:
#     """
#     Count ROOT files in a shared output directory that belong to one playlist.
#     Expected filename patterns:
#       kin_dist_mcle1_p6_...
#       kin_dist_datale1_p6_...
#     """
#     counts = {"mc": 0, "data": 0, "other": 0, "ambiguous": 0, "total": 0}

#     mc_prefix = f"kin_dist_mc{playlist_raw}_"
#     data_prefix = f"kin_dist_data{playlist_raw}_"

#     try:
#         for fn in os.listdir(path):
#             if not fn.endswith(".root"):
#                 continue

#             is_mc = fn.startswith(mc_prefix)
#             is_data = fn.startswith(data_prefix)

#             if is_mc or is_data:
#                 counts["total"] += 1

#             if is_mc and not is_data:
#                 counts["mc"] += 1
#             elif is_data and not is_mc:
#                 counts["data"] += 1
#             elif is_mc and is_data:
#                 counts["ambiguous"] += 1
#             else:
#                 counts["other"] += 1

#     except FileNotFoundError:
#         pass

#     return counts


# def format_block_log(block: dict, file_counts: Dict[str, int]) -> List[str]:
#     playlist_raw = block["playlist"]
#     playlist = strip_p6(playlist_raw) if playlist_raw else None
#     tag = block["selection_tag"]
#     out_tag = block["out_tag"]
#     config = block["config"]

#     exp_mc = block["expected_mc_jobs"]
#     exp_data = block["expected_data_jobs"]
#     exp_total = exp_mc + exp_data

#     in_dir = f"/pnfs/minerva/scratch/users/qvuong/{config}/{out_tag}_{tag}_hists/" if config and out_tag else "(unknown)"

#     seen_mc = file_counts["mc"]
#     seen_data = file_counts["data"]
#     seen_total = file_counts["total"]

#     lines = []
#     lines.append("==================================================")
#     lines.append(f"Config              : {config}")
#     lines.append(f"Playlist            : {playlist} (raw: {playlist_raw})")
#     lines.append(f"Selection tag       : {tag}")
#     lines.append(f"Out tag             : {out_tag}")
#     lines.append(f"Input dir           : {in_dir}")
#     lines.append(f"Expected MC jobs    : {exp_mc}")
#     lines.append(f"Expected DATA jobs  : {exp_data}")
#     lines.append(f"Expected TOTAL      : {exp_total}")
#     lines.append(f"Seen MC ROOT files  : {seen_mc}")
#     lines.append(f"Seen DATA ROOT files: {seen_data}")
#     lines.append(f"Seen TOTAL ROOT files: {seen_total}")
#     if file_counts["ambiguous"] or file_counts["other"]:
#         lines.append(f"Seen ambiguous ROOT files: {file_counts['ambiguous']}")
#         lines.append(f"Seen other ROOT files    : {file_counts['other']}")

#     if seen_mc != exp_mc:
#         lines.append(f"[WARN] MC mismatch (expected {exp_mc}, saw {seen_mc})")
#     else:
#         lines.append("[OK] MC count matches")

#     if seen_data != exp_data:
#         lines.append(f"[WARN] DATA mismatch (expected {exp_data}, saw {seen_data})")
#     else:
#         lines.append("[OK] DATA count matches")

#     if seen_total != exp_total:
#         lines.append(f"[WARN] TOTAL mismatch (expected {exp_total}, saw {seen_total})")
#     else:
#         lines.append("[OK] TOTAL count matches")

#     lines.append("==================================================")
#     lines.append("")
#     return lines


# def block_ok(block: dict, file_counts: Dict[str, int]) -> Tuple[bool, bool]:
#     ok_mc = (file_counts["mc"] == block["expected_mc_jobs"])
#     ok_data = (file_counts["data"] == block["expected_data_jobs"])
#     return ok_mc, ok_data


# def handle_block(block: dict, dry_run: bool, log_sink: Optional[List[str]] = None):
#     if not block:
#         return

#     playlist_raw = block["playlist"]
#     tag = block["selection_tag"]
#     out_tag = block["out_tag"]
#     config = block["config"]

#     if not playlist_raw or not tag or not out_tag or not config:
#         msg = f"[WARN] Incomplete block, skipping: {block}"
#         print(msg, file=sys.stderr)
#         if log_sink is not None:
#             log_sink.append(msg)
#             log_sink.append("")
#         return

#     playlist = strip_p6(playlist_raw)
#     in_dir = f"/pnfs/minerva/scratch/users/qvuong/{config}/{out_tag}_{tag}_hists/"
#     file_counts = count_root_files_for_playlist(in_dir, playlist_raw)
#     ok_mc, ok_data = block_ok(block, file_counts)
#     should_run = ok_mc and ok_data

#     print("\n".join(format_block_log(block, file_counts)))
#     if log_sink is not None:
#         log_sink.extend(format_block_log(block, file_counts))

#     cmd = [
#         "python", "combine_file_specialSampleModification.py",
#         "--playlist", playlist,
#         "--i", in_dir,
#         "--cal_POT",
#         "--ntuple_tag", "MAD",
#         "--selection_tag", tag,
#     ]

#     print("Command:")
#     print(" ", " ".join(cmd))
#     if log_sink is not None:
#         log_sink.append("Command:")
#         log_sink.append("  " + " ".join(cmd))
#         log_sink.append("")

#     if dry_run:
#         if should_run:
#             print("[DRY-RUN] would execute combine_file_specialSampleModification.py (counts OK for MC and DATA)")
#             if log_sink is not None:
#                 log_sink.append("[DRY-RUN] would execute combine_file_specialSampleModification.py (counts OK for MC and DATA)")
#                 log_sink.append("")
#         else:
#             print("[DRY-RUN] would SKIP combine_file_specialSampleModification.py (count mismatch)")
#             if log_sink is not None:
#                 log_sink.append("[DRY-RUN] would SKIP combine_file_specialSampleModification.py (count mismatch)")
#                 log_sink.append("")
#     else:
#         if should_run:
#             print("[RUNNING] counts OK for MC and DATA")
#             if log_sink is not None:
#                 log_sink.append("[RUNNING] counts OK for MC and DATA")
#                 log_sink.append("")
#             subprocess.run(cmd, input="\n", text=True)
#         else:
#             print("[SKIP] Not running combine_file_specialSampleModification.py because counts did not match for MC and/or DATA")
#             if log_sink is not None:
#                 log_sink.append("[SKIP] Not running combine_file_specialSampleModification.py because counts did not match for MC and/or DATA")
#                 log_sink.append("")


# def combined_output_path(playlist: str, selection_tag: str, is_data: bool, ntuple_tag: str = "MAD") -> str:
#     old_playlist = getattr(AnalysisConfig, "playlist", None)
#     old_selection_tag = getattr(AnalysisConfig, "selection_tag", None)
#     old_ntuple_tag = getattr(AnalysisConfig, "ntuple_tag", None)

#     try:
#         AnalysisConfig.playlist = playlist
#         AnalysisConfig.selection_tag = selection_tag
#         AnalysisConfig.ntuple_tag = ntuple_tag
#         return AnalysisConfig.SelectionHistoPath(playlist, is_data)
#     finally:
#         AnalysisConfig.playlist = old_playlist
#         AnalysisConfig.selection_tag = old_selection_tag
#         AnalysisConfig.ntuple_tag = old_ntuple_tag


# def maybe_hadd_fhc(selection_tag: str, dry_run: bool, log_sink=None):
#     mc_files = [
#         combined_output_path(pl, selection_tag, False)
#         for pl in EXPECTED_MC_PLAYLISTS
#     ]
#     data_files = [
#         combined_output_path(pl, selection_tag, True)
#         for pl in EXPECTED_DATA_PLAYLISTS
#     ]

#     mc_dir = os.path.dirname(mc_files[0]) if mc_files else "."
#     data_dir = os.path.dirname(data_files[0]) if data_files else "."

#     out_mc = os.path.join(mc_dir, f"kin_dist_mcleFHC_{selection_tag}_MAD.root")
#     out_data = os.path.join(data_dir, f"kin_dist_dataleFHC_{selection_tag}_MAD.root")

#     missing_mc = [f for f in mc_files if not os.path.exists(f)]
#     missing_data = [f for f in data_files if not os.path.exists(f)]

#     lines = []
#     lines.append("################ FHC MERGE CHECK ################")
#     lines.append(f"Selection tag: {selection_tag}")
#     lines.append(f"MC output dir : {mc_dir}")
#     lines.append(f"DATA output dir: {data_dir}")
#     lines.append(f"Found MC      : {len(mc_files) - len(missing_mc)}/{len(mc_files)}")
#     lines.append(f"Found DATA    : {len(data_files) - len(missing_data)}/{len(data_files)}")

#     if missing_mc:
#         lines.append("[WARN] Missing MC files:")
#         lines.extend([f"  {x}" for x in missing_mc])
#     else:
#         lines.append("[OK] All MC playlist outputs are present")

#     if missing_data:
#         lines.append("[WARN] Missing DATA files:")
#         lines.extend([f"  {x}" for x in missing_data])
#     else:
#         lines.append("[OK] All DATA playlist outputs are present")

#     print("\n".join(lines))
#     if log_sink is not None:
#         log_sink.extend(lines)
#         log_sink.append("")

#     if missing_mc or missing_data:
#         msg = "[SKIP] Not running hadd because not all expected merged playlist outputs exist"
#         print(msg)
#         if log_sink is not None:
#             log_sink.append(msg)
#             log_sink.append("")
#         return

#     cmd_mc = ["hadd", "-f", out_mc] + mc_files
#     cmd_data = ["hadd", "-f", out_data] + data_files

#     print("FHC MC hadd command:")
#     print(" ", " ".join(cmd_mc))
#     print("FHC DATA hadd command:")
#     print(" ", " ".join(cmd_data))

#     if log_sink is not None:
#         log_sink.append("FHC MC hadd command:")
#         log_sink.append("  " + " ".join(cmd_mc))
#         log_sink.append("FHC DATA hadd command:")
#         log_sink.append("  " + " ".join(cmd_data))
#         log_sink.append("")

#     if dry_run:
#         print("[DRY-RUN] would run FHC MC/DATA hadd")
#         if log_sink is not None:
#             log_sink.append("[DRY-RUN] would run FHC MC/DATA hadd")
#             log_sink.append("")
#     else:
#         print("[RUNNING] hadd for FHC MC")
#         subprocess.run(cmd_mc, check=False)
#         print("[RUNNING] hadd for FHC DATA")
#         subprocess.run(cmd_data, check=False)


# def parse_jobsub_command(line: str) -> Dict[str, Optional[str]]:
#     """
#     Parse the printed jobsub_submit command line and recover:
#       config, playlist, sample, count, selection_tag, out_tag
#     Expected tail shape:
#       ... file://.../testWrapper.sh CONFIG PLAYLIST SAMPLE COUNT SELECTION_TAG OUT_TAG [extra args...]
#     """
#     result = {
#         "config": None,
#         "playlist": None,
#         "sample": None,
#         "selection_tag": None,
#         "out_tag": None,
#     }

#     try:
#         toks = shlex.split(line.strip())
#     except ValueError:
#         return result

#     wrapper_idx = None
#     for i, tok in enumerate(toks):
#         if tok.startswith("file://") and tok.endswith(".sh"):
#             wrapper_idx = i
#             break

#     if wrapper_idx is None:
#         return result

#     tail = toks[wrapper_idx + 1:]
#     if len(tail) < 6:
#         return result

#     result["config"] = tail[0]
#     result["playlist"] = tail[1]
#     result["sample"] = tail[2]
#     result["selection_tag"] = tail[4]
#     result["out_tag"] = tail[5]
#     return result


# def main():
#     parser = argparse.ArgumentParser()
#     parser.add_argument("runningNotes", help="runningNotes.txt file")
#     parser.add_argument(
#         "--dry-run",
#         action="store_true",
#         help="Print commands only, do not execute combine_file_specialSampleModification.py",
#     )
#     parser.add_argument(
#         "--max-blocks",
#         type=int,
#         default=None,
#         help="Maximum number of playlist groups to process (default: all)",
#     )
#     args = parser.parse_args()

#     notes_basename = os.path.basename(args.runningNotes)
#     m_notes = RE_NOTES_STAMP.match(notes_basename)
#     notes_stamp = m_notes.group(1) if m_notes else "UNKNOWN"

#     grouped: Dict[Tuple[str, str, str, str], dict] = {}
#     blocks_processed = 0

#     log_lines: List[str] = []
#     tags_seen = set()

#     current_playlist = None
#     current_sample = None
#     current_njobs = None
#     current_out_tag = None

#     with open(args.runningNotes, "r", encoding="utf-8", errors="replace") as f:
#         for line in f:
#             m = RE_PLAYLIST_LINE.match(line)
#             if m:
#                 current_playlist = m.group(1)
#                 continue

#             m = RE_SAMPLE_LINE.match(line)
#             if m:
#                 current_sample = m.group(1).lower()
#                 continue

#             m = RE_NJOBS_LINE.match(line)
#             if m:
#                 current_njobs = int(m.group(1))
#                 continue

#             m = RE_OUTTAG_LINE.match(line)
#             if m:
#                 current_out_tag = m.group(1)
#                 continue

#             if RE_RUNNING_COMMAND.match(line):
#                 parsed = parse_jobsub_command(line)
#                 config = parsed["config"]
#                 playlist = parsed["playlist"] or current_playlist
#                 sample = parsed["sample"] or current_sample
#                 selection_tag = parsed["selection_tag"]
#                 out_tag = parsed["out_tag"] or current_out_tag

#                 if not all([config, playlist, sample, selection_tag, out_tag]):
#                     msg = f"[WARN] Could not fully parse command line, skipping:\n{line.rstrip()}"
#                     print(msg, file=sys.stderr)
#                     log_lines.append(msg)
#                     log_lines.append("")
#                     continue

#                 key = (config, playlist, selection_tag, out_tag)
#                 if key not in grouped:
#                     grouped[key] = {
#                         "config": config,
#                         "playlist": playlist,
#                         "selection_tag": selection_tag,
#                         "out_tag": out_tag,
#                         "expected_mc_jobs": 0,
#                         "expected_data_jobs": 0,
#                     }

#                 if sample == "mc":
#                     grouped[key]["expected_mc_jobs"] += (current_njobs or 0)
#                 elif sample == "data":
#                     grouped[key]["expected_data_jobs"] += (current_njobs or 0)

#                 tags_seen.add(selection_tag)

#                 current_playlist = None
#                 current_sample = None
#                 current_njobs = None
#                 current_out_tag = None

#     for block in grouped.values():
#         handle_block(block, args.dry_run, log_lines)
#         blocks_processed += 1
#         if args.max_blocks is not None and blocks_processed >= args.max_blocks:
#             print(f"[INFO] Reached --max-blocks={args.max_blocks}, stopping.")
#             break

#     for tag in sorted(tags_seen):
#         maybe_hadd_fhc(tag, args.dry_run, log_lines)

#     log_dir = "/exp/minerva/data/users/qvuong/runningNotes"
#     os.makedirs(log_dir, exist_ok=True)
#     log_name = f"Log_{notes_stamp}.txt"
#     log_path = os.path.join(log_dir, log_name)

#     header = [
#         f"Diagnostic log generated from runningNotes stamp: {notes_stamp}",
#         f"Mode: {'DRY-RUN' if args.dry_run else 'RUN'}",
#         f"Input runningNotes: {os.path.abspath(args.runningNotes)}",
#         "",
#     ]
#     try:
#         with open(log_path, "w", encoding="utf-8") as lf:
#             lf.write("\n".join(header + log_lines))
#         print(f"[INFO] Wrote diagnostic log: {log_path}")
#     except OSError as e:
#         print(f"[ERROR] Failed to write log file {log_path}: {e}", file=sys.stderr)


# if __name__ == "__main__":
#     main()

