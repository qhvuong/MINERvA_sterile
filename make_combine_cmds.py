#!/usr/bin/env python3
import os
import re
import sys
import argparse
import subprocess

RE_RUN_CMD = re.compile(r'^\s*Running command:\s*$')
RE_SEL_CMD = re.compile(r'^\s*python\s+selection/gridSelection\.py\b(.*)$')
RE_PLAYLIST = re.compile(r'--playlist\s+(\S+)')
RE_TAG = re.compile(r'--selection_tag\s+(\S+)')
RE_TAR = re.compile(r'qvuong-CCNUE_selection_(\d{4}-\d{2}-\d{2}-\d{6})\.tar\.gz')
RE_JOBS = re.compile(r'^\s*(\d+)\s+job\(s\)\s+submitted\b')

def strip_p6(s: str) -> str:
    return s[:-3] if s.endswith("_p6") else s

def count_root_files(path: str) -> int:
    try:
        return sum(1 for f in os.listdir(path) if f.endswith(".root"))
    except FileNotFoundError:
        return 0

def handle_block(block: dict, dry_run: bool):
    if not block:
        return

    playlist_raw = block["playlist"]
    tag = block["selection_tag"]
    stamp = block["stamp"]
    expected = block["expected_jobs"]

    if not playlist_raw or not tag or not stamp:
        print(f"[WARN] Incomplete block, skipping: {block}", file=sys.stderr)
        return

    playlist = strip_p6(playlist_raw)
    in_dir = f"/pnfs/minerva/persistent/users/qvuong/CCNUE_selection_{stamp}_hists/"
    seen = count_root_files(in_dir)

    print("==================================================")
    print(f"Playlist         : {playlist} (raw: {playlist_raw})")
    print(f"Selection tag    : {tag}")
    print(f"Stamp            : {stamp}")
    print(f"Input dir        : {in_dir}")
    print(f"Expected outputs : {expected}")
    print(f"ROOT files seen  : {seen}")

    if seen != expected:
        print(f"[WARN] count mismatch (expected {expected}, saw {seen})")
    else:
        print("[OK] file count matches")

    cmd = [
        "python", "combine_file.py",
        "--playlist", playlist,
        "--i", in_dir,
        "--cal_POT",
        "--ntuple_tag", "MAD",
        "--selection_tag", tag,
    ]

    print("Command:")
    print(" ", " ".join(cmd))

    if dry_run:
        print("[DRY-RUN] not executing")
    else:
        print("[RUNNING]")
        subprocess.run(cmd, input="\n", text=True)

    print("==================================================\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("runningNotes", help="runningNotes.txt file")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands only, do not execute combine_file.py",
    )
    parser.add_argument(
        "--max-blocks",
        type=int,
        default=None,
        help="Maximum number of submission blocks to process (default: all)",
    )
    args = parser.parse_args()

    current = None
    waiting_for_cmd = False
    blocks_processed = 0

    with open(args.runningNotes, "r", encoding="utf-8", errors="replace") as f:
        for line in f:

            # Case A: old format marker line
            if RE_RUN_CMD.match(line):
                if current is not None:
                    handle_block(current, args.dry_run)
                    blocks_processed += 1
                    if args.max_blocks is not None and blocks_processed >= args.max_blocks:
                        print(f"[INFO] Reached --max-blocks={args.max_blocks}, stopping.")
                        return

                current = {
                    "playlist": None,
                    "selection_tag": None,
                    "stamp": None,
                    "expected_jobs": 0,
                }
                waiting_for_cmd = True
                continue

            # Case B: old format, command comes right after marker
            if waiting_for_cmd:
                waiting_for_cmd = False
                m = RE_SEL_CMD.match(line)
                if m:
                    args_line = m.group(1)
                    mp = RE_PLAYLIST.search(args_line)
                    mt = RE_TAG.search(args_line)
                    current["playlist"] = mp.group(1) if mp else None
                    current["selection_tag"] = mt.group(1) if mt else None
                continue

            # Case C: NEW format — command appears directly (no "Running command:" line)
            m_direct = RE_SEL_CMD.match(line)
            if m_direct:
                # finalize previous block
                if current is not None:
                    handle_block(current, args.dry_run)
                    blocks_processed += 1
                    if args.max_blocks is not None and blocks_processed >= args.max_blocks:
                        print(f"[INFO] Reached --max-blocks={args.max_blocks}, stopping.")
                        return

                # start new block from this command line
                current = {
                    "playlist": None,
                    "selection_tag": None,
                    "stamp": None,
                    "expected_jobs": 0,
                }
                args_line = m_direct.group(1)
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

                mj = RE_JOBS.search(line)
                if mj:
                    current["expected_jobs"] += int(mj.group(1))

    # current = None
    # waiting_for_cmd = False
    # blocks_processed = 0

    # with open(args.runningNotes, "r", encoding="utf-8", errors="replace") as f:
    #     for line in f:
    #         if RE_RUN_CMD.match(line):
    #             if current is not None:
    #                 handle_block(current, args.dry_run)
    #                 blocks_processed += 1
    #                 if args.max_blocks is not None and blocks_processed >= args.max_blocks:
    #                     print(f"[INFO] Reached --max-blocks={args.max_blocks}, stopping.")
    #                     return

    #             current = {
    #                 "playlist": None,
    #                 "selection_tag": None,
    #                 "stamp": None,
    #                 "expected_jobs": 0,
    #             }
    #             waiting_for_cmd = True
    #             continue

    #         if waiting_for_cmd:
    #             waiting_for_cmd = False
    #             m = RE_SEL_CMD.match(line)
    #             if m:
    #                 args_line = m.group(1)
    #                 mp = RE_PLAYLIST.search(args_line)
    #                 mt = RE_TAG.search(args_line)
    #                 current["playlist"] = mp.group(1) if mp else None
    #                 current["selection_tag"] = mt.group(1) if mt else None
    #             continue

    #         if current is not None:
    #             mtar = RE_TAR.search(line)
    #             if mtar and current["stamp"] is None:
    #                 current["stamp"] = mtar.group(1)

    #             mj = RE_JOBS.search(line)
    #             if mj:
    #                 current["expected_jobs"] += int(mj.group(1))

    # Handle last block
    if current is not None:
        if args.max_blocks is None or blocks_processed < args.max_blocks:
            handle_block(current, args.dry_run)

if __name__ == "__main__":
    main()
