import ROOT

filelist = "playlist_mcle1p6.txt"

counts = {
    "POT_Used": 0,
    "POT_Total_only": 0,
    "No_Meta": 0,
    "No_POT_branch": 0,
    "Open_fail": 0,
}

with open(filelist) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        tf = None
        try:
            tf = ROOT.TFile.Open(line)
        except Exception as e:
            print("[OPEN FAIL]", line)
            print("   ", e)
            counts["Open_fail"] += 1
            continue

        if not tf or tf.IsZombie():
            print("[BAD OPEN]", line)
            counts["Open_fail"] += 1
            continue

        meta = tf.Get("Meta")
        if not meta:
            print("[NO META]", line)
            counts["No_Meta"] += 1
            tf.Close()
            continue

        has_used = bool(meta.GetBranch("POT_Used"))
        has_total = bool(meta.GetBranch("POT_Total"))

        if has_used:
            counts["POT_Used"] += 1
        elif has_total:
            print("[POT_Total ONLY]", line)
            counts["POT_Total_only"] += 1
        else:
            print("[NO POT BRANCH]", line)
            counts["No_POT_branch"] += 1

        tf.Close()

print("\nSummary:")
for k, v in counts.items():
    print(f"  {k:15s} = {v}")