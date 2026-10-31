
import ROOT
import math
import os

ROOT.gROOT.SetBatch(True)

directory = "/exp/minerva/data/users/qvuong/antinu_mu"

tags = [
    "ValidateCuts_CCnumubar_Inclusive_v1",
    "ValidateCuts_CCnumubar_FV_HasMINOSMatch_v1",
    "ValidateCuts_CCnumubar_FV_InverseHasNoBackExitingTracks_v1",
    "ValidateCuts_CCnumubar_FV_InverseHasNoBackExitingTracks_DeadTime_LeptonAngle_v1",
]

results = {}

for tag in tags:
    print("\n==========", tag, "==========")

    yields = {}

    for sample in ["data", "mc"]:
        path = os.path.join(
            directory,
            f"kin_dist_{sample}le5_{tag}_MAD.root"
        )

        f = ROOT.TFile.Open(path)

        if not f or f.IsZombie():
            print(sample, "FILE MISSING")
            continue

        h = f.Get("Eel")

        if not h:
            print(sample, "Eel MISSING")
            f.Close()
            continue

        integral = h.Integral(0, h.GetNbinsX() + 1)
        entries = h.GetEntries()

        print(
            sample,
            "entries =", entries,
            "raw integral =", integral
        )

        yields[sample] = integral
        f.Close()

    if "data" in yields and "mc" in yields:
        if yields["mc"] != 0 and math.isfinite(yields["mc"]):
            ratio = yields["data"] / yields["mc"]
            print("Raw data/MC =", ratio)
        else:
            ratio = float("nan")

        results[tag] = {
            "data": yields["data"],
            "mc": yields["mc"],
            "ratio": ratio,
        }

print("\n========== SUMMARY ==========")

for tag, values in results.items():
    print(
        tag,
        "data =", round(values["data"], 3),
        "mc =", round(values["mc"], 3),
        "raw D/MC =", round(values["ratio"], 5),
    )
