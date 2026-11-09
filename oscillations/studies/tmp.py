import ROOT

ROOT.gROOT.SetBatch(True)

mc_path = (
    "/exp/minerva/data/users/qvuong/nu_e/"
    "kin_dist_mcleFHC_CCnue_FullCuts_v4p3p1_MAD.root"
)

scale = 0.150312908217

f = ROOT.TFile.Open(mc_path)

categories = [
    "CCNuEQE", "CCNuEDelta", "CCNuEDIS",
    "CCNuE2p2h", "CCNuE", "CCNuEWrongSign",
    "NCPi0", "NCCohPi0", "CCPi0",
    "NCPi", "CCPi", "CCOther",
    "NCOther", "NuEElastic", "Other",
]

for region in ("EN4", "EN4_dEdX"):
    print(f"\n===== {region} =====")

    total = f.Get(region)
    for bin_idx in (1, 2):
        lo = total.GetXaxis().GetBinLowEdge(bin_idx)
        hi = total.GetXaxis().GetBinUpEdge(bin_idx)

        print(f"\nEnergy bin {bin_idx}: {lo:.2f}–{hi:.2f} GeV")
        print(f"{'Category':<22} {'MC at data POT':>15}")

        components_sum = 0.0

        for category in categories:
            h = f.Get(f"{region}_{category}")
            if not h:
                continue

            value = h.GetBinContent(bin_idx) * scale
            components_sum += value

            print(f"{category:<22} {value:15.3f}")

        print(f"{'Sum of categories':<22} {components_sum:15.3f}")
        print(
            f"{'Total histogram':<22} "
            f"{total.GetBinContent(bin_idx) * scale:15.3f}"
        )

f.Close()