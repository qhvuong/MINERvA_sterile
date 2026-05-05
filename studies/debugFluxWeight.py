# import ROOT

# flux_dir = "/exp/minerva/app/users/qvuong/MAT_AL9/CC-NuE-XSec/custom_plotutils/data/flux"

# def check_flux_ratios(pdg, playlist, energies):
#     cv_path = f"{flux_dir}/flux-gen2thin-pdg{pdg}-{playlist}.root"
#     gen_path = f"{flux_dir}/flux-g4numiv5-pdg{pdg}-{playlist}.root"

#     fcv = ROOT.TFile.Open(cv_path)
#     fgen = ROOT.TFile.Open(gen_path)

#     hcv = fcv.Get("flux_E_cvweighted")
#     hgen = fgen.Get("flux_E_unweighted")

#     if not hcv:
#         raise RuntimeError(f"Missing flux_E_cvweighted in {cv_path}")
#     if not hgen:
#         raise RuntimeError(f"Missing flux_E_unweighted in {gen_path}")
#     if not hcv.HasVertErrorBand("Flux"):
#         raise RuntimeError(f"No Flux band in {cv_path}")

#     eb = hcv.GetVertErrorBand("Flux")
#     n = eb.GetNHists()

#     print("\n====================================================")
#     print("pdg =", pdg, "playlist =", playlist)
#     print("CV file  =", cv_path)
#     print("GEN file =", gen_path)
#     print("n Flux universes =", n)
#     print("CV hist integral =", hcv.Integral())
#     print("Flux band CV integral =", eb.Integral())
#     print("bandCV/mainCV =", eb.Integral() / hcv.Integral() if hcv.Integral() else "NA")

#     for E in energies:
#         cv_val = hcv.Interpolate(E)
#         gen_val = hgen.Interpolate(E)
#         cv_weight = cv_val / gen_val if gen_val else float("nan")

#         ratios = []
#         full_weights = []

#         for i in range(n):
#             uval = eb.GetHist(i).Interpolate(E)
#             ratio = uval / cv_val if cv_val else float("nan")
#             ratios.append(ratio)
#             full_weights.append(ratio * cv_weight)

#         print("\nE =", E)
#         print("  cvweighted flux =", cv_val)
#         print("  unweighted flux  =", gen_val)
#         print("  CV flux weight cvweighted/unweighted =", cv_weight)
#         print("  universe/CV ratios min mean max = {:.8g} {:.8g} {:.8g}".format(
#             min(ratios), sum(ratios) / len(ratios), max(ratios)
#         ))
#         print("  full universe weights min mean max = {:.8g} {:.8g} {:.8g}".format(
#             min(full_weights), sum(full_weights) / len(full_weights), max(full_weights)
#         ))
#         print("  first 10 universe/CV ratios:")
#         for i, r in enumerate(ratios[:10]):
#             print("    {:3d}: {:.8g}".format(i, r))

# check_flux_ratios(
#     pdg=12,
#     playlist="minerva1",
#     energies=[19.723, 8.47426, 4.51719, 3.19317],
# )

# check_flux_ratios(
#     pdg=14,
#     playlist="minerva1",
#     energies=[3.47577, 8.6602, 17.1428],
# )



import ROOT

path = "/pnfs/minerva/persistent/users/qvuong/custom_plotutils/data/flux/LE13_nue.root"

def inspect_file(path):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError("Could not open " + path)

    print("\nKeys in file:")
    for key in f.GetListOfKeys():
        print(" ", key.GetName(), key.GetClassName())

    # Likely names from compute_flux.py:
    # flux_E_unweighted
    # flux_E_cvweighted
    for hname in ["flux_E_unweighted", "flux_E_cvweighted"]:
        h = f.Get(hname)
        if not h:
            print("\nMissing", hname)
            continue

        print("\n====================================================")
        print("Histogram:", hname)
        print("Title:", h.GetTitle())
        print("Main CV integral:", h.Integral())

        if not h.HasVertErrorBand("Flux"):
            print("No Flux band")
            continue

        eb = h.GetVertErrorBand("Flux")
        print("Flux band CV integral:", eb.Integral())
        print("bandCV/mainCV:", eb.Integral() / h.Integral() if h.Integral() else "NA")
        print("n Flux universes:", eb.GetNHists())

        energies = [3.19317, 3.47577, 4.51719, 8.47426, 8.6602, 17.1428, 19.723]

        for E in energies:
            cv = h.Interpolate(E)
            ratios = []

            for i in range(eb.GetNHists()):
                u = eb.GetHist(i).Interpolate(E)
                ratios.append(u / cv if cv else float("nan"))

            print("\nE =", E)
            print("  CV =", cv)
            print("  universe/CV min mean max = {:.8g} {:.8g} {:.8g}".format(
                min(ratios),
                sum(ratios) / len(ratios),
                max(ratios),
            ))
            print("  first 10 universe/CV ratios:")
            for i, r in enumerate(ratios[:10]):
                print("    {:3d}: {:.8g}".format(i, r))

    f.Close()

inspect_file(path)