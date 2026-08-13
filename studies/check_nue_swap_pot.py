#!/usr/bin/env python3

import ROOT


FILES = {
    "fhc_nue": (
        "/exp/minerva/data/users/qvuong/nu_e/"
        "kin_dist_mcleFHC_p8Tuples_CCnue_updatedFluxes_MAD.root"
    ),
    "fhc_swap": (
        "/exp/minerva/data/users/qvuong/nu_e_swapped/"
        "kin_dist_mcleFHC_p8Tuples_CCnueswap_updatedFluxes_MAD.root"
    ),
    "rhc_nue": (
        "/exp/minerva/data/users/qvuong/antinu_e/"
        "kin_dist_mcle5_p8Tuples_CCnuebar_updatedFluxes_MAD.root"
    ),
    "rhc_swap": (
        "/exp/minerva/data/users/qvuong/antinu_e_swapped/"
        "kin_dist_mcle5_p8Tuples_CCnuebarswap_updatedFluxes_MAD.root"
    ),
}


def inspect_file(label, path):
    print("\n" + "=" * 80)
    print(label)
    print(path)

    root_file = ROOT.TFile.Open(path)

    if not root_file or root_file.IsZombie():
        print("ERROR: could not open file")
        return

    print("\nTop-level keys:")
    for key in root_file.GetListOfKeys():
        print(
            "  {:40s}  {}".format(
                key.GetName(),
                key.GetClassName(),
            )
        )

    meta = root_file.Get("Meta")

    if not meta:
        print("\nNo Meta tree found.")
        root_file.Close()
        return

    print("\nMeta entries:", meta.GetEntries())

    branches = [
        branch.GetName()
        for branch in meta.GetListOfBranches()
    ]

    print("Meta branches:")
    for branch in branches:
        print("  ", branch)

    pot_branch = None

    for candidate in [
        "POT_Used",
        "POT_Total",
        "POT",
        "pot",
    ]:
        if candidate in branches:
            pot_branch = candidate
            break

    if pot_branch is None:
        print("\nNo recognized POT branch found.")
        root_file.Close()
        return

    total_pot = 0.0

    for entry in meta:
        total_pot += float(getattr(entry, pot_branch))

    print("\nPOT branch:", pot_branch)
    print("Summed POT: {:.16e}".format(total_pot))

    root_file.Close()


def get_pot(path):
    root_file = ROOT.TFile.Open(path)

    if not root_file or root_file.IsZombie():
        raise RuntimeError("Could not open {}".format(path))

    meta = root_file.Get("Meta")

    if not meta:
        raise RuntimeError("No Meta tree in {}".format(path))

    total = 0.0

    for entry in meta:
        total += float(entry.POT_Used)

    root_file.Close()
    return total


def main():
    pots = {}

    for label, path in FILES.items():
        inspect_file(label, path)
        pots[label] = get_pot(path)

    print("\n" + "=" * 80)
    print("POT COMPARISON")

    print(
        "FHC nominal/swap POT ratio = {:.12f}".format(
            pots["fhc_nue"] / pots["fhc_swap"]
        )
    )

    print(
        "RHC nominal/swap POT ratio = {:.12f}".format(
            pots["rhc_nue"] / pots["rhc_swap"]
        )
    )


if __name__ == "__main__":
    main()