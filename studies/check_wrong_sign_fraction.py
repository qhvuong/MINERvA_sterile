"""
Compute CCNuE wrong-sign fraction in MC truth.
"""

from config.AnalysisConfig import AnalysisConfig
from tools import Utilities
from tools.EventClassification import EventClassifier
from tools.KinematicsCalculator import KinematicsCalculator
from tools.SystematicsUniverse import GetAllSystematicsUniverses


WRONG_SIGN_CATEGORIES = ("CCNuEWrongSignQE", "CCNuEWrongSign")


def count_wrong_sign_fraction(chainwrapper):
    kin_cal = KinematicsCalculator(
        correct_beam_angle=True,
        correct_MC_energy_scale=False,
        calc_true=True,
        calc_reco=False,
    )
    eventClassifier = EventClassifier(
        classifiers=["Truth"],
        use_kin_cuts=True,
        use_sideband=[],
    )

    universes = GetAllSystematicsUniverses(chainwrapper, False)
    cv_universe = next(iter(universes.values()))[0]
    cv_universe.LoadTools(kin_cal, eventClassifier)

    nEvents = chainwrapper.GetEntries()
    if AnalysisConfig.testing and nEvents > 5000:
        nEvents = 5000

    counts = {
        "total_truth_signal": 0,
        "wrong_sign_qe": 0,
        "wrong_sign_nonqe": 0,
    }

    for counter in range(nEvents):
        if counter % 100000 == 0:
            print(counter)

        cv_universe.SetEntry(counter)

        if not cv_universe.IsVerticalOnly():
            kin_cal.CalculateKinematics(cv_universe)
            eventClassifier.Classify(cv_universe)

        if eventClassifier.is_true_signal:
            counts["total_truth_signal"] += 1

        if eventClassifier.truth_class == "CCNuEWrongSignQE":
            counts["wrong_sign_qe"] += 1
        elif eventClassifier.truth_class == "CCNuEWrongSign":
            counts["wrong_sign_nonqe"] += 1

    wrong_sign_total = counts["wrong_sign_qe"] + counts["wrong_sign_nonqe"]
    denom = counts["total_truth_signal"]
    frac = (wrong_sign_total / float(denom)) if denom else 0.0

    print("Total truth signal:", denom)
    print("Wrong-sign QE:", counts["wrong_sign_qe"])
    print("Wrong-sign non-QE:", counts["wrong_sign_nonqe"])
    print("Wrong-sign total:", wrong_sign_total)
    print("Wrong-sign / truth signal:", frac)
    print("Wrong-sign / truth signal (%):", frac * 100.0)


if __name__ == "__main__":
    chainwrapper = Utilities.fileChain(
        AnalysisConfig.playlist,
        "mc",
        AnalysisConfig.ntuple_tag,
        "Truth",
        AnalysisConfig.count[0],
        AnalysisConfig.count[1],
    )
    count_wrong_sign_fraction(chainwrapper)
