"""
  EventClassification.py:
   Divide events in the ntuple into nu_e CCQE signal sample
   and various background categories.
   Importable.
  
   Original author: J. Wolcott (jwolcott@fnal.gov)
                    February 2013
"""

import ROOT
from config.CutConfig import SAMPLE_CUTS,KINEMATICS_CUTS
# import tools.CutLibrary as CL
# print(">>> CutLibrary loaded from:", CL.__file__)
# print(">>> Has CalcApothem?", hasattr(CL, "CalcApothem"))
# print(">>> Has CalcApothem1?", hasattr(CL, "CalcApothem1"))
from tools.CutLibrary import CUTS
from config.SignalDef import SIGNAL_DEFINITION,TRUTH_CATEGORIES,EXTRA_OTHER



class EventClassifier(object):

    """ Class doing event classification (applies cuts, determines if event is true or reco signal, etc.). """

    # dictionaries are ordered alphabetically by key...
    # so we need to keep track of the order separately

    def __init__(self, classifiers=["Reco", "Truth"], is_PC=False, use_kin_cuts=False, use_sideband= None):

        # copy the cuts.
        # we use several classifiers simultaneously
        # for the background selection (each with
        # slightly modified cuts), so each classifier
        # needs to work with its own copy internally.
        self.cuts = CUTS
        self.use_kin_cuts = use_kin_cuts
        if not self.use_kin_cuts:
            del TRUTH_CATEGORIES["NonPhaseSpace"]

        self.samples = {"Signal"}

        if use_sideband is None:
            self.samples.update(list(SAMPLE_CUTS.keys()))
        else:
            self.samples.update(use_sideband)

        self.used_cuts = set().union(*[v for k,v in SAMPLE_CUTS.items() if k in self.samples])
        #remove inverse cuts, we use not in the future
        

        if self.use_kin_cuts:
            self.used_cuts.update(["Reco{}".format(cut_name) for cut_name in KINEMATICS_CUTS])

        for j in  [i for i in self.used_cuts if i.startswith("Inverse")]:
            if j[7:] not in self.used_cuts:
                self.used_cuts.add(j[7:])
            self.used_cuts.remove(j)

        #self.is_PC = is_PC

        #self._warned_r2 = False

        self.classifiers_to_use = []
        for classifier in classifiers:
            cl = "_Classify%s" % classifier.capitalize()
            if hasattr(self, cl):
                self.classifiers_to_use.append(getattr(self, cl))


        self.counter = [0,0]
        #self.RegisterBranches()

        #instantiate centralized cuts:
        #self.ZRange = ROOT.reco.ZRange(ROOT.PythonMinervaUniverse,ROOT.PlotUtils.detail.empty)("Tracker",5990,8340)
        #self.Apothem = ROOT.reco.Apothem(ROOT.PythonMinervaUniverse,ROOT.PlotUtils.detail.empty)(850)
        #self.HasMINOSMatch = ROOT.reco.HasMINOSMatch(ROOT.PythonMinervaUniverse,ROOT.PlotUtils.detail.empty)()


    def Classify(self,entry):
        """ Actually do the classification.  Should be called only once per event.

            Note that only the classifiers that were registered when the object was
            initialized will be used.  (By default they are all enabled.)"""

        self.tree = entry
        # Classify new event, clean up
        self.reco_cuts_passed = dict.fromkeys(self.used_cuts,False)
        self.is_reco_signal = False
        self.side_band =None
        try:
            self.new_truth = self.tree.ShortName() == "cv"
        except AttributeError:
            self.new_truth=True

        if self.new_truth:
            self.truth_class = None
            self.is_true_signal = False


        sc = True
        for cl in self.classifiers_to_use:
            sc = cl() and sc

        if self.new_truth:
            if self.side_band is not None:
                self.counter[0] += 1
            if self.is_true_signal:
                self.counter[1] += 1
        return sc

    def _ClassifyReco(self):
        """ Determine if the event passes the cuts on reconstructed quantities. """

        # Check all cuts for simplicity. A better way is only checking cuts that will be used.
        for cut_name in self.reco_cuts_passed:
            cut = self.cuts[cut_name]
            try :
                # I sorted the prong  so only check 0.
                self.reco_cuts_passed[cut_name] = cut.DoesEventPass(self.tree, 0)
            except (IndexError,AttributeError,TypeError):
                # this value wasn't defined because the reco didn't get this far.
                # don't crash, just assume this cut failed.
                continue
        # return classification result by selected cuts.
        for sample in self.samples:
            passed = all(self.reco_cuts_passed[cut_name] for cut_name in SAMPLE_CUTS[sample]
                         if not cut_name.startswith("Inverse"))
            passed = passed and not any(self.reco_cuts_passed[cut_name[7:]] for cut_name in SAMPLE_CUTS[sample]
                                        if cut_name.startswith("Inverse"))
            if self.use_kin_cuts:
                passed = passed and all(self.reco_cuts_passed["Reco{}".format(cut_name)] for cut_name in KINEMATICS_CUTS)
            if sample == "Signal":
                self.is_reco_signal = passed
            if passed:
                self.side_band = sample
                break #samples are exclusive
        return True

    def _ClassifyTruth(self):
        if not self.new_truth:
            return True

        self.truth_class = "Other"

        for sig_def,func in TRUTH_CATEGORIES.items():
            if func(self.tree):
                self.truth_class = sig_def
                break

        if self.truth_class in EXTRA_OTHER:
            self.truth_class = "Other"

        self.is_true_signal = self.truth_class in SIGNAL_DEFINITION
        return True


### No longer used
    def RegisterBranches(self):
        """ Excerpted here from SelectorPlotter
            so that we can always be sure we get
            the branches we need for this object. """

        branches_to_enable = [
            "wgt",

            "michel_*",

            "phys_n_dead*",

            "prong_part_score",
            "prong_CalCorrectedEnergy",
            "prong_DataFrac",
            "prong_dEdXMeanFrontTracker",
            "prong_ECALVisE",
            "prong_FirstFireFraction",
            "prong_HCALVisE",
            "prong_NonMIPClusFrac",
            "prong_ODVisE",
            "prong_SideECALVisE",
            "prong_TransverseGapScore",
            "prong_TotalVisE",
            "prong_binned_energy_bin_indices",
            "prong_binned_energy_bin_contents",
            "prong_projection_bin_width",

            "EMLikeTrackMultiplicity",
            "HasFiducialVertex",
            "HasNoBackExitingTracks",
            "HasNoNonEMExitingTracks",
            "HasNoVertexMismatch",
            "HasTracks",
            "NuECCQE_E",
            "NuECCQE_Q2",
            "Psi",
            "StartPointVertexMultiplicity",
#           "UpstreamODEnergy",
            "VertexTrackMultiplicity",

            "vtx",
        ]
        
        if any( [e in self.cut_sequence for e in EventClassifier.KINEMATICS_CUT_SEQUENCE] ):
            branches_to_enable += CutConfig.kin_calculator.RequiredBranches()
        
#       print "using classifiers:", self.classifiers_to_use
        if self._ClassifyTruth in self.classifiers_to_use:
#           print "got Truth."
            branches_to_enable.extend([
                "mc_incoming",   # incoming PDG
                "mc_incomingE",
                "mc_current",
                "mc_intType",
                "mc_nFSPart",
                "mc_FSPartPDG",
                "truth_pass_NuECCQE",
                "truth_IsFiducial",
            ])

        for branch in branches_to_enable:
#           print 'enabling branch:', branch
            self.tree.SetBranchStatus(branch, 1)
