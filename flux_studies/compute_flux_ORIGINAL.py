#!/bin/env python
"""
  compute_flux.py:
   Calculate the flux for a variety of beam configurations.
   
   For an example of how to extend this work
   and make other flux-related plots, see
   Ana/CCNuE/macros/flux/nue_flux_LE.py.
  
   Original author: J. Wolcott (jwolcott@fnal.gov)
                    November 2012
"""
import sys
import array
import numbers
import os.path
import optparse
import ROOT
ROOT.gSystem.Load("libFluxLoop.so")
# import PyCintex
import PlotUtils.LoadPlotUtilsLib
# import PlotUtils.MinervaStyle
class EnumeratorType(object):
        """ A helper class that makes closed enumerations easy """
        def __init__(self, enumeration={}, capitalization="standard"):
                self._capitalization = capitalization
                self._attrs = {}
                self._attrs_reversed = {}
        
                # dictionary-like things
                if hasattr(enumeration, "keys"):
                        self._attrs = enumeration.copy()
                # list-like things (values auto-assigned)
                elif hasattr(enumeration, "__iter__"):
                        self._attrs = dict( [ (k, i) for i, k in enumerate(enumeration) ] )
                        
                # enforce capitalization, key-value type consistency
                for key, val in self._attrs.items():
                        if self._capitalize_key(key) != key:
                                del self._attrs[key]
                                key = self._capitalize_key(key)
                                self._attrs[key] = val
                        if isinstance(key, int):
                                raise AttributeError("Can't use integers for enumeration keys!")
                        elif not isinstance(val, int):
                                raise AttributeError("Invalid enumeration value: '%s'.  Values must be integers." % val)
                
                        if val in self._attrs_reversed:
                                raise AttributeError("Enumeration values must be unique: got value %d multiple times." % val)
                        self._attrs_reversed[val] = key
                
        
        def __getattr__(self, key):
                key = self._capitalize_key(key)
                
                if key in self._attrs:
                        return self._attrs[key]
                elif key in self._attrs_reversed:
                        return self._attrs_reversed[key]
                
                raise AttributeError("Key '%s' not found in enumeration" % key)
                
        def __getitem__(self, key):
                return self.__getattr__(key)
                
        def __contains__(self, key):
                key = self._capitalize_key(key)
                
                return (key in self._attrs) or (key in self._attrs_reversed)
        
        def _capitalize_key(self, key):
                if not hasattr(key, "capitalize"):
                        return key
        
                if self._capitalization == "standard":
                        return key.capitalize()
                elif self._capitalization == "all_upper":
                        return key.upper()
                elif self._capitalization == "all_lower":
                        return key.lower()
                else:
                        return key
                
        def keys(self):
                return list(self._attrs.keys())
        
        def values(self):
                return list(self._attrs.values())
                
        def iterkeys(self):
                for key in self._attrs:
                        yield key
        
        def itervalues(self):
                for val in self._attrs.values():
                        yield val
        
        def iteritems(self):
                for key, val in self._attrs.items():
                        yield key, val
HornCurrent = EnumeratorType( ['FHC', 'RHC'], capitalization="all_upper" )
# note that DIS and Res are disabled because GENIE doesn't spit out
# cross-sections per nucleon for those channels...
InteractionChannel = EnumeratorType( ['Tot', 'QEl' ], capitalization=None ) 
#InteractionChannel = EnumeratorType( ['Tot', 'QEl', 'Res', 'DIS', ], capitalization=None )
InteractionCurrent = EnumeratorType( {'CC':1, 'NC':2}, capitalization="all" )
TargetNucleon = EnumeratorType( {'n': 2112, 'p': 2212}, capitalization="all_lower" )
# note that for now, since the flux calculation depends on the number of scattering centers
# (i.e., nuclei), and we only have a calculation of this for C12, C12 is currently
# the only option available...
TargetNucleus = EnumeratorType( ['C12',], capitalization=None )
TargetZ = EnumeratorType( {'C': 6,} )
#TargetNucleus = EnumeratorType( ['n', 'H1', 'H2', 'He4', 'C12', 'N14', 'O16', 'Al27', 'Si28', 'Cl35', 'Ar40', 'Ti48', 'Fe56', 'Pb207'], capitalization=None )
#TargetZ = EnumeratorType( {'H': 1, 'He': 2, 'C': 6, 'N': 7, 'O': 8, 'Al': 13, 'Si': 14, 'Cl': 17, 'Ar': 18, 'Ti': 22, 'Fe': 26, 'Pb': 82} )
NeutrinoFlavor = EnumeratorType( {'e': 12, 'mu': 14, 'tau': 16}, capitalization="all_lower" )
NeutrinoHelicity = EnumeratorType( {'PARTICLE': -1, 'ANTIPARTICLE': 1}, capitalization="all_upper" )
# use most recent official CCInclusive Ana processing
DefaultFileList = { \
    'FHC' : "/minerva/data/users/minervapro/New_CCInclusiveReco_POTbugFixed/playlist_summary/MC_minerva1_CCInclusiveReco.txt", \
    'RHC' : "/minerva/data/users/minervapro/New_CCInclusiveReco_POTbugFixed/playlist_summary/MC_minerva5_CCInclusiveReco.txt" \
    }
DefaultMeanPOTPerFile = { "FHC" : 2.192000E+20 / 440 , "RHC" : 8.975000E+20 / 1800 }
class FluxCalculator(object):
        N_E_BINS = 1000
        E_MIN = 0   # in GeV
        E_MAX = 100 # in GeV
        
        def __init__(self,
                hc=HornCurrent.FHC,
                curr=InteractionCurrent.CC,
                chan=InteractionChannel.Tot,
                nucleus=TargetNucleus.C12,
                nucleon=None,
                nu_flavor=NeutrinoFlavor.mu,
                nu_helicity=None,
                qgsp=False,
                testing=False,
                xsec_file="$MPARAMFILES/GENIE/spline_files/gxspl-nuclear-MINERVA_Full_v2126.root",
                filelist=None,
                pot_per_file=None,
                use_meta_tree=False,
                max_pot_testing=1e19,
                calc_errors=False,
                n_universes=None,
                use_ppfx=False
):
                self._horn_current = hc
                self._int_current = curr
                self._int_channel = chan
                self._target_nucleus = nucleus
                self._target_nucleon = nucleon
                self._nu_flavor = nu_flavor
                self._nu_helicity = nu_helicity
                
                self._qgsp = qgsp
                self._testing = testing
                self._xsec_filename = xsec_file
                self._xsec_file = None
                
                self._ntuple_filelist = filelist
                self._ntuple_chain = ROOT.TChain("Truth", "Truth")
                self._meta_chain = ROOT.TChain("Meta", "Meta")
                self._chain_is_loaded = False
                self._pot_per_file = pot_per_file
                self._use_meta_tree = use_meta_tree
                self._max_pot_testing = max_pot_testing
                
                self._calc_errors = calc_errors
                self._n_universes = n_universes
                if self._calc_errors and self._n_universes is None:
                        self._n_universes = 500
                self._use_ppfx = use_ppfx
                
                self._params = {}
                
                self.histos = {}
                
                self._additional_branches = []  # additional MC branches to enable
                
                # the event loop is written in C++ to make it faster...
                self._loop_obj = ROOT.FluxCalculatorLoop()
                
                self.Validate()
        
        def CalculateFlux(self):
                """ Actually calculates the flux.
                
                There are basically three steps:
                 (1) Convert the GENIE cross-section TGraph to a binned histogram
                     with the same binning as the event rate plots we'll make in (2).
                 (2) Plot the number of events (both unweighted and using CV reweighting)
                     using the true event classification from the MC truth.
                     Scale the number of events by the fiducial volume, the target
                     nucleus density, and the POT to get an event rate.
                 (3) Divide the event rates from (2) by the cross-section from (1)
                     to get both unweighted and CV-reweighted fluxes.
                
                """
                
                if not self.LoadChain():
                        return
                        
                if not self.LoadGENIEXsec():
                        return
                        
                self.CalculateParameters()
                # step 1: convert the relevant GENIE TGraph
                # to a histogram with the correct binning.
                # (this is just for the convenience of division
                # later on.)
                self.histos["x_sect_E"] = ROOT.TH1D("x_sect_E", self.histos["xsec"].GetTitle().split(";")[0] + ";E_{#nu} (GeV);#frac{d#sigma}{dE} (cm^{2} / %s)" % TargetNucleus[self._target_nucleus], self._params["n_E_bins"], self._params["E_min"], self._params["E_max"])
                for i in range(1, self._params["n_E_bins"]+1):
                        # get the value of the cross-section graph at the histogram bin center
                        self.histos["x_sect_E"].SetBinContent( i, self.histos["xsec"].Eval( self._params["E_bin_width"]*(i-0.5)) )
                        
                        # note that the cross-section graph contains no errors.
                        # rather than letting ROOT compute default (Gaussian) errors
                        # (which will be crazy, since these values are so small),
                        # we just don't allow errors on the cross-section to
                        # contribute to the flux error budget via this graph.
                        # we're dividing out exactly the same flux that was used
                        # to do the generation of these events anyway, so
                        # any uncertainties here are irrelevant.
                        self.histos["x_sect_E"].SetBinError( i, 0 )
                        
                # GENIE reports cross-sections in units of (x10^{-38})
                # which means we need to add the 10^{-38} here manually.
                self.histos["x_sect_E"].Scale(1e-38)
                # get the selection cuts and store them in the parameters
                # in case later review is necessary
                selection_cuts = self.GetSelectionCuts()
                self._params.update(selection_cuts)
                
                # step 2: loop over the events to calculate the event rates and error bands.
                # preparatory work: create histograms and error bands.
                # (use the MnvH1D machinery so we can do uncertainties properly.)
                error_bands = set([ "ppfx1_Total", "Flux_BeamFocus"]) if self._use_ppfx else set([ "Flux_Tertiary", "Flux_BeamFocus", "Flux_NA49", ])
#                for err in error_bands:
#                        self._ntuple_chain.SetBranchStatus("mc_wgt_%s" % err, 1)
                        
                for cut_name, cut in selection_cuts.items():
                        event_histo_name = "eventcount_E_%s" % cut_name
                        histogram = PlotUtils.MnvH1D(event_histo_name, "Event count (%s);E_{#nu} (GeV);Events / GeV" % cut_name, self._params["n_E_bins"], self._params["E_min"], self._params["E_max"])
                        self.histos[event_histo_name] = histogram
                        
                        hist_error_bands = set(histogram.GetErrorBandNames())
                        for band in error_bands - hist_error_bands:
                                histogram.AddVertErrorBand(band, self._n_universes)
                        # now, count the raw number of true events
                        # that are found within the fiducial volume...
                        print(("Calculating the flux for selection '%s'..." % cut_name))
                        print((" using selection cut:", cut))
                        self._ntuple_chain.Draw(">>evt_list_%s" % cut_name, cut)
                        evt_list = ROOT.gDirectory.Get("evt_list_%s" % cut_name)
                        evt_list.SetReapplyCut(True)
                        # call the event loop function from the C++ library
                        # (assume CV weighted unless specifically told otherwise?)
                        cvweighted = not("unweighted" in cut_name)
                        self._loop_obj.EventLoop(self._ntuple_chain, evt_list, histogram, "mc_incomingE", 0.001, cvweighted,)
                        # ... we'd like the event rate to be a true histogram,
                        # in which the bin size is irrelevant.
                        # to do that, we divide out the bin size.
                        histogram.Scale(1./self._params["E_bin_width"])
                                                
                        # ... then clone the event count histograms and scale them to get event rate
                        rate_histo_name = event_histo_name.replace("eventcount", "rate")
                        self.histos[rate_histo_name] = self.histos[event_histo_name].Clone(rate_histo_name)
                        self.histos[rate_histo_name].SetTitle("Event rate (%s);E_{#nu} (GeV);Events / %s / P.O.T. / GeV" % (cut_name,TargetNucleus[self._target_nucleus]))
                        
                        # the scale factor is (n_planes x targets per plane = # of target nuclei in fiducial volume) x (P.O.T.) x 1/1e4 (we store flux in units of m^2)
                        scale_factor = self._params["n_planes"] * self._params["n_scattering_centers"] * self._params["total_POT"] * 1e-4
                        self.histos[rate_histo_name].Scale( 1./(scale_factor) )
                
                        # step 3: divide the event rate by the cross-section to get the flux.
                        flux_histo_name = rate_histo_name.replace("rate", "flux")
                        flux_histo = self.histos[rate_histo_name].Clone(flux_histo_name)
                        if hasattr(flux_histo, "DivideSingle"):
#                                print "calling 'DivideSingle' on histogram", flux_histo.GetName()
                                flux_histo.DivideSingle(self.histos[rate_histo_name], self.histos["x_sect_E"])
                        else:
                                flux_histo.Divide(self.histos["x_sect_E"])
                        
                        nu_text = "#nu_{%s}" % NeutrinoFlavor[self._nu_flavor]
                        if self._nu_helicity == NeutrinoHelicity.ANTIPARTICLE:
                                nu_text = "#overbar{%s}" % nu_text
                        # note that the answer is ALWAYS in "/ GeV",
                        # no matter what the bin size is,
                        # because we already divided it out when creating
                        # the event rate histogram.
                        # (of course, if you rebin the result histogram, you'll still have to
                        #  scale it by whatever factor you changed the bins by.)
                        flux_histo.SetTitle("Flux (%s);Energy (GeV);%ss / m^{2} / P.O.T. / GeV" % (cut_name, nu_text))
                        
                        self.histos["flux_E_%s" % cut_name] = flux_histo
                        
#                        self.histos["cv_weights_%s" % cut_name] = ROOT.TH1D("cv_weights_%s" % cut_name, "CV weights", 40, 0.8, 1.2)
#                        self.histos["cv_weights_%s" % cut_name].SetDirectory(ROOT.gDirectory)
#                        self._ntuple_chain.Draw("mc_cvweight_total>>cv_weights_%s" % cut_name)
#                        self.histos["cv_weights_%s" % cut_name].SetDirectory(None)
        def FillErrorBands(self, histogram):
                """ Computes the flux/GENIE error bands on histogram 'histogram'. 
                
                Pretty slow since the event loop has to be done on the Python side.  
                This is alleviated as much as possible by leaving unnecessary branches
                disabled and using an event list.
                
                Please note that the errors are only correct for CV-reweighted samples! """
                
                if evt_list.GetN() < 1:
                        return
                        
                # only works for MnvH1Ds...
                if not hasattr(histogram, "GetErrorBandNames"):
                        return
                
                
        def _FillDefaultOptions(self):
                """ Applies the correct default options if the user hasn't specified them. """
                if self._nu_helicity is None:
                        if self._nu_flavor in (NeutrinoFlavor.mu, NeutrinoFlavor.e):
                                self._nu_helicity = NeutrinoHelicity.PARTICLE if self._horn_current == HornCurrent.FHC else NeutrinoHelicity.ANTIPARTICLE
                        else:
                                self._nu_helicity = NeutrinoHelicity.PARTICLE
                
                if self._target_nucleon is None and self._int_channel != InteractionChannel.Tot:
                        if self._int_current == InteractionCurrent.CC:
                                self._target_nucleon = TargetNucleon.n if self._nu_helicity == NeutrinoHelicity.PARTICLE else TargetNucleon.p
                        else:
                                self._target_nucleon = TargetNucleon.p
                
        def GetFiducialCut(self):
                """ Fiducial volume cut.  You can override this in a subclass if you like.
                    This version corresponds to the frozen detector CCQE one. """
                upstream_z_cut = "mc_vtx[2]>6117"    #this is just US of module 30
                downstream_z_cut = "mc_vtx[2]<8193"  #this is just DS of module 75
                hexagon_cut = "(TMath::Abs(mc_vtx[0])<850.0) && ((TMath::Abs(mc_vtx[1])<850.0/TMath::Sqrt(3.0))||TMath::Abs(mc_vtx[1]) < 850.0*2/TMath::Sqrt(3.0) - TMath::Abs(mc_vtx[0]/TMath::Sqrt(3.0)))"
                
                return "(%s)" % ") && (".join([upstream_z_cut, downstream_z_cut, hexagon_cut])
        def GetSelectionCuts(self):
                """ Event selection cuts.  You can override this in a subclass if you like.
                
                    By default, you get:
                      - fiducial volume cut (see GetFiducialCut())
                      - incoming neutrino cut (flavor and helicity)
                      - interaction channel and current cuts
                      - target nucleus (and, if applicable, nucleon) cuts
                    
                    The default returns two selections: CV-reweighted and unweighted.
                """
                #nan_cut = "mc_cvweight_total == mc_cvweight_total"
                
                interaction_cut = "mc_current == %d && mc_intType!=8" % self._int_current
                if self._int_channel != InteractionChannel.Tot:
                        interaction_cut += " && mc_intType == %d" % self._int_channel
                
                particle_cut = "mc_incoming == %d" % ( (-1 if self._nu_helicity == NeutrinoHelicity.ANTIPARTICLE else 1) * self._nu_flavor )
                
                # nucleus_cut = "mc_targetZ == %d" % TargetZ[[x for x in TargetNucleus[self._target_nucleus] if not x.isdigit()]]
                nucleus = ''.join([x for x in TargetNucleus[self._target_nucleus] if not x.isdigit()])
                nucleus_cut = "mc_targetZ == %d" % TargetZ[nucleus]

                if self._target_nucleon is not None:
                        nucleon_cut = "mc_targetNucleon == %d" % self._target_nucleon
#                else:
                nucleon_cut = None
                
                #cuts = [nan_cut, self.GetFiducialCut(), interaction_cut, particle_cut, nucleus_cut]
                cuts = [self.GetFiducialCut(), interaction_cut, particle_cut, nucleus_cut]
                if nucleon_cut is not None:
                        cuts.append(nucleon_cut)
                basic_cut = "(%s)" % ") && (".join(cuts)
                
                return { "unweighted": ROOT.TCut(basic_cut), "cvweighted": ROOT.TCut("mc_cvweight_total") * ROOT.TCut(basic_cut) }
                
        def CalculateParameters(self):
                """ Calculate some parameters needed in the plotting. """
                
                # note that self._params["total_POT"] is calculated in LoadChain()...
                
                self._params["n_E_bins"] = FluxCalculator.N_E_BINS # number of bins in plots vs. energy
                self._params["E_min"] = FluxCalculator.E_MIN # minimum energy in plots vs. energy
                self._params["E_max"] = FluxCalculator.E_MAX # maximum energy in plots vs. energy
    #these parameters must match the fiducial region used in GetFiducialCut
#                self._params["area"] = 2.503   # area (in m^2) of hexagon apothem 0.85 m
                self._params["area"] = 2.562   # area (in m^2) of hexagon apothem 0.86 m
                self._params["n_scattering_centers"] = 2.20722e+27  # this is only for Carbon, I'm afraid.  (it's scattering centers/(cm^2) for 1 plane)
                self._params["n_planes"] = (75-30+1)*2  # number of planes of interest (module 30-75 inclusive)
                self._params["E_bin_width"] = (self._params["E_max"]-self._params["E_min"]) / float(self._params["n_E_bins"])
                
                        
        def LoadChain(self, filelist=None, testing=None):
                """ Loads the ntuples to be used for calculating the flux
                from a file containing the list of ntuples. 
                
                Returns True if it does so successfully; otherwise, returns False. """
                if self._chain_is_loaded:
                        return
                filelist = filelist or self._ntuple_filelist
                testing = testing if testing is not None else self._testing
                
                self._ntuple_chain.Reset()
                
                # by default (if nothing else is specified),
                # we attempt to find the appropriate file list
                # for the secondary interaction model and
                # horn current in the configuration.
                if filelist is None:
                        
                        if self._qgsp:
                                filelist = "$FLUXROOT/file_lists/%s_%s_filelist.txt" % ("QGSP" if self._qgsp else "FTFP", HornCurrent[self._horn_current])
                                # these are due to the GENIE ray-tracing bug.
                                # the FHC sample is full detector and the RHC one is frozen
                                # (so they have different masses & volumes and the 
                                self._pot_per_file = 9.809e16 if self._horn_current == HornCurrent.FHC else 9.488e16
                        else:
                                self._pot_per_file = DefaultMeanPOTPerFile[ HornCurrent[self._horn_current] ]
                                filelist = DefaultFileList[ HornCurrent[self._horn_current] ]
                        print(("Selected default file list: '%s'" % filelist))
                        print(("POT per file:", self._pot_per_file))
                filelist = os.path.expandvars(filelist)
                try:
                        with open(filelist) as f:
                                i = 0
                                for line in f:
                                        line = line.strip()
                                        if len(line) == 0 or line[0] == '#':
                                                continue
                                        self._ntuple_chain.Add(line)
                                        if self._use_meta_tree:
                                                self._meta_chain.Add(line)
                                
                                        i += 1
                                        if self._testing:
                                                # If using meta tree, defer POT check until later
                                                if not self._use_meta_tree and self._pot_per_file is not None:
                                                        if i * self._pot_per_file >= self._max_pot_testing:
                                                                break
                                        # if self._testing and i*self._pot_per_file >= self._max_pot_testing:
                                        #         break
                except IOError:
                        print("Couldn't open filelist: '%s'.  No flux will be made." % filelist, file=sys.stderr)
                        return False
                        
                if self._use_meta_tree:
                        self._params["total_POT"] = self.GetPOTFromMetaTree()
                elif self._pot_per_file is not None:
                        self._params["total_POT"] = i * self._pot_per_file
                else:
                        raise ValueError("POT per file is None. Either set --mean_pot_per_file or use --use_meta_tree.")
                # else:
                #         self._params["total_POT"] = i*self._pot_per_file
                print(("Loaded %d ntuples (for a total of %g POT)." % (i, self._params["total_POT"])))
                
                # disable all but the branches of interest...
                self._ntuple_chain.SetBranchStatus("*", 0)
                
                branches_to_enable = set([ "mc_incomingE", "mc_intType", "mc_current", "mc_incoming", "mc_targetZ", "mc_targetNucleon", "mc_vtx*", "mc_cvweight_total", "mc_ppfx1_cvweight", "mc_hornCurrent_cvweight"]) | set(self._additional_branches)
                list(map(self._ntuple_chain.SetBranchStatus, branches_to_enable, [1,]*len(branches_to_enable)))
                
                self._chain_is_loaded = True
                return True
                
        def LoadGENIEXsec(self, xsec_file=None):
                """ Loads the GENIE cross-section from the specified ROOT file. """
                xsec_filename = xsec_file or self._xsec_filename
                
                self._xsec_file = ROOT.TFile(os.path.expandvars(xsec_filename))
                if not self._xsec_file.IsOpen():
                        raise OSError("Can't open file: '%s'" % xsec_filename)
                        
                # the GENIE cross-sections file is structured as follows:
                # + nutype_target  (<--- folder)
                #     +----> chan_curr_nucleon  (<--- histogram)
                #     +----> chan_curr_nucleon
                # + nutype_target
                #     +----> chan_curr_nucleon
                #     +----> chan_curr_nucleon
                # etc.
                xsec_folder_name = "nu_%s%s_%s".lower() % (NeutrinoFlavor[self._nu_flavor], "_bar" if self._nu_helicity == NeutrinoHelicity.ANTIPARTICLE else "", TargetNucleus[self._target_nucleus])
                xsec_folder = self._xsec_file.Get(xsec_folder_name)
                if not xsec_folder:
                        print("Couldn't find any cross-sections for this combination of neutrino flavor+helicity and target nucleus.  I looked for: '%s'" % xsec_folder_name, file=sys.stderr)
                        return False
                
                # for total cross-section, we don't want to restrict
                # to one particular nucleon (p or n).
                # the others are sorted that way in the GENIE file.
                if self._int_channel != InteractionChannel.Tot:
                        xsec_name = "%s_%s_%s" % (InteractionChannel[self._int_channel], InteractionCurrent[self._int_current], TargetNucleon[self._target_nucleon])
                else:
                        xsec_name = "%s_%s" % (InteractionChannel[self._int_channel], InteractionCurrent[self._int_current])
                xsec_name = xsec_name.lower()
                
                self.histos["xsec"] = xsec_folder.Get(xsec_name)
                if not self.histos["xsec"]:
                        print("Couldn't find the cross-sections for this process and nucleon.  I looked for: '%s' (in folder: '%s')" % (xsec_name, xsec_folder_name.GetName()), file=sys.stderr)
                        return False
                else:
                        print(("Loaded cross-section: %s / %s" % (xsec_folder.GetName(), xsec_name)))
                
                # reformat to give a bit more information
                self.histos["xsec"].SetTitle( "GENIE cross-section: %s / %s;E_{#nu} (GeV);#frac{d#sigma}{dE} (10^{-38} cm^{2} / GeV / %s" % (xsec_folder.GetName(), xsec_name, TargetNucleus[self._target_nucleus]) )
                self.histos["xsec"].SetName("GENIE_xsec")
                
                return True
        def GetPOTFromMetaTree(self):
                total_pot = 0.0
                for oneFile in self._meta_chain:
                        total_pot += oneFile.POT_Total
                return total_pot
        def PrintConfig(self):
                print(("  Horn current mode:", HornCurrent[self._horn_current]))
                print(("  Interaction current:", InteractionCurrent[self._int_current]))
                print(("  Interaction channel:", InteractionChannel[self._int_channel]))
                print(("  Target nucleus:", TargetNucleus[self._target_nucleus]))
                print(("  Target nucleon:", TargetNucleon[self._target_nucleon] if self._target_nucleon is not None else "auto-selected"))
                print(("  Neutrino flavor:", NeutrinoFlavor[self._nu_flavor]))
                print(("  Neutrino helicity:", NeutrinoHelicity[self._nu_helicity] if self._nu_helicity is not None else "auto-selected"))
                print(("  Secondary interaction model:", ("QGSP" if self._qgsp else "FTFP")))
                print(("  GENIE cross-section file:", self._xsec_filename))
                print(("  Ntuple file list:", self._ntuple_filelist if self._ntuple_filelist is not None else "auto-selected"))
                print(("  POT per file:", self._pot_per_file if self._pot_per_file is not None else "auto-selected"))
                print(("  Calculate errors:", self._calc_errors))
                if self._calc_errors:
                        print(("  Number of universes for 'many universes' flux errors calculation:", self._n_universes))
                print()
        def SaveHistosToDisk(self, filename, save_mnvh1ds=False):
                """ Saves the output flux histograms to the specified file.
                
                Be sure to call CalculateFlux() (or put something else in
                the 'histos' attribute of the FluxCalculator) before calling this
                or you'll wind up with an empty ROOT file..."""
                outfile = ROOT.TFile(filename, "RECREATE")
                if not outfile.IsOpen():
                        print("ERROR: couldn't open the output file ('%s').  Histograms won't be saved..." % filename, file=sys.stderr)
                        return
                print ("Saving plots...")
                for i, histo in enumerate(self.histos.values()):
                        print(("  saving histogram %d/%d with name:" % (i+1, len(self.histos)), histo.GetName()))
                        if hasattr(histo, "GetCVHistoWithError"):
                                if save_mnvh1ds:
                                        histo.Write()
                        
                                have_bands = len(histo.GetErrorBandNames()) > 0
                                if have_bands:
                                        print(("   (computing the errors... "), end=' ')
                                        sys.stdout.flush()
                                histo = histo.GetCVHistoWithError()
                                ROOT.SetOwnership(histo, False)        # the MnvH1D puts this into ROOT's garbage collector.  don't let it get deleted when my reference to it expires or we'll get double deletes.
                                if have_bands:
                                        print (" ... done.)")
                        histo.Write()
                
                for p, val in self._params.items():
                        if isinstance(val, numbers.Number):
                                param = ROOT.TParameter("double")(p, val)
                        else:
                                param = ROOT.TObjString(str(val))
                        param.Write()
                
                outfile.Write()
                outfile.Close()
                
                print ("... saving done.")
        def Validate(self):
                self._FillDefaultOptions()
        
                attrs_to_validate = {
                        "_horn_current": HornCurrent,
                        "_int_current": InteractionCurrent,
                        "_int_channel": InteractionChannel,
                        "_target_nucleus": TargetNucleus,
                        "_nu_flavor": NeutrinoFlavor,
                }
                
                for attr, enum in attrs_to_validate.items():
                        val = getattr(self, attr)
                        if val not in enum:
#                                print enum._attrs
#                                print enum._attrs_reversed
                                raise TypeError("Invalid value for attribute '%s': %s" % (attr, val))
                # for now, we can only handle C12...
                if self._target_nucleus != TargetNucleus.C12:
                        raise ValueError("For now, this tool can only calculate the flux using C12 as the target nucleus...")
                        
                if self._ntuple_filelist is not None and self._pot_per_file is None and not self._use_meta_tree:
                        raise ValueError("If you manually specify the ntuple filelist, you must also specify the POT per file!")
                        
                if self._nu_helicity is not None and self._nu_helicity not in NeutrinoHelicity:
                        raise TypeError("Invalid value for attribute '_nu_helicity': %s" % self._nu_helicity)
                if self._xsec_file is None:
                        self.LoadGENIEXsec()
                        
# if this is invoked from the command line,
# try to parse the arguments
def Bootstrap():
        parser = optparse.OptionParser(usage="usage: %prog [options]")
        parser.add_option("--RHC",
                dest="use_rhc",
                action="store_true",
                help="Use reverse horn current mode (RHC) instead of forward horn current (FHC).  FHC is the default.  If you supply both --FHC and --RHC, the trailing argument is the one that takes precedence.", 
                default=False
        )
        
        parser.add_option("--FHC",
                dest="use_rhc",
                action="store_false",
                help="Use forward horn current (FHC) mode.  This is on by default; the option is supported in case you wish to be explicit.  If you supply both --FHC and --RHC, the trailing argument is the one that takes precedence.")
        
        parser.add_option("--qgsp_flux",
                dest="use_qgsp",
                action="store_true",
                help="Use the older QGSP flux instead of FTFP.  FTFP is the default.",
                default=False
        )
        
        parser.add_option(
                "-t",
                "--test",
                dest="test",
                action="store_true",
                help="Use a smaller sample for testing purposes.  Default: %default.",
                default=False
        )
        
        parser.add_option("--channel",
                dest="channel",
                help="The interaction channel to use.  Choices: " + ", ".join(list(InteractionChannel.keys())) + ".  Note that not every channel has a cross-section for every flux.  Default: %default.  (Note that you will rarely want to use anything besides 'Tot' unless you know what you're doing.)",
                default="Tot"
        )
        
        parser.add_option("--current",
                dest="current",
                help="The interaction current to use.  Choices: " + ", ".join(list(InteractionCurrent.keys())) + ".  Default: %default.",
                default='CC'
        )
        parser.add_option("--nu_flavor",
                dest="flavor",
                help="Flavor of the neutrino you want the flux for.  Choices: " + ", ".join(list(NeutrinoFlavor.keys())) + ".  Default: '%default'.",
                default="mu",
        )
        parser.add_option("--nu_helicity",
                dest="helicity",
                action="store",
                type=int,
                help="Neutrino helicity.  Choices: +1 = right-handed = 'antineutrino', -1 = left-handed = 'neutrino'.  Default: 'right sign' according to horn current for nu_mu, nu_e; -1 for nu_tau.",
                default=None,
        )
        
        parser.add_option("--target_nucleus",
                dest="target_nucleus",
                help="Target nucleus to use.  Choices: " + ", ".join(list(TargetNucleus.keys())) + ".  (Note that not every nucleus has a cross-section for every flux.)  Default: '%default'.",
                default="C12",
        )
        
        parser.add_option("--target_nucleon",
                dest="target_nucleon",
                help="Target nucleon within the nucleus.  Choices: " + ", ".join(list(TargetNucleon.keys())) + ".  (Note that both choices will not make sense for certain combinations of horn current, process, and sign selection: e.g., neutrino CCQE on p is nonexistent.)  The default is to use both types of nucleons for the inclusive channel, and to select the 'typical' charge-current nucleon for CC events (e.g., 'n' for CCQE neutrino) and p for NC events when using a specific channel.  (See --channel, --current, --nu_helicity.)",
                default=None,
        )
        
        parser.add_option(
                "-f",
                "--xsec_file",
                dest="xsec_file",
                help="The input ROOT file containing GENIE's cross-sections.  Default: '%default'.",
                default="$MPARAMFILES/GENIE/spline_files/gxspl-nuclear-MINERVA_Full_v2126.root"
        )
        
        parser.add_option(
                "--filelist",
                dest="filelist",
                help="The file containing the list of MINERvA MC analysis ntuples you wish to use to get the flux.  The default is to choose the appropriate one from the file_list/ directory based on the horn current and flux model you gave as options.",
                default=None,
        )
        
        parser.add_option(
                "--mean_pot_per_file",
                dest="pot",
                help="Mean number of POT in one file. It is better to get this from the Meta tree than to specify it here. If you use an automatically selected filelist, this value will be automatically calculated.",
                type=float,
                default=None
        )
        parser.add_option(
                "--use_meta_tree",
                dest="meta",
                help="Use the Meta TTree to calculate the total POT based on a file list. Will ignore the --mean_pot_per_file",
                action="store_true",
                default=False
        )
        parser.add_option(
                "--max_pot_testing",
                dest="max_pot_testing",
                help="Maximum POT exposure to use in testing mode (if --pot_per_file is more than this, only 1 file will be used).  Default: %default.",
                type=float,
                default=1e19
        )
        
        parser.add_option(
                "--calc_errors",
                dest="do_errors",
                help="Calculate the flux errors.  THIS WILL TAKE A LONG TIME (though the time scales with '--n_universes').  Default: %default.",
                action="store_true",
                default=False,
        )
        
        parser.add_option(
                "--n_universes",
                dest="n_universes",
                help="Number of universes to use in flux error calculation (see '--calc_errors').  Default: %default.",
                type=int,
                default=500,
        )
        parser.add_option(
                "--use_ppfx",
                dest="use_ppfx",
                help="Use ppfx flux instead of NA49, Tertiary. BeamFocus used in either case.  Default: %default.",
                action="store_true",
                default=True,
        )
        
        parser.add_option(
                "-o",
                "--output_file",
                dest="output_file",
                help="The path where you want to save the ROOT file containing the output flux histograms.  Default: '%default'.",
                default="$PWD/flux_histos.root",
        )
        
        parser.add_option(
                "-n",
                "--no_action",
                dest="nop",
                action="store_true",
                help="Don't actually do anything; just print a summary of the configuration you chose.",
                default=False
        )
        (options, arguments) = parser.parse_args()
        
        # the FluxCalculator instance actually does the work
        try:
                calculator = FluxCalculator(
                        hc=HornCurrent.RHC if options.use_rhc else HornCurrent.FHC,
                        curr=InteractionCurrent[options.current],
                        chan=InteractionChannel[options.channel],
                        nucleus=TargetNucleus[options.target_nucleus],
                        nucleon=options.target_nucleon,
                        nu_flavor=NeutrinoFlavor[options.flavor],
                        nu_helicity=options.helicity,
                        qgsp=options.use_qgsp,
                        testing=options.test,
                        xsec_file=options.xsec_file,
                        filelist=options.filelist,
                        pot_per_file=options.pot,
                        use_meta_tree=options.meta,
                        max_pot_testing=options.max_pot_testing,                        
                        calc_errors=options.do_errors,
                        n_universes=options.n_universes,
                        use_ppfx=options.use_ppfx
                )
        except OSError as e:
                print("Couldn't open your cross-section file: '%s'" % options.xsec_file, file=sys.stderr)
                sys.exit(1)
        except (AttributeError, TypeError, ValueError) as e:
                print("There was an error trying to process your options: '%s'" % e, file=sys.stderr)
                parser.print_help()
                sys.exit(1)
                
        if options.nop:
                calculator.PrintConfig()
                sys.exit(0)
                
        return calculator, options.output_file
                
if __name__ == "__main__":
        print("test")
        calculator, output_file = Bootstrap()
        ROOT.gROOT.SetBatch(True)
        print ("\nCalculating the flux with the following options:")
        calculator.PrintConfig()
        if calculator._testing:
                print ("Note: using testing mode (max of 1e19 POT or 1 file with >1e19 POT)...")
        calculator.CalculateFlux()
        calculator.SaveHistosToDisk(filename=output_file, save_mnvh1ds=True)
        print ("\nAll done.  Bye.")