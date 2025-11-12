"""
    SystematicsUniverse.py:
    Define systematics shift universe. 
"""

import math
import os,sys
import random
import itertools
import operator
import random
import re

import ROOT
import PlotUtils
import ctypes

from config import SystematicsConfig
from tools import Utilities
from collections import OrderedDict
from tools import pcweight
#from .pcweight import GetModelWeight,MyWeighter
from tools import AvailableECorrection


OneSigmaShift = [-1.0,1.0]
M_e = 0.511 # MeV
M_mu =  105.658 #MeV
M_e_sqr = M_e**2
M_mu_sqr = M_mu**2

#Alias:
#ROOT.PythonMinervaUniverse = PlotUtils.DefaultCVUniverse

#MethodName: UpperCaseCamel
#universeShortName: lower_case_underscore

# The base universe, define functions to be used for all universes.
class CVUniverse(ROOT.PythonMinervaUniverse):
    is_pc = False
    def __init__(self, chain, nsigma = None): #nsigma is None for data because we don't shift data
        #ROOT.PythonMinervaUniverse.__init__(self, chain, 0 if nsigma is None else nsigma)
        super(CVUniverse,self).__init__(chain,0 if nsigma is None else nsigma)
        self.weight = None
        self.tuning_weight = None
        self.InitWithoutSuper(chain,nsigma)

    def InitWithoutSuper(self,chain,nsigma):
        self.nEntries = chain.GetEntries()
        self.mc = nsigma is not None
        self.chain=chain
        self.pcweighter = pcweight.MyWeighter

    # magic function that allow direct access to TBranch by Universe.TBranch, but prohibit any python magic variable being generated.
    def __getattr__(self,attrName):
        def makelist(x):
            return list(makelist(_) for _ in x) if hasattr(x, "__getitem__") else x

        if attrName == "kin_cal" or attrName == "classifier":
            #the two attrs can't be in chainwrapper 
            raise AttributeError()

        if attrName ==  "pcweighter":
            print (self.ShortName())

        # determine the dimension of leaf.
        try:
            return makelist(self.LeafGetters[attrName](attrName))
        except KeyError:
            #havn't see this function. Got to figure out how to get it.
            i = ctypes.c_int(1)
            leaf = self.chain.GetLeaf(attrName,-1)
            type = leaf.GetTypeName()
            if leaf.GetLeafCounter(i) or i.value > 1:
                dimension=1
            else:
                dimension=type.count("vector")
            # Assuming type names are Double_t,Int_t,Bool_t,vector<vector<Double_t>>, use regex to match the type
            regex="(<|^)([A-Za-z]+)(_t|>|$)"
            type = re.search(regex,type).group(2).capitalize()

            self.LeafGetters[attrName] = getattr(self, "Get{}{}".format("" if dimension==0 else ("Vec" if dimension==1 else "VecOfVec"),type))

        return makelist(self.LeafGetters[attrName](attrName))

    # return text when you print this class object. must reimplement or __getattr__ will hijack the method resolution
    def __repr__(self):
        return self.LatexName()+" Universe at sigma: "+ str(self.GetSigma())

    @property
    def n_prongs(self):
        return self.GetInt("n_prongs")

    # keep a reference to kinematicsCalulator and eventClassifier such that other function can pull info from universe
    def LoadTools(self,kinematicCalculator,eventClassifier):
        self.kin_cal = kinematicCalculator
        self.classifier = eventClassifier

    def SetCVEntry(self,n_entry):
        #go to another event, discard weight calculated.
        self.weight = None
        self.tuning_weight = None
        CVUniverse.LLR = None
        super(CVUniverse,self).SetEntry(n_entry)

    def SetLeptonType(self):
        if abs(SystematicsConfig.AnaNuPDG) == 12: 
            self.LeptonTheta = self.ElectronTheta
            self.LeptonEnergy = self.ElectronEnergy
            self.LeptonP3D = self.ElectronP3D
            self.M_lep_sqr =  M_e_sqr
            self.GetCorrection = self.GetLeakageCorrection
            return True
        else:
            self.LeptonTheta = self.GetThetamu
            self.LeptonEnergy  = self.GetEmu
            self.LeptonP3D = self.GetPmu
            self.M_lep_sqr =  M_mu_sqr
            self.GetCorrection = self.GetNuEFuzz
            return False

    @property
    def nsigma(self):
        return self.GetSigma() if self.mc else None

    def GetWeight(self,bkgtuning=True):
        if self.weight is None:
            if self.nsigma is None:
                self.weight = 1
            else:
                self.weight = self.GetModelWeight() if self.is_pc else self.GetStandardWeight()
        if bkgtuning and self.nsigma is None:
            if self.tuning_weight is None:
                self.tuning_weight = self.pcweighter.GetWeight(self)
            return self.weight *self.tuning_weight
        else:
            return self.weight

    def GetStandardWeight(self):
        weight = 1.0
        weight *= self.GetGenieWeight()

        ### take care of flavor swapped sample below ###
        pdg = self.mc_incoming
        newpdg = pdg
        if SystematicsConfig.USE_SWAPPED:
            if pdg == 12:
                newpdg = 14
            elif pdg == -12:
                newpdg = -14
            elif pdg == 14:
                newpdg = 12
            elif pdg == -14:
                newpdg = -12
            else:
                print("Not neutrino PDG? {}".format(pdg))
        pdg = newpdg
        weight *= self.GetFluxAndCVWeight(self.mc_incomingE*1e-3,pdg)

        weight *= self.GetLowRecoil2p2hWeight()
        weight *= self.GetRPAWeight()
        weight *= self.GetMyLowQ2PiWeight() # using MENU1PI for Aaron's result
        weight *= self.GetMyMinosEfficiencyWeight()

        ### MnvTune v4.3.1 block below ###
        #weight *= self.GetTargetMassWeight() # 
        weight *= self.GetCOHPionWeight() # Alejandro Ramirez's Result
        if self.mc_intType == 4:
            weight *= 1.4368 # diffractive coherent reweight from Aaron Mislivec

        return weight


    #def GetPCWeight(self):
    #    return 1.0
    def GetMyMinosEfficiencyWeight(self):
         if self.HasNoBackExitingTracks:
             return 1
         else:
             return self.GetMinosEfficiencyWeight()

    def GetMyLowQ2PiWeight(self, channel = SystematicsConfig.LowQ2PiWeightChannel):
        if channel is None:
            return 1
        else:
            return super(CVUniverse,self).GetLowQ2PiWeight(channel.upper())

    def GetModelWeight(self):
        w = GetModelWeight(self, 'Eel',"Pi0") # (Ee, Theta, ETh), (PCElectron, PCPhoton, PCPi0) 
        return w

    def GetCOHPionWeight(self):
        def GetTrueHighEpi():
            nFSpi = self.mc_nFSPart
            pionE = -1.0
            pionKE = -1.0
            for i in range(nFSpi):
                try:
                    pdg = self.mc_FSPartPDG[i]
                except:
                    continue

                if pdg != 211:
                    continue
                energy = self.mc_FSPartE[i]
                mass = 139.569
                tpi = energy - mass
                if tpi >= pionKE:
                    pionKE = tpi
                    pionE = energy
            return pionE

        def GetTrueAngleHighTpi():
            nFSpi = self.mc_nFSPart
            angle = -9999 #WRTbeam and in degrees
            pionKE = 0.0
            idk = -9999
            for i in range(nFSpi):
                try:
                    pdg = self.mc_FSPartPDG[i]
                except:
                    continue

                if pdg != 211:
                    continue
                energy = self.mc_FSPartE[i]
                mass = 139.569
                tpi = energy - mass
                if tpi >= pionKE:
                    pionKE = tpi
                    pimomentumvec = ROOT.TVector3(self.mc_FSPartPx[i], self.mc_FSPartPy[i],self.mc_FSPartPz[i])
                    r = ROOT.Math.RotationX(SystematicsConfig.BEAM_ANGLE)
                    s=r(pimomentumvec)
                    deg_wrtb = s.Theta()
                    angle = deg_wrtb

            #Making sure angle is only between 0 and pi
            if angle < 0.0:
                angle = -1.0*angle
            if angle > math.pi:
                angle = 2.0*math.pi - angle
            return angle*180./math.pi # Degrees
       
        weight = 1.0
        if self.mc_intType != 4:
           return 1.0
        elif self.mc_intType == 4:
           npi = self.mc_FSPartPDG.count(211)
           if npi == 0:
               return 1.0
           angle = GetTrueAngleHighTpi()
           epi = GetTrueHighEpi()/1000. # This is supposed to be the energy of the pion!!!! 
           if epi < 0:
               return 1.0
           else:
               weight *= self.GetCoherentPiWeight(angle, epi) #Inputs are in Degrees and GeVi
               return weight
        else:
            return 1.0

    #################### functions to be overide in systematics shift universes #################

    def ElectronEnergyRaw(self):
        return self.GetVecElem("prong_part_E",0,3)

    @Utilities.decorator_ReLU
    def ElectronEnergy(self):
        return self.ElectronEnergyRaw() + self.GetEMEnergyShift() #+ self.GetLeakageCorrection()

    def ElectronP3D(self):
        electronp = self.GetVecOfVecDouble("prong_part_E")
        #print([electronp[0][i] for i in range(4)])
        scale = self.GetEMEnergyShift()/electronp[0][3] if (electronp[0][3]>0) else 0
        p = ROOT.Math.XYZVector(*tuple(list(electronp[0])[:3]))*(1+scale)
        r = ROOT.Math.RotationX(SystematicsConfig.BEAM_ANGLE)
        s=r(p)
        if s is None:
            print (p,r,tuple(list(electronp[0])[:3]))
        return s

    def ElectronProtonAngle(self):
        electronp = self.GetVecOfVecDouble("prong_part_E")
        scale = self.GetEMEnergyShift()/electronp[0][3] if (electronp[0][3]>0) else 0
        electron_p = ROOT.TVector3(*tuple(list(electronp[0])[:3]))*(1+scale)
        
        protonX = self.MasterAnaDev_proton_Px_fromdEdx
        protonY = self.MasterAnaDev_proton_Py_fromdEdx
        protonZ = self.MasterAnaDev_proton_Pz_fromdEdx
        proton_p = ROOT.TVector3(*tuple([protonX,protonY,protonZ]))
        dotprod = proton_p.Dot(electron_p)
        mag = proton_p.Mag() * electron_p.Mag()
        if mag > 0:
            dotprod /= mag
        if protonX > -9998 and protonY > -9998 and protonZ > -9998:
            return math.degrees(math.acos(dotprod))
        else:
            return -5

    def ElectronTheta(self):
        return self.ElectronP3D().Theta()

    def PrimaryProtonTheta(self):
        theta = self.MasterAnaDev_proton_theta
        return math.degrees(theta) if theta>=0 else -1

    def ConeInsideE(self): 
        e = self.GetVecElem("prong_TrimmedCaloEnergy",0)/1e3
        psi = self.GetVecElem("Psi",0)
        #print e, 'inside', psi, 'psi'
        return e*psi

    def ConeOutsideE(self):
        localVec = list(self.GetVecDouble("ExtraEnergy"))
        e = sum(localVec)
	#print e, 'out'
        return e 

    @staticmethod
    def Scale4VectorMomentum(vet,scale):
        return vet


    def PrimaryProtonScore(self):
        score = self.MasterAnaDev_proton_score 
        return score if score>=0 else -1

    def PrimaryProtonTheta(self):
        theta = self.MasterAnaDev_proton_theta
        return math.degrees(theta) if theta>=0 else -1


    # =============== collcetion of all recoil energy stuff.========
    def GetEAvail(self):
        if self.GetAnaToolName() == "MasterAnaDev" or self.GetAnaToolName() == "CCQENu":
            return self.blob_recoil_E_tracker+self.blob_recoil_E_ecal 
        else:
            correction = AvailableECorrection.GetSmoothVisECorrection(self)*1e3 #Eavail correction 
            #return self.recoile_passive_tracker+self.recoile_passive_ecal
            return self.recoile_passive_tracker+self.recoile_passive_ecal+correction

    def GetNuEFuzz(self):
        try:
            fuzz = self.blob_nuefuzz_E_tracker+self.blob_nuefuzz_E_ecal
        except:
            fuzz = 0
        #return 0
        return max(0,fuzz) * SystematicsConfig.AVAILABLE_E_CORRECTION


    def GetLeakageCorrection(self): 
        return SystematicsConfig.LEAKAGE_CORRECTION(self.ElectronEnergyRaw()) + (SystematicsConfig.LEAKAGE_BIAS if self.nsigma is None else 0) 

    def Pi0_Additional_Leakage(self):
        return random.uniform(0,20) if self.mc and (self.mc_intType == 4 or (self.mc_current==2 and self.mc_intType==10)) else 0

    # A function wrapper that maps negative values to zero
    # returning MeV
    @Utilities.decorator_ReLU
    def AvailableEnergy(self):
        #return self.GetEAvail()* SystematicsConfig.AVAILABLE_E_CORRECTION
        if abs(SystematicsConfig.AnaNuPDG) == 12:
            return self.GetEAvail()* SystematicsConfig.AVAILABLE_E_CORRECTION -self.GetCorrection()
        else:
            return self.GetEAvail()

    @Utilities.decorator_ReLU
    def RecoilEnergy(self): 
        result = None
        if self.GetAnaToolName() == "MasterAnaDev" or self.GetAnaToolName() == "CCQENu":            
            return sum(self.__getattr__("blob_recoil_E_{}".format(i)) for i in ["tracker","ecal","hcal","od","nucl"])
        else:
            result = sum(self.__getattr__("recoile_passive_{}".format(i)) for i in ["tracker","ecal","hcal","odclus","nucl"])
        return result - self.GetCorrection() 

    def GetFuzzCorrection(self):
        return 0

    def GetEMEnergyShift(self):
        return self.prong_ECALCalibE[0]*SystematicsConfig.EM_ENERGY_SCALE_SHIFT_ECAL if self.nsigma is None else 0

    def dedxDelta(self):
        dedx = self.GetVecDouble("dEdXs")
        dedx_dz  = self.GetVecDouble("dEdXs_dz")
        dedx_proj  = self.GetVecDouble("dEdXs_proj")
        total = 0
        outer = 0
        if dedx.size() == 0:
            return 0.9999
        for i in range(0,dedx.size()):
            if dedx_proj[i]>400:
                outer+=dedx[i]
            total+=dedx[i]
        return outer/total 

class FluxUniverse:
    def __init__(self, chain, universe_number):
        self.universe = ROOT.PlotUtils.FluxUniverse[ROOT.PythonMinervaUniverse](chain, 1, universe_number)
        self.cv_universe = CVUniverse(chain, 1)
        if universe_number >= 100:
            raise ValueError(f"Universe index {universe_number} out of range. Max allowed for LE is 99.")

    def GetWeight(self,test=None):
        return self.GetStandardWeight()

    def GetStandardWeight(self):
        weight = self.cv_universe.GetStandardWeight()
        weight *= self.universe.GetWeightRatioToCV()
        return weight

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.universe, attr, getattr(self.cv_universe, attr, None))

    # @staticmethod
    # def GetSystematicsUniverses(chainwrapper):
    #     n = SystematicsConfig.NUM_FLUX_UNIVERSE
    #     universes = []
    #     for i in range(n):
    #         try:
    #             u = FluxUniverse(chainwrapper, i)
    #             # Smoke test with explicit values (avoid reading from the tree)
    #             _ = u.universe.GetFluxAndCVWeight(Enu=1.0, nu_pdg=14)  # or 12
    #             universes.append(u)
    #         except Exception as e:
    #             print(f"[WARNING] Skipping flux universe {i}: {e}")
    #             break
    #     return universes

    @staticmethod
    def GetSystematicsUniverses(chain):
        print(f"[FluxUniverse:GetSystematicsUniverses] building {SystematicsConfig.NUM_FLUX_UNIVERSE} universes")
        return [FluxUniverse(chain,  i ) for i in range(0,SystematicsConfig.NUM_FLUX_UNIVERSE)]

class GenieUniverse():
    def __init__(self,chain,nsigma,universe_name):
        self.universe = ROOT.PlotUtils.GenieUniverse[ROOT.PythonMinervaUniverse](chain, nsigma, universe_name)
        self.cv_universe = CVUniverse(chain, nsigma)  

    def GetWeight(self,test=None):
        return self.GetStandardWeight()

    def GetStandardWeight(self):
        weight = self.cv_universe.GetStandardWeight()
        weight *= self.universe.GetGenieWeight()
        return weight

    def __getattr__(self, attr):
        """Redirect attribute access to GenieUniverse first, then CVUniverse"""
        return getattr(self.universe, attr, getattr(self.cv_universe, attr, None))

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [GenieUniverse(chain,i,j) for j in SystematicsConfig.GENIE_UNIVERSES for i in OneSigmaShift]

class GenieRvx1piUniverse():
    def __init__(self,chain,nsigma,universe_name):
        self.universe = ROOT.PlotUtils.GenieRvx1piUniverse[ROOT.PythonMinervaUniverse](chain, nsigma, universe_name)
        self.cv_universe = CVUniverse(chain, nsigma)  

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.universe, attr, getattr(self.cv_universe, attr, None))

    def GetWeight(self,test=None):
        return self.GetStandardWeight()

    def GetStandardWeight(self):
        weight = self.cv_universe.GetStandardWeight()
        weight *= self.universe.GetGenieWeight()
        return weight

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [GenieRvx1piUniverse(chain,i,j) for j in ["Rvn1pi","Rvp1pi"] for i in OneSigmaShift]

class GenieFaCCQEUniverse():
    def __init__(self,chain,nsigma,universe_number):
        self.universe = ROOT.PlotUtils.GenieFaCCQEUniverse[ROOT.PythonMinervaUniverse](chain, nsigma, universe_number)
        self.cv_universe = CVUniverse(chain, nsigma)  

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.universe, attr, getattr(self.cv_universe, attr, None))

    def GetWeight(self,test=None):
        return self.GetStandardWeight()

    def GetStandardWeight(self):
        weight = self.cv_universe.GetStandardWeight()
        weight *= self.universe.GetGenieWeight()
        return weight

    @staticmethod
    def GetSystematicsUniverses(chain ):
        if ROOT.PythonMinervaUniverse.UseZExpansionFaReweight():
            return [GenieFaCCQEUniverse(chain,1.0,i) for i in range(0,SystematicsConfig.NumZExpansionUniverses)]
        else:
            return [GenieFaCCQEUniverse(chain,i,-1) for i in OneSigmaShift]

class GenieNormCCResUniverse():
    def __init__(self,chain,nsigma):
        self.universe = ROOT.PlotUtils.GenieNormCCResUniverse[ROOT.PythonMinervaUniverse](chain, nsigma)
        self.cv_universe = CVUniverse(chain, nsigma)  

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.universe, attr, getattr(self.cv_universe, attr, None))

    def GetWeight(self,test=None):
        return self.GetStandardWeight()

    def GetStandardWeight(self):
        weight = self.cv_universe.GetStandardWeight()
        weight *= self.universe.GetGenieWeight()
        return weight

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [GenieNormCCResUniverse(chain,i) for i in OneSigmaShift]

class GenieMaResUniverse():
    def __init__(self,chain,nsigma):
        self.universe = ROOT.PlotUtils.GenieMaResUniverse[ROOT.PythonMinervaUniverse](chain, nsigma)
        self.cv_universe = CVUniverse(chain, nsigma)  

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.universe, attr, getattr(self.cv_universe, attr, None))

    def GetWeight(self,test=None):
        return self.GetStandardWeight()

    def GetStandardWeight(self):
        weight = self.cv_universe.GetStandardWeight()
        weight *= self.universe.GetGenieWeight()
        return weight

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [GenieMaResUniverse(chain,i) for i in OneSigmaShift]

class GenieMvResUniverse():
    def __init__(self,chain,nsigma):
        self.universe = ROOT.PlotUtils.GenieMvResUniverse[ROOT.PythonMinervaUniverse](chain, nsigma)
        self.cv_universe = CVUniverse(chain, nsigma)  

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.universe, attr, getattr(self.cv_universe, attr, None))

    def GetWeight(self,test=None):
        return self.GetStandardWeight()

    def GetStandardWeight(self):
        weight = self.cv_universe.GetStandardWeight()
        weight *= self.universe.GetGenieWeight()
        return weight

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [GenieMvResUniverse(chain,i) for i in OneSigmaShift]

# I don't think nu_e analysis need minos shift
# class MinosUniverse(ROOT.PlotUtils.GenieUniverse(ROOT.PythonMinervaUniverse), CVUniverse):

class Universe2p2h():
    def __init__(self,chain,universe_number):
        self.universe = ROOT.PlotUtils.Universe2p2h[ROOT.PythonMinervaUniverse](chain, 1, universe_number)
        self.cv_universe = CVUniverse(chain, 1)  

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.universe, attr, getattr(self.cv_universe, attr, None))

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [Universe2p2h(chain, i) for i in SystematicsConfig.UNIVERSES_2P2H]

class RPAUniverse():
    def __init__(self,chain, universe_number,q2_region):
        self.universe = ROOT.PlotUtils.RPAUniverse[ROOT.PythonMinervaUniverse](chain, 1, universe_number, q2_region)
        self.cv_universe = CVUniverse(chain, 1)  

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.universe, attr, getattr(self.cv_universe, attr, None))

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [RPAUniverse(chain, i,j) for j in SystematicsConfig.RPA_UNIVERSES for i in SystematicsConfig.RPA_UNIVERSES[j] ]

class ResponseUniverse():
    def __init__(self,chain,nsigma,name):
        self.universe = ROOT.PlotUtils.ResponseUniverse[ROOT.PythonMinervaUniverse](chain, nsigma, name)
        self.cv_universe = CVUniverse(chain, 1)  
        self.re = "^blob_.*_E_(tracker|ecal|od|nucl|hcal)$"

    def __getattr__(self, attrName):
        """Redirect attribute access to ResponseUniverse, then CVUniverse"""
        if re.match(self.re, attrName) is None:
            return getattr(self.universe, attrName, getattr(self.cv_universe, attrName, None))
        else:
            # Fetching the attribute and performing the custom calculation
            r = getattr(self.universe, attrName, None)
            shift = getattr(self.universe, f"{attrName}_{self.m_name}", None) * self.m_frac_unc
            return r + self.nsigma * shift

    def RecoilClustersCalorimetry(self,splined):
        return super(ResponseUniverse,self).RecoilClustersCalorimetry(splined)+self.GetRecoilShift()

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [ResponseUniverse(chain,i,j) for j in SystematicsConfig.RESPONSE_BRANCHES for i in OneSigmaShift]

# class LowQ2PionUniverse():
#     def __init__(self,chain,nsigma,channel):
#         self.universe = ROOT.PlotUtils.LowQ2PionUniverse[ROOT.PythonMinervaUniverse](chain, nsigma, channel)
#         self.cv_universe = CVUniverse(chain, nsigma)  
#         self.channel = channel.upper() if channel is not None else None

#     def __getattr__(self, attr):
#         """Redirect attribute access to FluxUniverse first, then CVUniverse"""
#         return getattr(self.universe, attr, getattr(self.cv_universe, attr, None))

#     def GetMyLowQ2PiWeight(self):
#         if self.nsigma == 0 or self.channel is None:
#             return 1
#         else :
#             return super(LowQ2PionUniverse,self).GetLowQ2PiWeight(self.channel)

#     @staticmethod
#     def GetSystematicsUniverses(chain):
#         cvshifts = [LowQ2PionUniverse(chain,i,SystematicsConfig.LowQ2PiWeightChannel) for i in OneSigmaShift] if SystematicsConfig.LowQ2PiWeightChannel is not None else []
#         return cvshifts

class LowQ2PionUniverse:
    def __init__(self, chain, nsigma, channel):
        # Use the 2-arg ctor that exists in your library
        self.universe = ROOT.PlotUtils.LowQ2PionUniverse[ROOT.PythonMinervaUniverse](chain, nsigma)
        self.cv_universe = CVUniverse(chain, nsigma)
        self.nsigma = nsigma
        # Store channel for use when asking for the weight
        self.channel = (channel or "").lower()  # keep lowercase: GENIE/PlotUtils often expect this

    def __getattr__(self, attr):
        # delegate to the C++ syst universe first, then CV
        try:
            return getattr(self.universe, attr)
        except AttributeError:
            return getattr(self.cv_universe, attr)

    def GetMyLowQ2PiWeight(self):
        # nsigma==0 or no channel => no shift
        if not self.channel or self.nsigma == 0:
            return 1.0
        # Call the method on the wrapped C++ universe
        return self.universe.GetLowQ2PiWeight(self.channel)

    @staticmethod
    def GetSystematicsUniverses(chain):
        # Build ±1σ universes only if a channel is configured
        if SystematicsConfig.LowQ2PiWeightChannel is None:
            return []
        ch = SystematicsConfig.LowQ2PiWeightChannel
        return [LowQ2PionUniverse(chain, i, ch) for i in OneSigmaShift]


class MuonAngleXResolutionUniverse():
    def __init__(self,chain,nsigma):
        self.universe = ROOT.PlotUtils.MuonAngleXResolutionUniverse[ROOT.PythonMinervaUniverse](chain, nsigma)
        self.cv_universe = CVUniverse(chain, nsigma)  

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.universe, attr, getattr(self.cv_universe, attr, None))

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [MuonAngleXResolutionUniverse(chain,i) for i in OneSigmaShift]

class MuonResolutionUniverse():
    def __init__(self,chain,nsigma):
        self.universe = ROOT.PlotUtils.MuonResolutionUniverse[ROOT.PythonMinervaUniverse](chain, nsigma)
        self.cv_universe = CVUniverse(chain, nsigma)  

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.universe, attr, getattr(self.cv_universe, attr, None))

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [MuonResolutionUniverse(chain,i) for i in OneSigmaShift]

class MuonAngleYResolutionUniverse():
    def __init__(self,chain,nsigma):
        self.universe = ROOT.PlotUtils.MuonAngleYResolutionUniverse[ROOT.PythonMinervaUniverse](chain, nsigma)
        self.cv_universe = CVUniverse(chain, nsigma)  

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.universe, attr, getattr(self.cv_universe, attr, None))

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [MuonAngleYResolutionUniverse(chain,i) for i in OneSigmaShift]

class MuonUniverseMinerva():
    def __init__(self,chain,nsigma):
        self.universe = ROOT.PlotUtils.MuonUniverseMinerva[ROOT.PythonMinervaUniverse](chain, nsigma)
        self.cv_universe = CVUniverse(chain, nsigma)

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.universe, attr, getattr(self.cv_universe, attr, None))

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [MuonUniverseMinerva(chain,i) for i in OneSigmaShift]

class MuonUniverseMinos():
    def __init__(self,chain,nsigma):
        self.universe = ROOT.PlotUtils.MuonUniverseMinos[ROOT.PythonMinervaUniverse](chain, nsigma)
        self.cv_universe = CVUniverse(chain, nsigma)

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.universe, attr, getattr(self.cv_universe, attr, None))

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [MuonUniverseMinos(chain,i) for i in OneSigmaShift]

class MinosEfficiencyUniverse():
    def __init__(self,chain,nsigma):
        self.universe = ROOT.PlotUtils.MinosEfficiencyUniverse[ROOT.PythonMinervaUniverse](chain, nsigma)
        self.cv_universe = CVUniverse(chain, nsigma)

    def GetMyMinosEfficiencyWeight(self):
         if self.HasNoBackExitingTracks:
             return 1
         else:
             return super(MinosEfficiencyUniverse,self).GetMinosEfficiencyWeight()

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.universe, attr, getattr(self.cv_universe, attr, None))

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [MinosEfficiencyUniverse(chain,i) for i in OneSigmaShift]

###########################################################################
class ElectronEnergyShiftUniverse():
    def __init__(self,chain, nsigma,region):
        self.cv_universe = CVUniverse(chain,nsigma)
        self.region = region

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.cv_universe, attr, None)

    def GetEMEnergyShift(self):
        return self.nsigma*SystematicsConfig.EM_ENERGY_SCALE_UNCERTAINTY[self.region]*self.GetVecElem("prong_"+self.region+"CalibE",0)

    def ShortName(self):
        return "elE_"+self.region

    def LatexName(self):
        return "EM Energy Scale in "+self.region

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [ElectronEnergyShiftUniverse(chain, i,region) for region in SystematicsConfig.EM_ENERGY_SCALE_UNCERTAINTY for i in OneSigmaShift]


###########################################################################
class ElectronAngleShiftUniverse():
    def __init__(self,chain, nsigma):
        self.cv_universe = CVUniverse(chain,nsigma)
        self.axis_angle = random.random()*math.pi
        self.shift_angle = random.gauss(0,self.nsigma*SystematicsConfig.LEPTON_ANGLE_UNCERTAINTY)

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.cv_universe, attr, None)

    def ShortName(self):
        return "eltheta"

    def LatexName(self):
        return "Election Candidate Angle"

    def ElectronP3D(self):
        p = super(ElectronAngleShiftUniverse,self).ElectronP3D()
        unit_3_vector = p/p.R()
        normal_vector1 = p.Cross(ROOT.Math.XYZVector(1.0,0,0))
        if normal_vector1.R() != 0:
            normal_vector1 *= 1.0/normal_vector1.R()
        r1 = ROOT.Math.AxisAngle(unit_3_vector,self.axis_angle)
        rotation_axis = r1(normal_vector1)
        r2 = ROOT.Math.AxisAngle(rotation_axis,self.shift_angle)
        return r2(p)

    def ElectronTheta2(self):
        theta = super(ElectronAngleShiftUniverse,self).ElectronTheta()
        return theta+self.shift_angle

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [ElectronAngleShiftUniverse(chain,  1) for i in range(SystematicsConfig.NUM_UNIVERSE)]

    

###########################################################################
class BirksShiftUniverse():
    def __init__(self,chain, nsigma):
        self.cv_universe = CVUniverse(chain,nsigma)

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.cv_universe, attr, None)

    @property
    def prong_part_score(self):
        return self.prong_part_score_BirkUP if self.nsigma>0 else self.prong_part_score_BirkDOWN

    @property
    def prong_dEdXMeanFrontTracker(self):
        return self.prong_dEdXMeanFrontTracker_BirkUP if self.nsigma>0 else self.prong_dEdXMeanFrontTracker_BirkDOWN

    @property
    def prong_NonMIPClusFrac(self):
        return self.prong_NonMIPClusFrac_BirkUP if self.nsigma>0 else self.prong_NonMIPClusFrac_BirkDOWN

    def ShortName(self):
        return "birks"

    def LatexName(self):
        return "Birk's Constant"

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [BirksShiftUniverse(chain, i) for i in OneSigmaShift]

###########################################################################

class BeamAngleShiftUniverse():
    def __init__(self,chain, nsigma, x):
        self.cv_universe = CVUniverse(chain,nsigma)
        self.rotation =  ROOT.Math.RotationX(SystematicsConfig.BEAM_XANGLE_UNCERTAINTY*nsigma) if x else ROOT.Math.RotationY(SystematicsConfig.BEAM_YANGLE_UNCERTAINTY*nsigma)

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.cv_universe, attr, None)

    def ShortName(self):
        return "beam_angle"

    def LatexName(self):
        return "Beam Angle"

    def ElectronP3D(self):
        p = super(BeamAngleShiftUniverse,self).ElectronP3D()
        return self.rotation(p)

    @staticmethod
    def GetSystematicsUniverses(chain ):
        return [BeamAngleShiftUniverse(chain, i, j) for j in [True,False] for i in OneSigmaShift]


class GeantHadronUniverse():
    def __init__(self,chain,nsigma,pdg):
        self.universe = ROOT.PlotUtils.GeantHadronUniverse[ROOT.PythonMinervaUniverse](chain, nsigma, pdg)
        self.cv_universe = CVUniverse(chain, nsigma)  

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.universe, attr, getattr(self.cv_universe, attr, None))

    def GetWeight(self,test=None):
        return self.GetStandardWeight()

    def GetStandardWeight(self):
        weight = self.cv_universe.GetStandardWeight()
        weight *= self.universe.GetGeantHadronWeight()
        return weight

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [GeantHadronUniverse(chain,i,j) for j in SystematicsConfig.GEANT_PARTICLES for i in OneSigmaShift]

class TargetMassUniverse():
    def __init__(self,chain,nsigma):
        self.universe = ROOT.PlotUtils.TargetMassScintillatorUniverse[ROOT.PythonMinervaUniverse](chain, nsigma)
        self.cv_universe = CVUniverse(chain, nsigma)

    def GetWeight(self, test=None):
        return self.GetStandardWeight()

    def GetStandardWeight(self):
        weight = self.cv_universe.GetStandardWeight()
        weight*= self.GetWeightRatioToCV()
        return weight

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.universe, attr, getattr(self.cv_universe, attr, None))

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [TargetMassUniverse(chain, i) for i in OneSigmaShift]

class MKModelUniverse():
    def __init__(self,chain,nsigma):
        self.cv_universe = CVUniverse(chain, nsigma)  
        self.reweighter = ROOT.PlotUtils.MKReweighter(ROOT.PythonMinervaUniverse,ROOT.PlotUtils.detail.empty)()

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.cv_universe, attr, None)

    def IsVerticalOnly(self):
        return True

    def GetStandardWeight(self):
        weight = self.cv_universe.GetStandardWeight()
        weight*= self.reweighter.GetWeight(self.cv_universe,ROOT.PlotUtils.detail.empty())
        return weight

    def GetWeight(self, test=None):
        return self.GetStandardWeight()

    def ShortName(self):
        return "MK_model"

    def LatexName(self):
        return "MK Model"

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [MKModelUniverse(chain,1)]


class FSIWeightUniverse():
    def __init__(self,chain,nsigma,iweight):
        self.cv_universe = CVUniverse(chain,nsigma)
        self.reweighter = ROOT.PlotUtils.FSIReweighter(ROOT.PythonMinervaUniverse,ROOT.PlotUtils.detail.empty)((iweight+1)//2,(iweight+1)%2)

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.cv_universe, attr, None)

    def IsVerticalOnly(self):
        return True

    def GetWeight(self, test=None):
        return self.GetStandardWeight()

    def GetStandardWeight(self):
        weight = self.cv_universe.GetStandardWeight()
        weight*= self.reweighter.GetWeight(self.cv_universe,ROOT.PlotUtils.detail.empty())
        return weight

    def ShortName(self):
        return "fsi_weight"

    def LatexName(self):
        return "FSI weight"

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [FSIWeightUniverse(chain,1,i) for i in range(3)]

class SusaValenciaUniverse():
    def __init__(self,chain,nsigma):
        self.cv_universe = CVUniverse(chain,nsigma)
        self.reweighter = ROOT.PlotUtils.SuSAFromValencia2p2hReweighter(ROOT.PythonMinervaUniverse,ROOT.PlotUtils.detail.empty)()

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.cv_universe, attr, None)

    def IsVerticalOnly(self):
        return True

    def GetLowRecoil2p2hWeight(self):
        return 1.0

    def GetRPAWeight(self):
        return 1.0

    def ShortName(self):
        return "SuSA_Valencia_Weight"

    def LatexName(self):
        return "SuSA Valencia Weight"

    def GetWeight(self, test=None):
        return self.GetStandardWeight()

    def GetStandardWeight(self):
        weight = self.cv_universe.GetStandardWeight()
        weight*= self.reweighter.GetWeight(self.cv_universe,ROOT.PlotUtils.detail.empty())
        return weight

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [SusaValenciaUniverse(chain,1)]

class LeakageUniverse():
    def __init__(self,chain,nsigma):
        self.cv_universe = CVUniverse(chain,nsigma)

    def __getattr__(self, attr):
        """Redirect attribute access to FluxUniverse first, then CVUniverse"""
        return getattr(self.cv_universe, attr, None)

    def GetLeakageCorrection(self): 
        return super(LeakageUniverse,self).GetLeakageCorrection() + ( self.nsigma*SystematicsConfig.LEAKAGE_SYSTEMATICS if abs(self.mc_primaryLepton) ==11 else 0)

    def ShortName(self):
        return "Leakage_Uncertainty"

    def LatexName(self):
        return "Leakage Uncertainty"

    @staticmethod
    def GetSystematicsUniverses(chain):
        return [LeakageUniverse(chain,i) for i in OneSigmaShift]


def GetAllSystematicsUniverses(chain,is_data,is_pc =False,exclude=None,playlist=None):
    #use list because I want to control the order, process cv first, then vertical shifts, then lateral shifts
    universes = []
    CVUniverse.is_pc = is_pc
    CVUniverse.LeafGetters= {}
    if is_data:
        #data has only cv universe
        universes.append(CVUniverse(chain, None))
    else:
        #append cv universe
        universes.append(CVUniverse(chain,0))
        
        #Set Playlist, only MC cares about playlist
        if playlist is None:
            #attempt guess playlist from the chain
            playlist = Utilities.PlaylistLookup(chain.GetValue("mc_run",0))

        if not is_pc:
            CVUniverse.SetPlaylist(playlist)
            #Set NuE constraint
            CVUniverse.SetNuEConstraint(SystematicsConfig.USE_NUE_CONSTRAINT)
        else:
            CVUniverse.SetPlaylist("dummy")
            CVUniverse.SetNuEConstraint(False)

        # Set analysis nu PDG if provided as int
        # CVUniverse.SetAnalysisNuPDG(SystematicsConfig.AnaNuPDG)
        ana_pdg = getattr(SystematicsConfig, "AnaNuPDG", None)
        if isinstance(ana_pdg, int):
            CVUniverse.SetAnalysisNuPDG(int(ana_pdg))

        # #Set NonResPi weight
        # CVUniverse.SetNonResPiReweight(False)
        # CVUniverse.SetDeuteriumGeniePiTune(False)
        # CVUniverse.SetZExpansionFaReweight(False) 
        # CVUniverse.RPAMaterials(False)

        # Tunes / flux / materials
        CVUniverse.SetNonResPiReweight(True)
        CVUniverse.SetDeuteriumGeniePiTune(True)
        CVUniverse.SetZExpansionFaReweight(bool(SystematicsConfig.NumZExpansionUniverses)) 
        CVUniverse.SetNFluxUniverses(SystematicsConfig.NUM_FLUX_UNIVERSE)
        CVUniverse.RPAMaterials(True)

        # Truth tree toggle
        if chain.GetTree().GetName() == "Truth":
            CVUniverse.SetTruth(True)

        if exclude is None or "all" not in exclude:
            elastic_mode = bool(getattr(SystematicsConfig, "ELASTIC_MODE", False))
            abs_pdg = abs(ana_pdg) if isinstance(ana_pdg, int) else None

            # ===== Lepton reconstruction universes =====
            # For ν–e elastic (ELASTIC_MODE) or CC νe (|PDG|=12), use electron-side universes.
            if elastic_mode or abs_pdg == 12:
                universes.extend(ElectronEnergyShiftUniverse.GetSystematicsUniverses(chain))
                universes.extend(ElectronAngleShiftUniverse.GetSystematicsUniverses(chain))
            elif abs_pdg == 14:
                universes.extend(MuonUniverseMinerva.GetSystematicsUniverses(chain))
                universes.extend(MuonUniverseMinos.GetSystematicsUniverses(chain))
                universes.extend(MuonResolutionUniverse.GetSystematicsUniverses(chain))
                universes.extend(MuonAngleXResolutionUniverse.GetSystematicsUniverses(chain))
                universes.extend(MuonAngleYResolutionUniverse.GetSystematicsUniverses(chain))
                universes.extend(MinosEfficiencyUniverse.GetSystematicsUniverses(chain))
            else:
                # If PDG is weird and not explicitly in elastic mode
                if not elastic_mode:
                    raise ValueError(
                        "AnaNuPDG should be ±12 or ±14 (or set ELASTIC_MODE=True for ν–e). Got: {}".format(ana_pdg)
                    )

            # ===== Beam angle (affects θ reconstruction) =====
            universes.extend(BeamAngleShiftUniverse.GetSystematicsUniverses(chain))

            # ===== Detector response (ensure EM-only branches for ν–e via SystematicsConfig.RESPONSE_BRANCHES) =====
            universes.extend(ResponseUniverse.GetSystematicsUniverses(chain))

            print("Calc Flux Systematics...")
            # ===== Flux universes (per-event species weighting) =====
            universes.extend(FluxUniverse.GetSystematicsUniverses(chain))
            print("Done Flux Systematics...")

            # ===== Generator/tune/hadronic universes =====
            # Skip for ν–e elastic signal; include for CC analyses.
            if not elastic_mode:
                universes.extend(GenieUniverse.GetSystematicsUniverses(chain))
                universes.extend(GenieRvx1piUniverse.GetSystematicsUniverses(chain))
                universes.extend(GenieFaCCQEUniverse.GetSystematicsUniverses(chain))
                universes.extend(GenieMaResUniverse.GetSystematicsUniverses(chain))
                universes.extend(GenieMvResUniverse.GetSystematicsUniverses(chain))
                universes.extend(GenieNormCCResUniverse.GetSystematicsUniverses(chain))
                print("done Genie")
                universes.extend(Universe2p2h.GetSystematicsUniverses(chain))
                print("done 2p2h")
                universes.extend(RPAUniverse.GetSystematicsUniverses(chain))
                print("done RPA")
                #universes.extend(NonResonantPionUniverse.GetSystematicsUniverses(chain ))
                universes.extend(LowQ2PionUniverse.GetSystematicsUniverses(chain))
                print("done Low Q2")
                #universes.extend(LowQ2PionUniverseAlt.GetSystematicsUniverses(chain )) used for warping study variant
                #universes.extend(BirksShiftUniverse.GetSystematicsUniverses(chain ))
                universes.extend(MKModelUniverse.GetSystematicsUniverses(chain))
                universes.extend(FSIWeightUniverse.GetSystematicsUniverses(chain))
                universes.extend(SusaValenciaUniverse.GetSystematicsUniverses(chain))
                # universes.extend(GeantHadronUniverse.GetSystematicsUniverses(chain))  # keep off unless needed

            # ===== Leakage & target normalization =====
            universes.extend(LeakageUniverse.GetSystematicsUniverses(chain))
            universes.extend(TargetMassUniverse.GetSystematicsUniverses(chain))



    # Group universes in dict.
    univ_dict = OrderedDict()
    for univ in universes:
        univ_dict.setdefault(univ.ShortName(), []).append(univ)

    if exclude is not None and len(univ_dict)>1:
        for i in exclude:
            if i in univ_dict:
                print("deleting universe {} from {}".format(i,univ_dict[i]))
                del univ_dict[i]

    return univ_dict

#code for test=Noneing new cv universe.
if __name__ == "__main__":

    #path = "/exp/minerva/data/users/hsu/merged_files/merged-1558038217.root"
    #path = "/pnfs/minerva/persistent/users/jyhan/NuECCQE/v21r1p1_recoilE_final/mcme1A/grid/central_value/minerva/ana/v21r1p1/00/11/00/17/SIM_minerva_00110017_Subruns_0325_NuECCQE_Ana_Tuple_v21r1p1.root"
    path = "/pnfs/minerva/persistent/DataPreservation/p4/FullDetector/Merged_mc_ana_me1A_DualVertex_p4/MasterAnaDev_mc_AnaTuple_run00110011_Playlist.root"
    
    chainWrapper = PlotUtils.ChainWrapper("MasterAnaDev")
    chainWrapper.Add(path)
    chainWrapper.GetValue("prong_part_score",0)

    univs = GetAllSystematicsUniverses(chainWrapper,False)

    for u in univs.keys():
        print(u)
        for univ in univs[u]:
            univ.SetEntry(5)
            weight = univ.GetWeight(False)
            print(weight)
        print("")
