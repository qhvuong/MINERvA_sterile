"""
  KinematicsCalculator.py:
   Tools for calculating derived kinematics quantities.
  
   Original author: J. Wolcott (jwolcott@fnal.gov)
                    January 2015
"""

import math
import ROOT
import os
from tools import Utilities
from numpy import interp
from config import SystematicsConfig

# constants needed in calculations.
M_e = 0.511 # MeV
M_p = 938.272 # MeV
M_n = 939.565 # MeV
M_pion = 139.5701 # MeV
M_mu =  105.658 #MeV
QE_binding_E = 34 # MeV  (same as nu_mu PRL for FHC)

# these will be used repeatedly in calculations.
# better just do the arithmetic once here.
M_n_star = (M_n - QE_binding_E)
M_e_sqr = M_e**2
M_mu_sqr = M_mu**2
M_n_star_sqr = M_n_star**2
M_p_sqr = M_p**2
BEAM_ANGLE = math.radians(-3.3)
 
sin_ee_sqr = 0.4
sin_ue_sqr = 0.002
delta_m_sqr = 7.34 # eV^2

L = 1
E = 1

class KinematicsCalculator(object):
    DEFAULTS = {
        "correct_beam_angle": True,
        "calc_reco": True,
        "calc_true": True,
        "include_hadron_momenta": False,
        "is_pc" : False
    }
    def __init__(self, **args):
        # copy any user-specified values out of the args.
        # use defaults if any are not specified.
        for kw, val in KinematicsCalculator.DEFAULTS.items():
            if kw in args:
                val = args[kw]
                del args[kw]
            setattr(self, kw, val)

        self._clear = False
        self.new_truth = False
        self.Clear()
        self.Cals = [self.Clear]
        if self.calc_reco:
            self.Cals.append(self.CalculateRecoKinematics)
        if self.calc_true:
            self.Cals.append(self.CalculatePCTrueKinematics if self.is_pc else self.CalculateTrueKinematics)
        self.spliner = NuMuSpline()

    @staticmethod
    def calcq3(Q2, Enu, Elep):
        return math.sqrt(Q2+(Enu-Elep)**2)

    def AngleBetween(self, v1, v2):
        dotprod = v1.Dot(v2)
        mag = v1.Mag() * v2.Mag()
        if mag > 0:
            dotprod /= mag
#       print dotprod
        return math.degrees(math.acos(dotprod))

    def Clear(self):
        if self._clear:
            return
#       print "clearing kinematics"
        self.reco_q2_cal = None
        self.reco_theta_lep = None
        self.reco_theta_lep_rad_det = None
        self.reco_theta_lep_rad = None
        self.reco_thetaX_lep_rad_det = None
        self.reco_thetaX_lep_rad = None
        self.reco_thetaY_lep_rad_det = None
        self.reco_thetaY_lep_rad = None
        # self.reco_theta2D_lep_rad = None
        # self.reco_phi_lep_rad = None
        self.reco_cos_theta_lep = None

        self.reco_E_lep = None
        self.reco_E_nu_cal = None
        self.reco_W_DFR = None
        self.reco_visE = None
        self.reco_q0 = None
        self.reco_q3 = None
        self.LLR = None
        self.Euv = None
        self.Exuv = None
        self.reco_Etheta2 = None
        if self.new_truth :
            self.true_theta_lep_rad = None
            self.true_thetaX_lep_rad = None
            self.true_thetaY_lep_rad = None
            self.true_theta2D_lep_rad = None
            self.true_phi_lep_rad = None
            self.true_theta_lep = None
            self.true_E_lep = None
            self.true_P_lep = None
            self.true_q2_cal = None
            self.true_enu_genie = None
            self.true_q0 = None
            self.true_q3 = None
            self.true_q2 = None
            self.true_visE = None
            self.true_LeptonP3D = None
            self.true_LeptonP3D_det = None
            if self.include_hadron_momenta:
                self.proton_momenta = []
                self.proton_thetas = []
                self.pion_momenta = []
                self.pion_thetas = []
                self.pi0_momenta = []
                self.pi0_thetas = []
        self._clear = True

        return True

    def CalculateKinematics(self,event):
        self.event = event
        self.new_truth = event.ShortName() == "cv"
        sc = True
        for _ in self.Cals:
            try:
                sc = _() and sc
            except Exception as e:
                print(e)
                #raise(e)
                sc = False
        return sc

    def CalculateRecoKinematics(self):
        self._clear = False
        if self.event.SetLeptonType():
            nParticles = self.event.GetVecInt("prong_nParticles")
            if len(nParticles) < 1 or nParticles[0] < 1:
                return False
        event = self.event

        #decide if whether it is muon or electron event, and get lepton kinematics
        # self.reco_theta_lep_rad = event.LeptonTheta()
        # self.reco_thetaX_lep_rad = event.LeptonThetaX()
        # self.reco_thetaY_lep_rad = event.LeptonThetaY()
        # self.reco_theta2D_lep_rad = event.LeptonTheta2D()
        # self.reco_phi_lep_rad = event.LeptonPhi()
        # self.reco_vertex3D_beam      = event.Vertex3D_beam()
        self.reco_LeptonP3D_det      = event.LeptonP3D_det()
        self.reco_LeptonP3D          = event.LeptonP3D()
        # print(self.reco_LeptonP3D_det.Unit())
        # print(type(self.reco_LeptonP3D_det))
        # print([m for m in dir(self.reco_LeptonP3D_det) if m.lower() in ("mag","mag2","r","r2","rho","rho2","perp","perp2","pt","pt2","unit")])


        self.reco_theta_lep_rad_det  = self.reco_LeptonP3D_det.Theta()
        self.reco_theta_lep_rad      = self.reco_LeptonP3D.Theta()
        self.reco_thetaX_lep_rad_det = math.atan2(self.reco_LeptonP3D_det.X(), self.reco_LeptonP3D_det.Z())
        self.reco_thetaX_lep_rad     = math.atan2(self.reco_LeptonP3D.X(), self.reco_LeptonP3D.Z())
        self.reco_thetaY_lep_rad_det = math.atan2(self.reco_LeptonP3D_det.Y(), self.reco_LeptonP3D_det.Z())
        self.reco_thetaY_lep_rad     = math.atan2(self.reco_LeptonP3D.Y(), self.reco_LeptonP3D.Z())

        # self.reco_thetaY_lep_rad_det = event.LeptonThetaY_det()
        # self.reco_thetaY_lep_rad     = event.LeptonThetaY_beam()
        self.reco_E_lep = event.LeptonEnergy()/1e3
        self.M_lep_sqr = event.M_lep_sqr/1e6
        self.reco_P_lep = math.sqrt(max(0,self.reco_E_lep**2-self.M_lep_sqr))
        self.reco_cos_theta_lep = math.cos(self.reco_theta_lep_rad)
        self.reco_theta_lep = math.degrees(self.reco_theta_lep_rad)
        self.reco_Pt_lep = self.reco_P_lep*math.sin(self.reco_theta_lep_rad)

        #get appearance or disappearance probabilities
        #self.reco_appearance = sin_ue_sqr*np.sin(1.27*delta_m_sqr*L/E)**2
        #self.reco_dissapearance = sin_ee_sqr*np.sin(1.27*delta_m_sqr*L/E)**2\

        #now hadronic calorimetry
        self.reco_q0 = self.spliner.get_q0(self.event.RecoilEnergy()/1e3)
        self.reco_visE = self.event.AvailableEnergy()/1e3
        
        # calc q2, q3
        self.reco_E_nu_cal = self.reco_E_lep + self.reco_q0
        self.reco_E_nu_QE = KinematicsCalculator.Enu_QE(self.reco_E_lep, self.reco_P_lep, self.reco_theta_lep_rad, self.M_lep_sqr)
        self.reco_q2_cal = KinematicsCalculator.Q2_cal(self.reco_E_lep, self.reco_theta_lep_rad, self.reco_P_lep, self.M_lep_sqr, self.reco_E_nu_cal)

        self.reco_q3 = self.calcq3(self.reco_q2_cal,self.reco_E_nu_cal,self.reco_E_lep)
        self.pseudo_q3 = self.reco_E_nu_cal - self.reco_P_lep*self.reco_cos_theta_lep
        self.reco_W2 = (M_n/1e3)**2+2*M_n/1e3*self.reco_q0-self.reco_q2_cal
        self.reco_W = math.sqrt(self.reco_W2) if self.reco_W2>=0 else -1
        self.reco_Etheta2 = self.reco_E_lep * self.reco_theta_lep_rad**2
        #print(self.event.GetLowRecoilQ3(),self.reco_q3)
        viewE = (self.event.prong_XViewE[0], self.event.prong_VViewE[0],self.event.prong_UViewE[0])
        self.Ex = (viewE[0])/sum(viewE)
        self.Eu = (viewE[1])/sum(viewE)
        self.Ev = (viewE[2])/sum(viewE)
        self.Exuv = self.Ex-self.Eu-self.Ev
        self.Euv = self.Eu-self.Ev
            #self.LLR = self.event.GetLLR()
        # except RuntimeError:
            #this ntuple does not have these variables
            # pass

        # self.leftE = sum(self.event.ExtraEnergyClusters_energy[i] for i in range(self.event.ExtraEnergy_nClusters) if self.event.ExtraEnergyClusters_X[i]>self.event.vtx[0])
        # self.leftE += sum(self.event.eProngClusters_energy[i] for i in range(self.event.eProng_nClusters) if self.event.eProngClusters_X[i]>self.event.vtx[0])
        # self.rightE = sum(self.event.ExtraEnergyClusters_energy[i] for i in range(self.event.ExtraEnergy_nClusters) if self.event.ExtraEnergyClusters_X[i]<self.event.vtx[0])
        # self.rightE += sum(self.event.eProngClusters_energy[i] for i in range(self.event.eProng_nClusters) if self.event.eProngClusters_X[i]<self.event.vtx[0])
        return True

    def CalculatePCTrueKinematics(self):
        assert self.event.mc_processType == 5
        if not self.new_truth:
            return True
        self._clear = False
        # type 5 is mc

        #find electron index. 0 if not found for any reason.
        e_index = [i for i,x in enumerate(self.event.mc_FSPartPDG) if abs(x) == 11]
        e_index = e_index[0] if len(e_index)>0 else 0
        mom = ROOT.Math.XYZTVector(self.event.mc_FSPartPx[e_index], self.event.mc_FSPartPy[e_index], self.event.mc_FSPartPz[e_index], self.event.mc_FSPartE[e_index])
        # if self.correct_beam_angle:
        #     mom = ROOT.Math.RotationX(BEAM_ANGLE)*mom

        self.true_theta_lep = math.degrees(mom.Theta())
        self.true_E_lep = self.event.mc_FSPartE[0] / 1e3
        #self.true_lep_momentum = mom.Vect()
        self.true_visE = self.CalculateVisibleE(True,11)

        return True

    def CalculateTrueKinematics(self):
        assert self.event.mc_processType >= 0 and self.event.mc_processType != 5
        if not self.new_truth:
            return True
        self._clear = False

        self.true_enu_genie = self.event.GetEnuTrue()/1e3
        self.true_E_lep = self.event.GetElepTrue()/1e3
        self.true_P_lep = self.event.GetPlepTrue()/1e3
        self.true_q0 = self.event.Getq0True()/1e3
        self.true_q3 = self.event.Getq3True()/1e3
        self.true_q2 = self.event.GetQ2True()/1e6
        self.true_theta_lep_rad = self.event.GetThetalepTrue()
        self.true_theta_lep = math.degrees(self.true_theta_lep_rad)
        self.true_visE = self.CalculateVisibleE()
        self.true_Pt_lep = self.true_P_lep*math.sin(self.true_theta_lep_rad)

        # --- NEW: true thetaX, thetaY in beam coordinates ---
        # mc_primFSLepton = [px, py, pz, E] in MeV in detector coords.
        # TruthFunctions.h rotates by MinervaUnits::numi_beam_angle_rad before taking angles.
        px = self.event.GetVecElem("mc_primFSLepton", 0)
        py = self.event.GetVecElem("mc_primFSLepton", 1)
        pz = self.event.GetVecElem("mc_primFSLepton", 2)

        p3 = ROOT.TVector3(px, py, pz)              # beam
        self.true_LeptonP3D = ROOT.TVector3(p3)     # store a copy (beam)

        p_det = ROOT.TVector3(p3)                   # det = inverse rotation
        p_det.RotateX(-SystematicsConfig.BEAM_ANGLE)
        self.true_LeptonP3D_det = p_det

        # print(type(p3))
        # print([m for m in dir(p3) if m.lower() in ("mag","mag2","r","r2","rho","rho2","perp","perp2","pt","pt2","unit")])

        # Use the SAME convention as TruthFunctions.h:
        # p3lep.RotateX(MinervaUnits::numi_beam_angle_rad);
        # p3.RotateX(ROOT.PlotUtils.MinervaUnits.numi_beam_angle_rad)
        # det angles
        self.true_theta_lep_rad_det  = p_det.Theta()
        self.true_thetaX_lep_rad_det = math.atan2(p_det.X(), p_det.Z())
        self.true_thetaY_lep_rad_det = math.atan2(p_det.Y(), p_det.Z())

        # beam angles
        self.true_thetaX_lep_rad = math.atan2(p3.X(), p3.Z())
        self.true_thetaY_lep_rad = math.atan2(p3.Y(), p3.Z())
        self.true_thetaX_lep = math.degrees(self.true_thetaX_lep_rad)
        self.true_thetaY_lep = math.degrees(self.true_thetaY_lep_rad)
        # self.true_theta2D_lep_rad = math.sqrt(self.true_thetaX_lep_rad**2 + self.true_thetaY_lep_rad**2)
        # self.true_phi_lep_rad = math.atan2(p3.Y(), p3.X())
    
    
        vtxX = self.event.GetVecElem("mc_vtx", 0)
        vtxY = self.event.GetVecElem("mc_vtx", 1)
        vtxZ = self.event.GetVecElem("mc_vtx", 2)
        vtx3D = ROOT.TVector3(vtxX, vtxY, vtxZ)
        vtx3D.RotateX(SystematicsConfig.BEAM_ANGLE)
        self.true_vertex3D_beam = vtx3D
        # -----------------------------------------------

        if self.include_hadron_momenta:
            for part_idx in range(len(self.event.mc_FSPartPDG)):
                pdg = abs(self.event.mc_FSPartPDG[part_idx])
                if pdg not in (111, 211, 2212):
                    continue

                mom = ROOT.TVector3(self.event.mc_FSPartPx[part_idx], self.event.mc_FSPartPy[part_idx], self.event.mc_FSPartPz[part_idx])
                mom_mag = mom.Mag()/1e3
                theta = math.degrees(mom.Theta())
#               print pdg, mom_mag, theta

                if pdg == 2212:
                    self.proton_momenta.append(mom_mag)
                    self.proton_thetas.append(theta)
                elif pdg == 111:
                    self.pi0_momenta.append(mom_mag)
                    self.pi0_thetas.append(theta)
                else:
                    self.pion_momenta.append(mom_mag)
                    self.pion_thetas.append(theta)

        return True
        #       print E_e, self.E_nu_QE, self.Q2_QE

    @staticmethod
    def Enu_QE(E_lep, p_lep, theta_lep,m_lep_sqr):
        numerator = M_p_sqr/1e6 - M_n_star_sqr/1e6 - m_lep_sqr + 2 * M_n_star/1e3 * E_lep
        denominator = 2 * (M_n_star/1e3 - E_lep + p_lep * math.cos(theta_lep))
        
        if denominator == 0:
            return None
    
        return numerator/denominator

    @staticmethod
    def Q2_cal(E_lep, theta_lep, p_lep, M_lep_sqr, Enu):
        return max(0, 2 * Enu * (E_lep - p_lep * math.cos(theta_lep)) - M_lep_sqr)

    def CalculateVisibleE(self):
        totalE = lambda E:E
        PDG_energy_map ={
            2212: (lambda E: E-M_p),
            2112: (lambda E: 0),
            211: (lambda E: E-M_pion),
            -211: (lambda E: E-M_pion),
            111 : totalE,
            22 : totalE,
            11 : (lambda E:0),
            -11: (lambda E:0),
            13:  (lambda E:0),
            -13: (lambda E:0),
        }
        Eavail = 0.0
        for i in range(len(self.event.mc_FSPartPDG)):
            E = self.event.mc_FSPartE[i]
            pdg  = self.event.mc_FSPartPDG[i]
            #print (pdg,E)
            if PDG_energy_map.get(pdg) is not None:
                Eavail += PDG_energy_map.get(pdg)(E)
            elif pdg>1000000000 :
                continue
            elif pdg>=2000:
                Eavail+=E-M_p
            elif pdg<=-2000:
                Eavail+=E+M_p
            else:
                Eavail+=E

        #a = max(0,Eavail - E_primary_lepton)
        return max(0,Eavail)/1000

        # protonE = 0.0
        # neutronE = 0.0
        # pionE = 0.0
        # pizeroE = 0.0
        # muonE = 0.0
        # gammaE = 0.0
        # electron = 0.0
        # otherE = 0.0
        # missingE = 0.0
        # nucleonsE = 0.0

        # protonN = 0
        # neutronN = 0
        # pionN = 0
        # pizeroN = 0
        # muonN = 0
        # electronN = 0
        # gammaN = 0
        # otherN = 0

        # for i in range(len(self.event.mc_FSPartPDG)):
        #     E = self.event.mc_FSPartE[i]
        #     pdg  = self.event.mc_FSPartPDG[i]
        #     if pdg==2212:
        #         protonE += E - M_p
        #         if E - M_p > 5.0:
        #             protonN += 1
        #     elif pdg==2112:
        #         neutronE += E - M_n;
        #         if E - M_n > 5.0:
        #             neutronN += 1
        #     elif abs(pdg)==211:
        #         pionE += E
        #         pionN += 1
        #     elif pdg==111:
        #         pizeroE += E
        #         pizeroN += 1
        #     elif abs(pdg)==13:
        #         muonE += E
        #         muonN += 1
        #     elif abs(pdg)==22:
        #         gammaE += E
        #         gammaN += 1
        #     elif abs(pdg) != 11:
        #         #The Other Category.
        #         #cout << "other " <<  pdg << " energy " << E << endl;
        #         if pdg>=2000000000:
        #             missingE += E
        #         elif pdg>1000000000:
        #             nucleonsE += E
        #             # do nothing, assume negligible energy for nucleons
        #             # save me the trouble of asking what that nucleon's mass was.
        #         elif pdg>=2000:
        #             otherE += E - M_p
        #             otherN += 1
        #             #primarily strange baryons 3122s
        #         elif pdg<=-2000:
        #             otherE += E + M_p
        #             otherN += 1
        #             #primarily anti-baryons -2112 and -2212s,
        #         else:
        #             otherE += E
        #             otherN += 1
        #             #primarily kaons and eta mesons.

                # True visible energy is the kinetic energy of protons and pions, plus the total energies of everything else

        visibleE = protonE + pionE - pionN*M_pion + pizeroE + muonE + gammaE + otherE
        return max(0,visibleE)/1000.0
        # else :
        #     return max(0,self.event.recoile_passive)/1e3

class q3_spline(object):
    def __init__(self,s,xp,fp):
        self.xp = xp
        self.fp = fp
        self.s = s

    @staticmethod
    def initByFile(path,name):
        scale = 1.0
        x =[]
        y =[]
        reading = False
        with open(path) as f:
            lines = f.readlines()
            for line in lines:
                line = line.split();
                if len(line)==0:
                    continue
                if line[0]=="BEGIN" and line[1] == name:
                    reading = True
                    continue
                if not reading:
                    continue
                if line[0] == "END":
                    break
                if line[0]=="SCALE":
                    scale = float(line[1])
                    continue
                if line[0]=="POLYPOINT":
                    x.append(float(line[1]))
                    y.append(float(line[2]))
                    continue
        return scale,x,y


    def get_q0(self,eavail):
        return self.s * interp(eavail,self.xp,self.fp)
       

class DefaultSpline(q3_spline):
    def __init__(self):
        super(DefaultSpline,self).__init__(1,[0,50],[0.50])


class NuESpline(q3_spline):
    def __init__(self):
        super(NuESpline, self).__init__(*q3_spline.initByFile( "{}/Calibrations/energy_calib/CCNueCalorimetryTunning.txt".format(os.environ["MPARAMFILES"]),"CCNuE_Nu_Tracker"))

class NuMuSpline(q3_spline):
    def __init__(self):
         super(NuMuSpline, self).__init__(*q3_spline.initByFile( "{}/Calibrations/energy_calib/NukeCCInclusiveTunings.txt".format(os.environ["MPARAMFILES"]),"NukeCC_Nu_Tracker"))
