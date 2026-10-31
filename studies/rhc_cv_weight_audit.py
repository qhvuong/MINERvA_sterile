"""Read-only CV-weight audit helper for selection/eventSelection.py.

Use only for MC CV events after Classify(). Writes selected-event rows to CSV.
See integration snippet in chat. No production weights are modified.
"""
import csv
import math
from collections import defaultdict

class CVWeightAudit:
    fields = ["entry","pdg","intType","Etrue_GeV","Ereco_GeV","theta_deg",
              "truth_class","nonrespi","flux","mec","rpa","geant","minos",
              "product","cached","relative_diff","invalid_mec_kin"]
    def __init__(self, filename="rhc_cv_weight_audit.csv", max_rows=20000):
        self.f=open(filename,"w",newline="")
        self.writer=csv.DictWriter(self.f,fieldnames=self.fields)
        self.writer.writeheader()
        self.max_rows=max_rows
        self.rows=0
        self.counts=defaultdict(int)
        self.max_error=0.

    def record(self,u,entry):
        if self.rows>=self.max_rows or u.ShortName()!="cv":
            return
        # This helper matches the active V1 configuration (no deuterium,
        # no low-Q2-pion/coherent/diffractive terms). Confirm on GPVM.
        import tools.SystematicsUniverse as su
        if not su.USE_MNVTUNE_V1 or su.V1_ADD_DEUTERIUM or any(
            [su.V1_ADD_LOWQ2PI,su.V1_ADD_COH,su.V1_ADD_LE_COH,su.V1_ADD_DIFFCOH]
        ):
            raise RuntimeError("Audit requires current simplified V1 switches; update component product first")
        pdg=int(u.mc_incoming)
        if su.SystematicsConfig.USE_SWAPPED:
            pdg={12:14,-12:-14,14:12,-14:-12}.get(pdg,pdg)
        rvn=u.GetVecElem("truth_genie_wgt_Rvn1pi",2)
        rvp=u.GetVecElem("truth_genie_wgt_Rvp1pi",2)
        nonrespi=0.43 if rvn<1 or rvp<1 else 1.0
        flux=u.GetFluxAndCVWeight(u.mc_incomingE*1e-3,pdg)
        mec=1.0
        invalid=False
        if int(u.mc_intType)==8 and int(u.mc_targetZ)>=2:
            q0=u.Getq0True()*1e-3
            q3=u.Getq3True()*1e-3
            if math.isfinite(q0) and math.isfinite(q3):
                mec=u.GetLowRecoil2p2hWeight()
            else:
                invalid=True
        rpa=u.GetRPAWeight()
        geant=u.GetGeantHadronWeight()
        minos=u.GetMyMinosEfficiencyWeight()
        product=nonrespi*flux*mec*rpa*geant*minos
        cached=u.GetWeight(False)
        vals=[nonrespi,flux,mec,rpa,geant,minos,product,cached]
        if not all(math.isfinite(float(v)) for v in vals):
            raise RuntimeError(f"NONFINITE CV weight entry {entry}: {vals}")
        relative_diff=abs(product-cached)/max(abs(cached),1e-12)
        self.max_error=max(self.max_error,relative_diff)
        self.counts["events"]+=1
        self.counts["invalid_mec_kin"]+=int(invalid)
        if relative_diff>1e-8:
            self.counts["closure_fail"]+=1
            if self.counts["closure_fail"]<11:
                print("[CV WEIGHT MISMATCH]",entry,relative_diff,vals,flush=True)
        self.writer.writerow(dict(
            entry=entry,pdg=int(u.mc_incoming),intType=int(u.mc_intType),
            Etrue_GeV=float(u.mc_incomingE)*1e-3,
            Ereco_GeV=u.kin_cal.reco_E_lep,
            theta_deg=u.kin_cal.reco_theta_lep,
            truth_class=u.classifier.truth_class,
            nonrespi=nonrespi,flux=flux,mec=mec,rpa=rpa,geant=geant,
            minos=minos,product=product,cached=cached,
            relative_diff=relative_diff,invalid_mec_kin=int(invalid)))
        self.rows+=1

    def close(self):
        self.f.close()
        print("[CV AUDIT]",dict(self.counts),
              "max relative closure error",self.max_error,flush=True)
