#!/usr/bin/env python3
"""ME nu-e prediction using the existing LE Tomalak cross-section implementation.

Supply one MC playlist and an explicit exposure for a diagnostic, or supply
data POT associated with each MC playlist to form a full exposure mixture.
The default output is TRUE electron energy. Published data are efficiency-
corrected RECONSTRUCTED energy; a response is needed for a final comparison.
"""
import argparse
from array import array
import csv
import importlib.util
import json
import math
from pathlib import Path

FLAVORS = {"nue":12, "nuebar":-12, "numu":14, "numubar":-14}
HERE = Path(__file__).resolve().parent


def load_paper(path, np):
    p = json.loads(path.read_text())
    d = np.asarray(p["data_13"], dtype=float)
    c = np.asarray(p["covariance_13"], dtype=float)
    if d.shape != (13,) or c.shape != (13,13):
        raise ValueError("Expected Table IV's 13-entry data and 13x13 covariance")
    np.testing.assert_allclose(c,c.T,rtol=0,atol=1e-12)
    if np.linalg.eigvalsh(c).min() <= 0:
        raise ValueError("Published covariance is not positive definite")
    return p,d,c


def kernel(enu_edges, electron_edges, dsigma, me, ne_order, ee_order, np):
    """Integral dEnu dEe dsigma/dTe for unit flux density in each input bin."""
    xe,we = np.polynomial.legendre.leggauss(ne_order)
    xt,wt = np.polynomial.legendre.leggauss(ee_order)
    result = np.zeros((len(electron_edges)-1,len(enu_edges)-1))
    for j,(lo,hi) in enumerate(zip(enu_edges[:-1],enu_edges[1:])):
        # Split where the electron-energy integration limits change branch.
        # This avoids integrating threshold kinks with one Gaussian rule.
        thresholds=[]
        for edge in electron_edges:
            t=max(0.,float(edge)-me)
            threshold=(t+math.sqrt(t*t+2*me*t))/2
            if lo < threshold < hi:
                thresholds.append(threshold)
        breaks=sorted(set([float(lo),float(hi)]+thresholds))
        nodes=[]
        for left,right in zip(breaks[:-1],breaks[1:]):
            nodes.extend(zip((left+right)/2+(right-left)/2*xe,we*(right-left)/2))
        for enu,weight in nodes:
            if enu <= 0:
                continue
            tmax = 2*enu*enu/(me+2*enu)
            for i,(elow,ehigh) in enumerate(zip(electron_edges[:-1],electron_edges[1:])):
                tlow,thigh = max(0.,elow-me),min(tmax,ehigh-me)
                if thigh <= tlow:
                    continue
                ts=(tlow+thigh)/2+(thigh-tlow)/2*xt
                integral=sum(float(w)*dsigma(float(enu),float(t)) for w,t in zip(wt,ts))*(thigh-tlow)/2
                if not math.isfinite(integral) or integral < 0:
                    raise ValueError("Invalid integrated cross section at Enu={}".format(enu))
                result[i,j] += weight*integral
    return result


def self_test(np, paper):
    _,d,c=load_paper(paper,np)
    # For dsigma/dT=1, me=0, 1<Enu<2: these three areas are exactly 1/2.
    k=kernel(np.array([1.,2.]),np.array([0.,.5,1.,3.]),lambda e,t:1.,0.,8,16,np)
    np.testing.assert_allclose(k[:,0],[.5,.5,.5],rtol=1e-13,atol=1e-13)
    # A constant flux remains unchanged when its energy bin is subdivided.
    fine=kernel(np.array([1.,1.5,2.]),np.array([0.,.5,1.,3.]),lambda e,t:1.,0.,8,16,np)
    np.testing.assert_allclose(fine@np.ones(2),k@np.ones(1),atol=1e-13)
    # Linear exposure weighting and universe propagation.
    cv=k@np.array([2.])
    universes=np.array([[1.],[2.],[3.]])@k.T
    np.testing.assert_allclose(universes[1],cv)
    np.testing.assert_allclose(3*cv+5*cv,8*cv)
    np.testing.assert_allclose([d[:6].sum(),d[6:12].sum()],[1188.47,798.36],atol=1e-10)
    assert c[:12,:12].shape==(12,12) and np.any(c[:6,6:12]!=0)
    print("PASS: analytic cross-section integration, flux bin subdivision, exposure linearity, universe propagation, published data totals, covariance symmetry/positivity and cross-beam block.")


def read_flux(ROOT,np,path,hname,band_name,nuniv):
    f=ROOT.TFile.Open(str(path),"READ")
    if not f or f.IsZombie() or f.TestBit(ROOT.TFile.kRecovered):
        raise RuntimeError("Invalid/recovered flux file: "+str(path))
    try:
        h=f.Get(hname)
        if not h or not h.HasVertErrorBand(band_name):
            raise RuntimeError("Missing flux histogram/band in "+str(path))
        band=h.GetVertErrorBand(band_name)
        if band.GetNHists()!=nuniv:
            raise RuntimeError("{} has {} universes, expected {}".format(path,band.GetNHists(),nuniv))
        n=h.GetNbinsX()
        edges=np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1,n+2)])
        cv=np.array([h.GetBinContent(i) for i in range(1,n+1)])
        u=np.array([[band.GetHist(v).GetBinContent(i) for i in range(1,n+1)] for v in range(nuniv)])
        if not np.all(np.isfinite(edges)) or not np.all(np.diff(edges)>0) or not np.all(np.isfinite(cv)) or not np.all(np.isfinite(u)):
            raise RuntimeError("Invalid flux contents/binning")
        return edges,cv,u
    finally:
        f.Close()


def mnv_hist(ROOT,name,edges,cv,universes,band_name):
    h=ROOT.PlotUtils.MnvH1D(name,name,len(edges)-1,array("d",edges))
    h.SetDirectory(0)
    h.AddVertErrorBand(band_name,len(universes))
    band=h.GetVertErrorBand(band_name)
    for i in range(h.GetNbinsX()+2):
        value=float(cv[i-1]) if 1<=i<=len(cv) else 0.
        h.SetBinContent(i,value);h.SetBinError(i,0.)
        band.SetBinContent(i,value);band.SetBinError(i,0.)
        for u in range(len(universes)):
            hu=band.GetHist(u)
            hu.SetBinContent(i,float(universes[u,i-1]) if 1<=i<=len(cv) else 0.)
            hu.SetBinError(i,0.)
    return h


def write_matrix(ROOT,name,values):
    m=ROOT.TMatrixD(len(values),len(values[0]))
    for i,row in enumerate(values):
        for j,value in enumerate(row):
            m[i][j]=float(value)
    m.Write(name)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--beam",choices=["FHC","RHC"])
    ap.add_argument("--playlist",help="Single MC source, e.g. ME1D or ME5A")
    ap.add_argument("--pot",type=float,help="Explicit data POT for --playlist; never inferred from MC POT")
    ap.add_argument("--playlist-pot",action="append",default=[],metavar="ME1D=DATA_POT",
                    help="Repeat for the data exposures represented by distinct MC sources; total POT is their sum")
    ap.add_argument("--flux-dir",type=Path,default=Path("/exp/minerva/data/users/qvuong/flux_studies/producedFluxes_p8_ME_offsetPPFX"))
    ap.add_argument("--filename-pattern",default="{playlist}_{flavor}.root")
    ap.add_argument("--hist-name",default="flux_E_cvweighted")
    ap.add_argument("--flux-band-name",default="Flux")
    ap.add_argument("--n-universes",type=int,default=1000)
    ap.add_argument("--target-electrons",type=float,default=1.99e30)
    ap.add_argument("--flux-area-unit",choices=["m2","cm2"],default="m2")
    ap.add_argument("--tomalak-dir",type=Path,default=Path("/exp/minerva/data/users/qvuong/nueel_prediction_studies/tomalak_inputs"))
    ap.add_argument("--paper-inputs",type=Path,default=HERE/"me_nue_elastic_paper_inputs.json")
    ap.add_argument("--response-json",type=Path,help="Optional JSON: true_edges_GeV, reco_edges_GeV, matrix_reco_by_true. Pure energy migration, not selection efficiency.")
    ap.add_argument("--enu-order",type=int,default=8)
    ap.add_argument("--electron-order",type=int,default=32)
    ap.add_argument("--check-quadrature",action="store_true")
    ap.add_argument("--plots",action="store_true")
    ap.add_argument("--output",type=Path)
    ap.add_argument("--self-test",action="store_true")
    args=ap.parse_args()
    if args.self_test:
        import numpy as np
        self_test(np,args.paper_inputs)
        return 0
    if not args.beam:
        ap.error("--beam is required")
    if not math.isfinite(args.target_electrons) or args.target_electrons<=0 or args.n_universes<=0 or min(args.enu_order,args.electron_order)<2:
        ap.error("Require positive target/universe counts and quadrature orders >=2")
    if args.playlist:
        if args.playlist_pot or args.pot is None or not math.isfinite(args.pot) or args.pot<=0:
            ap.error("Use --playlist with --pot, or repeated --playlist-pot, not both")
        exposures={args.playlist:args.pot}
    else:
        if not args.playlist_pot or args.pot is not None:
            ap.error("Supply --playlist and --pot, or repeated --playlist-pot")
        exposures={}
        for item in args.playlist_pot:
            try:
                key,text_value=item.split("=",1)
                value=float(text_value)
            except ValueError:
                ap.error("Expected PLAYLIST=DATA_POT: "+item)
            if key in exposures or not math.isfinite(value) or value<=0:
                ap.error("Duplicate playlist or nonpositive data POT: "+item)
            exposures[key]=value
    for playlist in exposures:
        if not playlist.startswith(("ME1",) if args.beam=="FHC" else ("ME5","ME6")):
            ap.error("Playlist does not match beam: "+playlist)
    if args.beam=="RHC":
        print("RHC normalization is explicit: paper quotes both 1.2e21 and 1.12e21 POT; verify the intended exposure. Target electrons: {}".format(args.target_electrons),flush=True)
    import ROOT
    import PlotUtils.LoadPlotUtilsLib
    import numpy as np
    ROOT.gROOT.SetBatch(True);ROOT.TH1.AddDirectory(False)
    paper,data,cov=load_paper(args.paper_inputs,np)
    response=None
    true_edges=np.array(paper["numerical_edges_GeV"],dtype=float)
    reco_edges=true_edges.copy()
    if args.response_json:
        payload=json.loads(args.response_json.read_text())
        true_edges=np.array(payload["true_edges_GeV"],dtype=float)
        reco_edges=np.array(payload["reco_edges_GeV"],dtype=float)
        response=np.array(payload["matrix_reco_by_true"],dtype=float)
        np.testing.assert_allclose(reco_edges,paper["numerical_edges_GeV"],rtol=0,atol=1e-10)
        if response.shape!=(len(reco_edges)-1,len(true_edges)-1) or not np.all(np.isfinite(response)) or np.any(response<0) or np.any(response.sum(axis=0)>1+1e-8):
            raise ValueError("Invalid energy-migration matrix")
    if not np.all(np.isfinite(true_edges)) or not np.all(np.diff(true_edges)>0) or true_edges[0]<0:
        raise ValueError("Invalid true electron energy edges")
    if args.plots and response is None and not np.array_equal(true_edges,np.array(paper["numerical_edges_GeV"])):
        raise ValueError("Plotting truth against paper data requires matching energy bins")
    base_path=HERE/"make_nue_elastic_from_flux_TomalakXS_fixedTotalNorm_mnv.py"
    spec=importlib.util.spec_from_file_location("legacy_tomalak_xsec_for_ME",base_path)
    xs=importlib.util.module_from_spec(spec);spec.loader.exec_module(xs)
    csv_names={"nue":"electronNeutrinos","nuebar":"electronAntineutrinos","numu":"muonNeutrinos","numubar":"muonAntineutrinos"}
    cv_by_flavor={};u_by_flavor={};inputs=[];quad_checks=[]
    for flavor,pdg in FLAVORS.items():
        suffix=csv_names[flavor]
        rad=xs.RadCorrTomalak(str(args.tomalak_dir/("f1_"+suffix+".csv")),str(args.tomalak_dir/("f10_"+suffix+".csv")),flavor)
        dsigma=lambda e,t,flavor=flavor,rad=rad:xs.dsigma_dTe_tomalak(e,t,flavor,rad)
        cv_out=np.zeros(len(true_edges)-1);u_out=np.zeros((args.n_universes,len(true_edges)-1))
        kernels={}
        for playlist,pot in exposures.items():
            path=args.flux_dir/args.filename_pattern.format(playlist=playlist,tag=playlist[2:],flavor=flavor,pdg=pdg)
            edges,cv,u=read_flux(ROOT,np,path,args.hist_name,args.flux_band_name,args.n_universes)
            cache_key=tuple(edges)
            if cache_key not in kernels:
                kernels[cache_key]=kernel(edges,true_edges,dsigma,xs.ME,args.enu_order,args.electron_order,np)
            k=kernels[cache_key]
            factor=pot*args.target_electrons*(1e-4 if args.flux_area_unit=="m2" else 1.)
            contribution=factor*(k@cv)
            cv_out+=contribution
            u_out+=factor*(u@k.T)
            inputs.append(dict(flavor=flavor,playlist=playlist,path=str(path),data_POT=pot,
                               enu_min=float(edges[0]),enu_max=float(edges[-1]),input_bins=len(cv),predicted_events=float(contribution.sum())))
            if args.check_quadrature:
                finer=kernel(edges,true_edges,dsigma,xs.ME,2*args.enu_order,2*args.electron_order,np)
                refined=factor*(finer@cv)
                quad_checks.append(dict(flavor=flavor,playlist=playlist,bin_absolute_difference=(refined-contribution).tolist(),
                    max_absolute_difference=float(np.max(np.abs(refined-contribution))),
                    total_fractional_difference=float((refined.sum()-contribution.sum())/contribution.sum()) if contribution.sum() else None))
                print("  doubled-quadrature maximum bin change: {:.6g} events".format(quad_checks[-1]["max_absolute_difference"]),flush=True)
            print("{} {}: {:.6f} predicted events".format(playlist,flavor,contribution.sum()),flush=True)
        cv_by_flavor[flavor]=cv_out;u_by_flavor[flavor]=u_out
    cv_total=sum(cv_by_flavor.values());u_total=sum(u_by_flavor.values())
    if not np.all(np.isfinite(u_total)):
        raise ValueError("Nonfinite prediction")
    output=args.output or Path("me_nue_elastic_"+args.beam+".root")
    output.parent.mkdir(parents=True,exist_ok=True)
    metadata=dict(beam=args.beam,exposures=exposures,total_POT=sum(exposures.values()),target_electrons=args.target_electrons,
        input_flux_area_unit=args.flux_area_unit,input_flux_content="differential per GeV",inputs=inputs,
        xsec_module=str(base_path),xsec_total_scales=xs.TOMALAK_TOTAL_XSEC_SCALE,
        xsec_note="Same existing LE scalar Tomalak implementation and temporary total-normalization factors; not an independently updated NLO calculation.",
        observable="true total electron energy" if response is None else "energy-migration-folded prediction",
        true_edges_GeV=true_edges.tolist(),response_json=str(args.response_json) if args.response_json else None,
        model_limits=["Published data are efficiency-corrected reconstructed energy, not unfolded true energy.",
            "No detector migration is applied without --response-json. Do not use the default truth prediction as a final reco fit input.",
            "Incoming fluxes must be unconstrained by these same nu-e data to avoid double counting.",
            "Input MC statistical errors are not propagated; stored prediction bin errors are zero. Flux systematics are carried by the Flux universes.",
            "No oscillations or flux posterior constraint are applied here. Universe indices are preserved across inputs.",
            "Exact data exposure mapping to MC sources and RHC normalization require confirmation; MC POT is not substituted for data exposure."],
        paper_source=paper["source"],paper_normalization_notes=paper["normalization"],quadrature_checks=quad_checks)
    metadata["nonflux_covariance_role"]="Published measurement covariance supplies the data/non-flux term. Do not add another data statistical term; the prediction Flux band is separate."
    fout=ROOT.TFile(str(output),"RECREATE")
    if fout.IsZombie():
        raise RuntimeError("Cannot create "+str(output))
    keep=[]
    for flavor in ["total"]+list(FLAVORS):
        cv=cv_total if flavor=="total" else cv_by_flavor[flavor]
        u=u_total if flavor=="total" else u_by_flavor[flavor]
        hist=mnv_hist(ROOT,"h_nue_elastic_"+flavor,true_edges,cv,u,args.flux_band_name)
        hist.GetXaxis().SetTitle("True total electron energy (GeV)")
        fout.cd();hist.Write();keep.append(hist)
        if response is not None:
            folded=mnv_hist(ROOT,"h_nue_elastic_reco_"+flavor,reco_edges,response@cv,u@response.T,args.flux_band_name)
            folded.GetXaxis().SetTitle("Reconstructed electron energy (GeV)")
            fout.cd();folded.Write();keep.append(folded)
    paper_edges=paper["numerical_edges_GeV"]
    for beam,start in [("FHC",0),("RHC",6)]:
        h=ROOT.TH1D("paper_data_"+beam,"Published data: total measurement errors",6,array("d",paper_edges))
        h.SetDirectory(0)
        for i in range(6):
            h.SetBinContent(i+1,float(data[start+i]));h.SetBinError(i+1,float(math.sqrt(cov[start+i,start+i])))
        fout.cd();h.Write();keep.append(h)
    write_matrix(ROOT,"paper_measurement_covariance_12",cov[:12,:12])
    write_matrix(ROOT,"paper_measurement_covariance_FHC",cov[:6,:6])
    write_matrix(ROOT,"paper_measurement_covariance_RHC",cov[6:12,6:12])
    # Explicit aliases for the role these measurement matrices have in the fit.
    write_matrix(ROOT,"cov_nonflux_ME_FHC_RHC",cov[:12,:12])
    write_matrix(ROOT,"cov_nonflux_ME_FHC",cov[:6,:6])
    write_matrix(ROOT,"cov_nonflux_ME_RHC",cov[6:12,6:12])
    ROOT.TObjString(json.dumps(metadata)).Write("prediction_metadata_json")
    fout.Close()
    output.with_suffix(".json").write_text(json.dumps(metadata,indent=2))
    with output.with_suffix(".csv").open("w",newline="") as stream:
        writer=csv.writer(stream);writer.writerow(["energy_low_GeV","energy_high_GeV","true_total"]+list(FLAVORS))
        for i in range(len(cv_total)):
            writer.writerow([true_edges[i],true_edges[i+1],cv_total[i]]+[cv_by_flavor[f][i] for f in FLAVORS])
    np.savetxt(output.with_name(output.stem+"_paper_covariance_12.csv"),cov[:12,:12],delimiter=",")
    if args.plots:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig,(ax,ratio)=plt.subplots(2,1,figsize=(9,7),sharex=True,gridspec_kw={"height_ratios":[3,1]},layout="constrained")
        pred=cv_total if response is None else response@cv_total
        pu=u_total if response is None else u_total@response.T
        pedges=true_edges if response is None else reco_edges
        ax.stairs(pred,pedges,color="#007f8b",label="True-energy prediction" if response is None else "Response-folded prediction")
        ax.stairs(pred+np.std(pu,axis=0),pedges,color="#007f8b",ls=":",label="Flux ensemble standard deviation")
        ax.stairs(np.maximum(0,pred-np.std(pu,axis=0)),pedges,color="#007f8b",ls=":")
        sl=slice(0,6) if args.beam=="FHC" else slice(6,12)
        centers=(np.array(paper_edges[:-1])+np.array(paper_edges[1:]))/2
        centers[-1]=14.5
        errors=np.sqrt(np.diag(cov)[sl]);observed=data[sl]
        ax.errorbar(centers,observed,yerr=errors,fmt="ko",label="Published data (total measurement error)")
        if len(pred)==6:
            ratio.errorbar(centers,np.divide(observed,pred,out=np.full(6,np.nan),where=pred>0),
                yerr=np.divide(errors,pred,out=np.full(6,np.nan),where=pred>0),fmt="ko")
        ratio.axhline(1,color="gray",ls=":")
        ratio.set(xlabel="Electron energy (GeV); final bin includes all energies above 9 GeV",ylabel="Data / prediction",xlim=(.8,20))
        ax.set(ylabel="Events per published bin",title="ME {}: diagnostic nu-e prediction, {:.4g} POT".format(args.beam,sum(exposures.values())))
        ax.legend(fontsize=9)
        if response is None:
            ax.text(.02,.02,"Diagnostic only: true-energy model vs efficiency-corrected reco data; no migration applied.",transform=ax.transAxes,fontsize=8)
        fig.savefig(output.with_suffix(".png"),dpi=170);plt.close(fig)
    print("Wrote {} | total prediction {:.6f} | no fitted normalization".format(output,cv_total.sum()))
    if args.beam=="RHC":
        print("NOTE: publication quotes both 1.2e21 and 1.12e21 POT; confirm exact exposure and target normalization before final use.")
    return 0


if __name__=="__main__":
    raise SystemExit(main())
