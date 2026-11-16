#!/usr/bin/env python3
"""Plot native-tuple-derived vs legacy ME BeamFocus response, without ROOT.

Uses the exact common-bin comparison CSV, with optional direct tuple-audit
points. Universe ensembles are compared by their marginal spread, not by
assuming the legacy and tuple universe indices are paired.
"""
import argparse
import csv
import json
from pathlib import Path


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def make_plots(comparison, audits, playlists, emin, emax, out):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    names={12:r"$\nu_e$",-12:r"$\bar\nu_e$",14:r"$\nu_\mu$",-14:r"$\bar\nu_\mu$"}
    for playlist in playlists:
        fig,axs=plt.subplots(4,2,figsize=(13,12),sharex=True,layout="constrained")
        for row,(pdg,name) in enumerate(names.items()):
            key="gen2thin-pdg{}-minervame{}".format(pdg,playlist[2:])
            rows=sorted([r for r in comparison if r["case"]==key and r["band"]=="Flux_BeamFocus"
                         and float(r["energy_low"])>=emin-1e-8 and float(r["energy_high"])<=emax+1e-8],
                        key=lambda r:float(r["energy_low"]))
            if not rows:
                raise ValueError("Missing BeamFocus comparison: "+key)
            # Draw each bin independently so gaps never acquire invented values.
            for side,color,label in [("produced","#007f8b","Produced (native tuple response)"),
                                     ("reference","#c54b32","Legacy ME flux")]:
                for j,r in enumerate(rows):
                    edges=[float(r["energy_low"]),float(r["energy_high"])]
                    cv=float(r[side+"_cv"])
                    if cv<=0:
                        continue
                    rms=float(r[side+"_weight_rms"])
                    if not np.isfinite(rms):
                        continue
                    axs[row,0].plot(edges,[100*rms]*2,color=color,lw=1.6,label=label if j==0 else None)
                    q16=float(r[side+"_weight_q16"]);q84=float(r[side+"_weight_q84"])
                    median=float(r[side+"_weight_median"])
                    axs[row,1].fill_between(edges,[q16]*2,[q84]*2,color=color,alpha=.18,
                                            label=label+": 16–84%" if j==0 else None)
                    axs[row,1].plot(edges,[median]*2,color=color,lw=1)
            for audit_index,(label,records) in enumerate(audits):
                selected=[r for r in records if r["playlist"]==playlist and int(r["pdg"])==pdg
                          and float(r["energy_low"])>=emin-1e-8 and float(r["energy_high"])<=emax+1e-8]
                if not selected:
                    continue
                x=np.array([(float(r["energy_low"])+float(r["energy_high"]))/2 for r in selected])
                y=np.array([float(r["tuple_RMS_percent"]) for r in selected])
                marker=["o","x","+"][audit_index%3]
                color=["black","#794ca4","#c89300"][audit_index%3]
                axs[row,0].plot(x,y,ls="none",marker=marker,ms=3.5,color=color,label=label)
                if all("tuple_q16" in r for r in selected):
                    low=np.array([float(r["tuple_q16"]) for r in selected])
                    mid=np.array([float(r["tuple_median"]) for r in selected])
                    high=np.array([float(r["tuple_q84"]) for r in selected])
                    axs[row,1].errorbar(x,mid,yerr=np.array([mid-low,high-mid]),fmt=marker,
                                       ms=3,color=color,lw=.65,label=label+": median / 16–84%")
            axs[row,0].set_ylabel(name+"\nRMS about weight 1 (%)")
            axs[row,1].set_ylabel(name+"\nUniverse / CV weight")
            axs[row,1].axhline(1,color="gray",ls=":",lw=.8)
            for ax in axs[row]:
                ax.set_xlim(emin,emax);ax.grid(alpha=.15)
                ax.legend(fontsize=7,framealpha=.7)
        for ax in axs[-1]:
            ax.set_xlabel("True neutrino energy (GeV)")
        axs[0,0].set_title("Beam-focusing uncertainty size")
        axs[0,1].set_title("Median line and central 68% of universe weights")
        fig.suptitle(playlist+": BeamFocus response only\nWhole-bin comparison; no pairing of legacy and tuple universe indices",fontsize=13)
        for extension in ("png","pdf"):
            fig.savefig(out/(playlist+"_beamfocus_response."+extension),dpi=170)
        plt.close(fig)
        print("Saved "+playlist+" BeamFocus response comparison")


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--comparison-dir",type=Path,default=Path("me_flux_common_bin_comparison"))
    ap.add_argument("--audit-dirs",type=Path,nargs="*",default=[])
    ap.add_argument("--playlists",nargs="+",default=["ME1D","ME1M","ME1N","ME5A","ME6A"])
    ap.add_argument("--emin",type=float,default=.5)
    ap.add_argument("--emax",type=float,default=20.)
    ap.add_argument("--out",type=Path,default=Path("me_beamfocus_response_plots"))
    args=ap.parse_args()
    if args.emax<=args.emin:
        ap.error("Require emax > emin")
    comparison=read_csv(args.comparison_dir/"weights_by_energy.csv")
    audits=[]
    for directory in args.audit_dirs:
        meta=json.loads((directory/"inputs.json").read_text())
        selection=meta.get("selection","")
        label="Direct tuple: FV only" if selection.startswith("fv-only") else "Direct tuple: producer selection"
        audits.append((label,read_csv(directory/"summary.csv")))
    args.out.mkdir(parents=True,exist_ok=True)
    make_plots(comparison,audits,args.playlists,args.emin,args.emax,args.out)
    (args.out/"plot_inputs.json").write_text(json.dumps(dict(comparison_dir=str(args.comparison_dir),
        audit_dirs=[str(p) for p in args.audit_dirs],playlists=args.playlists,
        note="Produced curves use the production carbon-CC FV selection; they are not labelled FV-only. Direct audit points retain their own selection. Low-CV bins remain visible; use audit event counts to assess support."),indent=2))


if __name__=="__main__":
    main()
