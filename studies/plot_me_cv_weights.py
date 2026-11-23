#!/usr/bin/env python3
"""Plot saved ME CV diagnostics using NumPy/Matplotlib only; no ROOT or config.

Usage: python plot_me_cv_weights.py /path/to/ME1D_ryan
"""
import argparse
import csv
import gzip
import json
from pathlib import Path
import numpy as np

COMPONENTS=['genie','flux','2p2h','RPA','lowQ2pi','MEcoherent','diffractive','geant','MINOS','elasticFSI']


def plot_directory(directory,output=None):
    directory=Path(directory)
    output=Path(output) if output else directory
    info=json.loads((directory/'run.json').read_text())
    if info.get('status')!='COMPLETED':
        raise RuntimeError('Diagnostic run did not complete; inspect run.json')
    arrays=[];categories=[]
    with gzip.open(directory/'selected_events.csv.gz','rt') as f:
        for row in csv.DictReader(f):
            arrays.append([float(row[k]) for k in COMPONENTS])
            categories.append(row['category'])
    values=np.asarray(arrays,dtype=float)
    cats=np.asarray(categories)
    if len(values)!=info['selected']:
        raise RuntimeError('Event CSV does not match the recorded selected-event count')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    output.mkdir(parents=True,exist_ok=True)
    shown=[c for c in ['CCNuMuQE','CCNuMu2p2h','CCNuMuDelta','CCNuMuDIS','CCNuMuWrongSign'] if np.any(cats==c)]
    fig,axes=plt.subplots(2,5,figsize=(17,7),constrained_layout=True)
    for j,(ax,name) in enumerate(zip(axes.flat,COMPONENTS)):
        for cat in shown:
            x=values[cats==cat,j]
            ax.hist(x,bins=np.linspace(0,5,101),histtype='step',weights=np.ones(len(x))/len(x),label=cat.replace('CCNuMu',''))
        ax.axvline(1,color='gray',ls=':')
        ax.set(title=name,xlabel='Component weight',ylabel='Fraction / 0.05',yscale='log',xlim=(0,5))
    axes.flat[0].legend(fontsize=7)
    fig.suptitle('ME '+info['playlist']+' — selected-event CV components; tails above 5 listed in CSV')
    fig.savefig(output/'component_distributions.png',dpi=160);plt.close(fig)
    with np.load(directory/'spectra.npz',allow_pickle=False) as saved:
        edges=saved['edges']
        fig,ax=plt.subplots(figsize=(10,5),constrained_layout=True)
        for label in ['All',*shown]:
            n=saved[label+'_w431'];d=saved[label+'_w430']
            ax.stairs(np.divide(n,d,out=np.full_like(n,np.nan),where=d!=0),edges,label=label.replace('CCNuMu',''))
        ax.axhline(1,color='gray',ls='--')
        ax.set(xlabel='Neutrino Energy Estimator (GeV)',ylabel='v4.3.1 / v4.3.0',title='ME '+info['playlist']+' — same events, FSI effect',xlim=(.8,20))
        ax.legend();ax.grid(alpha=.2)
        fig.savefig(output/'FSI_effect_by_category.png',dpi=180);plt.close(fig)
    print('Plots:',output)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory',type=Path)
    parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args()
    plot_directory(args.directory,args.output_dir)
