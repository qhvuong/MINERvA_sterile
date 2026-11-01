#!/usr/bin/env python3
"""Read-only neutron CV diagnostic from existing event-level audit CSV.
Does not change production, does not calculate efficiencies, and is not POT normalized.
Requires columns: pdg, intType, Ereco_GeV, geant, cached.
"""
import argparse
import csv
import math
from collections import defaultdict

EDGES = [0.8,2.0,2.5,3.0,3.8,4.6,5.5,6.5,8.0,10.0,12.5,16.0,20.0]

def binned(x):
    for lo,hi in zip(EDGES,EDGES[1:]):
        if lo <= x < hi or (hi == EDGES[-1] and x == hi):
            return f'{lo:g}-{hi:g}'
    return None

def bucket():
    return {'count':0,'cv':0.0,'off':0.0,'wmin':float('inf'),'wmax':-float('inf'),'changed':0}

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--audit',required=True)
    ap.add_argument('--output',default='neutron_cv_impact.csv')
    args=ap.parse_args()
    grouped=defaultdict(bucket)
    invalid=defaultdict(int)
    n=0
    with open(args.audit,newline='') as f:
        for r in csv.DictReader(f):
            n+=1
            try:
                pdg=int(float(r['pdg']))
                typ=int(float(r['intType']))
                e=float(r['Ereco_GeV'])
                w=float(r['geant'])
                cv=float(r['cached'])
            except (KeyError,TypeError,ValueError):
                invalid['parse']+=1;continue
            if not all(math.isfinite(v) for v in (e,w,cv)) or w<=0:
                invalid['nonfinite_or_nonpositive']+=1;continue
            species='numubar' if pdg==-14 else 'numu' if pdg==14 else 'other'
            keys=[('all','all'),(species,'all'),(species,f'intType={typ}')]
            b=binned(e)
            if b:keys.append((species,b))
            else:invalid['outside_E_mu_bins']+=1
            for k in keys:
                d=grouped[k]
                d['count']+=1;d['cv']+=cv;d['off']+=cv/w
                d['wmin']=min(d['wmin'],w);d['wmax']=max(d['wmax'],w)
                d['changed']+=int(abs(w-1)>1.e-8)
    fields=['species','region','n','cv_sum','no_neutron_cv_sum','cv_over_no_cv','impact_percent','changed_count','min_w','max_w']
    with open(args.output,'w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        for (species,region),d in sorted(grouped.items()):
            ratio=d['cv']/d['off'] if d['off']!=0 else float('nan')
            row=dict(species=species,region=region,n=d['count'],cv_sum=d['cv'],no_neutron_cv_sum=d['off'],cv_over_no_cv=ratio,impact_percent=(ratio-1)*100,changed_count=d['changed'],min_w=d['wmin'],max_w=d['wmax'])
            writer.writerow(row)
            print(f'{species:9} {region:13} n={d["count"]:5d} CV={d["cv"]:10.2f} noCV={d["off"]:10.2f} impact={(ratio-1)*100:+7.2f}% changed={d["changed"]:5d} geant_range=[{d["wmin"]:.3g},{d["wmax"]:.3g}]')
    print(f'\nRead {n} rows. Issues: {dict(invalid)}. Output: {args.output}')
    print('IMPORTANT: Counterfactual weight removal only for audited selected events, not an efficiency or a complete production comparison.')

if __name__=='__main__':main()
