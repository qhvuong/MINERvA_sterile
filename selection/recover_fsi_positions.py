#!/usr/bin/env python3
"""Recover final tables after the v1 TOTAL-row validation bug; no ROOT required."""
import argparse
import csv
import gzip
import itertools
import json
import math
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from tools.FSIProbabilityBalance import BalanceAccumulator
from tools.FSIPositions import PositionAccumulator


def write_csv(path,rows):
    with path.open('w',newline='') as f:
        if rows:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def recover(source,destination):
    source=Path(source).resolve();destination=Path(destination).resolve()
    audit=json.loads((source/'audit.json').read_text())
    if audit.get('status')!='failed' or audit.get('error')!='Position event partition does not close':
        raise RuntimeError('Recovery only supports the exact v1 final partition failure')
    destination.mkdir(parents=True,exist_ok=False)
    position=PositionAccumulator();balance=BalanceAccumulator();count=selected=0
    variant_names=('standard','pure_all','pure_retained','pure_retained_delta','without_0p70','without_1p03','extend_three_rule_to_ge4')
    budget_names=('single_hadron_balance','product_correlation','four_slot_cap','resonance_adjustment','event_combination_rules')
    with gzip.open(source/'events.csv.gz','rt',newline='') as ef,gzip.open(source/'hadrons.csv.gz','rt',newline='') as hf:
        groups=iter(itertools.groupby(csv.DictReader(hf),key=lambda h:int(h['chain_entry'])))
        pending=next(groups,None);previous=-1
        for event in csv.DictReader(ef):
            entry=int(event['chain_entry'])
            if entry<=previous:raise RuntimeError('Events not strictly ordered')
            previous=entry;records=[]
            if pending is not None and pending[0]<entry:raise RuntimeError('Unmatched hadron event')
            if pending is not None and pending[0]==entry:
                for row in pending[1]:
                    h=dict(row)
                    for k in ('pdg','fate','record_index','first_daughter','last_daughter'):h[k]=int(h[k])
                    for k in ('ke_MeV','raw1','adjusted1','effective1','start_r_fm','probe_displacement_fm','tracking_boundary_fm'):h[k]=float(h[k])
                    for k in ('eligible','retained'):
                        if h[k] not in ('True','False'):raise RuntimeError('Invalid boolean')
                        h[k]=h[k]=='True'
                    if h['category']!=event['category'] or h['stages']!=event['stages'] or h['targetA']!=event['targetA']:raise RuntimeError('Event/hadron context mismatch')
                    if h['position_group'] not in ('inside','outside'):raise RuntimeError('Invalid position group')
                    records.append(h)
                pending=next(groups,None)
            if len(records)!=int(event['preFSI_hadrons']):raise RuntimeError('Hadron count mismatch')
            stages=event['stages'].split(';');w=float(event['generator_baseline'])
            variants={k:float(event[k]) for k in variant_names};budget={k:float(event[k]) for k in budget_names}
            if not all(math.isfinite(v) for v in [w]+list(variants.values())+list(budget.values())):raise RuntimeError('Nonfinite event data')
            if not math.isclose(sum(budget.values()),variants['standard']-1,rel_tol=1e-9,abs_tol=1e-7):raise RuntimeError('Event budget mismatch')
            args=(event['category'],stages,w)
            position.add(*args,variants['standard'],records,int(event['targetA']))
            balance.add(*args,variants,budget,records,int(event['ntot']),int(event['targetA']),float(event['true_Enu_GeV']))
            count+=1;selected+=int('selected_all_EN4' in stages)
        if pending is not None:raise RuntimeError('Trailing unmatched hadron rows')
    rows=balance.event_rows();position.validate(rows)
    for name,rs in [('event_balance',rows),('hadron_balance',balance.hadron_rows()),('position_event_balance',position.event_rows()),('position_hadron_balance',position.hadron_rows())]:write_csv(destination/(name+'.csv'),rs)
    (destination/'summary.json').write_text(json.dumps(dict(totals=[r for r in rows if r['split']=='all']),indent=2)+'\n')
    audit.update(status='complete',recovered_from=str(source),recovery_reason=audit.pop('error'),cohort_events=count,selected_events=selected,position_partition_validated=True,indexed_classifier_validated=True,
                 recovery_note='Reaggregated saved records after the final v1 TOTAL-row bug. Per-event ROOT checks were performed by the original loop; not rerun here. Raw gzip files remain in recovered_from.')
    (destination/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print('Recovered',count,'events into',destination)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source');p.add_argument('destination');a=p.parse_args();recover(a.source,a.destination)
