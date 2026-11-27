#!/usr/bin/env python3
"""Print event partitions and carbon-neutron fate balances from one completed run."""
import csv
import json
import sys
from pathlib import Path
p=Path(sys.argv[1])
a=json.loads((p/'audit.json').read_text())
if a.get('status')!='complete':raise SystemExit('Audit is not complete')
print('Playlist/tree/PDG:',a['playlist'],a['tree'],a['AnaNuPDG'])
print('Assumed geometry:',a['position_geometry'])
rows=list(csv.DictReader((p/'position_event_balance.csv').open()))
for basis in ('unit','generator_flux'):
 for scope in ('all_species','neutrons'):
  selected=[r for r in rows if r['category']=='DIS' and r['stage']=='before' and r['basis']==basis and r['scope']==scope]
  baseline=sum(float(r['baseline']) for r in selected);delta=sum(float(r['delta']) for r in selected)
  print('\nDIS',basis,scope,'total change %',100*delta/baseline if baseline else None)
  for group in ('with_outside','without_outside'):
   subset=[r for r in selected if r['event_group']==group]
   d=sum(float(r['delta']) for r in subset)
   print(group,'events',sum(int(r['events']) for r in subset),'delta',d,'percentage points',100*d/baseline if baseline else None)
print('\nCarbon neutron single-hadron balance: category KE position N_noScatter N_elastic gain loss balance')
for r in csv.DictReader((p/'position_hadron_balance.csv').open()):
 if r['category'] in ('QE','DIS') and r['stage']=='before' and r['basis']=='unit' and r['targetA']=='12' and r['pdg']=='2112':
  print(*(r[k] for k in ['category','KE_MeV','position_group','no_scatter','elastic','gain','loss','balance']))
