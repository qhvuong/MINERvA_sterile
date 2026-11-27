"""Diagnostic-only inspection of the installed FSI code and probability accounting."""
import ast
import math
import re
import struct
from pathlib import Path

SMALL_BOOST = struct.unpack('f',struct.pack('f',1.03))[0]
KE_EDGES = (0.,25.,30.,35.,40.,50.,75.,100.,125.,150.,200.,300.,500.,1000.,2000.,3800.)


def cpp_class(header):
    # Ignore comments/strings when finding class braces, preserving offsets.
    scrub = re.sub(r'//[^\n]*|/\*[\s\S]*?\*/|"(?:\\.|[^"\\])*"',
                   lambda m:' '*len(m.group()),header)
    match = re.search(r'\bclass\s+weight_fsi\b',scrub)
    if not match:raise RuntimeError('weight_fsi class not found')
    opening = scrub.index('{',match.end());level=0
    for i in range(opening,len(scrub)):
        if scrub[i]=='{':level+=1
        elif scrub[i]=='}':
            level-=1
            if level==0:
                return header[match.start():i+1]+';'
    raise RuntimeError('Unterminated FSI class')


def diagnostic_cpp(header,source,installed_header):
    """Clone the actual installed calculator under a new class name.

    Only access control and the class name change. The loaded production class
    and source files remain untouched. The clone is validated event by event.
    """
    declaration=cpp_class(header).replace('private:','public:')
    match=re.search(r'PlotUtils::weight_fsi::weight_fsi\s*\(',source)
    if not match:raise RuntimeError('FSI constructor not found')
    body=source[match.start():]
    declaration=re.sub(r'\bweight_fsi\b','FSIBalanceCalculator',declaration)
    body=re.sub(r'\bweight_fsi\b','FSIBalanceCalculator',body)
    import json
    return '\n'.join([
        '#include <iostream>','#include <cmath>','#include <algorithm>',
        '#include <vector>','#include <TMath.h>','#include "PlotUtils/ChainWrapper.h"',
        '#include '+json.dumps(str(Path(installed_header).resolve())),
        'using std::sqrt;', 'namespace PlotUtils { '+declaration+' }',body,
        '''namespace PlotUtils {
        std::vector<double> FSIBalanceRawWeights(FSIBalanceCalculator& calc,
                  int fate, double ke, int A, int pdg) {
            double w1=-1., w2=-1.;
            calc.GetQEElasticWeight(w1,w2,fate,ke,A,pdg);
            return {w1,w2};
        }}'''])


def generator_function(source,namespace):
    """Production V430 with only final GEANT/MINOS blocks removed."""
    tree=ast.parse(source)
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='CVPythonUniverse')
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='GetStandardWeightV4p3p0')
    expected=[('V4P3P0_ADD_GEANT_HADRON','GetGeantHadronWeight'),
              ('V4P3P0_ADD_MINOS_EFFICIENCY','GetMyMinosEfficiencyWeight')]
    for node,(flag,method) in zip(fn.body[-3:-1],expected):
        exemplar=ast.parse('if {}:\n    weight *= self.{}()'.format(flag,method)).body[0]
        if ast.dump(node,include_attributes=False)!=ast.dump(exemplar,include_attributes=False):
            raise RuntimeError('Production detector-weight blocks changed')
    if len(fn.body[-3:-1])!=2 or not isinstance(fn.body[-1],ast.Return):
        raise RuntimeError('Production return changed')
    fn.body=fn.body[:-3]+fn.body[-1:];fn.name='_fsi_balance_generator_weight'
    module=ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[]))
    private=dict(namespace)
    exec(compile(module,'<production generator+flux baseline>','exec'),private)
    return private[fn.name]


def combine(pions,nucleons,int_type,factor=.70,boost=SMALL_BOOST,extend=False):
    """Independent transcription of the standard event combination; checked against C++."""
    pi=(list(pions)+[-1.]*4)[:4];nu=(list(nucleons)+[-1.]*4)[:4]
    ntot=sum(x>=-.001 for x in pi+nu)
    w=1.
    if ntot==1:
        w=nu[0] if nu[0]>=-.001 else pi[0]
    elif ntot==2:
        for i in range(4):
            if pi[i]>=-.001:w*=pi[i]
            if int_type==2 and abs(pi[i]-1.)<=.0001:w*=1.1
            if abs(nu[i])<.001:w*=0.
            elif nu[i]>=.999:w*=1.+(nu[i]-1.)*factor
    elif ntot==3 or extend and ntot>=4:
        for i in range(4):
            if abs(pi[i])<.001:w*=0.
            elif pi[i]>=1.005:w*=pi[i]
            elif pi[i]>=.999:w*=boost
            if abs(nu[i])<.001:w*=0.
            elif nu[i]>=1.005:w*=1.+(nu[i]-1.)*factor
            elif nu[i]>=.999:w*=boost
    return w


def inspect_hadrons(hadrons,int_type,res_id,charm,target_a):
    """Input hadrons sorted by descending pre-FSI KE, with exact C++ raw weights."""
    pi=[];nu=[];records=[];pure_all=1.;pure_retained=1.;pure_delta=1.;linear=0.
    enabled=not charm and target_a>=4
    for h in hadrons:
        pdg,fate,ke,raw=h['pdg'],h['fate'],h['ke_MeV'],h['raw1']
        if not all(math.isfinite(x) for x in (ke,raw)) or ke<0:
            raise RuntimeError('Invalid intranuclear hadron: {}'.format(h))
        adjusted=raw
        if int_type==2 and res_id==0 and 100.<ke<300. and raw>1.05:
            adjusted=1.+(raw-1.)*.90
        candidate=enabled and fate in (1,3) and not(pdg==2212 and ke<125.)
        usable=candidate and raw>=-.001
        effective=raw if usable else 1.
        pure_all*=effective;linear+=effective-1.
        retained=False
        if candidate:
            weights=nu if pdg in (2212,2112) else pi
            if len(weights)<4:
                weights.append(adjusted);retained=True
        if retained and raw>=-.001:
            pure_retained*=raw;pure_delta*=adjusted
        records.append(dict(h,adjusted1=adjusted,effective1=effective,
                            retained=retained,eligible=usable))
    ntot=sum(x>=-.001 for x in pi+nu)
    standard=combine(pi,nu,int_type)
    budget=dict(single_hadron_balance=linear,
                product_correlation=pure_all-1.-linear,
                four_slot_cap=pure_retained-pure_all,
                resonance_adjustment=pure_delta-pure_retained,
                event_combination_rules=standard-pure_delta)
    if not math.isclose(sum(budget.values()),standard-1.,rel_tol=1e-10,abs_tol=1e-10):
        raise RuntimeError('Event decomposition failed')
    variants=dict(standard=standard,pure_all=pure_all,pure_retained=pure_retained,
                  pure_retained_delta=pure_delta,
                  without_0p70=combine(pi,nu,int_type,factor=1.),
                  without_1p03=combine(pi,nu,int_type,boost=1.),
                  extend_three_rule_to_ge4=combine(pi,nu,int_type,extend=True))
    if not all(math.isfinite(x) for x in list(variants.values())+list(budget.values())):
        raise RuntimeError('Nonfinite diagnostic product/budget')
    return records,variants,budget,ntot


def energy_bin(value):
    import bisect
    if not math.isfinite(value) or value<0:raise ValueError('Invalid kinetic energy')
    i=bisect.bisect_right(KE_EDGES,value)
    return '{}_to_{}'.format(KE_EDGES[i-1],KE_EDGES[i] if i<len(KE_EDGES) else 'inf')


class BalanceAccumulator:
    def __init__(self):
        self.events={}
        self.hadrons={}

    def add(self,category,stages,baseline,variants,budget,records,ntot,target_a,true_energy):
        if not math.isfinite(baseline):raise ValueError('Nonfinite baseline')
        # Event variance accounting, including correlated hadrons within an event.
        dimensions=[('all','all'),('contributors',str(ntot)),
                    ('preFSI_hadrons',str(len(records))),('targetA',str(target_a)),
                    ('true_Enu_GeV', '{}_to_{}'.format(int(true_energy//2)*2,int(true_energy//2)*2+2))]
        for cat in ('TOTAL',category):
            for stage in stages:
                for basis,w in (('unit',1.),('generator_flux',baseline)):
                    for split,value in dimensions:
                        key=(cat,stage,basis,split,value)
                        row=self.events.setdefault(key,dict(events=0,sum_baseline=0.))
                        row['events']+=1;row['sum_baseline']+=w
                        for name,v in variants.items():
                            x=w*v
                            row['sum_'+name]=row.get('sum_'+name,0.)+x
                            row['sum_delta2_'+name]=row.get('sum_delta2_'+name,0.)+(w*(v-1.))**2
                        for name,v in budget.items():
                            row['budget_'+name]=row.get('budget_'+name,0.)+w*v
                    temporary={}
                    for h in records:
                        # KE and multiplicity refer to hadrons before FSI, independent of fate.
                        hkey=(cat,stage,basis,str(target_a),str(h['pdg']),energy_bin(h['ke_MeV']),str(len(records)))
                        r=temporary.setdefault(hkey,dict(events_with_hadrons=1,hadrons=0,
                            no_scatter=0,elastic=0,other_fates=0,eligible=0,retained=0,
                            excluded_negative_or_sentinel=0,gain=0.,loss=0.,balance=0.,sum_hadron_baseline=0.))
                        r['hadrons']+=1
                        r['no_scatter' if h['fate']==1 else 'elastic' if h['fate']==3 else 'other_fates']+=1
                        r['eligible']+=int(h['eligible']);r['retained']+=int(h['retained'])
                        r['excluded_negative_or_sentinel']+=int(h['raw1']<-.001)
                        delta=w*(h['effective1']-1.)
                        r['gain']+=max(0.,delta);r['loss']+=min(0.,delta);r['balance']+=delta
                        r['sum_hadron_baseline']+=w
                    for key,r in temporary.items():
                        dest=self.hadrons.setdefault(key,{name:0 for name in r})
                        for name,x in r.items():dest[name]+=x
                        dest['sum_balance2_by_event']=dest.get('sum_balance2_by_event',0.)+r['balance']**2

    def event_rows(self):
        rows=[]
        marginal={}
        for key,row in self.hadrons.items():
            marginal[key[:3]]=marginal.get(key[:3],0.)+row['balance']
        for key,r in sorted(self.events.items()):
            row=dict(zip(('category','stage','basis','split','value'),key));row.update(r)
            baseline=r['sum_baseline']
            for name in list(r):
                if name.startswith('sum_') and not name.startswith('sum_delta2_') and name!='sum_baseline':
                    row['ratio_'+name[4:]]=r[name]/baseline if baseline else None
            budget=sum(v for k,v in r.items() if k.startswith('budget_'))
            expected=r['sum_standard']-baseline
            if not math.isclose(budget,expected,rel_tol=1e-9,abs_tol=1e-7):
                raise RuntimeError('Summed event budget failed')
            if key[3]=='all' and not math.isclose(r['budget_single_hadron_balance'],
                    marginal.get(key[:3],0.),rel_tol=1e-9,abs_tol=1e-7):
                raise RuntimeError('Hadron-table/event-table accounting mismatch')
            rows.append(row)
        return rows

    def hadron_rows(self):
        rows=[]
        for key,r in sorted(self.hadrons.items()):
            row=dict(zip(('category','stage','basis','targetA','pdg','KE_MeV','preFSI_hadrons'),key))
            row.update(r)
            row['mean_effective_weight']=1.+r['balance']/r['sum_hadron_baseline'] if r['sum_hadron_baseline'] else None
            row['MC_balance_scale']=math.sqrt(r['sum_balance2_by_event'])
            rows.append(row)
        return rows
