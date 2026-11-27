"""Read-only FSI position diagnostics; geometry assumes explicitly recorded defaults."""
import math
import re
from collections import defaultdict
from tools.FSIProbabilityBalance import energy_bin


def extract_class(header, name):
    start = re.search(r'\bclass\s+' + name + r'\s*\{', header)
    if not start:
        raise RuntimeError('Missing class ' + name)
    opening = header.index('{', start.start()); depth = 0
    for i in range(opening, len(header)):
        if header[i] == '{': depth += 1
        elif header[i] == '}':
            depth -= 1
            if depth == 0: return header[start.start():i+1] + ';'
    raise RuntimeError('Unterminated class ' + name)


def indexed_fate_cpp(header, source):
    """Clone the installed classifier, carrying the original index through its sort."""
    particle = extract_class(header, 'intranuke_particle')
    particle = particle.replace('int __pdg;', 'int __record_index;\nint __pdg;')
    if '__record_index' not in particle: raise RuntimeError('Particle declaration changed')
    fate = extract_class(header, 'fate_code')
    start = source.index('PlotUtils::fate_code::fate_code()')
    end = source.index('PlotUtils::weight_fsi::weight_fsi()')
    body = source[start:end]
    needle = 'fIntranukeParticles.push_back(particle);'
    if body.count(needle) != 1: raise RuntimeError('Classifier insertion changed')
    body = body.replace(needle, 'particle.__record_index = i;\n' + needle)
    if body.count('decreasingKE()') != 1: raise RuntimeError('Classifier sorting changed')
    body = body.replace('decreasingKE()', '[](const intranuke_particle& a, const intranuke_particle& b) { return a.__ke > b.__ke; }')
    text = 'namespace PlotUtils {\n' + particle + '\n' + fate + '\n}\n' + body
    for old, new in [('intranuke_particle','PositionIntranukeParticle'),('fate_code','PositionFateCode')]:
        text = re.sub(r'\b' + old + r'\b',new,text)
    return '#include <algorithm>\n#include <cmath>\n#include <cstdio>\nusing std::sqrt;\n' + text


POSITION_FIELDS = ['record_index','first_daughter','last_daughter','start_x_fm','start_y_fm','start_z_fm',
 'start_r_fm','probe_displacement_fm','tracking_boundary_fm','position_group','daughter_r_fm']


def geometry(index, ids, statuses, fd, ld, positions, incoming, target_a, r0, nr, step):
    n=len(ids)
    if any(len(v)!=n for v in [statuses,fd,ld]+list(positions)):
        raise RuntimeError('Generator vector length mismatch')
    if not 0<=index<n: raise RuntimeError('Invalid original hadron index')
    probes=[i for i in range(n) if ids[i]==incoming and statuses[i]==0]
    if len(probes)!=1: raise RuntimeError('Expected one initial-state probe')
    xyz=[float(v[index]) for v in positions]
    probe=[float(v[probes[0]]) for v in positions]
    if not all(math.isfinite(x) for x in xyz+probe): raise RuntimeError('Nonfinite position')
    radius=math.sqrt(sum(x*x for x in xyz))
    displacement=math.sqrt(sum((x-y)**2 for x,y in zip(xyz,probe)))
    boundary=r0*target_a**(1./3.)*nr+step
    daughter=int(fd[index]); daughter_r=''
    if 0<=daughter<n:
        coords=[float(v[daughter]) for v in positions]
        if not all(math.isfinite(x) for x in coords): raise RuntimeError('Nonfinite daughter position')
        daughter_r=math.sqrt(sum(x*x for x in coords))
    return dict(record_index=index,first_daughter=daughter,last_daughter=int(ld[index]),
      start_x_fm=xyz[0],start_y_fm=xyz[1],start_z_fm=xyz[2],start_r_fm=radius,
      probe_displacement_fm=displacement,tracking_boundary_fm=boundary,
      position_group='outside' if radius>=boundary else 'inside',daughter_r_fm=daughter_r)


class PositionAccumulator:
    def __init__(self): self.hadrons={}; self.events={}

    def add(self,category,stages,baseline,standard,records,target_a):
        for stage in stages:
            for basis,w in [('unit',1.),('generator_flux',baseline)]:
                for scope in ['all_species','neutrons']:
                    outside=any(h['eligible'] and h['position_group']=='outside' and
                                (scope=='all_species' or h['pdg']==2112) for h in records)
                    key=(category,stage,basis,str(target_a),scope,'with_outside' if outside else 'without_outside')
                    r=self.events.setdefault(key,dict(events=0,baseline=0.,corrected=0.,delta=0.))
                    r['events']+=1;r['baseline']+=w;r['corrected']+=w*standard;r['delta']+=w*(standard-1.)
                for h in records:
                    key=(category,stage,basis,str(target_a),str(h['pdg']),energy_bin(h['ke_MeV']),h['position_group'])
                    r=self.hadrons.setdefault(key,dict(hadrons=0,no_scatter=0,elastic=0,other_fates=0,
                      eligible=0,gain=0.,loss=0.,balance=0.,sum_displacement_fm=0.))
                    r['hadrons']+=1;r['eligible']+=int(h['eligible'])
                    r['no_scatter' if h['fate']==1 else 'elastic' if h['fate']==3 else 'other_fates']+=1
                    delta=w*(h['effective1']-1.)
                    r['gain']+=max(0.,delta);r['loss']+=min(0.,delta);r['balance']+=delta
                    r['sum_displacement_fm']+=h['probe_displacement_fm']

    def hadron_rows(self):
        keys=('category','stage','basis','targetA','pdg','KE_MeV','position_group')
        return [dict(zip(keys,k),**r) for k,r in sorted(self.hadrons.items())]

    def event_rows(self):
        keys=('category','stage','basis','targetA','scope','event_group')
        return [dict(zip(keys,k),**r) for k,r in sorted(self.events.items())]

    def validate(self, reference):
        # Each scope partitions entire events once, including events with zero eligible hadrons.
        for row in reference:
            if row['split']!='all':continue
            for scope in ['all_species','neutrons']:
                selected=[r for k,r in self.events.items() if (row['category']=='TOTAL' or k[0]==row['category']) and k[1:3]==(row['stage'],row['basis']) and k[4]==scope]
                for dest,src in [('baseline','sum_baseline'),('corrected','sum_standard')]:
                    if not math.isclose(sum(r[dest] for r in selected),row[src],rel_tol=1e-9,abs_tol=1e-7):
                        raise RuntimeError('Position event partition does not close: {} {} {} {}'.format(row['category'],row['stage'],row['basis'],scope))
                if sum(r['events'] for r in selected)!=row['events']:raise RuntimeError('Position event count mismatch')
                htotal=sum(r['balance'] for k,r in self.hadrons.items() if (row['category']=='TOTAL' or k[0]==row['category']) and k[1:3]==(row['stage'],row['basis']))
                if not math.isclose(htotal,row['budget_single_hadron_balance'],rel_tol=1e-9,abs_tol=1e-7):
                    raise RuntimeError('Position hadron accounting does not close')
