"""Straight-line geometry only: no nuclear-density or interaction-probability model."""
import math
from tools.FSIProbabilityBalance import energy_bin

DIRECTION_FIELDS=['probe_x_fm','probe_y_fm','probe_z_fm','px_MeV','py_MeV','pz_MeV',
 'momentum_MeV','direction_valid','cos_start_radial','cos_probe_radial',
 'displacement_parallel_fm','displacement_transverse_fm','path_to_boundary_fm',
 'probe_path_to_boundary_fm','path_difference_fm']


def ray_length(x,u,radius):
    """Length of forward ray inside a sphere, including outside/inward starts."""
    dot=sum(a*b for a,b in zip(x,u));r2=sum(a*a for a in x)
    disc=dot*dot+radius*radius-r2
    if disc<=0:return 0.
    root=math.sqrt(disc)
    return max(0.,-dot+root-max(0.,-dot-root))


def direction_geometry(record,probe,momentum):
    x=[record['start_'+a+'_fm'] for a in ('x','y','z')]
    p=list(map(float,momentum));probe=list(map(float,probe))
    if not all(math.isfinite(v) for v in x+p+probe):raise RuntimeError('Nonfinite direction input')
    r=record['tracking_boundary_fm'];norm=math.sqrt(sum(v*v for v in p))
    out=dict(zip(['probe_x_fm','probe_y_fm','probe_z_fm'],probe))
    out.update(zip(['px_MeV','py_MeV','pz_MeV'],p));out.update(momentum_MeV=norm,direction_valid=norm>0)
    for key in DIRECTION_FIELDS[8:]:out[key]=''
    if norm==0:return out
    u=[v/norm for v in p]
    def cosine(v):
        mag=math.sqrt(sum(a*a for a in v))
        return max(-1.,min(1.,sum(a*b for a,b in zip(v,u))/mag)) if mag else ''
    delta=[a-b for a,b in zip(x,probe)];parallel=sum(a*b for a,b in zip(delta,u))
    transverse=math.sqrt(sum((a-parallel*b)**2 for a,b in zip(delta,u)))
    # Outside starts do not enter the shown transport loop, even if directed inward.
    path=ray_length(x,u,r) if record['position_group']=='inside' else 0.
    reference=ray_length(probe,u,r) if sum(v*v for v in probe)<r*r else 0.
    out.update(cos_start_radial=cosine(x),cos_probe_radial=cosine(probe),
      displacement_parallel_fm=parallel,displacement_transverse_fm=transverse,
      path_to_boundary_fm=path,probe_path_to_boundary_fm=reference,path_difference_fm=reference-path)
    return out


def interval(value,edges):
    if value=='':return 'undefined'
    import bisect
    i=bisect.bisect_right(edges,value)
    return '{}_to_{}'.format(edges[i-1] if i else '-inf',edges[i] if i<len(edges) else 'inf')


class DirectionAccumulator:
    def __init__(self):self.data={}
    def add(self,category,stages,baseline,records,target_a):
        if category not in ('QE','DIS') or target_a!=12:return
        for h in records:
            if h['pdg']!=2112:continue
            dimensions=[('path',interval(h['path_to_boundary_fm'],[0.,2.,4.,6.,8.,10.,12.,16.,20.])),
              ('radius_direction',interval(h['start_r_fm'],[0.,1.,2.,3.,4.,6.,10.,15.])+'|'+interval(h['cos_start_radial'],[-1.,-.5,0.,.5,1.000001]))]
            for stage in stages:
                for basis,w in [('unit',1.),('generator_flux',baseline)]:
                    for split,value in dimensions:
                        key=(category,stage,basis,energy_bin(h['ke_MeV']),split,value)
                        r=self.data.setdefault(key,dict(hadrons=0,no_scatter=0,elastic=0,other_fates=0,eligible=0,gain=0.,loss=0.,balance=0.))
                        r['hadrons']+=1;r['eligible']+=int(h['eligible'])
                        r['no_scatter' if h['fate']==1 else 'elastic' if h['fate']==3 else 'other_fates']+=1
                        delta=w*(h['effective1']-1.)
                        r['gain']+=max(0.,delta);r['loss']+=min(0.,delta);r['balance']+=delta
    def rows(self):
        return [dict(zip(('category','stage','basis','KE_MeV','split','value'),k),targetA=12,pdg=2112,**r) for k,r in sorted(self.data.items())]

    def validate(self,reference):
        expected={}
        fields=('hadrons','no_scatter','elastic','other_fates','eligible','gain','loss','balance')
        for row in reference:
            if row['category'] not in ('QE','DIS') or int(row['targetA'])!=12 or int(row['pdg'])!=2112:continue
            key=tuple(row[k] for k in ('category','stage','basis','KE_MeV'))
            r=expected.setdefault(key,{f:0. for f in fields})
            for f in fields:r[f]+=row[f]
        for split in ('path','radius_direction'):
            actual={}
            for key,row in self.data.items():
                if key[4]!=split:continue
                r=actual.setdefault(key[:4],{f:0. for f in fields})
                for f in fields:r[f]+=row[f]
            if actual.keys()!=expected.keys():raise RuntimeError('Direction partition keys mismatch')
            for key,r in expected.items():
                for f in fields:
                    if not math.isclose(actual[key][f],r[f],rel_tol=1e-9,abs_tol=1e-7):
                        raise RuntimeError('Direction partition does not close: {} {} {}'.format(split,key,f))
