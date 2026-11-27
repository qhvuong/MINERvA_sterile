"""Trace the production V430 factors and compare independent tune runs."""
import ast
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path


def close(a, b):
    return math.isclose(float(a), float(b), rel_tol=1e-10, abs_tol=1e-12)


class WeightTrace:
    def __init__(self, source, namespace):
        tree = ast.parse(source)
        cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'CVPythonUniverse')
        method = next(n for n in cls.body if isinstance(n, ast.FunctionDef)
                      and n.name == 'GetStandardWeightV4p3p0')
        self.method_sha256 = hashlib.sha256(ast.dump(method, include_attributes=False).encode()).hexdigest()
        self.labels = []
        self.expressions = {}
        self.values = {}
        owner = self

        class Instrument(ast.NodeTransformer):
            def visit_AugAssign(self, node):
                if isinstance(node.target, ast.Name) and node.target.id == 'weight':
                    if not isinstance(node.op, ast.Mult):
                        raise RuntimeError('Unexpected nonmultiplicative production weight operation')
                    value = node.value
                    if isinstance(value, ast.Call):
                        name = value.func.attr if isinstance(value.func, ast.Attribute) else getattr(value.func, 'id', 'call')
                    elif isinstance(value, ast.Name):
                        name = value.id
                    elif isinstance(value, ast.Constant):
                        name = 'constant_{}'.format(value.value)
                    else:
                        name = 'expression'
                    label = '{:02d}_{}'.format(len(owner.labels), name)
                    owner.labels.append(label)
                    owner.expressions[label] = ast.dump(value, include_attributes=False)
                    node.value = ast.Call(func=ast.Name(id='_pair_capture', ctx=ast.Load()),
                                          args=[ast.Constant(label), value], keywords=[])
                return node

        method = Instrument().visit(method)
        if not self.labels:
            raise RuntimeError('No production factors found')
        method.name = '_pair_traced_v430'
        module = ast.fix_missing_locations(ast.Module(body=[method], type_ignores=[]))
        private = dict(namespace, _pair_capture=self.capture)
        exec(compile(module, '<traced production V430>', 'exec'), private)
        self.function = private[method.name]

    def capture(self, label, value):
        value = float(value)
        if not math.isfinite(value):
            raise RuntimeError('Nonfinite component {}'.format(label))
        if label in self.values:
            raise RuntimeError('Factor executed more than once: {}'.format(label))
        self.values[label] = value
        return value

    def evaluate(self, universe):
        self.values = {}
        total = float(self.function(universe))
        components = {label: self.values.get(label, 1.) for label in self.labels}
        product = 1.
        # Runtime order, including conditional factors, matches the production expression.
        for value in self.values.values():
            product *= value
        if not math.isfinite(total) or not close(total, product):
            raise RuntimeError('Traced factors do not reconstruct V430 total')
        return total, components, set(self.values)


def load_run(directory):
    directory = Path(directory)
    audit = json.loads((directory/'audit.json').read_text())
    if audit['status'] != 'complete':
        raise ValueError('Run is not complete: {}'.format(directory))
    events = {}
    with (directory/'events.csv').open() as stream:
        for row in csv.DictReader(stream):
            key = int(row['chain_entry'])
            if key in events:
                raise ValueError('Duplicate event key {}'.format(key))
            events[key] = row
    if len(events) != audit['selected_DIS_events']:
        raise ValueError('Event count does not match audit')
    return audit, events


def compare_runs(directory430, directory431, output):
    a, events430 = load_run(directory430)
    b, events431 = load_run(directory431)
    if (a['tune'], b['tune']) != ('v4p3p0', 'v4p3p1'):
        raise ValueError('Supply V430 first and V431 second')
    # These define the event keys, sample, and shared baseline calculation.
    required = ('input_files','processed_entries','playlist','count','AnaNuPDG',
                'signal_cuts','kinematics_cuts','EN4_edges_GeV','method_sha256',
                'component_labels','component_expressions','shared_tune_flags')
    differences = [name for name in required if a.get(name) != b.get(name)]
    if differences:
        raise ValueError('Incompatible runs: {}'.format(', '.join(differences)))
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    common = sorted(events430.keys() & events431.keys())
    missing430 = sorted(events431.keys()-events430.keys())
    missing431 = sorted(events430.keys()-events431.keys())
    labels = a['component_labels']
    fields = ['chain_entry','truth_class','EN4_GeV','region','w430','w431','delta',
              'actual_ratio','fsi_evaluated_430','fsi_evaluated_431','expected_w431',
              'fsi_only_closure','changed_components','changed_applied_flags']
    for label in labels:
        fields += [label+'_430',label+'_431']
    changes = {label: dict(component=label, different_events=0, changed_applied_flags=0,
                           max_abs_difference=0.) for label in labels}
    sums = defaultdict(lambda: [0,0.,0.,0.])
    invalid = []
    influential = []
    with (out/'matched_events.csv').open('w',newline='') as f:
        writer = csv.DictWriter(f,fieldnames=fields)
        writer.writeheader()
        for key in common:
            old, new = events430[key], events431[key]
            context_changed = [name for name in ('truth_class','region','mc_incoming','mc_intType')
                               if old[name] != new[name]]
            for name in ('EN4_GeV','true_Enu_GeV'):
                if not close(old[name],new[name]):
                    context_changed.append(name)
            if context_changed:
                invalid.append(dict(chain_entry=key,reason='event context changed',fields=context_changed))
            w430, w431 = float(old['active_total']), float(new['active_total'])
            fsi430, fsi431 = float(old['fsi_evaluated']), float(new['fsi_evaluated'])
            expected = w430*fsi431
            closure = close(w431,expected) and close(fsi430,fsi431)
            changed, changed_applied = [], []
            row = dict(chain_entry=key, truth_class=old['truth_class'], EN4_GeV=old['EN4_GeV'],
                       region=old['region'],w430=w430,w431=w431,delta=w431-w430,
                       actual_ratio=w431/w430 if w430 else None,
                       fsi_evaluated_430=fsi430,fsi_evaluated_431=fsi431,
                       expected_w431=expected,fsi_only_closure=closure)
            for label in labels:
                x,y = float(old[label]),float(new[label])
                row[label+'_430'],row[label+'_431'] = x,y
                changes[label]['max_abs_difference'] = max(changes[label]['max_abs_difference'],abs(x-y))
                if not close(x,y):
                    changes[label]['different_events'] += 1
                    changed.append(label)
                if old[label+'_applied'] != new[label+'_applied']:
                    changes[label]['changed_applied_flags'] += 1
                    changed_applied.append(label)
            row.update(changed_components=';'.join(changed),changed_applied_flags=';'.join(changed_applied))
            if not closure:
                invalid.append(dict(chain_entry=key,reason='FSI-only closure failed'))
            writer.writerow(row)
            influential.append(dict(chain_entry=key,truth_class=old['truth_class'],
                                    region=old['region'],w430=w430,w431=w431,
                                    fsi=fsi431,delta=w431-w430))
            for category in ('TOTAL_selected_muon_DIS',old['truth_class']):
                s=sums[(category,old['region'])]
                for i,x in enumerate((1,w430,w431,expected)):
                    s[i] += x
    totals = [dict(category=c,region=r,events=v[0],sum_w430=v[1],sum_w431=v[2],
                   delta=v[2]-v[1],ratio=v[2]/v[1] if v[1] else None,
                   expected_sum_w431=v[3]) for (c,r),v in sorted(sums.items())]
    for filename, rows in [('component_checks.csv',list(changes.values())),('integrated_yields.csv',totals)]:
        with (out/filename).open('w',newline='') as f:
            if rows:
                writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    passed = bool(common) and not (missing430 or missing431 or invalid or
                                   any(r['different_events'] or r['changed_applied_flags'] for r in changes.values()))
    summary = dict(status='passed' if passed else 'failed',matched_events=len(common),
                   missing_in_430=missing430,missing_in_431=missing431,failures=invalid,
                   totals=totals,top_regular_events=sorted(
                       [e for e in influential if e['region']=='regular'],key=lambda e:abs(e['delta']),reverse=True)[:20])
    (out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    return summary
