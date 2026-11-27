"""Common-weight selection accounting; no ROOT dependency."""
import ast
import bisect
import math
from collections import defaultdict


NEUTRON_EDGES = (0., 25., 30., 35., 40., 50., 75., 100., 125., 200., 500.)
TRUE_ENERGY_EDGES = (0., 2., 3., 4., 6., 8., 10., 15., 20.)
DETECTOR_TERMS = {
    'V4P3P0_ADD_GEANT_HADRON': 'GetGeantHadronWeight',
    'V4P3P0_ADD_MINOS_EFFICIENCY': 'GetMyMinosEfficiencyWeight',
}


def generator_flux_function(source, namespace):
    """Use the production V430 expression with only its two detector blocks removed.

    Reject source drift rather than guessing which operations are detector terms.
    The generated function has private globals; the production module is unchanged.
    """
    tree = ast.parse(source)
    classes = [n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'CVPythonUniverse']
    if len(classes) != 1:
        raise RuntimeError('Cannot identify CVPythonUniverse')
    methods = [n for n in classes[0].body if isinstance(n, ast.FunctionDef)
               and n.name == 'GetStandardWeightV4p3p0']
    if len(methods) != 1:
        raise RuntimeError('Cannot identify the production V430 method')
    method = methods[0]
    removed = []
    for statement in method.body:
        if isinstance(statement, ast.If) and isinstance(statement.test, ast.Name):
            flag = statement.test.id
            if flag in DETECTOR_TERMS:
                expected = ast.parse('if {}:\n    weight *= self.{}()'.format(flag, DETECTOR_TERMS[flag])).body[0]
                if ast.dump(statement, include_attributes=False) != ast.dump(expected, include_attributes=False):
                    raise RuntimeError('Detector block changed: {}'.format(flag))
                removed.append(flag)
    if removed != list(DETECTOR_TERMS):
        raise RuntimeError('Expected exactly the known GEANT and MINOS detector blocks')
    # They must remain the last operations before returning the weight.
    if not (isinstance(method.body[-1], ast.Return) and
            isinstance(method.body[-1].value, ast.Name) and method.body[-1].value.id == 'weight'):
        raise RuntimeError('V430 return changed')
    tail = method.body[-3:-1]
    if any(not isinstance(n, ast.If) or not isinstance(n.test, ast.Name) or
           n.test.id != flag for n, flag in zip(tail, DETECTOR_TERMS)):
        raise RuntimeError('Detector terms are no longer the final two operations')
    method.body = method.body[:-3] + method.body[-1:]
    method.name = '_diagnostic_generator_flux_v430'
    module = ast.fix_missing_locations(ast.Module(body=[method], type_ignores=[]))
    private = dict(namespace)
    exec(compile(module, '<diagnostic production V430 without detector terms>', 'exec'), private)
    return private[method.name], ast.dump(module, include_attributes=False)


def interval(value, edges):
    if value is None:
        return 'none'
    if not math.isfinite(value) or value < edges[0]:
        raise ValueError('Invalid diagnostic energy {}'.format(value))
    i = bisect.bisect_right(edges, value)
    return '{}_to_{}'.format(edges[i-1], edges[i] if i < len(edges) else 'inf')


def group(weight):
    if not math.isfinite(weight):
        raise ValueError('Nonfinite FSI weight')
    if weight < 0:
        return 'negative'
    if weight == 0:
        return 'zero'
    if abs(weight-1.) <= 1e-8:
        return 'near_one'
    if weight < 1:
        return 'suppressed'
    return 'one_to_two' if weight <= 2 else 'above_two'


class SelectionAccounting:
    def __init__(self):
        self.data = defaultdict(lambda: [0, 0., 0., 0., 0., 0.])

    def add(self, category, stages, baseline, elastic, npi, nnuc, min_neutron_ke, true_energy):
        corrected = baseline*elastic
        delta = corrected-baseline
        if not all(math.isfinite(x) for x in (baseline, corrected, elastic, delta, delta*delta, baseline*baseline)):
            raise ValueError('Nonfinite event weight')
        dimensions = (
            ('all', 'all'),
            ('weight_group', group(elastic)),
            ('contributors', 'npi={};nnuc={}'.format(npi, nnuc)),
            ('ntot', str(npi+nnuc)),
            ('min_fate1_neutron_KE_MeV', interval(min_neutron_ke, NEUTRON_EDGES)),
            ('contributors_and_neutron_KE', 'npi={};nnuc={};neutron={}'.format(
                npi, nnuc, interval(min_neutron_ke, NEUTRON_EDGES))),
            ('true_Enu_GeV', interval(true_energy, TRUE_ENERGY_EDGES)),
        )
        # Each split partitions events exactly once; do not sum different split types.
        for category_label in ('TOTAL', category):
            for stage in stages:
                for split, value in dimensions:
                    row = self.data[(category_label, stage, split, value)]
                    for i, number in enumerate((1, baseline, corrected, delta,
                                                baseline*baseline, delta*delta)):
                        row[i] += number

    def rows(self):
        result = []
        for (category, stage, split, value), v in sorted(self.data.items()):
            result.append(dict(category=category, stage=stage, split=split, value=value,
                               events=v[0], sum_w430=v[1], sum_w431=v[2], delta=v[3],
                               ratio=v[2]/v[1] if v[1] else None,
                               sum_w430_squared=v[4], sum_delta_squared=v[5]))
        return result

    def efficiencies(self, first_stage, last_stage):
        before = {(c, split, value): v for (c, stage, split, value), v in self.data.items()
                  if stage == first_stage}
        result = []
        for (c, split, value), b in sorted(before.items()):
            a = self.data.get((c, last_stage, split, value), [0]*6)
            eff430 = a[1]/b[1] if b[1] else None
            eff431 = a[2]/b[2] if b[2] else None
            result.append(dict(category=c, split=split, value=value,
                               before_events=b[0], after_events=a[0],
                               before_w430=b[1], before_w431=b[2],
                               after_w430=a[1], after_w431=a[2],
                               before_ratio=b[2]/b[1] if b[1] else None,
                               after_ratio=a[2]/a[1] if a[1] else None,
                               efficiency430=eff430, efficiency431=eff431,
                               efficiency_ratio=eff431/eff430 if eff430 and eff431 is not None else None))
        return result
