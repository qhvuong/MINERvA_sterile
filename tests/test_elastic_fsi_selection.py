import unittest
from tools.ElasticFSISelection import (SelectionAccounting, generator_flux_function,
                                      interval, NEUTRON_EDGES)


SOURCE = '''
class CVPythonUniverse:
    def GetStandardWeightV4p3p0(self):
        weight = self.physics()
        if V4P3P0_ADD_GEANT_HADRON:
            weight *= self.GetGeantHadronWeight()
        if V4P3P0_ADD_MINOS_EFFICIENCY:
            weight *= self.GetMyMinosEfficiencyWeight()
        return weight
'''


class SelectionTests(unittest.TestCase):
    def test_detector_removal_without_calling_detector_functions(self):
        class Event:
            def physics(self):
                return 3.25
            def GetGeantHadronWeight(self):
                raise AssertionError('Detector factor evaluated before selection')
            GetMyMinosEfficiencyWeight = GetGeantHadronWeight
        namespace = dict(V4P3P0_ADD_GEANT_HADRON=True, V4P3P0_ADD_MINOS_EFFICIENCY=True)
        fn, _ = generator_flux_function(SOURCE, namespace)
        self.assertEqual(fn(Event()), 3.25)
        self.assertNotIn('_diagnostic_generator_flux_v430', namespace)

    def test_source_drift_rejected(self):
        for changed in (
            SOURCE.replace('weight *= self.GetMyMinosEfficiencyWeight()', 'weight *= self.other()'),
            SOURCE.replace('        return weight', '        weight *= 2\n        return weight'),
            SOURCE.replace('V4P3P0_ADD_GEANT_HADRON', 'RENAMED_FLAG'),
        ):
            with self.assertRaises(RuntimeError):
                generator_flux_function(changed, {})

    def test_before_after_efficiency_identity_and_partitions(self):
        a = SelectionAccounting()
        a.add('DIS', ['before','after'], 2., 3., 2, 1, 33.66, 6.)
        a.add('DIS', ['before'], 3., 0., 1, 1, 70., 6.)
        rows = a.rows()
        for split in {r['split'] for r in rows}:
            before = [r for r in rows if r['category']=='DIS' and r['stage']=='before' and r['split']==split]
            self.assertEqual(sum(r['events'] for r in before), 2)
            self.assertEqual(sum(r['sum_w430'] for r in before), 5.)
            self.assertEqual(sum(r['sum_w431'] for r in before), 6.)
        efficiency = next(r for r in a.efficiencies('before','after')
                          if r['category']=='DIS' and r['split']=='all')
        self.assertEqual(efficiency['efficiency430'], .4)
        self.assertEqual(efficiency['efficiency431'], 1.)
        self.assertAlmostEqual(efficiency['efficiency_ratio'],
                               efficiency['after_ratio']/efficiency['before_ratio'])

    def test_zero_weight_and_empty_selected_strata(self):
        a = SelectionAccounting()
        a.add('2p2h', ['before'], 2., 0., 0, 1, None, 4.)
        r = next(r for r in a.efficiencies('before','after') if r['split']=='all')
        self.assertEqual(r['after_events'], 0)
        self.assertIsNone(r['efficiency431'])
        self.assertIsNone(r['efficiency_ratio'])

    def test_neutron_boundaries(self):
        self.assertEqual(interval(None,NEUTRON_EDGES), 'none')
        self.assertEqual(interval(29.999,NEUTRON_EDGES), '25.0_to_30.0')
        self.assertEqual(interval(30.,NEUTRON_EDGES), '30.0_to_35.0')
        self.assertEqual(interval(500.,NEUTRON_EDGES), '500.0_to_inf')


if __name__ == '__main__':
    unittest.main()
