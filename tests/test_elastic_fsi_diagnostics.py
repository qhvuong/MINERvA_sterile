import math
import unittest
from tools.ElasticFSIDiagnostics import Accounting, weight_group


class DiagnosticsTests(unittest.TestCase):
    def test_boundaries_and_nonfinite(self):
        cases = [(-0.1, 'negative'), (0., 'zero'), (0.5, 'suppressed'),
                 (1.-1e-10, 'near_one'), (1.+1e-10, 'near_one'),
                 (1.5, 'enhanced_to_2'), (2., 'enhanced_to_2'), (2.1, 'above_2')]
        for value, expected in cases:
            self.assertEqual(weight_group(value), expected)
        for value in (math.nan, math.inf, -math.inf):
            with self.assertRaises(ValueError):
                weight_group(value)

    def test_weighted_closure_and_root_flow_edges(self):
        a = Accounting([0.8, 2., 20.])
        for i, (energy, baseline, elastic) in enumerate([
            (0.7, 4., 3.), (0.8, 10., 0.), (2., 20., 1.5), (19., 5., 2.), (20., 100., 3.)
        ]):
            a.add('DIS', energy, baseline, elastic, {'chain_entry': i})
        totals = {(r['category'], r['scope']): r for r in a.totals()}
        regular = totals['TOTAL', 'regular']
        self.assertEqual(regular['sum_w430'], 35.)
        self.assertEqual(regular['sum_w431'], 40.)
        self.assertEqual(regular['delta'], 5.)
        self.assertEqual(regular['events'], 3)
        self.assertEqual(totals['DIS', 'regular'], dict(regular, category='DIS'))
        self.assertEqual(totals['TOTAL', 'overflow']['sum_w430'], 100.)
        self.assertEqual(totals['TOTAL', 'underflow']['sum_w430'], 4.)
        self.assertTrue(any(r['scope'] == 'bin_2' and r['events'] == 2 for r in a.rows()))

    def test_rank_by_contribution_separately_from_elastic_weight(self):
        a = Accounting([0., 1.], top_n=1)
        for i, baseline, elastic in [(0, 1., 10.), (1, 100., 1.5), (2, 1000., 0.)]:
            a.add('DIS', 0.5, baseline, elastic, {'chain_entry': i})
        top = {r['ranking']: r['events'][0] for r in a.top_events()}
        self.assertEqual(top['elastic']['chain_entry'], 0)
        self.assertEqual(top['absolute_delta']['chain_entry'], 2)
        self.assertEqual(len(a.zeros['DIS', 'regular']), 1)

    def test_zero_baseline_and_signed_weights(self):
        a = Accounting([0., 1.])
        a.add('QE', 0.5, 0., 1., {})
        self.assertIsNone(a.totals()[0]['ratio'])
        a.add('QE', 0.5, -2., 1.5, {})
        self.assertAlmostEqual(a.totals()[0]['delta'], -1.)
        with self.assertRaises(ValueError):
            a.add('QE', 0.5, float('nan'), 1., {})


if __name__ == '__main__':
    unittest.main()
