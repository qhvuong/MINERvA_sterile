import itertools
import math
import unittest
from tools.FSIProbabilityBalance import inspect_hadrons,combine,BalanceAccumulator,cpp_class,SMALL_BOOST


def hadron(pdg=2112,fate=1,ke=200.,raw=2.):
    return dict(pdg=pdg,fate=fate,ke_MeV=ke,raw1=raw)


class BalanceTests(unittest.TestCase):
    def test_one_hadron_compensation(self):
        total=0.
        for p,fate,weight in [(0.25,1,2.),(.25,3,0.),(.5,4,-1.)]:
            _,v,b,_=inspect_hadrons([hadron(fate=fate,raw=weight)],3,0,0,12)
            total+=p*v['standard']
            self.assertAlmostEqual(sum(b.values()),v['standard']-1.)
        self.assertAlmostEqual(total,1.)

    def test_normalized_single_weights_do_not_guarantee_standard_event_normalization(self):
        standard=pure=0.
        outcomes=[(.25,1,2.),(.25,3,0.),(.5,4,-1.)]
        for a,b in itertools.product(outcomes,repeat=2):
            _,v,_,_=inspect_hadrons([hadron(fate=a[1],raw=a[2]),hadron(fate=b[1],raw=b[2])],3,0,0,12)
            standard+=a[0]*b[0]*v['standard'];pure+=a[0]*b[0]*v['pure_all']
        self.assertAlmostEqual(pure,1.)
        self.assertAlmostEqual(standard,.930625)

    def test_proton_threshold_does_not_apply_to_neutrons(self):
        records,v,_,n=inspect_hadrons([hadron(pdg=2212,ke=100.),hadron(ke=40.,raw=9.)],8,0,0,12)
        self.assertFalse(records[0]['eligible'])
        self.assertTrue(records[1]['eligible'])
        self.assertEqual(v['standard'],9.)
        self.assertEqual(n,1)

    def test_sentinel_and_exact_zero_are_distinct(self):
        _,v,_,n=inspect_hadrons([hadron(ke=20.,raw=-1.)],3,0,0,12)
        self.assertEqual(v['standard'],1.)
        self.assertEqual(n,0)
        _,v,_,n=inspect_hadrons([hadron(fate=3,raw=0.)],3,0,0,12)
        self.assertEqual(v['standard'],0.)
        self.assertEqual(n,1)

    def test_cap_and_four_or_more_bypass_budget(self):
        _,v,b,n=inspect_hadrons([hadron(pdg=211) for _ in range(6)],3,0,0,12)
        self.assertEqual(v['pure_all'],64.)
        self.assertEqual(v['pure_retained'],16.)
        self.assertEqual(v['standard'],1.)
        self.assertEqual(n,4)
        self.assertEqual(b['four_slot_cap'],-48.)
        self.assertEqual(sum(b.values()),0.)

    def test_damping_reduces_example_and_pions_are_not_damped(self):
        _,v,_,_=inspect_hadrons([hadron(pdg=211,raw=3.),hadron(pdg=111,raw=2.),hadron(raw=4.)],3,0,0,12)
        self.assertAlmostEqual(v['pure_all'],24.)
        self.assertAlmostEqual(v['standard'],18.6)
        self.assertEqual(v['without_0p70'],24.)

    def test_float_small_boost_and_resonance_adjustment(self):
        self.assertAlmostEqual(combine([1.,1.],[1.],3),SMALL_BOOST**3)
        _,v,b,_=inspect_hadrons([hadron(raw=2.)],2,0,0,12)
        self.assertAlmostEqual(v['standard'],1.9)
        self.assertAlmostEqual(b['resonance_adjustment'],-.1)

    def test_charm_is_not_reweighted(self):
        _,v,b,n=inspect_hadrons([hadron(raw=20.)],3,0,1,12)
        self.assertEqual(v['standard'],1.)
        self.assertEqual(v['pure_all'],1.)
        self.assertEqual(n,0)

    def test_event_accounting_keeps_hadron_covariance(self):
        records,v,b,n=inspect_hadrons([hadron(),hadron()],3,0,0,12)
        acc=BalanceAccumulator();acc.add('DIS',['before'],2.,v,b,records,n,12,6.)
        row=next(r for r in acc.hadron_rows() if r['category']=='DIS' and r['basis']=='unit')
        self.assertEqual(row['balance'],2.)
        self.assertEqual(row['sum_balance2_by_event'],4.)
        event=next(r for r in acc.event_rows() if r['category']=='DIS' and r['basis']=='unit' and r['split']=='all')
        self.assertAlmostEqual(event['ratio_standard'],2.89)

    def test_cpp_class_extraction_ignores_comment_braces(self):
        header='namespace X { class weight_fsi { /* } */ private: int a; }; }'
        self.assertEqual(cpp_class(header),'class weight_fsi { /* } */ private: int a; };')


if __name__=='__main__':unittest.main()
