import math
import unittest
from pathlib import Path
from tools.FSIPositions import geometry, indexed_fate_cpp, PositionAccumulator

class PositionsTest(unittest.TestCase):
    def test_boundary_and_original_not_daughter(self):
        b=1.4*12**(1/3)*3+.05
        result=geometry(1,[ -14,2112,2112],[0,14,1],[-1,2,-1],[-1,2,-1],
                        [[0,b,50],[0,0,0],[0,0,0]],-14,12,1.4,3,.05)
        self.assertEqual(result['position_group'],'outside')
        self.assertAlmostEqual(result['start_r_fm'],b)
        self.assertEqual(result['daughter_r_fm'],50)
        self.assertEqual(result['record_index'],1)
    def test_displacement_uses_probe(self):
        result=geometry(1,[-14,2112],[0,14],[-1,-1],[-1,-1],[[2,5],[0,4],[0,0]],-14,12,1.4,3,.05)
        self.assertAlmostEqual(result['probe_displacement_fm'],5)
        self.assertEqual(result['position_group'],'inside')
    def test_invalid_probe_rejected(self):
        with self.assertRaises(RuntimeError):
            geometry(0,[2112],[14],[-1],[-1],[[0],[0],[0]],-14,12,1.4,3,.05)
    def test_event_partition_not_duplicated_by_hadrons(self):
        a=PositionAccumulator()
        h=dict(eligible=True,position_group='outside',pdg=2112,ke_MeV=200,fate=1,effective1=2,probe_displacement_fm=4)
        a.add('DIS',['before'],2,4,[h,h],12)
        a.add('DIS',['before'],3,1,[],12)
        for basis,b,c,linear in [('unit',2,5,2),('generator_flux',5,11,4)]:
            ref=dict(category='DIS',stage='before',basis=basis,split='all',events=2,sum_baseline=b,sum_standard=c,budget_single_hadron_balance=linear)
            a.validate([ref])
            a.validate([dict(ref,category='TOTAL')])
            with self.assertRaises(RuntimeError):
                a.validate([dict(ref,category='TOTAL',sum_standard=c+1)])
        for r in a.event_rows():
            if r['event_group']=='with_outside':self.assertEqual(r['events'],1)
    def test_instrument_installed_source_shape(self):
        root=Path(__file__).resolve().parents[2]/'fsi_balance_build'
        if not root.exists():self.skipTest('Original source fixture not distributed')
        s=indexed_fate_cpp((root/'weighters/weight_fsi.h').read_text(),(root/'weight_fsi.cxx').read_text())
        self.assertIn('particle.__record_index = i;',s)
        self.assertIn('class PositionFateCode',s)
        self.assertIn('class PositionIntranukeParticle',s)
        self.assertNotIn('decreasingKE()',s)

if __name__=='__main__':unittest.main()
