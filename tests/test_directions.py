import unittest
from tools.FSIDirections import direction_geometry,ray_length,DirectionAccumulator

def record(x,group='inside'):
    return dict(zip(('start_x_fm','start_y_fm','start_z_fm'),x),tracking_boundary_fm=10.,position_group=group,start_r_fm=sum(v*v for v in x)**.5)

class DirectionTest(unittest.TestCase):
    def test_outward_and_inward(self):
        self.assertAlmostEqual(direction_geometry(record([3,0,0]),[0,0,0],[100,0,0])['path_to_boundary_fm'],7)
        self.assertAlmostEqual(direction_geometry(record([3,0,0]),[0,0,0],[-100,0,0])['path_to_boundary_fm'],13)
    def test_displacement_shortens_path(self):
        r=direction_geometry(record([3,0,0]),[1,0,0],[100,0,0])
        self.assertEqual(r['path_difference_fm'],2)
        self.assertEqual(r['displacement_parallel_fm'],2)
        self.assertEqual(r['displacement_transverse_fm'],0)
    def test_outside_inward_does_not_enter_transport(self):
        self.assertEqual(ray_length([12,0,0],[-1,0,0],10),20)
        self.assertEqual(direction_geometry(record([12,0,0],'outside'),[0,0,0],[-100,0,0])['path_to_boundary_fm'],0)
    def test_zero_momentum_is_explicit(self):
        r=direction_geometry(record([1,0,0]),[0,0,0],[0,0,0])
        self.assertFalse(r['direction_valid']);self.assertEqual(r['path_to_boundary_fm'],'')
    def test_accumulator_partitions_each_split(self):
        h=dict(pdg=2112,ke_MeV=350,fate=1,eligible=True,effective1=2,**record([3,0,0]))
        h.update(direction_geometry(h,[1,0,0],[100,0,0]))
        a=DirectionAccumulator();a.add('DIS',['before'],2,[h],12)
        rows=a.rows();self.assertEqual(len(rows),4)
        ref=dict(category='DIS',stage='before',targetA=12,pdg=2112,KE_MeV='300.0_to_500.0',hadrons=1,no_scatter=1,elastic=0,other_fates=0,eligible=1,gain=1.,loss=0.,balance=1.)
        a.validate([dict(ref,basis='unit'),dict(ref,basis='generator_flux',gain=2.,balance=2.)])
        for row in rows:self.assertEqual(row['balance'],1 if row['basis']=='unit' else 2)

if __name__=='__main__':unittest.main()
