import csv
import json
import tempfile
import unittest
from pathlib import Path
from tools.FSIEventPair import WeightTrace,compare_runs


SOURCE='''
class CVPythonUniverse:
    def GetStandardWeightV4p3p0(self):
        weight = 1.0
        weight *= self.genie()
        if ENABLE:
            weight *= self.detector()
        return weight
'''


class Tests(unittest.TestCase):
    def test_trace_conditional_factor_and_runtime_product(self):
        class Event:
            def genie(self):return 2.
            def detector(self):return 3.
        t=WeightTrace(SOURCE,dict(ENABLE=True))
        total,components,applied=t.evaluate(Event())
        self.assertEqual(total,6.)
        self.assertEqual(list(components.values()),[2.,3.])
        self.assertEqual(len(applied),2)
        t=WeightTrace(SOURCE,dict(ENABLE=False))
        total,components,applied=t.evaluate(Event())
        self.assertEqual(total,2.)
        self.assertEqual(list(components.values()),[2.,1.])
        self.assertEqual(len(applied),1)

    def test_untraced_operations_fail_closure(self):
        class Event:
            def genie(self):return 2.
            def detector(self):return 3.
        t=WeightTrace(SOURCE.replace('return weight','return weight * 2'),dict(ENABLE=True))
        with self.assertRaises(RuntimeError):t.evaluate(Event())

    def make_run(self,path,tune,changed=False,extra=False):
        path.mkdir()
        audit=dict(status='complete',tune=tune,input_files=['file.root'],processed_entries=3,
                   playlist='me1A_p8',count=[0,1],AnaNuPDG=14,signal_cuts=['cut'],kinematics_cuts=[],
                   EN4_edges_GeV=[0.8,20.],method_sha256='same',component_labels=['factor'],
                   component_expressions={'factor':'expr'},shared_tune_flags={},selected_DIS_events=3 if extra else 2)
        (path/'audit.json').write_text(json.dumps(audit))
        rows=[]
        for entry,fsi in [(1,3.),(2,0.)]+([(3,1.)] if extra else []):
            factor=2.1 if changed and entry==1 else 2.
            rows.append(dict(chain_entry=entry,truth_class='CCNuMuDIS',region='regular',
                             mc_incoming=14,mc_intType=3,EN4_GeV=4.,true_Enu_GeV=4.,
                             active_total=factor*(fsi if tune=='v4p3p1' else 1.),
                             fsi_evaluated=fsi,factor=factor,factor_applied=1))
        with (path/'events.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

    def test_independent_runs_zero_fsi_and_exact_sum(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);self.make_run(p/'a','v4p3p0');self.make_run(p/'b','v4p3p1')
            summary=compare_runs(p/'a',p/'b',p/'out')
            self.assertEqual(summary['status'],'passed')
            r=next(r for r in summary['totals'] if r['category']=='CCNuMuDIS')
            self.assertEqual(r['sum_w430'],4.)
            self.assertEqual(r['sum_w431'],6.)
            self.assertEqual(r['ratio'],1.5)

    def test_changed_baseline_component_detected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);self.make_run(p/'a','v4p3p0');self.make_run(p/'b','v4p3p1',changed=True)
            s=compare_runs(p/'a',p/'b',p/'out')
            self.assertEqual(s['status'],'failed')
            with (p/'out'/'component_checks.csv').open() as stream:
                rows=list(csv.DictReader(stream))
            self.assertEqual(rows[0]['different_events'],'1')

    def test_unmatched_event_detected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);self.make_run(p/'a','v4p3p0');self.make_run(p/'b','v4p3p1',extra=True)
            s=compare_runs(p/'a',p/'b',p/'out')
            self.assertEqual(s['status'],'failed')
            self.assertEqual(s['missing_in_430'],[3])


if __name__=='__main__':unittest.main()
