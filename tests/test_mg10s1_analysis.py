"""MG046 controller and inference qualification with known synthetic distances."""
from pathlib import Path
import sys,json,unittest
import numpy as np
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import mg10s1 as m
import evaluate_mg10s1 as e
import bounded_pool as pool
from runner import seed_for


class Apparatus(unittest.TestCase):
    def test_CDF_and_labels(self):
        self.assertAlmostEqual(e.empirical_distance(np.array([0,0,1]),np.array([0,1,1])),1/3)
        self.assertEqual(e.empirical_distance(np.array([2,2]),np.array([3,3])),1)
        self.assertEqual(e.empirical_distance(np.array([0,1,3]),np.array([0,1,3])),0)
        self.assertEqual(e.classify(True,False,[(0,.05)]*96),"QUALIFIED ON TESTED SMALL SYSTEMS")
        self.assertEqual(e.classify(True,False,[(.051,.08)]*96),"NOT QUALIFIED ON TESTED CASES")
        self.assertEqual(e.classify(True,False,[(.04,.06)]*96),"UNRESOLVED")
        self.assertEqual(e.classify(False,False,[(0,.01)]*96),"UNRESOLVED")
        self.assertEqual(e.classify(True,False,[(0,.01)]*95),"UNRESOLVED")
        self.assertEqual(e.classify(True,True,[(0,.01)]*96),"QUALIFICATION FAILED")

    def test_synthetic_band_coverage_and_power(self):
        rng=np.random.default_rng(seed_for(10,94,0,4601,0))
        rows=[];saved={};trials=256
        for label,gap in (("null",0.),("inside",.01),("outside",.10)):
            a=rng.binomial(m.R,.5,(trials,96));b=rng.binomial(m.R,.5+gap,(trials,96))
            distance=abs(a/m.R-b/m.R)
            low=np.maximum(0,distance-e.half_error());high=np.minimum(1,distance+e.half_error())
            coverage=float(np.mean(np.all((low<=gap)&(high>=gap),axis=1)))
            qualified=float(np.mean(np.all(high<=.05,axis=1)))
            self.assertGreaterEqual(coverage,.98)
            if gap==0:self.assertGreaterEqual(qualified,.90)
            if gap==.10:self.assertEqual(qualified,0)
            rows.append(dict(label=label,true_distance=gap,coverage=coverage,qualification_frequency=qualified,trials=trials))
            saved[label+"_a"]=a;saved[label+"_b"]=b
        folder=m.BASE/"qualification"
        np.savez_compressed(folder/"synthetic_families.npz",**saved)
        (folder/"synthetic_summary.json").write_text(json.dumps(rows,indent=2)+"\n",encoding="utf-8")

    def test_owned_process_stops_and_resume(self):
        folder=m.BASE/"qualification/controller"
        worker=m.ROOT/"src/run_mg10s1.py"
        job=dict(id="development_stop",case=0,method=3,stage=4601,reps=512,chunk=512)
        for name,ceiling,timeout in (("timeout",30,.05),("budget",.05,30)):
            result=pool.run([job],folder/name,worker,m.ROOT,ceiling,timeout=timeout)
            self.assertEqual(len(result["records"]),1)
            self.assertEqual(result["records"][0]["status"],"TIMEOUT" if name=="timeout" else "BUDGET-SHORT")
            resumed=pool.run([job],folder/name,worker,m.ROOT,ceiling,timeout=timeout)
            self.assertEqual(resumed["charged_worker_seconds"],result["charged_worker_seconds"])
        with self.assertRaises(KeyboardInterrupt):
            with patch.object(pool.time,"sleep",side_effect=KeyboardInterrupt):pool.run([job],folder/"interruption",worker,m.ROOT,30)
        result=json.loads((folder/"interruption/ledger.json").read_text())
        self.assertEqual(result["records"][0]["status"],"INTERRUPTED")


if __name__=="__main__":unittest.main(verbosity=2)
