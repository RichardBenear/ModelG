"""MG049 confidence-family and owned-worker adapter qualification."""
from pathlib import Path
import sys,json,unittest
import numpy as np
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import evaluate_mg10s2 as e
import bounded_pool as pool
from runner import seed_for

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10S2"


class Apparatus(unittest.TestCase):
    def test_signed_CDF_and_decisions(self):
        self.assertEqual(e.distance(np.array([-2.,-1.]),np.array([1.,2.])),1.)
        self.assertAlmostEqual(e.distance(np.array([-1.,0,0]),np.array([-1.,0,1])),1/3)
        good=[dict(lower=0.,upper=.08)]*192
        self.assertEqual(e.classify(good,192,"OK"),"OK")
        self.assertEqual(e.classify(good[:-1],192,"OK"),"UNRESOLVED")
        self.assertEqual(e.classify([dict(lower=.081,upper=.12)]*192,192,"OK"),"NOT QUALIFIED ON TESTED CASES")
        self.assertEqual(e.classify([dict(lower=.07,upper=.09)]*192,192,"OK"),"UNRESOLVED")
        self.assertEqual(e.classify(good,192,"OK",True),"QUALIFICATION FAILED")

    def test_synthetic_families(self):
        rng=np.random.default_rng(seed_for(10,94,0,4901,0))
        arrays={};rows=[]
        for m in (192,96):
            self.assertLessEqual(e.width(m),.04)
            for label,gap in (("null",0.),("inside",.01),("outside",.15)):
                a=rng.binomial(e.R,.5,(256,m));b=rng.binomial(e.R,.5+gap,(256,m))
                d=abs(a-b)/e.R
                lower=np.maximum(0,d-e.width(m));upper=np.minimum(1,d+e.width(m))
                coverage=float(np.mean(np.all((lower<=gap)&(upper>=gap),axis=1)))
                qualified=float(np.mean(np.all(upper<=.08,axis=1)))
                self.assertGreaterEqual(coverage,.98)
                if gap==0:self.assertGreaterEqual(qualified,.70)
                if gap==.15:self.assertEqual(qualified,0)
                key=f"M{m}_{label}"
                arrays[key+"_a"]=a;arrays[key+"_b"]=b
                rows.append(dict(M=m,label=label,true_distance=gap,coverage=coverage,qualification_frequency=qualified,trials=256))
        np.savez_compressed(BASE/"qualification/synthetic_families.npz",**arrays)
        (BASE/"qualification/synthetic_summary.json").write_text(json.dumps(rows,indent=2)+"\n",encoding="utf-8")

    def test_owned_stops_and_resume(self):
        folder=BASE/"qualification/controller"
        job=dict(id="development_stop",kind="snapshot",case=0,method=0,block=0,reps=512,chunk=512,offset=0,stage=4901,warmup=False)
        script=ROOT/"src/run_mg10s2.py"
        for name,cap,timeout in (("budget",.05,30),("timeout",30,.05)):
            result=pool.run([job],folder/name,script,ROOT,cap,timeout=timeout,poll_seconds=.25)
            self.assertEqual(result["records"][0]["status"],"BUDGET-SHORT" if name=="budget" else "TIMEOUT")
            second=pool.run([job],folder/name,script,ROOT,cap,timeout=timeout,poll_seconds=.25)
            self.assertEqual(result["charged_worker_seconds"],second["charged_worker_seconds"])
        with self.assertRaises(KeyboardInterrupt):
            with patch.object(pool.time,"sleep",side_effect=KeyboardInterrupt):pool.run([job],folder/"interruption",script,ROOT,30,poll_seconds=.25)
        data=json.loads((folder/"interruption/ledger.json").read_text())
        self.assertEqual(data["records"][0]["status"],"INTERRUPTED")


if __name__=="__main__":unittest.main(verbosity=2)
