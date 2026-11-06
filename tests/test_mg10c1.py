"""MG044 development qualification, no scientific trajectories."""
from pathlib import Path
import sys,unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import mg10c1 as m
from runner import reference
from measure import smooth_periodic,detect


class Followup(unittest.TestCase):
    def test_controller_stops_and_resume(self):
        import run_mg10c1 as driver
        from unittest.mock import patch
        import json
        base=m.BASE/"qualification/controller"
        job=dict(id="development_stop",kind="full",method="E2",epsilon=0.,mode=7,T=5)
        for name,cap,timeout in (("budget",.05,30),("timeout",30,.05)):
            result=driver.pool([job],base/name,cap,timeout=timeout)
            self.assertEqual(len(result["records"]),1)
            self.assertEqual(result["records"][0]["status"],"BUDGET-SHORT" if name=="budget" else "TIMEOUT")
            resumed=driver.pool([job],base/name,cap,timeout=timeout)
            self.assertEqual(resumed["charged_worker_seconds"],result["charged_worker_seconds"])
            self.assertEqual(len(resumed["records"]),1)
        with self.assertRaises(KeyboardInterrupt):
            with patch.object(driver.time,"sleep",side_effect=KeyboardInterrupt):
                driver.pool([job],base/"interrupt",30)
        result=json.loads((base/"interrupt/ledger.json").read_text())
        self.assertEqual(result["records"][0]["status"],"INTERRUPTED")

    def test_decision_boundaries(self):
        from evaluate_mg10c1 import classify_reference,classify_amplification
        self.assertEqual(classify_reference(True,1e-9,1e-5,0,True),"HOMOGENEOUS REFERENCE RESOLVED")
        for complete,agreement,error,span,ledger in [(False,0,0,0,True),(True,1e-6,0,0,True),(True,0,1e-3,0,True),(True,0,0,1e-10,True),(True,0,0,0,False)]:
            self.assertEqual(classify_reference(complete,agreement,error,span,ledger),"UNRESOLVED")
        good=[dict(eligible_count=3,relative_error_max=.05)]*4
        self.assertEqual(classify_amplification(True,100,True,good),"CONTROLLED AMPLIFICATION SUPPORTED")
        self.assertEqual(classify_amplification(True,99,True,good),"NO >=100x AMPLIFICATION IN THE TESTED FAMILY/HORIZON")
        self.assertEqual(classify_amplification(True,100,False,good),"UNRESOLVED")
        self.assertEqual(classify_amplification(True,100,True,[dict(eligible_count=2,relative_error_max=0)]*4),"UNRESOLVED")

    def test_detector_fixtures(self):
        ref=reference(29.,.2)
        fixtures=m.fixtures(ref)
        folder=m.BASE/"qualification";folder.mkdir(parents=True,exist_ok=True)
        states={};rows=[]
        for name,field,expected in fixtures:
            actual=len(m.resolved_minima(field,ref))
            legacy=len(detect(smooth_periodic(field-m.SCALE[1],.2),ref["x"],.2,ref["thr"],ref["lam"]))
            states[name]=field
            rows.append(dict(name=name,expected=expected,actual=actual,legacy=legacy))
        import json
        np.savez_compressed(folder/"detector_fixtures.npz",**states)
        (folder/"detector_fixtures.json").write_text(json.dumps(rows,indent=2)+"\n",encoding="utf-8")
        self.assertTrue(all(r["actual"]==r["expected"] for r in rows),str(rows))

    def test_reduced_independent_rhs(self):
        for background in ([15,16,2,.1,.2],[4,7,1,.6,.4],[140/9,140/9,29/(140/9),0,0]):
            np.testing.assert_allclose(m.background_rhs(0,np.array(background)),m.background_ref_rhs(0,np.array(background)),rtol=0,atol=3e-13)
        primary=m.reduced(dict(kind="reduced_primary",T=2))
        reference=m.reduced(dict(kind="reduced_reference",T=2))
        np.testing.assert_allclose(primary["state"],reference["state"],rtol=0,atol=2e-9)
        self.assertLess(max(primary["interpolation_error"],reference["interpolation_error"]),1e-10)

    def test_mode_matrix_against_full_drift(self):
        rng=np.random.default_rng(104499)
        z=m.old.initial("H");z[3:5]=np.array([.2,.1])[:,None]
        background=z[:5,0]
        direction=rng.normal(size=4)
        scales=m.SCALE[[0,1,2,4]]
        eps=1e-6
        for mode in (1,16,64,250):
            cosine=np.cos(2*np.pi*mode*np.arange(500)/500)
            variation=np.zeros_like(z)
            variation[[0,1,2,4]]=direction[:,None]*scales[:,None]*cosine
            plus=m.old.ode_rhs(0,(z+eps*variation).ravel(),1,True).reshape(15,500)
            minus=m.old.ode_rhs(0,(z-eps*variation).ravel(),1,True).reshape(15,500)
            numeric=((plus-minus)/(2*eps))[[0,1,2,4]]/scales[:,None]
            factor=1 if mode==250 else 2
            coefficient=factor*np.mean(numeric*cosine,axis=1)
            predicted=m.mode_matrix(background,[mode])[0]@direction
            np.testing.assert_allclose(coefficient,predicted,rtol=2e-7,atol=2e-7)

    def test_full_adapter_and_accounting(self):
        for method in ("E2","Oplus"):
            for epsilon in (0.,1e-6):
                job=dict(method=method,epsilon=epsilon,mode=7,T=1)
                raw=m.full(job)
                self.assertTrue(np.array_equal(raw["times"],[0,1]))
                audit=m.old.accounting(raw["times"],raw["state"],dict(mode="FINITE",capacity=1,method="E2" if method=="E2" else "O"))
                self.assertTrue(audit["pass"],str(audit))
                self.assertTrue(m.old.valid_state(raw["state"][-1]))
                if method=="E2" and epsilon==0:
                    np.testing.assert_array_equal(np.ptp(raw["state"][:,:3],axis=2),0)


if __name__=="__main__":unittest.main(verbosity=2)
