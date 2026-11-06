"""MG042 mandatory deterministic qualification; development data only."""
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import mg10_reservoir as m
from lattice import P17, euler_mean_map, tau_leap_euler_1d
from runner import seed_for
from evaluate_mg10 import classify, complete


class Qualification(unittest.TestCase):
    def test_integer_ledger(self):
        self.assertTrue(np.array_equal(m.STOICH.sum(axis=0),np.zeros(7,dtype=int)))
        events = np.array([11,3,2,5,4,2,1],dtype=np.int64)
        initial = np.array([30]*7,dtype=np.int64)
        final = initial.copy()
        for k in range(7):
            for _ in range(events[k]):
                final += m.STOICH[:,k]
        np.testing.assert_array_equal(final-initial,m.STOICH@events)
        self.assertEqual(int(final.sum()),int(initial.sum()))
        # Explicit B consumption and Z creation; R4 conserves total count.
        np.testing.assert_array_equal(m.STOICH[:,3],[0,-1,1,0,-1,1,0])
        np.testing.assert_array_equal(m.STOICH[:,4],[0,1,-1,0,0,0,0])

    def test_ideal_mean_map(self):
        for nonuniform in (False,True):
            for steps in m.STEPS.values():
                z=m.initial("H",12)
                if nonuniform:
                    z[:3] += .01*np.sin(np.arange(12))[None,:]
                baseline=np.array(euler_mean_map(P17,.2,1/steps,1.,*z[:3].copy()))
                m.advance_euler(z,1,False,steps,steps)
                np.testing.assert_allclose(z[:3],baseline,rtol=0,atol=2e-11)
                np.testing.assert_array_equal(z[3:5],0)

    def test_pointwise_rhs_and_jacobian(self):
        rng=np.random.default_rng(104299)
        for finite in (False,True):
            z=m.initial("H",9)
            z[:3] *= rng.uniform(.95,1.05,(3,9))
            if finite: z[3:5]=rng.uniform(.01,.2,(2,9))
            rhs=m.ode_rhs(0,z.ravel(),10,finite).reshape(z.shape)
            if not finite:
                next_ideal=np.array(euler_mean_map(P17,.2,1e-5,1e-5,*z[:3].copy()))
                np.testing.assert_allclose(rhs[:3],(next_ideal-z[:3])/1e-5,rtol=0,atol=4e-10)
            direction=rng.normal(size=z.size)
            eps=1e-5
            numeric=(m.ode_rhs(0,z.ravel()+eps*direction,10,finite)-m.ode_rhs(0,z.ravel()-eps*direction,10,finite))/(2*eps)
            analytic=m.ode_jac(0,z.ravel(),10,finite)@direction
            np.testing.assert_allclose(analytic,numeric,rtol=2e-7,atol=1e-6)

    def test_zero_stocks_and_mirror(self):
        z=m.initial("H",8); z[3:5]=1
        rhs=m.ode_rhs(0,z.ravel(),1,True).reshape(z.shape)
        np.testing.assert_array_equal(rhs[[5,8]],0)
        z=m.initial("H",8); z[:3]+=.01*np.sin(np.arange(8))
        mirror=z[:,::-1].copy()
        m.advance_euler(z,1,True,1500,150)
        m.advance_euler(mirror,1,True,1500,150)
        np.testing.assert_allclose(z[:,::-1],mirror,rtol=0,atol=1e-12)

    def test_accounting_and_corrected_oracles(self):
        for cap in (1,1000000):
            for mode in ("FINITE","CHEMOSTAT"):
                for method in ("E0","E1","E2","O"):
                    job=dict(capacity=cap,mode=mode,method=method,start="H",n=8,T=2)
                    times,z=m.trajectory(job)
                    self.assertTrue(m.valid_state(z[-1]))
                    audit=m.accounting(times,z,job)
                    self.assertTrue(audit["pass"],str((job,audit)))
                    self.assertEqual(audit["external_A_total"]==0,mode=="FINITE")
                    self.assertEqual(audit["external_B_total"]==0,mode=="FINITE")

    def test_unchanged_ideal_stochastic(self):
        import mg8
        from mg8 import J
        z=m.initial("H",16)
        a=[np.rint(1000*c).astype(np.int64) for c in z[:3]]
        b=[v.copy() for v in a]
        args=(14.,29.,.1,1.,1.,1.,12.)
        seed=seed_for(10,99,0,0,0)
        first=tau_leap_euler_1d(*a,1000.,*args,.2,1/1500,40,0.,0.,np.arange(16)*.2,20,seed)
        second=mg8.run_stoch(*b,J,1000.,*args,.2,1/1500,40,20,seed)
        for x,y in zip(a,b): np.testing.assert_array_equal(x,y)
        for k in range(3): np.testing.assert_array_equal(first[k],second[k])

    def test_labels_and_invalid_records(self):
        from run_mg10 import budget_state
        self.assertEqual(budget_state(10,[3,4],20,1,10),"RUNNING")
        self.assertEqual(budget_state(10,[3,7],20,1,10),"BUDGET-SHORT")
        self.assertEqual(budget_state(10,[3],20,10,10),"TIMEOUT")
        self.assertEqual(classify(True,True,True,True,True,.005,.006),"RECOVERS BASELINE")
        self.assertEqual(classify(True,True,True,True,True,.02,.03),"BASELINE NOT RECOVERED")
        self.assertEqual(classify(True,True,True,False,False,.001,.001),"BASELINE NOT RECOVERED")
        self.assertEqual(classify(True,True,True,True,True,.009,.011),"UNRESOLVED AT MARGIN")
        self.assertEqual(classify(True,False,True,True,True,.001,.001),"NUMERICALLY UNRESOLVED")
        self.assertEqual(classify(False,True,True,True,True,.001,.001),"NUMERICALLY UNRESOLVED")
        times=np.arange(301);z=np.zeros((301,15,500))
        self.assertTrue(complete(times,z))
        self.assertFalse(complete(times[:-1],z[:-1]))
        self.assertFalse(complete(times[1:],z[1:]))
        for status in ("BUDGET-SHORT","TIMEOUT","INTERRUPTED","FAILED"):
            self.assertFalse(complete(times,z,status))
        z[70,0,0]=np.nan
        self.assertFalse(complete(times,z))
        z[70,0,0]=-1e-9
        self.assertFalse(complete(times,z))


if __name__=="__main__": unittest.main(verbosity=2)
