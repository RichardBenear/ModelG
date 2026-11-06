"""MG046 development-only finite stochastic qualification, no production keys."""
import sys,json,unittest
from pathlib import Path
import numpy as np
import numba as nb
from scipy.stats import binom
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import mg10s1 as m
from runner import seed_for


@nb.njit(cache=True)
def steps(initial,a0,b0,omega,diff,reps,seed):
    np.random.seed(seed)
    counts=np.zeros((reps,7),dtype=np.int64);delta=np.zeros((reps,7),dtype=np.int64)
    for r in range(reps):
        z=initial.copy();f=np.zeros_like(z);h=np.zeros((6,z.shape[1]),dtype=np.int64)
        dt,prob=m._leap(z,f,h,a0,b0,omega,diff,np.ones(7),.00025)
        if dt!=.00025 or prob>.25+1e-14:raise ValueError("Step check")
        counts[r]=f[:,0];delta[r]=z[:,0]-initial[:,0]
    return counts,delta


class SamplerQualification(unittest.TestCase):
    def test_integer_ledger_and_bounds(self):
        rows=[];folder=m.BASE/"qualification";folder.mkdir(parents=True,exist_ok=True)
        for case_index in range(4):
            case=m.fixture(case_index)
            for method in (0,1,3):
                data,stats=m.simulate(case,method,128,seed_for(10,case_index,method,4601,0))
                audit=m.ledger(data,case["state"])
                self.assertTrue(audit["pass_check"],str(audit))
                self.assertTrue(np.all(data[:,:,10]<=case["state"][4]))
                self.assertTrue(np.all(data[:,:,7]<=case["state"][3]))
                if method:self.assertLessEqual(stats[:,2].max(),.25+1e-14)
                if case_index==2:
                    np.testing.assert_array_equal(data[:,:,[3,4,7,10]],0)
                rows.append(dict(case=case_index,method=method,audit=audit,max_probability=float(stats[:,2].max())))
                np.savez_compressed(folder/f"development_case{case_index}_m{method}.npz",data=data,stats=stats)
        (folder/"bounds.json").write_text(json.dumps(rows,indent=2)+"\n",encoding="utf-8")

    def test_single_channel_A_laws(self):
        case=m.fixture(0);case["state"]=np.array([[0],[0],[0],[5],[0],[0],[0]],dtype=np.int64);case["a0"]=np.array([5.]);case["b0"]=np.array([1.])
        mask=np.array([1.,0,0,0,0,0,0])
        rows=[];replicates=16384
        bound=np.sqrt(np.log(12/1e-6)/(2*replicates))
        for method in (0,1):
            data,stats=m.simulate(case,method,replicates,seed_for(10,90,method,4601,0),mask)
            self.assertTrue(m.ledger(data,case["state"])["pass_check"])
            for k,t in enumerate(m.OBS):
                probability=np.exp(-14/5*t) if method==0 else np.exp(round(t/.001)*np.log1p(-14/5*.001))
                empirical=np.array([(data[:,k,3,0]<=count).mean() for count in range(6)])
                expected=binom.cdf(np.arange(6),5,probability)
                error=float(np.max(abs(empirical-expected)))
                self.assertLessEqual(error,bound)
                rows.append(dict(method=method,time=float(t),cdf_error=error,bound=float(bound)))
            np.savez_compressed(m.BASE/f"qualification/A_only_m{method}.npz",data=data,stats=stats)
        (m.BASE/"qualification/A_laws.json").write_text(json.dumps(rows,indent=2)+"\n",encoding="utf-8")

    def test_joint_R3_law_and_drift(self):
        case=m.fixture(1);reps=65536
        counts,delta=steps(case["state"],case["a0"],case["b0"],case["omega"],case["diff"],reps,seed_for(10,91,1,4601,0))
        dt=.00025
        expected=np.array([14,2,.1*12,29*12,12*11*3,.9*12,0])*dt
        deviation=abs(counts.mean(axis=0)-expected)
        tolerance=7*np.sqrt((expected+expected**2)/reps)+1e-12
        self.assertTrue(np.all(deviation<=tolerance),str((deviation,tolerance)))
        self.assertTrue(np.all(counts[:,3]<=1))
        # K3 is exactly Bernoulli(lambda3*dt), not a clipped proposal.
        p=29*12*dt
        self.assertLess(abs(np.mean(counts[:,3])-p),7*np.sqrt(p*(1-p)/reps))
        target=np.array([expected[0]-expected[1]+expected[2],expected[1]-expected[2]-expected[3]+expected[4]-expected[5],expected[3]-expected[4],-expected[0],-expected[3],expected[3],expected[5]])
        self.assertTrue(np.all(abs(delta.mean(axis=0)-target)<=7*np.sqrt((np.var(delta,axis=0)+1e-12)/reps)+1e-12))
        result=dict(expected_flux=expected.tolist(),observed_flux=counts.mean(axis=0).tolist(),tolerance=tolerance.tolist(),R3_max=int(counts[:,3].max()),replicates=reps)
        np.savez_compressed(m.BASE/"qualification/one_step_samples.npz",counts=counts,delta=delta)
        (m.BASE/"qualification/one_step_drift.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")

    def test_unchanged_IDEAL(self):
        from lattice import tau_leap_euler_1d,P17,steady_state
        import mg8
        initial=[np.full(16,round(1000*x),dtype=np.int64) for x in steady_state(P17)]
        a=[x.copy() for x in initial];b=[x.copy() for x in initial]
        args=(14.,29.,.1,1.,1.,1.,12.)
        seed=seed_for(10,92,0,4601,0)
        old=tau_leap_euler_1d(*a,1000.,*args,.2,1/1500,40,0.,0.,np.arange(16)*.2,20,seed)
        control=mg8.run_stoch(*b,mg8.J,1000.,*args,.2,1/1500,40,20,seed)
        for i in range(3):np.testing.assert_array_equal(old[i],control[i]);np.testing.assert_array_equal(a[i],b[i])

    def test_replay_and_invalid_inputs(self):
        case=m.fixture(3)
        for method in (0,3):
            first=m.simulate(case,method,16,seed_for(10,93,method,4601,0))
            second=m.simulate(case,method,16,seed_for(10,93,method,4601,0))
            for a,b in zip(first,second):np.testing.assert_array_equal(a,b)
        broken=m.fixture(0);broken["a0"][0]=0
        with self.assertRaises(ValueError):m.simulate(broken,1,1,99)
        with self.assertRaises(ValueError):m._exits4(5,-.1,0.,0.,0.)
        with self.assertRaises(ValueError):m._exits4(5,.9,.2,0.,0.)
        self.assertEqual(m._exits4(5,1.,0.,0.,0.),(5,0,0,0))


if __name__=="__main__":unittest.main(verbosity=2)
