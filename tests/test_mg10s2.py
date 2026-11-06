"""MG049 development-only actual-lattice and oracle qualification."""
from pathlib import Path
import sys,json,unittest
import numpy as np
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import mg10s2 as m
from mg10s2_measure import measure
from runner import seed_for


class LatticeQualification(unittest.TestCase):
    @classmethod
    def setUpClass(cls):(m.BASE/"qualification").mkdir(parents=True,exist_ok=True)

    def test_snapshot_oracle_packing(self):
        rows=[]
        for case in range(4):
            fixture=m.fixture(case)
            np.testing.assert_array_equal(fixture["diff"],[25,25,300])
            self.assertTrue(np.all(fixture["initial"][4]>=fixture["initial"][1]))
            for step in range(3):
                for sampler in (0,1):
                    seed=seed_for(10,case,2*step+sampler,4901,0)
                    event,pop,stats=m.snapshot(fixture,step,sampler,8,seed)
                    self.assertEqual(event.dtype,np.uint32)
                    np.testing.assert_array_equal(pop,m.reconstruct(fixture["initial"],event))
                    self.assertTrue(np.all(pop>=0))
                    self.assertTrue(np.all(pop.sum(axis=(1,2))==fixture["initial"].sum()))
                    if sampler==0:np.testing.assert_array_equal(stats[:,0],1);np.testing.assert_array_equal(stats[:,1],0)
                    changed=event.copy();changed[0,0,0]+=1
                    self.assertNotEqual(m.population_hash(pop),m.population_hash(m.reconstruct(fixture["initial"],changed)))
                    rows.append(dict(case=case,step=step,sampler=sampler,population_sha256=m.population_hash(pop)))
        (m.BASE/"qualification/snapshot_checks.json").write_text(json.dumps(rows,indent=2)+"\n",encoding="utf-8")

    def test_unchanged_kernels(self):
        fixture=m.fixture(0);seed=seed_for(10,90,0,4901,0)
        event,pop,stats=m.snapshot(fixture,0,0,4,seed)
        raw,other=m.s1.leap_chunk(fixture["initial"],fixture["a0"],fixture["b0"],fixture["omega"],fixture["diff"],np.ones(7),m.TAUS[0],np.array([m.TAUS[0]]),4,seed)
        np.testing.assert_array_equal(event,raw[:,0,7:20]);np.testing.assert_array_equal(pop,raw[:,0,:7]);np.testing.assert_array_equal(stats,other)
        fixture=m.fixture(4)
        first=m.short(fixture,0,4,seed)
        second=m.s1.ssa_chunk(fixture["initial"],fixture["a0"],fixture["b0"],fixture["omega"],fixture["diff"],np.ones(7),m.OBS,4,seed)
        for a,b in zip(first,second):np.testing.assert_array_equal(a,b)

    def test_FFT_scaling_phase_and_Nyquist(self):
        n=16;omega=100.;t=.0005;scale=np.sqrt(n*omega*t)
        i=np.arange(n);pop=np.zeros((1,3,n));pop[0,0]=2*scale/n*np.cos(2*np.pi*i/n);pop[0,1]=2*scale/n*np.sin(2*np.pi*i/n)
        values,z=measure(pop,np.zeros((3,n)),omega,t,[1,8])
        np.testing.assert_allclose(values[0,0],[1,1,0,0,0,0,0,-1],rtol=0,atol=1e-8)
        shifted,_=measure(np.roll(pop,3,axis=-1),np.zeros((3,n)),omega,t,[1,8])
        np.testing.assert_allclose(values[...,:6],shifted[...,:6],rtol=0,atol=1e-8)
        for s in range(3):pop[0,s]=(s+1)*scale/n*(-1.)**i
        values,z=measure(pop,np.zeros((3,n)),omega,t,[1,8])
        np.testing.assert_allclose(values[0,1],[1,4,9,2,3,6,2,0],rtol=0,atol=1e-8)
        self.assertEqual(z[0,1,1].imag,0)

    def test_covariance_event_enumeration(self):
        n=4;t=1/24000;omega=1.
        initial=np.array([[10,12,14,16],[8,11,14,17],[1,2,3,4],[40]*4,[40]*4,[0]*4,[0]*4],dtype=np.int64)
        fixture=dict(initial=initial,a0=np.full(n,40.),b0=np.full(n,40.),omega=omega,diff=m.DIFF)
        def enumeration(mode,instantaneous):
            phase=np.exp(-2j*np.pi*mode*np.arange(n)/n);cov=np.zeros((3,3),complex)
            for node in range(n):
                G,X,Y,A=initial[:4,node]
                groups=[(A,[14*t/40],[(3,0)]),
                        (G,[t,25*t,25*t],[(0,1),(0,-1),(0,-2)]),
                        (X,[.1*t,29*t,.9*t,25*t,25*t],[(1,0),(1,2),(1,3),(1,-1),(1,-2)]),
                        (Y,[X*(X-1)*t,300*t,300*t],[(2,1),(2,-1),(2,-2)])]
                for count,probabilities,channels in groups:
                    vectors=[]
                    for source,target in channels:
                        delta=np.zeros((3,n))
                        if source<3:delta[source,node]-=1
                        if target in (-1,-2):delta[source,(node-1)%n if target==-1 else (node+1)%n]+=1
                        elif target<3:delta[target,node]+=1
                        vectors.append(delta@phase/np.sqrt(n*omega*t))
                    mean=sum(p*v for p,v in zip(probabilities,vectors))
                    matrix=sum(p*np.outer(v,v.conj()) for p,v in zip(probabilities,vectors))
                    if not instantaneous:matrix-=np.outer(mean,mean.conj())
                    cov+=count*matrix
            return cov
        for instant in (False,True):
            actual=m.covariance(fixture,t,[1,2],instant)
            for k,mode in enumerate((1,2)):np.testing.assert_allclose(actual[k],enumeration(mode,instant),rtol=1e-12,atol=1e-10)
        mean=m.conditional_mean(fixture,t)
        self.assertAlmostEqual(float((mean-initial[:3]).sum()),t*(14*n-.9*initial[1].sum()),places=10)

    def test_partial_shard_retention(self):
        import run_mg10s2 as driver
        folder=m.BASE/"qualification/partial";folder.mkdir(exist_ok=True)
        job=dict(id="development_partial",kind="snapshot",case=0,method=0,block=0,reps=8,chunk=4,offset=0,stage=4901,warmup=False)
        original=m.snapshot;count=[0]
        def interrupted(*args):
            count[0]+=1
            if count[0]==2:raise RuntimeError("Intentional development interruption")
            return original(*args)
        with patch.object(m,"snapshot",side_effect=interrupted):self.assertEqual(driver.worker(job,folder),2)
        record=json.loads((folder/"development_partial.json").read_text())
        self.assertEqual(record["replicates"],4);self.assertEqual(len(record["shards"]),1)
        self.assertTrue((folder/record["shards"][0]["file"]).exists())


if __name__=="__main__":unittest.main(verbosity=2)
