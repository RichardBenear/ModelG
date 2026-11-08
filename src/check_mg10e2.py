"""MG074 separate projection/paired-statistical and complete raw/fork audit."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):os.environ[key]="1"
from pathlib import Path
import json,hashlib,math
import numpy as np
from scipy.special import stdtr
from scipy.optimize import brentq
from check_mg10d1_development import check_raw
from mg10d1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10E2"
def main(base=BASE,synthetic=False):
    base=Path(base);r=json.loads((base/"registration.json").read_text());e=json.loads((base/"evaluation.json").read_text())
    if not synthetic:
        for name,sha in r["source_hashes"].items():
            path=ROOT/name;content=path.read_bytes() if name.endswith(".npz") else path.read_text(encoding="utf-8").replace("\r\n","\n").encode()
            assert hashlib.sha256(content).hexdigest()==sha,name
    with np.load(ROOT/"results/MG003/ref_profiles.npz") as raw:ref=raw["h0.2_b29.0"][1:]
    kernel=np.fft.ifft(np.exp(-.5*(.4*2*np.pi*np.fft.fftfreq(500,d=.2))**2)).real
    dense=kernel[(np.arange(500)[:,None]-np.arange(500)[None,:])%500]
    template=np.array([dense@(f-f.mean()) for f in ref]);norm=np.sum(template**2,axis=-1)
    # Gaussian convolution is symmetric: move it to the template for an
    # independent linear-functional calculation on raw unsmoothed counts.
    weights=np.array([dense@f for f in template])/norm[:,None]
    score_base=dense@template[1];score_weights=np.array([np.roll(score_base,j) for j in range(500)])
    late={};checks=[];parents={}
    for job in r["parent_jobs"]:
        path=base/"parents"/(job["id"]+".npz");meta=path.with_suffix(".json")
        if not path.exists():continue
        checks.append(check_raw(path))
        with np.load(path) as raw:
            if meta.exists() and json.loads(meta.read_text()).get("status")=="COMPLETE":
                assert np.array_equal(raw["times"],np.arange(21)/4);assert raw["seed"]==job["seed"];assert hashlib.sha256(path.read_bytes()).hexdigest()==json.loads(meta.read_text())["raw_sha256"];parents[job["id"]]=(raw["data"][0,-1,:7].copy(),raw["a0"].copy(),raw["b0"].copy(),hashlib.sha256(path.read_bytes()).hexdigest())
    for job in r["branch_jobs"]:
        name=job["id"];path=base/"branches"/(name+".npz")
        if not path.exists():continue
        checks.append(check_raw(path))
        v=e["records"][name];assert hashlib.sha256(path.read_bytes()).hexdigest()==v["raw_sha256"]
        with np.load(path) as raw:
            data=raw["data"];times=raw["times"];fields=data[0,:,:3]/1000.
            arm=int(raw["arm"]);switch=float(raw["switch"]);assert arm==job["arm"] and switch==5.
            before=np.ones(500);after=np.ones(500)
            if arm:
                center=0 if arm==3 else 250
                for n in range(500):
                    if min((n-center)%500,(center-n)%500)<=25:before[n]=0
            if arm==1:after=before.copy()
            np.testing.assert_array_equal(raw["activity_before"],before);np.testing.assert_array_equal(raw["activity_after"],after)
            assert np.array_equal(times,np.arange(len(times))/4) and len(times)<=101
            assert raw["seed"]==job["seed"]
            for i in range(1,len(times)):
                off=(before if times[i-1]<switch else after)==0
                assert not np.any(data[0,i,7,off]-data[0,i-1,7,off]) and not np.any(data[0,i,20,off]-data[0,i-1,20,off])
            if not synthetic:
                parent=parents[job["parent"]];np.testing.assert_array_equal(raw["initial"],parent[0]);np.testing.assert_array_equal(data[0,0,:7],parent[0]);assert not np.any(data[0,0,7:])
                np.testing.assert_array_equal(raw["a0"],parent[1]);np.testing.assert_array_equal(raw["b0"],parent[2])
                assert raw["clamp_a"] and raw["clamp_b"]
                assert str(raw["parent_sha256"])==parent[3]
                if v["complete"]:
                    meta=json.loads(path.with_suffix(".json").read_text());assert meta["parent_sha256"]==parent[3] and meta["raw_sha256"]==v["raw_sha256"]
            centered=fields-fields.mean(axis=-1,keepdims=True);projection=np.sum(centered*weights[None],axis=-1)
            np.testing.assert_allclose(projection,v["projection"],rtol=1e-7,atol=1e-8)
            np.testing.assert_allclose(fields.mean(axis=-1),v["field_mean"],rtol=1e-12,atol=1e-12)
            # Dense smoothing and direct adjoint scores at all500 grid shifts;
            # no FFT correlation or production phase helper used here.
            smooth=centered@dense.T
            smooth_rms=np.sqrt(np.mean(smooth**2,axis=-1));moving=v["movement_metrics"]
            np.testing.assert_allclose(smooth_rms,moving["smooth_rms"],atol=1e-9,rtol=1e-7)
            scores=centered[:,1]@score_weights.T
            best=[next(j for j,x in enumerate(row) if x>=row.max()-1e-12*max(1.,abs(row.max()))) for row in scores]
            phases=[float((j*.2+50)%100-50) for j in best];assert phases==moving["x_phase_grid"]
            shifted=np.array([np.roll(template,j,axis=-1) for j in best]);amp=np.sum(smooth*shifted,axis=-1)/norm[None]
            np.testing.assert_allclose(amp,moving["x_aligned_projection"],atol=1e-8,rtol=1e-7)
            error=np.sqrt(np.mean((smooth-amp[:,:,None]*shifted)**2,axis=-1));np.testing.assert_allclose(error,moving["x_aligned_residual_rms"],atol=1e-8,rtol=1e-7)
            if len(times)>=21:
                early=(projection[0]/2+projection[20]/2+projection[1:20].sum(axis=0))/20;np.testing.assert_allclose(early,v["early_projection"],atol=1e-8)
            else:assert v["early_projection"] is None
            assert data[0,:,7].sum(axis=-1).tolist()==v["source_counts"] and data[0,:,20:22].sum(axis=-1).tolist()==v["external_inputs"]
            assert v["complete"]==bool(path.with_suffix(".json").exists() and json.loads(path.with_suffix(".json").read_text()).get("status")=="COMPLETE" and np.array_equal(times,np.arange(101)/4))
            if v["complete"]:
                value=(projection[80]/2+projection[100]/2+projection[81:100].sum(axis=0))/20;late[name]=value
                np.testing.assert_allclose(value,v["late_projection"],rtol=1e-7,atol=1e-8)
                profile=(fields[80]/2+fields[100]/2+fields[81:100].sum(axis=0))/20;np.testing.assert_allclose(profile,v["late_physical_profile"],atol=1e-10)
    for k in range(2):
        values=[]
        for b in range(8):
            ids=[f"a{a}_s{s}_k{k}_b{b}" for s in range(2) for a in range(4)]
            if not all(name in late for name in ids):assert f"k{k}_b{b}" not in e["block_effects"];continue
            z=lambda s,a:late[f"a{a}_s{s}_k{k}_b{b}"]
            for name,a,c in (("primary_restore_continuous",2,1),("restore_sham",2,0),("remote_restore_sham",3,0),("restore_remote",2,3)):
                effect=z(1,a)-z(1,c)-z(0,a)+z(0,c);np.testing.assert_allclose(effect,e["block_effects"][f"k{k}_b{b}"][name],atol=1e-8)
                if name=="primary_restore_continuous":values.append(float(effect[1]))
            np.testing.assert_allclose(np.array([z(1,a)-z(0,a) for a in range(4)]),e["block_effects"][f"k{k}_b{b}"]["seed_excess_by_arm"],atol=1e-8)
        summary=e["primary_by_step"][str(k)];n=len(values);assert summary["n"]==n
        if n:
            mean=math.fsum(values)/n;np.testing.assert_allclose(mean,summary["mean"],atol=1e-8)
            if n>1:
                variance=math.fsum((x-mean)**2 for x in values)/(n-1);critical=brentq(lambda x:stdtr(n-1,x)-.9875,0,100);half=critical*math.sqrt(variance/n)
                np.testing.assert_allclose([mean-half,mean+half],summary["ci_nominal"],atol=1e-8)
            else:assert summary["ci_nominal"] is None
        else:assert summary["mean"] is None and summary["ci_nominal"] is None
    assert len(parents)==e["parents_complete"] and len(late)==e["branches_complete"]
    expected="COMPLETE RESTORATION ESTIMATION" if len(parents)==32 and len(late)==128 else "INCOMPLETE RESTORATION ESTIMATION";assert e["status"]==expected
    charge=None
    if not synthetic:
        p=json.loads((base/"parents/ledger.json").read_text());assert p["ceiling"]==144000 and p["prior_charge"]==0
        charge=p["charged_worker_seconds"]
        if (base/"branches/ledger.json").exists():
            b=json.loads((base/"branches/ledger.json").read_text());assert b["ceiling"]==144000 and b["prior_charge"]==charge;charge=b["charged_worker_seconds"]
        if expected.startswith("COMPLETE "):assert not p["stop_reason"] and not b["stop_reason"] and len(p["records"])==32 and len(b["records"])==128
    atomic_json(base/"independent_check.json",dict(passed=True,status=expected,raw_checks=checks,charged_worker_seconds=charge,
      dependencies="Same author; shared constants/reference/NumPy, prior explicit integer checker. Dense Gaussian adjoint projection and all-shift grid-fit/norm/residual audit on unsmoothed raw fields, explicit temporal weights, manual contrasts and fsum/Student-CDF inversion; not external replication."))
if __name__=="__main__":main()
