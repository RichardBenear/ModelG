"""Separate raw reconstruction, Poisson recurrence, contrast and budget audit."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):os.environ[key]="1"
from pathlib import Path
import hashlib,json,math,time
import numpy as np
from scipy.optimize import brentq
from scipy.special import stdtr
from check_mg10d1_development import check_raw

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10E1/preparation"
def read(path):return json.loads(path.read_text(encoding="utf-8"))
def main():
    started=time.perf_counter();plan=read(BASE/"preparation_plan.json")
    for name,expected in plan["source_hashes"].items():
        path=ROOT/name;content=path.read_bytes() if name.endswith(".npz") else path.read_text(encoding="utf-8").replace("\r\n","\n").encode()
        assert hashlib.sha256(content).hexdigest()==expected,name
    masks=[]
    for arm in range(4):
        expected=np.ones(500)
        if arm in (1,2):
            center=250 if arm==1 else 0
            for node in range(500):
                if min((node-center)%500,(center-node)%500)<=25:expected[node]=0
        if arm==3:expected.fill(.898)
        np.testing.assert_array_equal(expected,plan["masks"][str(arm)]);masks.append(expected)
    q=read(BASE/"qualification/qualification_results.json");assert q["passed"] and all(v["passed"] for v in q["tests"])
    raw_checks=[]
    for item in q["raw_checks"]:
        v=check_raw(ROOT/item["file"]);assert v["sha256"]==item["sha256"];raw_checks.append(v)
    for arm in range(4):
        with np.load(BASE/f"qualification/mask_a{arm}.npz") as raw:
            np.testing.assert_array_equal(raw["a_activity"],masks[arm])
            assert not np.any(raw["data"][:,:,7,masks[arm]==0]) and not np.any(raw["data"][:,:,20,masks[arm]==0])
    for test in read(BASE/"qualification/source_tests.json"):
        path=BASE/"qualification"/test["raw"];assert hashlib.sha256(path.read_bytes()).hexdigest()==test["sha256"]
        with np.load(path) as raw:counts=np.sort(raw["data"][:,0,7,0]);alpha=float(raw["a_activity"][0])
        lam=14000/24000*alpha;mass=math.exp(-lam);cdf=0.;distance=0.
        for k in range(int(counts[-1])+1):
            if k:mass*=lam/k
            cdf+=mass;empirical=int(np.searchsorted(counts,k,side="right"))/len(counts)
            distance=max(distance,abs(empirical-cdf))
        epsilon=math.sqrt(math.log(6000)/(2*len(counts)))
        np.testing.assert_allclose([distance,epsilon,lam,float(counts.mean())],[test["distance"],test["epsilon"],test["expected_mean"],test["mean"]],atol=1e-13)
        assert distance+epsilon<=.03 and (alpha!=0 or counts[-1]==0)
    saved=read(BASE/"qualification/contrast_results.json")
    with np.load(BASE/"qualification/contrast_inputs.npz") as raw:
        late=raw["late"];shift=raw["shift"];sample=raw["sample"];linear=raw["linear"]
    average=(linear[60]/2+linear[80]/2+linear[61:80].sum(axis=0))/20
    np.testing.assert_allclose(average,saved["late_average"])
    for name,a,b in (("primary_local_remote",1,2),("local_sham",1,0),("remote_sham",2,0),("uniform_sham",3,0)):
        value=late[1,a]-late[1,b]-late[0,a]+late[0,b]
        shifted=late+shift[None]
        np.testing.assert_allclose(value,saved["effects"][name],atol=1e-12)
        np.testing.assert_allclose(value,shifted[1,a]-shifted[1,b]-shifted[0,a]+shifted[0,b],atol=1e-12)
    mean=math.fsum(sample)/len(sample);variance=math.fsum((float(x)-mean)**2 for x in sample)/(len(sample)-1)
    critical=brentq(lambda x:stdtr(len(sample)-1,x)-.9875,0,100)
    half=critical*math.sqrt(variance/len(sample));np.testing.assert_allclose([mean-half,mean+half],saved["primary_interval"]["ci_nominal"],rtol=1e-10)
    parent_forecast={};branch_forecast={}
    for s in range(2):
        for k in range(2):
            parent_path=BASE/f"devparents/parent_s{s}_k{k}.json";meta=read(parent_path);assert meta["status"]=="COMPLETE"
            checked=check_raw(parent_path.parent/meta["raw"]);assert checked["sha256"]==meta["sha256"];raw_checks.append(checked)
            with np.load(parent_path.parent/meta["raw"]) as raw:state=raw["data"][0,-1,:7].copy();a0=raw["a0"].copy();b0=raw["b0"].copy()
            old=[read(ROOT/f"results/MG10P3/preparation/timing/timing_s{s}_k{k}_r{r}.json")["compute_seconds"] for r in range(2)]
            parent_forecast[f"s{s}_k{k}"]=math.fsum(old)/2*5*1.25+15
            for a in range(4):
                path=BASE/f"timing/timing_a{a}_s{s}_k{k}.json";meta=read(path);assert meta["status"]=="COMPLETE"
                checked=check_raw(path.parent/meta["raw"]);assert checked["sha256"]==meta["sha256"];raw_checks.append(checked)
                with np.load(path.parent/meta["raw"]) as raw:
                    np.testing.assert_array_equal(raw["initial"],state);np.testing.assert_array_equal(raw["data"][0,0,:7],state)
                    assert not np.any(raw["data"][0,0,7:]);assert raw["clamp_a"] and raw["clamp_b"]
                    np.testing.assert_array_equal(raw["a0"],a0);np.testing.assert_array_equal(raw["b0"],b0);np.testing.assert_array_equal(raw["a_activity"],masks[a])
                    assert not np.any(raw["data"][:,:,7,masks[a]==0]) and not np.any(raw["data"][:,:,20,masks[a]==0])
                branch_forecast[f"a{a}_s{s}_k{k}"]=meta["compute_seconds"]*80*1.25+15
    ledger=read(BASE/"timing/ledger.json");assert not ledger["stop_reason"] and len(ledger["records"])==16
    total=8*(math.fsum(parent_forecast.values())+math.fsum(branch_forecast.values()))
    allowance=ledger["charged_worker_seconds"]+time.perf_counter()-started+5;assert allowance<900
    result=dict(passed=True,status="E1 PREPARATION PASS; PRODUCTION NOT APPROVED OR EXECUTED",raw_checks=raw_checks,
        parent_forecast_seconds=parent_forecast,branch_forecast_seconds=branch_forecast,forecast_total_seconds=total,
        forecast_core_hours=total/3600,eight_worker_wall_hours=total/28800,within_75_percent_reserve=total<=97200,
        controlled_development_seconds=ledger["charged_worker_seconds"],total_preparation_allowance_seconds=allowance,
        dependencies="Same author; shared reference/constants and earlier explicit stoichiometric reconstruction. Separate Poisson recurrence, empirical CDF, explicit trapezoidal weights, contrast arithmetic, fsum variance and Student CDF inversion. Not external replication.")
    (BASE/"independent_check.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k not in ("raw_checks","parent_forecast_seconds","branch_forecast_seconds")},indent=2))
if __name__=="__main__":main()
