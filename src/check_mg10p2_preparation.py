"""MG061 direct three-field metrics and binomial-CDF inversion checks."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import json,hashlib,math
import numpy as np
from scipy.stats import binom
from scipy.optimize import brentq
from check_mg10d1_development import check_raw

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10P2/preparation"


def main():
    plan=json.loads((BASE/"preparation_plan.json").read_text(encoding="utf-8"))
    for name,expected in plan["source_hashes"].items():
        p=ROOT/name;data=p.read_bytes() if name.endswith(".npz") else p.read_text(encoding="utf-8").replace("\r\n","\n").encode()
        assert hashlib.sha256(data).hexdigest()==expected,name
    q=json.loads((BASE/"qualification/qualification_results.json").read_text(encoding="utf-8"));assert q["passed"]
    with np.load(ROOT/"results/MG003/ref_profiles.npz") as raw:profile=raw["h0.2_b29.0"][1:]
    k=2*np.pi*np.fft.fftfreq(500,d=.2);kernel=np.fft.ifft(np.exp(-.5*(.4*k)**2)).real
    convolution=kernel[(np.arange(500)[:,None]-np.arange(500)[None,:])%500]
    def smooth(f):
        x=f-f[0];x-=x.mean();y=convolution@x;return y-y.mean()
    template=np.array([smooth(f) for f in profile])
    with np.load(BASE/"qualification/synthetic_metric_inputs.npz") as raw:fields=raw["data"][0,0,:3]/1000
    smoothed=np.array([smooth(f) for f in fields]);proj=np.sum(smoothed*template,axis=1)/np.sum(template**2,axis=1)
    saved=json.loads((BASE/"qualification/synthetic_metric_results.json").read_text(encoding="utf-8"))
    for name,values in dict(field_mean=fields.mean(axis=1),contrast_rms=fields.std(axis=1),smooth_rms=np.sqrt(np.mean(smoothed**2,axis=1)),projection=proj,residual_rms=np.sqrt(np.mean((smoothed-proj[:,None]*template)**2,axis=1))).items():
        np.testing.assert_allclose(values,np.array(saved["metrics"][name])[0,0],rtol=1e-9,atol=1e-10)
    bounds=json.loads((BASE/"qualification/binomial_bounds.json").read_text(encoding="utf-8"))
    for v in bounds:
        n,k=v["n"],v["k"]
        upper=1. if k==n else -math.expm1(math.log(.0125)/n) if k==0 else brentq(lambda p:binom.cdf(k,n,p)-.0125,0.,1.,xtol=1e-14)
        lower2=0. if k==0 else brentq(lambda p:binom.sf(k-1,n,p)-.025,0.,1.,xtol=1e-14)
        upper2=1. if k==n else brentq(lambda p:binom.cdf(k,n,p)-.025,0.,1.,xtol=1e-14)
        assert abs(upper-v["upper"])<1e-10
        np.testing.assert_allclose(v["interval"],[lower2,upper2],rtol=1e-9,atol=1e-10)
    raw_checks=[];timing={}
    for c in range(2):
        for s in range(2):
            for step in range(2):
                path=BASE/f"timing/timing_c{c}_s{s}_k{step}.json";d=json.loads(path.read_text(encoding="utf-8"));raw=path.parent/d["raw"]
                checked=check_raw(raw);assert checked["sha256"]==d["sha256"];raw_checks.append(checked)
                with np.load(raw) as a:
                    factor=1 if c==0 else .82;A=round(4200000*factor);B=round(135333333*factor)
                    ae=14*A/4200000;be=29*B/135333333;bg=np.array([ae/.9,ae/.9,be/(ae/.9)])
                    fields=np.repeat(bg[:,None],500,axis=1)+(profile if s else 0)
                    expected=np.zeros((7,500),dtype=np.int64);expected[:3]=np.rint(1000*fields).astype(np.int64);expected[3]=A;expected[4]=B
                    assert np.array_equal(a["initial"],expected)
                    assert np.all(a["a0"]==4200000) and np.all(a["b0"]==135333333)
                    assert a["clamp_a"]==(c==0) and a["clamp_b"]==(c==0)
                timing[f"c{c}_s{s}_k{step}"]=d["compute_seconds"]*10*1.25*128+16*5
    ledger=json.loads((BASE/"timing/ledger.json").read_text(encoding="utf-8"));assert not ledger["stop_reason"] and len(ledger["records"])==8
    projected=sum(timing.values())
    result=dict(status="PREPARATION AND INDEPENDENT CHECK PASS; PRODUCTION NOT REGISTERED",passed=True,
        reference_rms=np.sqrt(np.mean(template**2,axis=1)).tolist(),metric_fields_checked=3,binomial_count_cases_checked=len(bounds),prepared_raw_checks=raw_checks,
        candidate_paths=1024,candidate_jobs=128,forecast_by_cell_seconds=timing,forecast_seconds=projected,forecast_core_hours=projected/3600,
        eight_worker_wall_minutes=projected/8/60,production_candidate_ceiling_core_hours=4,forecast_within_75_percent_ceiling=projected<=10800,
        development_charged_seconds=ledger["charged_worker_seconds"],
        dependencies="Independent dense convolution/projection, explicit state construction and binomial-CDF inversion; shares NumPy/SciPy, frozen reference and earlier explicit integer reconstruction. Same author, not external replication.")
    (BASE/"independent_check.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k not in ("prepared_raw_checks","forecast_by_cell_seconds")},indent=2))


if __name__=="__main__":main()
