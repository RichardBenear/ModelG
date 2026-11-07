"""MG065 independent episodes, direct trigonometric alignment and averages."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import json,hashlib,time,math
import numpy as np
from check_mg10d1_development import check_raw

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10P3/preparation"


def independent_episodes(observations,times):
    n=len(times);out=[]
    def accepted(i):return observations[i]["status"]=="REFERENCE_MATCH"
    def linked(i,j):
        dt=times[j]-times[i];distance=abs((observations[j]["position"]-observations[i]["position"]+50)%100-50)
        return dt<=1+1e-10 and distance<=1.5*dt+1e-10
    starts=[i for i in range(n) if accepted(i) and (i==0 or not accepted(i-1) or not linked(i-1,i))]
    for i in starts:
        last=i
        for j in range(i+1,n):
            if not accepted(j) or not linked(j-1,j):break
            last=j
        if last==n-1:reason="END_OF_RECORD"
        elif not accepted(last+1):reason=observations[last+1]["status"]
        elif times[last+1]-times[last]>1+1e-10:reason="SAMPLING_GAP"
        else:reason="MOTION_SPLIT"
        duration=float(times[last]-times[i])
        out.append(dict(start=float(times[i]),last=float(times[last]),duration=duration,persistent=duration>=20-1e-10,termination=reason))
    return out


def main():
    started=time.perf_counter()
    plan=json.loads((BASE/"preparation_plan.json").read_text(encoding="utf-8"))
    for name,expected in plan["source_hashes"].items():
        p=ROOT/name;content=p.read_bytes() if name.endswith(".npz") else p.read_text(encoding="utf-8").replace("\r\n","\n").encode()
        assert hashlib.sha256(content).hexdigest()==expected,name
    q=json.loads((BASE/"qualification/qualification_results.json").read_text(encoding="utf-8"));assert q["passed"]
    episodes=json.loads((BASE/"qualification/episode_inputs.json").read_text(encoding="utf-8"))
    for case in episodes:
        actual=independent_episodes(case["observations"],case["times"])
        assert actual==case["result"],case["kind"]
        assert any(r["persistent"] for r in actual)==case["expected_persistent"]
    k=2*np.pi*np.fft.fftfreq(500,d=.2);kernel=np.fft.ifft(np.exp(-.5*(.4*k)**2)).real
    offsets=(np.arange(500)[:,None]-np.arange(500)[None,:])%500
    convolution=kernel[offsets]
    def smooth(f):
        z=f-f[0];z-=z.mean();y=convolution@z;return y-y.mean()
    def shifted(f,position):
        # Even-N periodic real trigonometric interpolant. This directly
        # evaluates the Dirichlet kernel rather than FFT phase rotation.
        d=(np.arange(500)+position/.2+250)%500-250
        weights=np.ones(500);mask=abs(d)>1e-12
        weights[mask]=np.sin(np.pi*d[mask])/(500*np.tan(np.pi*d[mask]/500))
        return weights[offsets]@f
    with np.load(ROOT/"results/MG003/ref_profiles.npz") as raw:profile=raw["h0.2_b29.0"][1:]
    template=np.array([smooth(f) for f in profile]);norm2=np.sum(template**2,axis=1)
    averages=json.loads((BASE/"qualification/average_results.json").read_text(encoding="utf-8"))
    for case in averages:
        path=BASE/"qualification"/case["raw"];assert hashlib.sha256(path.read_bytes()).hexdigest()==case["sha256"]
        with np.load(path) as raw:fields=raw["fields"];phases=raw["phases"]
        initial=np.array([[smooth(f) for f in frame] for frame in fields])
        for aligned,key in ((False,"fixed"),(True,"aligned")):
            values=np.array([[shifted(f,p) for f in frame] for frame,p in zip(initial,phases)]) if aligned else initial
            projection=np.sum(values*template[None],axis=-1)/norm2[None]
            residual=values-projection[:,:,None]*template[None]
            individual=np.sqrt(np.mean(residual**2,axis=-1));average=values.mean(axis=0)
            coeff=np.sum(average*template,axis=1)/norm2
            error=np.sqrt(np.mean((average-coeff[:,None]*template)**2,axis=1));mean=individual.mean(axis=0)
            saved=case[key]
            results=dict(profile=average,projection=coeff,residual_rms=error,contrast_rms=np.sqrt(np.mean(average**2,axis=1)),
                         mean_individual_residual_rms=mean,individual_projection=projection,individual_residual_rms=individual)
            for name,value in results.items():np.testing.assert_allclose(value,saved[name],rtol=1e-8,atol=1e-9)
            for j in range(3):
                if error[j]<=1e-12:assert saved["gain"][j] is None
                else:np.testing.assert_allclose(mean[j]/error[j],saved["gain"][j],rtol=1e-7,atol=1e-7)
    raw_checks=[]
    for item in q["raw_checks"]:
        checked=check_raw(ROOT/item["file"]);assert checked["sha256"]==item["sha256"];raw_checks.append(checked)
    forecast={};source_times={}
    for start in range(2):
        for step in range(2):
            seconds=[]
            for rep in range(2):
                path=BASE/f"timing/timing_s{start}_k{step}_r{rep}.json";v=json.loads(path.read_text(encoding="utf-8"))
                checked=check_raw(path.parent/v["raw"]);assert checked["sha256"]==v["sha256"];raw_checks.append(checked)
                seconds.append(v["compute_seconds"])
            mean=math.fsum(seconds)/2
            source_times[f"s{start}_k{step}"]=seconds
            forecast[f"s{start}_k{step}"]=mean*25*1.25+15
    total=sum(forecast.values())*16
    ledger=json.loads((BASE/"timing/ledger.json").read_text(encoding="utf-8"));assert not ledger["stop_reason"] and len(ledger["records"])==8
    allowance=ledger["charged_worker_seconds"]+time.perf_counter()-started+5
    assert allowance<900
    result=dict(status="P3 PREPARATION AND SEPARATE RECOMPUTATION PASS; SCIENTIFIC PRODUCTION NOT REGISTERED",passed=True,
        episode_cases=len(episodes),averaging_cases=len(averages),direct_interpolation="Periodic even-N real Dirichlet kernel, independent of FFT shift",
        raw_checks=raw_checks,timing_seconds_by_class=source_times,forecast_seconds_per_path=forecast,forecast_total_seconds=total,
        forecast_core_hours=total/3600,eight_worker_wall_hours=total/3600/8,production_candidate_ceiling_core_hours=24,
        within_75_percent_reserve=total<=64800,controlled_development_seconds=ledger["charged_worker_seconds"],total_preparation_allowance_seconds=allowance,
        dependencies="Same author; shares constants/reference/NumPy and earlier explicit reconstruction. Independent candidate-start episode scan, dense smoothing/direct trigonometric interpolation/averaging. Not external replication.")
    (BASE/"independent_check.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k not in ("raw_checks","timing_seconds_by_class","forecast_seconds_per_path")},indent=2))


if __name__=="__main__":main()
