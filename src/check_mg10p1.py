"""MG058 independent dense-convolution/direct-dot motif instrument oracle."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import json,hashlib
import numpy as np

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10P1"


def circular(a,b):return min((a-b)%100.,(b-a)%100.)


def audit():
    registration=json.loads((BASE/"registration.json").read_text(encoding="utf-8"))
    for name,expected in registration["source_hashes"].items():
        path=ROOT/name
        raw=path.read_bytes() if name.endswith(".npz") else path.read_text(encoding="utf-8").replace("\r\n","\n").encode()
        assert hashlib.sha256(raw).hexdigest()==expected,name
    p=registration["parameters"];n=p["N"];h=p["h"];length=p["L"]
    k=2*np.pi*np.fft.fftfreq(n,d=h)
    kernel=np.fft.ifft(np.exp(-.5*(p["sigma"]*k)**2)).real
    convolution=kernel[(np.arange(n)[:,None]-np.arange(n)[None,:])%n]
    def smooth(field):
        z=field-field[0];z=z-np.sum(z)/n
        y=convolution@z
        return y-np.sum(y)/n
    with np.load(ROOT/"results/MG003/ref_profiles.npz") as raw:template=smooth(raw["h0.2_b29.0"][2])
    basis=np.array([np.roll(template,shift) for shift in range(n)])
    norm2=np.sum(template**2);grid=(np.arange(n)-n//2)*h
    def assess(field):
        y=smooth(field);energy=float(np.sum(y**2))
        floor=1e-8*max(1.,140/9,float(np.max(np.abs(template))))
        if np.sqrt(energy/n)<=floor:return dict(status="FLAT",position=None,amplitude=0.,error=None,local_energy=None)
        scores=basis@y;peak=float(np.max(scores));best=int(np.where(scores>=peak-1e-12*max(1.,abs(peak)))[0][0]);a=float(scores[best]/norm2)
        result=dict(status="NO_REFERENCE_MATCH",position=None,amplitude=a,error=None,local_energy=None)
        if min(abs(a-p["amp_min"]),abs(a-p["amp_max"]))<=p["border"]:result["status"]="BORDERLINE";return result
        if a<=p["amp_min"] or a>=p["amp_max"]:return result
        far=np.array([circular(i*h,best*h)>p["ambiguity_radius"] for i in range(n)])
        if max(scores[far])>=p["ambiguity_ratio"]*scores[best]:result["status"]="AMBIGUOUS_MATCH";return result
        residual=y-a*basis[best];error=float(np.sqrt(np.sum(residual**2))/(a*np.sqrt(norm2)))
        position=(best*h+length/2)%length-length/2
        local=np.array([circular(v,position)<10 for v in grid]);fraction=float(np.sum(y[local]**2)/energy)
        result.update(error=error,local_energy=fraction)
        if abs(error-p["error_max"])<=p["border"] or abs(fraction-p["energy_min"])<=p["border"]:result["status"]="BORDERLINE";return result
        if error>p["error_max"] or fraction<p["energy_min"]:return result
        c0,c1,c2=scores[(best-1)%n],scores[best],scores[(best+1)%n];curvature=c0+c2-2*c1
        delta=.5*(c0-c2)/curvature if curvature<0 else 0.
        delta=min(.5,max(-.5,float(delta)))
        result.update(status="REFERENCE_MATCH",position=float(((best+delta)*h+length/2)%length-length/2))
        return result
    def compare(a,b):
        assert a["status"]==b["status"],(a,b)
        for key in ("amplitude","error","local_energy"):
            if a[key] is None:assert b[key] is None
            else:np.testing.assert_allclose(a[key],b[key],rtol=1e-9,atol=1e-9)
        if a["position"] is None:assert b["position"] is None
        else:assert circular(a["position"],b["position"])<=1e-7
    data=json.loads((BASE/"fixtures/evaluation.json").read_text(encoding="utf-8"))
    for item in data["raw"]:assert hashlib.sha256((BASE/"fixtures"/item["path"]).read_bytes()).hexdigest()==item["sha256"]
    checks=[];comparisons=0
    with np.load(BASE/"fixtures/frames.npz") as raw:fields=raw["fields"]
    assert len(fields)==len(registration["frames"])==len(data["frames"])
    for case,field,item in zip(registration["frames"],fields,data["frames"]):
        actual=assess(field);compare(actual,item["result"]);comparisons+=1
        expected=case["expected"]
        passed=actual["status"]=="REFERENCE_MATCH" and circular(actual["position"],case["position"])<=.1+1e-8 if expected=="match" else actual["status"]!="REFERENCE_MATCH" if expected=="not_match" else actual["status"]==expected
        assert passed==item["passed"]
        checks.append(dict(id=case["id"],passed=bool(passed),status=actual["status"]))
    for case,item in zip(registration["sequences"],data["sequences"]):
        with np.load(BASE/"fixtures"/(case["id"]+".npz")) as raw:fields=raw["fields"];times=raw["times"]
        results=[assess(f) for f in fields]
        for a,b in zip(results,item["result"]["observations"]):compare(a,b);comparisons+=1
        status="INITIAL_UNRESOLVED";duration=None;reason=None;persistent=False;matched=None
        if results[0]["status"]=="REFERENCE_MATCH" and circular(results[0]["position"],case["x0"])<=h/2+p["border"]:
            last=0
            for index in range(1,len(times)):
                dt=times[index]-times[index-1]
                if dt>p["max_gap"]+p["border"]:reason="SAMPLING_GAP";break
                if results[index]["status"]!="REFERENCE_MATCH":reason=results[index]["status"];break
                if circular(results[index]["position"],results[index-1]["position"])>p["max_speed"]*dt+p["border"]:reason="MOTION_AMBIGUITY";break
                last=index
            status="CONTINUOUS" if reason is None else "CENSORED";duration=float(times[last]-times[0]);matched=float(times[last]);persistent=duration>=p["persist"]-p["border"]
        result=item["result"]
        assert (status,duration,reason,persistent,matched)==(result["status"],result["visible_duration"],result["censor_reason"],result["persistent_reference_signal"],result["matched_through"])
        passed=(status,duration,reason,persistent)==(case["status"],case["duration"],case["reason"],case["persistent"])
        assert passed==item["passed"];checks.append(dict(id=case["id"],passed=passed,status=status))
    assert len(data["invalid"])==registration["invalid_checks"]
    for item in data["invalid"]:checks.append(dict(id="invalid"+str(item["id"]),passed=item["passed"]))
    passed=all(c["passed"] for c in checks);assert passed==data["passed"]
    result=dict(status="INDEPENDENT INSTRUMENT RECOMPUTATION PASS",qualified=passed,fixture_checks=checks,raw_frame_comparisons=comparisons,
        verdict="QUALIFIED ON REGISTERED SYNTHETIC FIXTURES" if passed else "QUALIFICATION FAILED",
        dependencies="Same author; shares reference, constants, NumPy and Gaussian Fourier multiplier definition; uses dense convolution/direct dots and separate decisions/tracking, not instrument functions. Invalid-API rejection outcomes are reported from the instrument's explicit tests, not separately reimplemented.")
    (BASE/"independent_check.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    failures=[c for c in checks if not c["passed"]]
    lines=["# MG10-P1 single-motif instrument qualification", "", "**"+result["verdict"]+"**. PI closure pending.","",
        f"{len(checks)} registered fixture/API assertions; {comparisons} frame assessments independently reconstructed with dense convolution/direct scores. Failed assertions: {len(failures)}.","",
        "This supports only a single isolated reference-like X motif and sampled field continuity. It does not establish physical particles, token lineage, multi-particle tracking, unique physical origin, death times, indefinite survival or equivalence of prior simulations. No D1/D2 array was reclassified and no stochastic trajectory generated.","",
        "Gaps, multiple motifs and excessive motion censor at the first unresolved observation. Identical replacement at the same location is observationally indistinguishable; the fixture intentionally demonstrates this limit. Flat/depleted fields, prescribed synthetic noise, extended waves and the registered broad impostors are controls, not a quantified false-positive rate for actual Model G noise.","",
        "Full input arrays, exact matrices/seeds/expectations, classification/continuity records, source/reference hashes and separate oracle are retained. Read execution.json for controller charges/stops. All frozen assertions, including failures if any, are preserved.", ""]
    if failures:lines += ["## Failed assertions", ""]+[f"- {c['id']}: {c.get('status','API check')}" for c in failures]
    (BASE/"MG10P1_results.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    print(f"{result['verdict']}: {len(checks)} assertions, {comparisons} independently reconstructed frames",flush=True)


if __name__=="__main__":audit()
