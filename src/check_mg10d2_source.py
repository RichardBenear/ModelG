"""MG055 independent scalar histogram/Decimal-CDF and raw-stream checker."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
from decimal import Decimal,localcontext
import json,hashlib,math
import numpy as np
from check_mg10d1_development import check_raw

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10D2"


def main():
    registration=json.loads((BASE/"diagnostic_registration.json").read_text(encoding="utf-8"))
    for name,expected in registration["source_hashes"].items():
        actual=hashlib.sha256((ROOT/name).read_text(encoding="utf-8").replace("\r\n","\n").encode()).hexdigest()
        assert actual==expected,name
    folder=BASE/"development/diagnostic"
    q=json.loads((folder/"qualification.json").read_text(encoding="utf-8"));assert q["passed"]
    checks=[]
    for item in q["primary"]+q["supplements"]:
        path=folder/item["raw"];assert hashlib.sha256(path.read_bytes()).hexdigest()==item["sha256"]
        with np.load(path) as raw:
            h=raw["hist"];overflow=raw["overflow"]
            assert int(h.sum())==item["job"]["reps"] and len(overflow)==h[-1]
            total=sum(int(h[k])*k for k in range(64))+sum(int(k) for k in overflow)
            squares=sum(int(h[k])*k*k for k in range(64))+sum(int(k)**2 for k in overflow)
            assert total==int(raw["total"]) and squares==int(raw["squares"])
        support=list(range(64))
        if len(overflow):support=list(range(max(64,int(overflow.max())+1)))
        counts=[int(h[k]) if k<64 else sum(int(v)==k for v in overflow) for k in support]
        n=item["job"]["n"];p=item["job"]["p"];R=item["job"]["reps"]
        with localcontext() as ctx:
            ctx.prec=80;dn=Decimal(n);dp=Decimal.from_float(p);dq=1-dp
            pmf=(dn*dq.ln()).exp();cdf=pmf;cumulative=0;distance=0.
            for k,count in enumerate(counts):
                if k:
                    pmf=pmf*(dn-k+1)*dp/(Decimal(k)*dq);cdf+=pmf
                cumulative+=count
                distance=max(distance,abs(cumulative/R-float(cdf)))
            distance=max(distance,1-float(cdf))
        comparisons=6 if item in q["primary"] else 1
        alpha=.01 if comparisons==6 else .001
        epsilon=math.sqrt(math.log(2*comparisons/alpha)/(2*R))
        assert abs(distance-item["cdf_distance"])<2e-12
        assert abs(distance+epsilon-item["upper"])<2e-12
        assert total/R==item["mean"]
        assert abs(squares/R-(total/R)**2-item["variance"])<1e-12
        checks.append(dict(id=item["job"]["id"],mean=total/R,expected_mean=n*p,lower=max(0,distance-epsilon),upper=distance+epsilon,passed=distance+epsilon<=item["margin"]))
    revised=[c for c in checks if c["id"].endswith("_m1")];assert len(revised)==2 and all(c["passed"] for c in revised)
    assert all(c["passed"] for c in checks if c["id"].startswith("supplement"))
    raw_checks=[]
    for item in q["raw"]:
        checked=check_raw(ROOT/item["file"]);assert checked["sha256"]==item["sha256"];raw_checks.append(checked)
    forecast={}
    for a in range(2):
        for s in range(2):
            for k in range(2):
                path=BASE/f"development/timing/timing_a{a}_s{s}_k{k}.json";v=json.loads(path.read_text(encoding="utf-8"))
                checked=check_raw(path.parent/v["raw"]);assert checked["sha256"]==v["sha256"];raw_checks.append(checked)
                forecast[f"a{a}_s{s}_k{k}"]=v["compute_seconds"]*600*1.25+5
    ledger=json.loads((BASE/"development/timing/ledger.json").read_text(encoding="utf-8"))
    assert not ledger["stop_reason"] and len(ledger["records"])==8 and ledger["charged_worker_seconds"]<900
    projected=sum(forecast.values())*4
    result=dict(status="INDEPENDENT SOURCE QUALIFICATION CHECK PASS",passed=True,source_comparisons=checks,raw_checks=raw_checks,
        law_reference="80-digit Decimal finite-binomial PMF/CDF from exact binary float p, independent of sampler and diagnostic float recurrence; shared NumPy for raw arrays",
        forecast_by_class_seconds=forecast,forecast_seconds=projected,forecast_core_hours=projected/3600,
        within_24_hour_75_percent_reserve=projected<=.75*86400,qualification_charged_seconds=ledger["charged_worker_seconds"],
        legacy_resolved_cdf_discrepancies=[c["id"] for c in checks if c["id"].endswith("_m0") and c["lower"]>0])
    (BASE/"development/independent_source_check.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    lines=["# MG055 source-law qualification", "", "Corrected source sampler QUALIFIED on the registered scalar tests; original D1 data remain unchanged.","",
        "|Case/sampler|Mean|Expected|CDF lower|CDF upper|Within diagnostic margin|", "|---|---:|---:|---:|---:|---|"]
    for c in checks:lines.append(f"|{c['id']}|{c['mean']:.9g}|{c['expected_mean']:.9g}|{c['lower']:.7g}|{c['upper']:.7g}|{c['passed']}|")
    lines += ["",f"Corrected32-path forecast: {projected/3600:.6f}core-hours;24-hour cap reserve check: {result['within_24_hour_75_percent_reserve']}.",
        f"Diagnostic/qualification charge: {ledger['charged_worker_seconds']:.6f}worker-seconds (<900).", "",
        "Raw histograms/moments/overflow values and seeds retained; simultaneous99% primary DKW bounds and independent80-digit finite-binomial CDF agree. Only corrected primary samplers and all supplemental fixtures are mandatory qualification gates; original-law failures are preserved diagnostic findings. Source-law qualification is not full spatial-trajectory equivalence. Original D1 labels are not rescored. All32 corrected trajectories require separate prospective registration before execution."]
    (BASE/"MG10D2_source_qualification.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k not in ("raw_checks","source_comparisons","forecast_by_class_seconds")},indent=2),flush=True)


if __name__=="__main__":main()
