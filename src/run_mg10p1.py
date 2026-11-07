"""MG058 prospective single-motif fixture registration and bounded execution."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import argparse,hashlib,json,platform,sys,time,traceback
from mg10d1_pool import run,atomic_json

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10P1"
FILES=("src/mg10p1.py","src/run_mg10p1.py","src/check_mg10p1.py","src/mg10d1_pool.py","src/runner.py","spec/MG10P1_proposal.md")


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def hashes():
    d={p:hashlib.sha256((ROOT/p).read_text(encoding="utf-8").replace("\r\n","\n").encode()).hexdigest() for p in FILES}
    d["results/MG003/ref_profiles.npz"]=digest(ROOT/"results/MG003/ref_profiles.npz")
    return d


def frame_matrix():
    from runner import seed_for
    cases=[]
    for position in (-49.8,-25.,-.1,0.,.1,12.4,49.8):cases.append(dict(kind="reference",position=position,background=140/9,amplitude=1.,expected="match"))
    for background in (12.,140/9,30.,1000.):cases.append(dict(kind="reference",position=0.,background=background,amplitude=1.,expected="match"))
    for amplitude in (.51,1.,1.49):cases.append(dict(kind="reference",position=0.,background=30.,amplitude=amplitude,expected="match"))
    for i in range(16):cases.append(dict(kind="noise_reference",position=(i-8)*.2,background=20.,amplitude=1.,sigma=(.15,.3)[i%2],seed=seed_for(10,i,0,5801,0),expected="match"))
    for background in (0.,1.,10.,140/9,1000.):cases.append(dict(kind="flat",background=background,expected="FLAT"))
    for i in range(16):cases.append(dict(kind="noise",background=20.,sigma=.3,seed=seed_for(10,i,1,5801,0),expected="not_match"))
    for mode in (1,8,32,100):
        for amplitude in (3.,10.):cases.append(dict(kind="wave",mode=mode,amplitude=amplitude,expected="not_match"))
    for width in (3.,6.,10.):cases.append(dict(kind="well",width=width,expected="not_match"))
    for positions in ((-15.,15.),(-30.,0.,30.)):cases.append(dict(kind="multiple",positions=positions,expected="not_match"))
    for amplitude in (.4,1.6):cases.append(dict(kind="reference",position=0.,background=30.,amplitude=amplitude,expected="not_match"))
    for amplitude in (.5,1.5):cases.append(dict(kind="reference",position=0.,background=30.,amplitude=amplitude,expected="BORDERLINE"))
    return [dict(c,id=f"f{i:03d}") for i,c in enumerate(cases)]


def sequence_matrix():
    cases=[dict(kind="stationary",end=30,x0=0.,status="CONTINUOUS",duration=30.,persistent=True,reason=None),
           dict(kind="moving",end=30,x0=0.,status="CONTINUOUS",duration=30.,persistent=True,reason=None),
           dict(kind="wrap",end=30,x0=49.,status="CONTINUOUS",duration=30.,persistent=True,reason=None),
           dict(kind="subcell",end=30,x0=0.,status="CONTINUOUS",duration=30.,persistent=True,reason=None),
           dict(kind="stationary",end=20,x0=0.,status="CONTINUOUS",duration=20.,persistent=True,reason=None),
           dict(kind="stationary",end=19,x0=0.,status="CONTINUOUS",duration=19.,persistent=False,reason=None)]
    for kind,reason,duration,persistent in (("permanent_gap","FLAT",9.,False),("temporary_gap","FLAT",9.,False),
        ("jump","MOTION_AMBIGUITY",9.,False),("sampling_gap","SAMPLING_GAP",2.,False),
        ("multiple_onset","AMBIGUOUS_MATCH",9.,False),("merger","AMBIGUOUS_MATCH",9.,False),("late_gap","FLAT",24.,True)):
        cases.append(dict(kind=kind,end=30,x0=0.,status="CENSORED",duration=duration,persistent=persistent,reason=reason))
    cases.append(dict(kind="wrong_initial",end=30,x0=0.,status="INITIAL_UNRESOLVED",duration=None,persistent=False,reason=None))
    cases.append(dict(kind="indistinguishable_replacement",end=30,x0=0.,status="CONTINUOUS",duration=30.,persistent=True,reason=None))
    return [dict(c,id=f"t{i:03d}") for i,c in enumerate(cases)]


def register():
    import numpy as np
    import mg10p1 as model
    file=BASE/"registration.json"
    if file.exists():raise ValueError("Already registered")
    raw,template=model.reference();x=(np.arange(500)-250)*.2
    atomic_json(file,dict(status="REGISTERED BEFORE FIXTURE GENERATION/EVALUATION",stage="MG058",authorization="proceed to next step after D2 closure; measurement qualification only",
        source_hashes=hashes(),frames=frame_matrix(),sequences=sequence_matrix(),invalid_checks=7,
        parameters=dict(N=500,h=.2,L=100.,sigma=.4,amp_min=.5,amp_max=1.5,error_max=.5,energy_min=.8,ambiguity_radius=2.,ambiguity_ratio=.9,border=1e-10,max_speed=1.5,max_gap=1.,persist=20.),
        static_reference=dict(rms=float(np.sqrt(np.mean(template**2))),energy_within10=float(np.sum(template[abs(x)<10]**2)/np.sum(template**2))),
        ceiling_seconds=300,fixture_timeout_seconds=240,oracle_timeout_seconds=60,workers=1,threads_per_library=1,
        environment=dict(python=sys.version,platform=platform.platform(),numpy=np.__version__),
        retention="All synthetic frame/sequence arrays, exact cases/seeds/expected labels, classifications, hashes, oracle and controller ledgers retained in results/MG10P1 and uploaded",
        scientific_scope="No stochastic trajectories or D1/D2 rescoring; sampled single-reference field continuity, not physical particle lineage/survival"))


def build_frame(case,profile):
    import numpy as np
    def shifted(position):return np.fft.ifft(np.fft.fft(profile)*np.exp(-2j*np.pi*np.fft.fftfreq(500)*position/.2)).real
    kind=case["kind"]
    if kind in ("reference","noise_reference"):
        f=case["background"]+case["amplitude"]*shifted(case["position"])
        if kind=="noise_reference":f+=np.random.default_rng(case["seed"]).normal(0,case["sigma"],500)
    elif kind=="flat":f=np.full(500,case["background"])
    elif kind=="noise":f=case["background"]+np.random.default_rng(case["seed"]).normal(0,case["sigma"],500)
    elif kind=="wave":f=20+case["amplitude"]*np.cos(2*np.pi*case["mode"]*np.arange(500)/500)
    elif kind=="well":f=20-8*np.exp(-.5*(((np.arange(500)-250)*.2)/case["width"])**2)
    else:f=30+sum((shifted(p) for p in case["positions"]),np.zeros(500))
    return np.asarray(f,dtype=float)


def build_sequence(case,profile):
    import numpy as np
    times=np.arange(case["end"]+1,dtype=float)
    if case["kind"]=="sampling_gap":times=times[times!=3]
    fields=[]
    for t in times:
        kind=case["kind"];position=t if kind=="moving" else 49+t if kind=="wrap" else .1*t if kind=="subcell" else 10. if kind=="wrong_initial" or kind=="jump" and t>=10 else 0.
        f=build_frame(dict(kind="reference",background=30.,amplitude=1.,position=position),profile)
        if kind=="permanent_gap" and t>=10 or kind=="temporary_gap" and t==10 or kind=="late_gap" and t>=25:f=np.full(500,30.)
        if kind=="multiple_onset" and t>=10:f=build_frame(dict(kind="multiple",positions=(0.,30.)),profile)
        if kind=="merger" and t>=10:
            gap=max(0.,10-2*(t-10));f=build_frame(dict(kind="multiple",positions=(-gap,gap)),profile)
        fields.append(f)
    return np.asarray(fields),times


def evaluate(directory):
    import numpy as np
    import mg10p1 as model
    registration=json.loads((BASE/"registration.json").read_text(encoding="utf-8"));profile,_=model.reference()
    frame_fields=np.array([build_frame(c,profile) for c in registration["frames"]])
    np.savez_compressed(directory/"frames.npz",fields=frame_fields)
    outputs=[];checks=[];manifest=[dict(path="frames.npz",sha256=digest(directory/"frames.npz"))]
    for case,field in zip(registration["frames"],frame_fields):
        result=model.assess(field);expected=case["expected"]
        ok=result["status"]=="REFERENCE_MATCH" and model.distance(result["position"],case["position"])<=.1+1e-8 if expected=="match" else result["status"]!="REFERENCE_MATCH" if expected=="not_match" else result["status"]==expected
        outputs.append(dict(case=case,result=result,passed=bool(ok)));checks.append(bool(ok))
    sequences=[]
    for case in registration["sequences"]:
        fields,times=build_sequence(case,profile);path=directory/(case["id"]+".npz");np.savez_compressed(path,fields=fields,times=times)
        manifest.append(dict(path=path.name,sha256=digest(path)))
        result=model.continuity(fields,times,case["x0"])
        ok=result["status"]==case["status"] and result["visible_duration"]==case["duration"] and result["persistent_reference_signal"]==case["persistent"] and result["censor_reason"]==case["reason"]
        sequences.append(dict(case=case,result=result,passed=bool(ok)));checks.append(bool(ok))
    valid=np.full(500,20.);nan=valid.copy();nan[0]=np.nan
    invalid=[lambda:model.assess(valid[:-1]),lambda:model.assess(nan),lambda:model.assess(-np.ones(500)),lambda:model.assess(np.ones((1,500))),
             lambda:model.continuity(np.array([valid,valid]),[0,0]),lambda:model.continuity(np.array([valid,valid]),[1,0]),lambda:model.continuity(np.array([valid]),[0],float("nan"))]
    invalid_results=[]
    for index,call in enumerate(invalid):
        try:call();ok=False
        except ValueError:ok=True
        invalid_results.append(dict(id=index,passed=ok));checks.append(ok)
    result=dict(status="FIXTURE CHECKS PASS" if all(checks) else "QUALIFICATION FAILED",passed=all(checks),frames=outputs,sequences=sequences,invalid=invalid_results,raw=manifest,checks=len(checks))
    atomic_json(directory/"evaluation.json",result)
    return result


def worker(job,directory):
    directory=Path(directory)
    if hashes()!=json.loads((BASE/"registration.json").read_text(encoding="utf-8"))["source_hashes"]:raise ValueError("Frozen source changed")
    if job["kind"]=="fixtures":
        result=evaluate(directory)
        atomic_json(directory/(job["id"]+".json"),dict(status="COMPLETE" if result["passed"] else "QUALIFICATION FAILED",job=job))
    else:
        from check_mg10p1 import audit
        audit();atomic_json(directory/(job["id"]+".json"),dict(status="COMPLETE",job=job))


def execute():
    directory=BASE/"fixtures"
    if (directory/"ledger.json").exists():raise ValueError("Already attempted; no automatic retry")
    first=run([dict(id="fixtures",kind="fixtures")],directory,Path(__file__).resolve(),ROOT,300,workers=1,timeout=240,poll_seconds=.25)
    if (directory/"evaluation.json").exists():
        oracle=run([dict(id="oracle",kind="oracle")],BASE/"oracle",Path(__file__).resolve(),ROOT,300,workers=1,prior=first["charged_worker_seconds"],timeout=60,poll_seconds=.25)
        atomic_json(BASE/"execution.json",dict(fixtures=first,oracle=oracle))
    else:atomic_json(BASE/"execution.json",dict(status="INCOMPLETE",fixtures=first))


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--register",action="store_true");parser.add_argument("--execute",action="store_true");parser.add_argument("--worker",nargs=2);args=parser.parse_args()
    if args.register:register()
    elif args.execute:execute()
    elif args.worker:
        j=json.loads(args.worker[0]);d=Path(args.worker[1])
        try:worker(j,d)
        except BaseException:atomic_json(d/(j["id"]+".json"),dict(status="FAILED",job=j,error=traceback.format_exc()));raise
    else:parser.error("Choose register/execute/worker")
