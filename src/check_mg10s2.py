"""MG049 raw checking without new sampler/evaluator imports.

Shared: canonical spectral transform/rounding, original seed helper, NumPy.
Initial means, event reconstruction, population hashes, covariance enumeration,
CDF merging, confidence bands and decision labels are independently expressed.
"""
from pathlib import Path
import argparse,hashlib,json,math
import numpy as np
from mg10s2_measure import measure,NAMES
from runner import seed_for

ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/"results/MG10S2"
TAUS=(1/6000,1/12000,1/24000);OBS=(1/8000,1/4000,1/2000)


def fixture(case):
    n=500 if case<4 else 16
    omega=(100,1000)[case//2] if case<4 else (100,1000)[case-4]
    c=np.broadcast_to(np.array([140/9,140/9,29/(140/9)])[:,None],(3,n)).copy()
    if case<4 and case%2:
        with np.load(ROOT/"results/MG003/ref_profiles.npz") as raw:c+=raw["h0.2_b29.0"][1:]
    initial=np.zeros((7,n),np.int64);initial[:3]=np.rint(omega*c).astype(np.int64)
    initial[3]=round(omega*4200);initial[4]=round(omega*29*(140/9)*300)
    return initial,omega


def replay(initial,events):
    e=events.astype(np.int64);n=initial.shape[1]
    vectors=np.array([[1,0,0,-1,0,0,0],[-1,1,0,0,0,0,0],[1,-1,0,0,0,0,0],
                      [0,-1,1,0,-1,1,0],[0,1,-1,0,0,0,0],[0,-1,0,0,0,0,1],[-1,0,0,1,0,0,0]])
    result=np.broadcast_to(initial,(len(e),7,n)).copy()
    for channel in range(7):result+=e[:,channel,None,:]*vectors[channel][None,:,None]
    incoming_left=[(i+1)%n for i in range(n)];incoming_right=[(i-1)%n for i in range(n)]
    for s in range(3):
        left,right=e[:,7+2*s],e[:,8+2*s]
        result[:,s]+=left[:,incoming_left]+right[:,incoming_right]-left-right
    assert np.all(e>=0) and np.all(result>=0) and np.all(result.sum(axis=(1,2))==initial.sum())
    assert np.all(e[:,6]==0)
    return result


def mean_state(initial,omega,t):
    n=initial.shape[1];out=initial[:3].astype(float).copy()
    for i in range(n):
        g,x,y=[float(v) for v in initial[:3,i]]
        f=[14*omega,g,x/10,29*x,x*(x-1)*y/omega**2,9*x/10]
        increments=[f[0]-f[1]+f[2],f[1]-f[2]-f[3]+f[4]-f[5],f[3]-f[4]]
        for s,d in enumerate((25,25,300)):
            increments[s]+=d*(float(initial[s,(i-1)%n])+float(initial[s,(i+1)%n])-2*float(initial[s,i]))
        out[:,i]+=t*np.array(increments)
    return out


def covariance_enumeration(initial,omega,t,modes):
    n=initial.shape[1];results=[];instant_results=[]
    for mode in modes:
        phase=np.exp(-2j*np.pi*mode*np.arange(n)/n)/np.sqrt(n*omega*t)
        covariance=np.zeros((3,3),complex);instant=np.zeros((3,3),complex)
        for i in range(n):
            G,X,Y,A=initial[:4,i];p=phase[i]
            eye=np.eye(3);left,right=phase[(i-1)%n]-p,phase[(i+1)%n]-p
            groups=[(A,[14*omega*t/A],[eye[0]*p]),
                    (G,[t,25*t,25*t],[(eye[1]-eye[0])*p,eye[0]*left,eye[0]*right]),
                    (X,[t/10,29*t,9*t/10,25*t,25*t],[(eye[0]-eye[1])*p,(eye[2]-eye[1])*p,-eye[1]*p,eye[1]*left,eye[1]*right]),
                    (Y,[X*(X-1)/omega**2*t,300*t,300*t],[(eye[1]-eye[2])*p,eye[2]*left,eye[2]*right])]
            for count,probabilities,vectors in groups:
                mu=sum(prob*v for prob,v in zip(probabilities,vectors))
                second=sum(prob*np.outer(v,v.conj()) for prob,v in zip(probabilities,vectors))
                covariance+=count*(second-np.outer(mu,mu.conj()));instant+=count*second
        results.append(covariance);instant_results.append(instant)
    return np.array(results),np.array(instant_results)


def merged_distance(first,second):
    a=np.sort(first);b=np.sort(second);i=j=0;maximum=0.
    while i<len(a) or j<len(b):
        value=min(a[i] if i<len(a) else np.inf,b[j] if j<len(b) else np.inf)
        while i<len(a) and a[i]==value:i+=1
        while j<len(b) and b[j]==value:j+=1
        maximum=max(maximum,abs(i/len(a)-j/len(b)))
    return maximum


def synthetic_check():
    folder=BASE/"qualification";summary=json.loads((folder/"synthetic_summary.json").read_text())
    with np.load(folder/"synthetic_families.npz") as raw:
        for row in summary:
            key=f"M{row['M']}_{row['label']}";d=abs(raw[key+"_a"]-raw[key+"_b"])/16384
            w=math.sqrt(2*math.log(4*row["M"]/.005)/16384)
            low=np.maximum(0,d-w);high=np.minimum(1,d+w)
            coverage=float(np.mean(np.all((low<=row["true_distance"])&(high>=row["true_distance"]),axis=1)))
            qualified=float(np.mean(np.all(high<=.08,axis=1)))
            assert coverage==row["coverage"] and qualified==row["qualification_frequency"]
            assert coverage>=.98
            if row["true_distance"]==0:assert qualified>=.70
            if row["true_distance"]==.15:assert qualified==0


def audit(development):
    folder=BASE/("timing" if development else "production")
    ledger=json.loads((folder/"ledger.json").read_text())
    registration=None if development else json.loads((BASE/"registration.json").read_text())
    cost=ledger["prior_charge"]+sum(r["charged_worker_seconds"] for r in ledger["records"])
    assert abs(cost-ledger["charged_worker_seconds"])<1e-8
    arrays={};indices={};cov_cache={};shards=0;complete=0;invalid=False
    for record in ledger["records"]:
        job=record["job"]
        if record["status"]!="COMPLETE":invalid|=record["status"]=="QUALIFICATION FAILED";continue
        case,method=job["case"],job["method"];initial,omega=fixture(case)
        local=[];stage=4902 if development else 4903
        for part in record["shards"]:
            file=folder/part["file"]
            assert hashlib.sha256(file.read_bytes()).hexdigest()==part["sha256"]
            seed=seed_for(10,case,method,stage,part["chunk"])
            assert seed==part["seed"]
            with np.load(file) as raw:
                assert int(raw["seed"])==seed and float(raw["omega"])==omega
                np.testing.assert_array_equal(raw["initial"],initial)
                np.testing.assert_array_equal(raw["a0"],initial[3]);np.testing.assert_array_equal(raw["b0"],initial[4])
                modes=[32,250] if case<4 else [1,8];np.testing.assert_array_equal(raw["modes"],modes)
                if case<4:
                    t=TAUS[method//2];assert float(raw["time"])==t and raw["events"].dtype==np.uint32
                    pop=replay(initial,raw["events"])
                    sha=hashlib.sha256(np.asarray(pop,dtype="<i8").tobytes(order="C")).hexdigest()
                    assert sha==part["population_sha256"]
                    np.testing.assert_allclose(raw["mean"],mean_state(initial,omega,t),rtol=0,atol=2e-10)
                    values,_=measure(pop[:,:3],raw["mean"],omega,t,modes)
                    np.testing.assert_array_equal(values,raw["measurements"])
                    if method%2==0:
                        np.testing.assert_array_equal(raw["stats"][:,0],1);np.testing.assert_array_equal(raw["stats"][:,1],0)
                    key=case,method//2
                    if key not in cov_cache:cov_cache[key]=covariance_enumeration(initial,omega,t,modes)
                    q,q0=cov_cache[key]
                    np.testing.assert_allclose(raw["q_snapshot"],q,rtol=1e-10,atol=1e-7)
                    np.testing.assert_allclose(raw["q_instant"],q0,rtol=1e-10,atol=1e-7)
                else:
                    data=raw["data"];np.testing.assert_array_equal(raw["observations"],OBS)
                    checked=[]
                    for k,t in enumerate(OBS):
                        pop=replay(initial,data[:,k,7:20]);np.testing.assert_array_equal(pop,data[:,k,:7])
                        values,_=measure(pop[:,:3],initial[:3],omega,t,modes);checked.append(values)
                    values=np.stack(checked,axis=1)
                    np.testing.assert_array_equal(values,raw["measurements"])
                assert np.all(np.isfinite(values)) and len(values)==part["replicates"]
                assert raw["stats"][:,2].max()<=.25+1e-14
                local.append(values);shards+=1
            indices.setdefault((case,method),[]).append(part["chunk"])
        joined=np.concatenate(local);assert len(joined)==job["reps"]
        arrays.setdefault((case,method),[]).append(joined);complete+=1
        if registration is not None:
            assert job==next(j for j in registration["jobs"] if j["id"]==job["id"])
            assert [p["seed"] for p in record["shards"]]==registration["seeds"][job["id"]]
    groups={}
    for key,parts in arrays.items():
        count=64 if key[0]<4 else 8
        if not development:count=16384
        combined=np.concatenate(parts)
        if len(combined)==count:
            groups[key]=combined
            expected=list(range(1 if development else (32 if key[0]<4 else 64)))
            assert sorted(indices[key])==expected
    return groups,cost,shards,complete,invalid


def check(development=False):
    synthetic_check();groups,cost,shards,complete,invalid=audit(development)
    folder=BASE/("timing" if development else "production")
    if development:
        assert complete==32 and len(groups)==32
        result=dict(status="PASS",jobs=complete,shards=shards,integer_and_population_hashes=True,analytic_covariances=True,charged_worker_seconds=cost)
    else:
        original=json.loads((folder/"evaluation.json").read_text())
        rows={"snapshot":[],"short":[]}
        for family,m in (("snapshot",192),("short",96)):
            w=math.sqrt(2*math.log(4*m/.005)/16384)
            for row in original[family]:
                case=row["case"];mode_index=([32,250] if family=="snapshot" else [1,8]).index(row["mode"])
                observable=NAMES.index(row["observable"])
                if family=="snapshot":
                    a,b=groups[case,2*row["method"]],groups[case,2*row["method"]+1]
                    a,b=a[:,mode_index,observable],b[:,mode_index,observable]
                else:
                    time_index=list(OBS).index(row["time"])
                    a,b=groups[case,3][:,time_index,mode_index,observable],groups[case,0][:,time_index,mode_index,observable]
                d=merged_distance(a,b);low,high=max(0.,d-w),min(1.,d+w)
                assert abs(d-row["distance"])<1e-14 and abs(low-row["lower"])<1e-14 and abs(high-row["upper"])<1e-14
                rows[family].append(dict(**row))
        def label(family,m,positive):
            r=rows[family]
            if invalid:return "QUALIFICATION FAILED"
            if len(r)!=m:return "UNRESOLVED"
            if all(x["upper"]<=.08 for x in r):return positive
            if any(x["lower"]>.08 for x in r):return "NOT QUALIFIED ON TESTED CASES"
            return "UNRESOLVED"
        a=label("snapshot",192,"QUALIFIED ON TESTED SNAPSHOT LAWS");b=label("short",96,"QUALIFIED ON TESTED SHORT-TIME SPATIAL LAWS")
        if invalid:overall="QUALIFICATION FAILED"
        elif a=="QUALIFIED ON TESTED SNAPSHOT LAWS" and b=="QUALIFIED ON TESTED SHORT-TIME SPATIAL LAWS" and complete==88:overall="QUALIFIED WITHIN REGISTERED S2 SCOPE"
        elif "NOT QUALIFIED ON TESTED CASES" in (a,b):overall="NOT QUALIFIED ON TESTED CASES"
        else:overall="UNRESOLVED"
        assert (a,b,overall)==(original["snapshot_verdict"],original["short_verdict"],original["verdict"])
        result=dict(status="AGREES",verdict=overall,snapshot_verdict=a,short_verdict=b,jobs=complete,shards=shards,charged_worker_seconds=cost,
                    comparisons_checked=len(rows["snapshot"])+len(rows["short"]),
                    dependencies="Shared canonical spectral transform/rounding and original seed helper; independent reconstruction/means/covariance enumeration/CDF merging. No sampler/evaluator import. Same author, not external replication.")
    (folder/"recomputed.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(result)


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--development",action="store_true");check(p.parse_args().development)
