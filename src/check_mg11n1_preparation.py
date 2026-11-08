"""MG079 separate integer/graph reconstruction and generator-action references."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,hashlib,math,time,shutil
from decimal import Decimal,localcontext
import numpy as np
from scipy.sparse.linalg import expm_multiply
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11N1/preparation'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def check_raw(path):
    with np.load(path) as z:
        data=z['data'];initial=z['initial'];ef=z['environment'];n=initial.shape[1];src=z['src'];dst=z['dst'];weight=z['weight'];times=z['times']
        assert data.dtype==np.int64 and ef.dtype==np.int64 and data.shape[2:]==(20,n) and ef.shape==(len(data),len(times),2,len(src))
        assert np.all(data>=0) and np.all(ef>=0) and np.all(np.diff(ef,axis=1)>=0) and np.all(np.diff(data[:,:,7:],axis=1)>=0)
        expected=np.broadcast_to(initial,data[:,:,:7].shape).copy();rx=data[:,:,7:14]
        for a,b,c in ((3,0,0),(0,1,1),(1,0,2),(4,5,3),(2,1,4),(1,6,5),(0,3,6)):expected[:,:,a]-=rx[:,:,c];expected[:,:,b]+=rx[:,:,c]
        expected[:,:,1]-=rx[:,:,3];expected[:,:,2]+=rx[:,:,3]
        for i in range(n):
            for s in range(3):
                left=data[:,:,14+2*s,i];right=data[:,:,15+2*s,i];expected[:,:,s,i]-=left+right;expected[:,:,s,(i-1)%n]+=left;expected[:,:,s,(i+1)%n]+=right
        edges={}
        for edge,(i,j,w) in enumerate(zip(src,dst,weight)):
            assert i!=j and w>0 and (int(i),int(j)) not in edges;edges[int(i),int(j)]=float(w)
            for s in range(2):expected[:,:,3+s,i]-=ef[:,:,s,edge];expected[:,:,3+s,j]+=ef[:,:,s,edge]
        assert all(edges.get((j,i))==w for (i,j),w in edges.items())
        np.testing.assert_array_equal(data[:,:,:7],expected)
        assert np.all(data[:,:,:7].sum(axis=(2,3))==initial.sum())
        if times[0]==0:np.testing.assert_array_equal(data[:,0,:7],np.broadcast_to(initial,data[:,0,:7].shape));assert not np.any(data[:,0,7:]) and not np.any(ef[:,0])
        assert np.array_equal(z['stats'][:,0],data[:,-1,7:14].sum(axis=(1,2))+data[:,-1,14:20].sum(axis=(1,2))+ef[:,-1].sum(axis=(1,2)))
        for s in range(2):
            if z['kappa'][s]==0:assert not np.any(ef[:,:,s])
        for c in range(7):
            if z['mask'][c]==0:assert not np.any(rx[:,:,c])
        for s in range(3):
            if z['diff'][s]==0:assert not np.any(data[:,:,14+2*s:16+2*s])
        # Components constrain environment transfers themselves, not active hops.
        return dict(file=str(Path(path).relative_to(ROOT)),sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),replicas=len(data),checkpoints=len(times),integer_reconstruction=True,closed_total=True)
def categorical(data,kind):
    p=data[:,:,:7];n=p.shape[-1];expected=np.zeros_like(p)
    if kind=='transport':
        assert np.all(p[:,:,3].sum(axis=-1)==1) and np.all(p[:,:,4].sum(axis=-1)==1);expected[:,:,3:5]=p[:,:,3:5]
        label=p[:,:,3].argmax(axis=-1)*n+p[:,:,4].argmax(axis=-1)
    elif kind=='A_source':
        assert np.all(p[:,:,3].sum(axis=-1)+p[:,:,0].sum(axis=-1)==1);expected[:,:,0]=p[:,:,0];expected[:,:,3]=p[:,:,3]
        label=np.where(p[:,:,3].sum(axis=-1)==1,p[:,:,3].argmax(axis=-1),n+p[:,:,0].argmax(axis=-1))
    else:
        converted=data[:,:,10].sum(axis=-1);assert np.all((converted==0)|(converted==1)) and np.all(p[:,:,4].sum(axis=-1)+converted==1)
        expected[:,:,4]=p[:,:,4];expected[:,:,1,0]=1-converted;expected[:,:,2,0]=converted;expected[:,:,5,0]=converted
        label=np.where(converted==1,n,p[:,:,4].argmax(axis=-1))
    np.testing.assert_array_equal(p,expected);return label
def probability(raw,kind,t):
    n=raw['initial'].shape[-1];size=n*n if kind=='transport' else 2*n if kind=='A_source' else n+1;Q=np.zeros((size,size));p=np.zeros(size)
    def edge(a,b,rate):Q[a,b]+=rate;Q[a,a]-=rate
    if kind=='transport':
        ia=int(raw['initial'][3].argmax());ib=int(raw['initial'][4].argmax());p[n*ia+ib]=1
        for a in range(n):
            for b in range(n):
                for i,j,w in zip(raw['src'],raw['dst'],raw['weight']):
                    if i==a:edge(n*a+b,n*j+b,raw['kappa'][0]*w)
                    if i==b:edge(n*a+b,n*a+j,raw['kappa'][1]*w)
    else:
        species=3 if kind=='A_source' else 4;p[int(raw['initial'][species].argmax())]=1
        for i,j,w in zip(raw['src'],raw['dst'],raw['weight']):edge(i,j,raw['kappa'][species-3]*w)
        for i in range(n):
            rate=14*float(raw['omega'])/raw['a0'][i] if kind=='A_source' else 29*raw['initial'][1,i]/raw['b0'][i]
            edge(i,n+i if kind=='A_source' else n,rate)
    result=expm_multiply(Q.T*t,p);assert abs(result.sum()-1)<1e-10 and result.min()>-1e-12;return result
def main():
    started=time.perf_counter();plan=read(BASE/'preparation_plan.json')
    for name,sha in plan['source_hashes'].items():
        p=ROOT/name;content=p.read_bytes() if name.endswith('.npz') else p.read_text(encoding='utf-8').replace('\r\n','\n').encode();assert hashlib.sha256(content).hexdigest()==sha
    q=read(BASE/'qualification/qualification_results.json');assert q['passed'] and all(t['passed'] for t in q['tests']);checks=[]
    for v in q['raw_checks']:
        x=check_raw(ROOT/v['file']);assert x['sha256']==v['sha256'];checks.append(x)
    with np.load(BASE/'qualification/continuum_inputs.npz') as z:
        state=z['state'];G,X,Y,A,B=state[:5];a0=float(z['a0']);b0=float(z['b0']);ref=np.zeros_like(state)
        for i in range(3):
            f1=14*A[i]/a0;f2=G[i];f2r=.1*X[i];f3=29*X[i]*B[i]/b0;f4=X[i]**2*Y[i];f5=.9*X[i]
            ref[:,i]=[f1-f2+f2r,f2-f2r-f3+f4-f5,f3-f4,-f1,-f3,f3,f5]
            for s,d in enumerate((25.,25.,300.)):ref[s,i]+=d*(state[s,(i-1)%3]+state[s,(i+1)%3]-2*state[s,i])
        np.testing.assert_allclose(ref,z['actual'],rtol=1e-12,atol=1e-12)
    limits=read(BASE/'qualification/linear_limits.json')
    # Independently use generator-action evolution instead of dense matrix exp.
    for record in limits['mixing']:
        initial=np.zeros((7,3),dtype=np.int64);initial[3,0]=1
        raw=dict(initial=initial,src=np.array([0,1,1,2]),dst=np.array([1,0,2,1]),weight=np.array([1.,1.,.5,.5]),kappa=np.array([record['kappa'],0.]),omega=1.,a0=np.array([20.,28.,35.]),b0=np.ones(3))
        law=probability(raw,'A_source',1.);np.testing.assert_allclose(law,record['law'],atol=1e-11)
    with localcontext() as ctx:
        ctx.prec=50
        for record in limits['capacity']:
            c=Decimal(str(record['capacity']));value=14*c*(1-(-1/c).exp());assert abs(float(value)-record['output'])<1e-11
    np.testing.assert_allclose(expm_multiply(1000*np.array(limits['disconnected_generator']).T,np.array([1.,0,0,0])),limits['disconnected_law'],atol=1e-11)
    timings={};descriptive=[]
    for cell in range(4):
        for method in range(2 if cell==3 else 1):
            m=read(BASE/f'timing/timing_c{cell}_m{method}.json');assert m['status']=='COMPLETE';path=BASE/'timing'/m['raw'];v=check_raw(path);assert v['sha256']==m['sha256'];checks.append(v)
            checks.append(check_raw(BASE/f'timing/timing_c{cell}_m{method}_warm.npz'));timings[cell,method]=m['compute_seconds']/m['reps']
            if cell<3:
                kind=('transport','A_source','B_contact')[cell]
                with np.load(path) as raw:
                    labels=categorical(raw['data'],kind)
                    for i,t in enumerate(raw['times'][1:],1):
                        p=probability(raw,kind,float(t));cdf=np.cumsum(np.bincount(labels[:,i],minlength=len(p)))/len(labels);descriptive.append(dict(cell=cell,time=float(t),distance=float(abs(cdf-np.cumsum(p)).max()),scope='Development descriptive only; not scored using production criteria'))
    forecasts={j['id']:timings[j['cell'],j['engine']]*4096*2+10 for j in plan['candidate_jobs']};total=math.fsum(forecasts.values())
    ledger=read(BASE/'timing/ledger.json');assert not ledger['stop_reason'] and len(ledger['records'])==5 and all(r['status']=='COMPLETE' and r['returncode']==0 for r in ledger['records'])
    allowance=ledger['charged_worker_seconds']+time.perf_counter()-started+5;assert allowance<900 and shutil.disk_usage(ROOT).free>=2*1024**3
    # Array bound for data20N, edge history2E, stats and initials/config separate.
    raw_bytes=0
    from mg11n1_cases import cell as factory
    for j in plan['candidate_jobs']:
        case,g,k,t,m=factory(j['cell']);raw_bytes+=j['reps']*len(t)*(20*g['n']+2*len(g['src']))*8
    result=dict(passed=True,status='N1 DEVELOPMENT PASS; PRODUCTION NOT REGISTERED',raw_checks=checks,formal_assertions=len(q['tests']),development_cdf_descriptive=descriptive,
       forecast_seconds_per_job=forecasts,forecast_worker_seconds=total,forecast_core_hours=total/3600,eight_worker_wall_minutes=total/480,within_75_percent_reserve=total<=2700,candidate_ceiling_core_hours=1,
       raw_array_bound_bytes=raw_bytes,controlled_development_seconds=ledger['charged_worker_seconds'],total_preparation_allowance_seconds=allowance,
       dependencies='Same author; shared laws/constants/NumPy/SciPy and qualified legacy disconnected clock. Separate explicit reaction/edge reconstruction, joint-state generators with expm_multiply, scalar continuum arithmetic and50-digit capacity reference. No external replication or nonlinear particle/full-lattice stochastic qualification.')
    (BASE/'independent_check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('raw_checks','forecast_seconds_per_job','development_cdf_descriptive')},indent=2))
if __name__=='__main__':main()
