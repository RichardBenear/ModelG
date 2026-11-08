"""MG073 independent schedule, dense metric, paired and budget recomputation."""
import os
for key in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS","NUMEXPR_NUM_THREADS"):os.environ[key]="1"
from pathlib import Path
import json,hashlib,math,time,shutil
import numpy as np
from scipy.optimize import brentq
from scipy.special import stdtr
from check_mg10d1_development import check_raw
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG10E2/preparation'
def read(path):return json.loads(path.read_text(encoding='utf-8'))
def expected_activity(arm,time,switch):
    out=np.ones(500)
    if arm==0 or (arm in (2,3) and time>=switch):return out
    center=0 if arm==3 else 250
    for n in range(500):
        if min((n-center)%500,(center-n)%500)<=25:out[n]=0
    return out
def check_schedule(path):
    checked=check_raw(path)
    with np.load(path) as raw:
        data=raw['data'];times=raw['times'];arm=int(raw['arm']);switch=float(raw['switch']);initial=raw['initial']
        assert times[0]==0 and not np.any(data[0,0,7:]) and np.array_equal(data[0,0,:7],initial)
        assert raw['clamp_a'] and raw['clamp_b']
        # Only constrain complete observation intervals lying inside the off
        # period; a restoration-crossing interval can contain later source events.
        for i in range(1,len(times)):
            if times[i]<=switch or arm==1:
                alpha=expected_activity(arm,float(times[i-1]),switch);off=alpha==0
                assert not np.any(data[0,i,7,off]-data[0,i-1,7,off])
                assert not np.any(data[0,i,20,off]-data[0,i-1,20,off])
    return checked
def interval(values):
    n=len(values);mean=math.fsum(values)/n if n else None;ci=None
    if n>1:
        var=math.fsum((float(x)-mean)**2 for x in values)/(n-1)
        t=brentq(lambda x:stdtr(n-1,x)-.9875,0,100);half=t*math.sqrt(var/n);ci=[mean-half,mean+half]
    return dict(n=n,mean=mean,ci_nominal=ci)
def main():
    started=time.perf_counter();p=read(BASE/'preparation_plan.json')
    for name,sha in p['source_hashes'].items():
        path=ROOT/name;content=path.read_bytes() if name.endswith('.npz') else path.read_text(encoding='utf-8').replace('\r\n','\n').encode();assert hashlib.sha256(content).hexdigest()==sha,name
    for a in range(4):
        np.testing.assert_array_equal(expected_activity(a,0,5),p['activities_before'][str(a)])
        np.testing.assert_array_equal(expected_activity(a,5,5),p['activities_after'][str(a)])
    q=read(BASE/'qualification/qualification_results.json');assert q['passed'] and all(v['passed'] for v in q['tests']);checks=[]
    for v in q['raw_checks']:
        result=check_schedule(ROOT/v['file']);assert result['sha256']==v['sha256'];checks.append(result)
    with np.load(ROOT/'results/MG003/ref_profiles.npz') as raw:ref=raw['h0.2_b29.0'][1:]
    kernel=np.fft.ifft(np.exp(-.5*(.4*2*np.pi*np.fft.fftfreq(500,d=.2))**2)).real
    dense=kernel[(np.arange(500)[:,None]-np.arange(500)[None,:])%500]
    template=np.array([dense@(f-f.mean()) for f in ref]);norm=np.sum(template**2,axis=-1)
    with np.load(BASE/'qualification/metric_inputs.npz') as raw:fields=raw['fields']
    saved=read(BASE/'qualification/metric_results.json')['metrics']
    smooth=np.array([[dense@(f-f.mean()) for f in frame] for frame in fields]);fixed=np.sum(smooth*template[None],axis=-1)/norm
    np.testing.assert_allclose(fixed,saved['projection'],atol=1e-9);np.testing.assert_allclose(fields.mean(axis=-1),saved['field_mean'],atol=1e-11)
    np.testing.assert_allclose(np.sqrt(np.mean(smooth**2,axis=-1)),saved['smooth_rms'],atol=1e-10)
    for i,frame in enumerate(smooth):
        scores=np.array([float(frame[1]@np.roll(template[1],j)) for j in range(500)])
        best=next(j for j,x in enumerate(scores) if x>=scores.max()-1e-12*max(1.,abs(scores.max())))
        shifted=np.stack([np.roll(f,best) for f in template]);amp=np.sum(frame*shifted,axis=-1)/norm
        assert saved['x_phase_grid'][i]==(best*.2+50)%100-50
        np.testing.assert_allclose(amp,saved['x_aligned_projection'][i],atol=1e-9)
        np.testing.assert_allclose(np.sqrt(np.mean((frame-amp[:,None]*shifted)**2,axis=-1)),saved['x_aligned_residual_rms'][i],atol=1e-9)
    with np.load(BASE/'qualification/contrast_inputs.npz') as raw:linear=raw['linear'];late=raw['late'];shift=raw['shift']
    c=read(BASE/'qualification/contrast_results.json');avg=(linear[80]/2+linear[100]/2+linear[81:100].sum(axis=0))/20;np.testing.assert_allclose(avg,c['late_average'])
    for name,a,b in (('primary_restore_continuous',2,1),('restore_sham',2,0),('remote_restore_sham',3,0),('restore_remote',2,3)):
        value=late[1,a]-late[1,b]-late[0,a]+late[0,b];np.testing.assert_allclose(value,c['effects'][name],atol=1e-12)
        z=late+shift[None];np.testing.assert_allclose(value,z[1,a]-z[1,b]-z[0,a]+z[0,b],atol=1e-12)
    for mode in ('full','missing','prefix'):
        for k in range(2):
            values=[]
            for b in range(8):
                if mode!='full' and k==0 and b==0:continue
                z=lambda s,a:np.array(c['records'][f'a{a}_s{s}_k{k}_b{b}']['late_projection'])
                values.append(float((z(1,2)-z(1,1)-z(0,2)+z(0,1))[1]))
            expected=interval(values);actual=c[mode][str(k)];assert expected['n']==actual['n'];np.testing.assert_allclose(expected['mean'],actual['mean'],atol=1e-12);np.testing.assert_allclose(expected['ci_nominal'],actual['ci_nominal'],atol=1e-10)
    parents={};branches={};timing=[]
    for s in range(2):
        for k in range(2):
            times=[read(ROOT/f'results/MG10E1/parents/parent_s{s}_k{k}_b{b}.json')['simulation_seconds'] for b in range(8)]
            parents[f's{s}_k{k}']=math.fsum(times)/8*1.25+15
            with np.load(ROOT/f'results/MG10E1/preparation/devparents/parent_s{s}_k{k}.npz') as raw:state=raw['data'][0,-1,:7].copy();a0=raw['a0'].copy();b0=raw['b0'].copy()
            for a in range(4):
                meta=read(BASE/f'timing/timing_a{a}_s{s}_k{k}.json');assert meta['status']=='COMPLETE';path=BASE/'timing'/meta['raw'];checked=check_schedule(path);assert checked['sha256']==meta['sha256'];checks.append(checked)
                with np.load(path) as raw:
                    np.testing.assert_array_equal(raw['initial'],state);np.testing.assert_array_equal(raw['a0'],a0);np.testing.assert_array_equal(raw['b0'],b0);assert raw['switch']==.125 and raw['arm']==a
                    np.testing.assert_array_equal(raw['times'],[0.,.125,.25]);assert raw['seed']==meta['seed']
                branches[f'a{a}_s{s}_k{k}']=meta['compute_seconds']*100*1.25+15;timing.append(meta['compute_seconds'])
    ledger=read(BASE/'timing/ledger.json');assert len(ledger['records'])==16 and not ledger['stop_reason'] and all(v['status']=='COMPLETE' and v['returncode']==0 for v in ledger['records'])
    total=8*(math.fsum(parents.values())+math.fsum(branches.values()));allowance=ledger['charged_worker_seconds']+time.perf_counter()-started+5;assert allowance<900
    assert shutil.disk_usage(ROOT).free>=4*1024**3
    result=dict(passed=True,status='E2 PREPARATION PASS; PRODUCTION NOT REGISTERED',raw_checks=checks,synthetic_fields=len(fields),qualification_assertions=len(q['tests']),parent_forecast_seconds=parents,branch_forecast_seconds=branches,
      forecast_total_seconds=total,forecast_core_hours=total/3600,eight_worker_wall_hours=total/28800,within_75_percent_reserve=total<=108000,proposed_cap_core_hours=40,
      controlled_development_seconds=ledger['charged_worker_seconds'],total_preparation_allowance_seconds=allowance,free_bytes=shutil.disk_usage(ROOT).free,
      dependencies='Same author, shared constants/reference/NumPy/SciPy and frozen event kernel. Separate per-leap schedule oracle, explicit integer reconstruction, dense smoothing/direct all-shift scores, trapezoidal sums/manual contrasts and fsum/Student-CDF inversion; not external replication.')
    (BASE/'independent_check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('raw_checks','parent_forecast_seconds','branch_forecast_seconds')},indent=2))
if __name__=='__main__':main()
