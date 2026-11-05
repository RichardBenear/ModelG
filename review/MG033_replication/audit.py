"""MG9d result audit: no production imports and no integration.
Usage: python audit.py /path/to/ModelG
"""
import json, sys, hashlib
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.signal import find_peaks, resample

root=Path(sys.argv[1]); p=root/'results/MG9d'; lam=100/34
ev=json.loads((p/'evaluation.json').read_text())
ledger=json.loads((p/'ledger.json').read_text())
manifest=json.loads((p/'snapshots_manifest.json').read_text())
z=np.load(p/'snapshots.npz')
keys=['R0','R1','R3','R2a','R2b','R2c','R4']
records={k:json.loads((p/'runs'/f'{k}.json').read_text()) for k in keys}
out={'revision':'d3bcfd35ad5b85a88ac59a6cdfd1e4cd78f06d98','runs':{}}
errors={k:[] for k in ['d','S','S_ref','C_edge','W','C_amp','C_phi2','minima','residual']}
def sha(y): return hashlib.sha256(np.ascontiguousarray(y,dtype=np.float64).tobytes()).hexdigest()
def shift(f,h,s):
    return np.fft.ifft(np.fft.fft(f)*np.exp(-2j*np.pi*np.fft.fftfreq(len(f),h)*s)).real
def rhs(y,h):
    g,x,y=y.reshape(3,-1); x0=140/9
    lap=lambda a:(np.roll(a,1)-2*a+np.roll(a,-1))/h**2
    q=29.3*x+x0*x0*y+29.3/x0*x*x+2*x0*x*y+x*x*y
    return np.concatenate([lap(g)-g+.1*x,lap(x)+g-x+q,12*lap(y)-q])
assert set(z.files)==set(manifest) and len(z.files)==77
for key in z.files: assert sha(z[key])==manifest[key]
out['verified_snapshot_hashes']=77
with np.load(root/'results/MG9c/snapshots.npz') as old:
    src=old['front_D02_29.30000_t30000'].reshape(3,500)
alt=(-1.)**np.arange(500)
y0=src-(src@alt/500)[:,None]*alt
starts={'R0':y0,'R1':y0,'R3':y0[:,(-np.arange(500))%500],
        'R4':resample(y0,1000,axis=1)}
for k,s in [('R2a',.05),('R2b',.1),('R2c',.15)]: starts[k]=np.array([shift(f,.2,s) for f in y0])
out['max_start_error']=max(float(abs(z[k+'_t0']-a.ravel()).max()) for k,a in starts.items())
assert out['max_start_error']<1e-12

for key,r in records.items():
    n,h=r['N'],r['h']; x=(np.arange(n)-n//2)*h; Epat=r['E_pat']
    assert r['complete'] and r['success'] and r['t_reached']==10000 and r['stop_reason'] is None
    assert r['seed_scale']==0 and r['b']==29.3 and Epat==6.60086115019687
    assert r['solver']==dict(rtol=1e-10 if key=='R1' else 1e-8,atol=1e-12 if key=='R1' else 1e-10,max_step=.1 if key=='R1' else .5)
    assert [s['t'] for s in r['diag']]==list(range(0,10001,50))
    assert [s['t'] for s in r['track']]==list(range(0,10001,5))
    assert len(r['segments'])==10 and all(s['complete'] and len(s['attempts'])==1 and s['attempts'][0]['success'] for s in r['segments'])
    assert sha(z[key+'_t0'])==r['init_sha256']
    assert all(sha(z[s['name']])==s['sha256'] for s in r['snapshots'])
    kernel=np.zeros(n); half=lam/(2*h)
    for j in range(-int(np.ceil(half))-1,int(np.ceil(half))+2):
        kernel[j%n]=max(0,min(j+.5,half)-max(j-.5,-half))
    kernel/=kernel.sum(); mult=np.fft.fft(kernel)**2
    smooth=lambda a:np.fft.ifft(np.fft.fft(a)*mult)
    def measures(f):
        E=2*abs(smooth((f-smooth(f).real)*np.exp(-2j*np.pi*x/lam)))
        above=E>=Epat/2; cross=np.flatnonzero(above[:-1]!=above[1:])
        assert len(cross)==2
        edges=np.array([x[j]+h*(Epat/2-E[j])/(E[j+1]-E[j]) for j in cross])
        pk,_=find_peaks(-np.tile(f,3),prominence=Epat/2,distance=int(np.ceil(lam/(2*h))))
        pk=pk[(pk>=n)&(pk<2*n)]-n
        a,b,c=f[(pk-1)%n],f[pk],f[(pk+1)%n]
        mins=x[pk]+h*(a-c)/(2*(a-2*b+c))
        return dict(C_edge=float(edges.mean()),W=float(np.diff(edges)[0]),
                    C_amp=float(x@(E*E)/sum(E*E)),C_phi2=float(x@(f*f)/sum(f*f)),minima=mins)
    f0=z[key+'_t0'][n:2*n]; ref=z[key+'_t2000'][n:2*n]
    m0=measures(f0); F=np.fft.fft(f0); freq=2*np.pi*np.fft.fftfreq(n,h)
    calc={}
    for t in range(0,10001,1000):
        f=z[f'{key}_t{t}'][n:2*n]
        def grad(d):
            A=F*np.exp(-1j*freq*d)
            return float(2*np.dot(np.fft.ifft(A).real-f,np.fft.ifft(-1j*freq*A).real))
        d=brentq(grad,-lam/4,lam/4,xtol=1e-14)
        q=measures(f); q['d']=d
        q['S']=float(np.linalg.norm(f-shift(f0,h,d))/np.linalg.norm(f0))
        q['minima']=np.array([q['minima'][np.argmin(abs(q['minima']-a))] for a in m0['minima'][1:-1]])
        for field in ('C_edge','C_amp','C_phi2'): q[field]-=m0[field]
        calc[t]=q
        if t>=2000:q['S_ref']=float(np.linalg.norm(f-shift(ref,h,d-calc[2000]['d']))/np.linalg.norm(ref))
        stored=r['diag'][t//50]
        for field in errors:
            if field in q:
                err=float(np.max(abs(np.asarray(q[field])-np.asarray(stored[field])))); errors[field].append(err)
                assert err<(2e-8 if field in ('d','S','S_ref') else 1e-10),(key,t,field,err)
    res=float(abs(rhs(z[key+'_t10000'],h)).max()); errors['residual'].append(abs(res-r['residual']))
    assert abs(res-r['residual'])<1e-10
    # Re-evaluate all six validity rules from the full saved diagnostic/track records.
    ds=r['diag']; ts=r['track']; failed=[]
    if not r['complete'] or len(ds)!=201:failed.append(1)
    if any(s['status']!='ok' or s['breakup'] or s['remote'] for s in ts):failed.append(2)
    if any(abs(s['n']-ts[0]['n'])>1 for s in ts) or abs(ds[-1]['W']-ds[0]['W'])>=lam/2:failed.append(3)
    if any(s['ambiguous'] for s in ds):failed.append(4)
    if ds[-1]['S_ref']>1e-3:failed.append(5)
    keep=[i for i in range(len(r['ref_minima'])) if all(s['minima'][i] is not None for s in ds)]
    dd=ds[-1]['d']-ds[40]['d']; band=.1*abs(dd)+5e-4
    changes={f:ds[-1][f]-ds[40][f] for f in ('C_edge','C_amp','C_phi2')}
    changes['m']=float(np.mean([ds[-1]['minima'][i]-ds[40]['minima'][i] for i in keep]))
    if any(abs(v-dd)>band for v in changes.values()):failed.append(6)
    assert failed==ev['validity'][key]['failed']
    speed={k:(calc[b]['d']-calc[a]['d'])/(b-a) for k,a,b in [('a',2000,6000),('b',6000,10000),('v',2000,10000)]}
    assert max(abs(speed[k]-ev['speeds_raw'][key][k]) for k in speed)<5e-13
    out['runs'][key]=dict(speeds=speed,failed=failed,final_residual=res,S_ref=calc[10000]['S_ref'],
        width_change=calc[10000]['W']-calc[0]['W'],minima_matched=len(keep),consistency_changes=changes,
        window_widths={t:calc[t]['W'] for t in (2000,6000,10000)})

# Independently check the registered rule branches reached by these actual data.
v={k:r['speeds'] for k,r in out['runs'].items()}; vf=5e-9
assert all(not out['runs'][k]['failed'] for k in keys if k!='R4')
assert out['runs']['R4']['failed']==[5]
assert all(abs(v['R1'][k]-v['R0'][k])<=.05*abs(v['R0'][k])+vf for k in ('a','b','v'))
assert abs(-v['R3']['v']-v['R0']['v'])<=.01*abs(v['R0']['v'])+vf
assert all(abs(v[k]['v']-v['R0']['v'])>.20*abs(v['R0']['v'])+3*vf for k in ('R2a','R2b','R2c'))
assert all(abs(v[k][f])>=3*vf for k in ('R0','R1') for f in ('a','b','v'))
assert all(0<v[k]['b']/v[k]['a']<=.9 for k in ('R0','R1'))
labels={'C1':'NUMERICS-AGREE','C3':'SYMMETRIC','C2':'ALIGNMENT-DEPENDENT','C4':'DECREASING','C5':'R4 INVALID / LOST-IDENTITY'}
assert labels==ev['components'] and ev['reading']=='MESH EFFECT SUPPORTED'
out.update(components=labels,reading='MESH EFFECT SUPPORTED',baseline_speed_ratio=v['R0']['b']/v['R0']['a'],
    max_numerics_speed_difference=max(abs(v['R0'][f]-v['R1'][f]) for f in ('a','b','v')),
    max_absolute_errors={k:max(a) for k,a in errors.items()})
assert [a['key'] for a in ledger['integrations']]==keys
assert all(a['status']=='done' and a['wall_s']==records[a['key']]['wall_s'] for a in ledger['integrations'])
out['integrations']=len(ledger['integrations']);out['core_h']=sum(a['wall_s'] for a in ledger['integrations'])/3600
assert out['core_h']==ev['core_h'] and out['integrations']==7
out['scope']='All 77 saved fields and hashes; all 1407 diagnostic and 14007 tracker records. No production imports or rerun trajectories. Fields exist only every 1000 units.'
Path(__file__).with_name('audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
