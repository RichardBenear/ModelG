"""Independent MG030 records/fields audit. No production imports or ODE solves.
Run: python review/MG030_replication/audit.py CHECKOUT
"""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks

base=Path(sys.argv[1]); root=base/'results/MG9c'
read=lambda p: json.loads(p.read_text())
z=np.load(root/'snapshots.npz'); manifest=read(root/'snapshots_manifest.json')
lam=100/34; x0=140/9
assert set(z.files)==set(manifest)
for k in z.files:
    assert hashlib.sha256(np.ascontiguousarray(z[k],dtype=np.float64).tobytes()).hexdigest()==manifest[k]

def rhs(y,b,h):
    g,x,y=y.reshape(3,-1)
    lap=lambda f:(np.roll(f,1)-2*f+np.roll(f,-1))/h**2
    # Expanded reaction polynomial, independent of production concentration subtraction.
    nonlin=b*x+x0*x0*y+(b/x0)*x*x+2*x0*x*y+x*x*y
    return np.array([lap(g)-g+.1*x,lap(x)+g-x+nonlin,12*lap(y)-nonlin])

# FFT circular convolution instead of the production rolling-sum implementation.
def envelope(f,h):
    n=len(f); kernel=np.zeros(n); half=lam/(2*h)
    for j in range(-int(np.ceil(half))-1,int(np.ceil(half))+2):
        kernel[j%n]+=max(0,min(j+.5,half)-max(j-.5,-half))
    kernel/=kernel.sum(); mult=np.fft.fft(kernel)**2
    smooth=lambda a:np.fft.ifft(np.fft.fft(a)*mult)
    x=(np.arange(n)-n//2)*h
    return 2*abs(smooth((f-smooth(f).real)*np.exp(-2j*np.pi*x/lam)))

def field_edges(y,h,E):
    n=len(y)//3; f=y.reshape(3,n)[1]; e=envelope(f,h); m=e>=E/2
    changes=np.flatnonzero(m!=np.roll(m,1))
    assert len(changes)==2 and m[n//2] and not m[0]
    vals=[]
    for j in changes:
        frac=(E/2-e[j-1])/(e[j]-e[j-1])
        vals.append(((j-1)-n//2+frac)*h)
    return np.array(vals)

def cells(y,h,E):
    f=y.reshape(3,-1)[1]; n=len(f)
    peaks,_=find_peaks(-np.tile(f,3),prominence=E/2,distance=int(np.ceil(lam/(2*h))))
    return int(sum((peaks>=n)&(peaks<2*n)))

out={'reviewed_revision':'810ad932de904d633c1411c5a5c6593e8cc4e47a','snapshot_hashes':len(z.files),'runs':{}}
max_edge_err=0; max_res_err=0
for p in sorted((root/'runs').glob('*.json')):
    r=read(p)
    if 'track' not in r: continue
    tr=r['track']; t=np.array([a['t'] for a in tr]); xl=np.array([a['xL'] for a in tr]); xr=np.array([a['xR'] for a in tr])
    assert r['complete'] and r['success'] and r['t_reached']==r['T']==t[-1]
    assert np.array_equal(t,np.arange(0,r['T']+1,5))
    assert all(a['status']=='ok' and not a['breakup'] and not a['remote'] for a in tr)
    # All observed tracks have less than half a wavelength of drawdown after t=100;
    # hence no registered up/down reversal of that magnitude.
    for a in [-xl[t>=100],xr[t>=100]]:
        assert np.max(np.maximum.accumulate(a)-a)<lam/2
    windows=([(0,1500),(1500,3000),(3000,6000)] if r['kind']=='pert' else
             [(1800,3600),(6000,15000),(15000,30000)] if r['extended'] else
             [(1800,3600),(3600,6000),(6000,12000)])
    ds=np.array([[sign*(np.interp(b,t,a)-np.interp(c,t,a)) for c,b in windows]
                 for a,sign in [(xl,-1),(xr,1)]])
    speeds=ds/np.array([b-a for a,b in windows]); late=t>=t[-1]-1000
    v=np.array([np.polyfit(t[late],-xl[late],1)[0],np.polyfit(t[late],xr[late],1)[0]])
    counts=np.array([a['n'] for a in tr]); nr=int(np.ptp(counts[t>=windows[0][0]]))
    final=z[r['snapshots'][-1]['name']]; F=rhs(final,r['b'],r['h']); residual=float(abs(F).max())
    expected_times=([1500,3000,6000] if r['kind']=='pert' else
                    [1800,3600,6000,12000]+([15000,30000] if r['extended'] else []))
    assert [s['t'] for s in r['snapshots']]==expected_times
    max_res_err=max(max_res_err,abs(residual-r['residual']))
    decreasing=(np.all(np.diff(abs(speeds),axis=1)<0,axis=1))|(abs(ds[:,2])<=lam/400)
    checks=dict(displacement=bool(np.all(abs(ds[:,2])<=lam/40)),velocity=bool(np.all(abs(v)<=lam/40000)),
                decreasing=bool(np.all(decreasing)),residual=residual<=1e-5,cells=nr<=1)
    assert np.max(abs(ds[:,2]))<lam/4
    label='SETTLED' if all(checks.values()) else 'UNRESOLVED-SLOW'
    assert label==r['label']
    for j,side in enumerate(['L','R']):
        np.testing.assert_allclose(ds[j],r['info']['D'][side],atol=1e-12)
    for snap in r['snapshots']:
        assert snap['sha256']==manifest[snap['name']]
        arr=z[snap['name']]; edges=field_edges(arr,r['h'],r['E_pat'])
        idx=int(snap['t']/5)
        max_edge_err=max(max_edge_err,float(max(abs(edges-[xl[idx],xr[idx]]))))
        assert cells(arr,r['h'],r['E_pat'])==counts[idx]
    # For extended runs, reconstruct the residual at the original horizon too.
    residuals={str(int(s['t'])):float(abs(rhs(z[s['name']],r['b'],r['h'])).max()) for s in r['snapshots']}
    replay=None
    if r['b'] in [29.25,29.30,29.35] and p.stem.startswith('front_D02_'):
        old=read(base/'results/MG9/runs'/f"front_grid_{r['b']:.4f}.json")
        assert tr[:361]==old['track']; assert r['E_pat']==old['E_pat']
        replay=True
    width=xr-xl; center=(xr+xl)/2
    a,b=windows[-1]; mask=t>=a
    stats=dict(label=label,checks=checks,residual=residual,registered_residual=r['residual'],D3=ds[:,2].tolist(),
               velocities=v.tolist(),cell_range=nr,counts=sorted(set(counts.tolist())),residual_history=residuals,replay_exact=replay,
               width_change_final_window=float(width[-1]-np.interp(a,t,width)),
               center_shift_final_window=float(center[-1]-np.interp(a,t,center)),final_width=float(width[-1]))
    if r['kind']=='main':
        stats['width_mean_rate_6000_to_end']=float((width[-1]-np.interp(6000,t,width))/(t[-1]-6000))
        # A constant-width envelope is not sufficient to show rigid field translation.
        f=final.reshape(3,-1); k=2*np.pi*np.fft.fftfreq(r['N'],d=r['h'])
        grad=np.fft.ifft(np.fft.fft(f,axis=1)*(1j*k),axis=1).real
        drift=-float(np.sum(F*grad)/np.sum(grad*grad))
        remainder=F+drift*grad
        stats['descriptive_translation_projection']=dict(speed=drift,
            relative_L2_remainder=float(np.linalg.norm(remainder)/np.linalg.norm(F)))
    out['runs'][p.stem]=stats

assert max_edge_err<1e-10 and max_res_err<1e-10
out['max_snapshot_edge_error']=max_edge_err;out['max_final_residual_error']=max_res_err
prep={}
for key,h,expected in [('prep01_29.22000',.1,34),('prepL_29.22000',.2,68)]:
    r=read(root/'runs'/f'{key}.json'); arr=z[key]; E=float(envelope(arr.reshape(3,-1)[1],h).mean())
    n=cells(arr,h,E); res=float(abs(rhs(arr,29.22,h)).max())
    assert abs(res-r['residual'])<1e-10 and res<1e-6
    assert (n==expected)==(r['status']=='PERSISTS')
    prep[key]=dict(cells=n,residual=res,status=r['status'])
out['preparations']=prep
main=read(root/'runs/front_D02_29.22000.json'); ref=np.array([main['track'][-1][a] for a in ['xL','xR']])
out['perturbation_edge_shifts']={}
for name in ['p','m']:
    r=read(root/'runs'/f'pert_{name}05_29.22000.json')
    edges=np.array([r['track'][-1][a] for a in ['xL','xR']]);delta=edges-ref
    assert max(abs(delta))<lam/10 and r['label']=='SETTLED'
    out['perturbation_edge_shifts'][name]=delta.tolist()
wide=read(root/'runs/front_D02w_29.22000.json')
out['width_trial_edge_difference']=(np.array([wide['track'][-1][a] for a in ['xL','xR']])-ref).tolist()
lg=read(root/'ledger.json')['integrations'];assert len(lg)==len(set(v['key'] for v in lg))==11
out['core_hours']=sum(v['wall_s'] for v in lg)/3600
out['integrations']=len(lg)
out['Q1']='UNRESOLVED';out['Q2']='NOT ELIGIBLE'
ev=read(root/'evaluation.json');assert ev['Q1']==out['Q1'] and ev['Q2']['eligible'] is False
assert abs(ev['core_h']-out['core_hours'])<1e-12
Path(__file__).with_name('audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
