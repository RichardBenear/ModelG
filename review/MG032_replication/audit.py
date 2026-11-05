"""Static MG9d implementation audit; never integrates an MG9d state.
Usage: python audit.py /path/to/ModelG
"""
import hashlib
import json
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
from scipy.signal import resample
from scipy.optimize import brentq

repo = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(repo / 'src'))
import run_mg9d as D

def shift(f, h, d):
    k = 2*np.pi*np.fft.fftfreq(len(f), h)
    return np.fft.ifft(np.fft.fft(f)*np.exp(-1j*k*d)).real

def rhs(y, h):
    g, x, z = y.reshape(3, -1)
    lap = lambda f: (np.roll(f, 1)-2*f+np.roll(f, -1))/h**2
    x0 = 140/9
    reaction = D.B*x + x0*x0*z + D.B/x0*x*x + 2*x0*x*z + x*x*z
    return np.concatenate([lap(g)-g+.1*x, lap(x)+g-x+reaction, 12*lap(z)-reaction])

with np.load(repo/'results/MG9c/snapshots.npz') as z:
    src = z[D.SRC_KEY].copy()
a = src.reshape(3, 500)
# Explicitly subtract the alternating mode, independently of the production rFFT filter.
alternating = (-1.)**np.arange(500)
y0 = (a - (a@alternating/500)[:, None]*alternating).ravel()
states = {'source': (src, .2), 'R0': (y0, .2), 'R1': (y0, .2),
          'R3': (y0.reshape(3, 500)[:, (-np.arange(500))%500].ravel(), .2),
          'R4': (resample(y0.reshape(3, 500), 1000, axis=1).ravel(), .1)}
for name, delta in [('R2a', .05), ('R2b', .1), ('R2c', .15)]:
    states[name] = (np.concatenate([shift(f, .2, delta) for f in y0.reshape(3,500)]), .2)
prod = D.start_states()
stored = json.loads((repo/'results/MG9d/transfer.json').read_text())
out = {'source_sha256': hashlib.sha256(src.tobytes()).hexdigest(),
       'y0_max_change': float(abs(y0-src).max()),
       'state_max_difference': {k: float(abs(y-prod[k]['y']).max()) for k,(y,h) in states.items() if k!='source'},
       'rhs_max': {k: float(abs(rhs(y,h)).max()) for k,(y,h) in states.items()}}
out['rhs_max_difference_from_stored'] = max(abs(v-stored['rhs_max'][k]) for k,v in out['rhs_max'].items())
assert max(out['state_max_difference'].values()) < 1e-12
assert out['rhs_max_difference_from_stored'] < 1e-10

# Analytic Jacobian versus independent centered directional finite differences.
rng = np.random.default_rng(32)
jac_checks = {}
for name in ('R0','R4'):
    y,h=states[name]; n=len(y)//3
    v=rng.normal(size=len(y)); v/=np.linalg.norm(v)
    _,f,j=D.make_rhs_jac(D.B,n,h)
    eps=1e-5
    fd=(rhs(y+eps*v,h)-rhs(y-eps*v,h))/(2*eps)
    jv=j(0,y)@v
    jac_checks[name]=float(np.linalg.norm(fd-jv)/np.linalg.norm(jv))
    assert jac_checks[name]<1e-7
out['jacobian_relative_directional_error']=jac_checks

# Rerun the production static calibration, then independently recover shifts by
# solving the analytic derivative of the FULL-FFT least-squares objective.
fl=D.floor_validation()
out['production_floor_rerun']={k:fl[k] for k in ('delta','v_f','max_err_per_measure','resolved_threshold')}
cases=[]
for name in ('R0','R4'):
    y,h=states[name]; n=len(y)//3; f=y[n:2*n]
    k=2*np.pi*np.fft.fftfreq(n,h); F=np.fft.fft(f)
    for delta in (.001,.004,.01,.05,.1,.15):
        for s in (delta,-delta):
            target=shift(f,h,s)
            def gradient(d):
                phase=F*np.exp(-1j*k*d)
                moved=np.fft.ifft(phase).real
                deriv=np.fft.ifft(-1j*k*phase).real
                return float(2*np.dot(moved-target,deriv))
            recovered=brentq(gradient,-.2,.2,xtol=1e-14)
            measured=D.recover(f,target,h,D.e_pat(),n)
            cases.append(dict(grid=name,true=s,independent=recovered,production=measured['d'],
                              ambiguous=measured['ambiguous'],matched=measured['matched']))
out['independent_alignment_max_error']=max(abs(c['independent']-c['true']) for c in cases)
out['independently_generated_shift_production_max_error']=max(abs(c['production']-c['true']) for c in cases)
out['shift_cases']=cases
assert out['independent_alignment_max_error']<1e-12
assert out['independently_generated_shift_production_max_error']<2e-9
assert all(not c['ambiguous'] and c['matched']==16 for c in cases)
assert fl['v_f']==5e-9

# Failure probes are mocks, not numerical trajectories. Exercise both actual
# solve_ivp empty-list output shape and a normalized empty-array equivalent.
failure=[]
for form in ('list','array'):
    empty_t=[] if form=='list' else np.array([])
    empty_y=[] if form=='list' else np.empty((1500,0))
    with patch.object(D,'solve_ivp',return_value=SimpleNamespace(success=False,t=empty_t,y=empty_y)) as mock:
        try:
            r=D.mg9d_ode(D.B,500,.2,1000,y0,seed_scale=0,t_eval=np.array([5.,10.]))
            result={'returned_complete':r['complete']}
        except Exception as e:
            result={'exception':type(e).__name__,'message':str(e)}
        failure.append(dict(form=form,attempts=mock.call_count,**result))
out['empty_solver_failure']=failure

# An interrupted attempt is currently treated as zero elapsed time.
lg={'integrations':[dict(key='R0',wall_s=None,started='2026-11-01T00:00:00')],
    'cap':10,'cutoff_core_s':3600}
with patch.object(D,'ledger',return_value=lg), patch.object(D,'_ledger_write'):
    out['unresolved_prior_attempt_allows_new_start']=D._ledger_start('R0')==1

# A nearly exhausted budget permits an arbitrarily long next run; there is no
# intra-run time check. This probe tests the entry guard only.
lg={'integrations':[dict(key='R0',wall_s=3599.)],'cap':10,'cutoff_core_s':3600}
with patch.object(D,'ledger',return_value=lg), patch.object(D,'_ledger_write'):
    out['3599_seconds_allows_new_start']=D._ledger_start('R1')==1

out['platform']=D.platform()
out['scope']='Static arrays, Jacobian differences, calibration, and mocked failure paths; no MG9d integration.'
dest=Path(__file__).with_name('audit.json')
dest.write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in out.items() if k!='shift_cases'},indent=2))
