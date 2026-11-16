"""MG190 preparation only: spatial forcing adapter, backgrounds and null benchmark."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import sys,time,hashlib,traceback
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
from lattice import P17,lap_periodic,lattice_ode_1d
from runner import reference
from run_mg176 import ROOT
from mg10e2c1_pool import atomic_json
PUBLIC=ROOT/'results/MG190_gravity_preparation';RAW=ROOT/'local_raw/MG190_gravity_preparation'
EPS=(0.,.00025,.0005,.001,.0025,.005,.01,.02)
FILES=('src/gravity_mg190.py','spec/MG190_gravity_preparation.md','src/lattice.py','src/runner.py','results/MG003/ref_profiles.npz')
def binding():return {v:hashlib.sha256((ROOT/v).read_bytes()).hexdigest() for v in FILES}

class SpatialODE:
 def __init__(self,h,a,b):
  self.n=round(100/h);self.h=h;self.x=(np.arange(self.n)-self.n//2)*h;self.lap=lap_periodic(self.n,h);self.a=np.broadcast_to(a,(self.n,)).copy();self.b=np.broadcast_to(b,(self.n,)).copy();self.I=sp.identity(self.n,format='csc')
 def rhs(self,t,y):
  G,X,Y=y.reshape(3,self.n);f=X*X*Y
  return np.concatenate((self.lap@G-G+.1*X+self.a,self.lap@X+G-(1+self.b)*X+f,12*(self.lap@Y)+self.b*X-f))
 def jac(self,t,y):
  G,X,Y=y.reshape(3,self.n);L=self.lap;I=self.I
  return sp.bmat([[L-I,.1*I,None],[I,L+sp.diags(-1-self.b+2*X*Y),sp.diags(X*X)],[None,sp.diags(self.b-2*X*Y),12*L-sp.diags(X*X)]],format='csc')
 def steady(self):
  X=self.a/.9;y=np.array([X,X,self.b/X]).ravel()
  for i in range(15):
   error=float(np.max(abs(self.rhs(0,y))))
   if error<1e-9:return y,error,i
   change=spsolve(self.jac(0,y),-self.rhs(0,y));step=1.
   while np.min(y+step*change)<=0 or np.max(abs(self.rhs(0,y+step*change)))>=error:
    step/=2
    if step<1e-8:raise RuntimeError('Newton line search')
   y+=step*change
  raise RuntimeError('Background Newton convergence')
 def evolve(self,y,T):
  start=time.perf_counter();sol=solve_ivp(self.rhs,(0,T),y,method='BDF',jac=self.jac,rtol=1e-8,atol=1e-10,max_step=.5,t_eval=[0,T]);assert sol.success and sol.y.min()>0
  return sol.y,time.perf_counter()-start

def bc(a):
 X0=a/.9
 f=lambda z:(1+np.exp(z)-.1/(1+np.exp(z)))*(X0**2+12*np.exp(z))/(12*np.exp(z))
 return float(minimize_scalar(f,bounds=(-15,15),method='bounded').fun)

def register():
 PUBLIC.mkdir(exist_ok=False);atomic_json(PUBLIC/'development_plan.json',dict(source=binding(),seconds_cap=180,workers=1,no_gradient_particle_runs=True,eps=EPS,reference_PDF_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (RAW/'references').glob('*.pdf')}))

def run():
 import json
 plan=json.loads((PUBLIC/'development_plan.json').read_text());assert plan['source']==binding() and not(PUBLIC/'attempt.json').exists();start=time.perf_counter();atomic_json(PUBLIC/'attempt.json',dict(status='RUNNING'));rows=[]
 try:
  model=SpatialODE(.2,14.,29.);r=reference(29.,.2,100.);y=(r['S0'][:,None]+r['phi']).ravel();direction=np.sin(np.arange(len(y))*1.3);fd=(model.rhs(0,y+1e-6*direction)-model.rhs(0,y-1e-6*direction))/(2e-6);np.testing.assert_allclose(fd,model.jac(0,y)@direction,atol=5e-6,rtol=1e-7)
  old=lattice_ode_1d(P17,500,.2,2.,seed_scale=0.,y0=r['phi'].ravel(),t_eval=[0.,2.]);new,seconds=model.evolve(y,2.);error=float(abs(new[:,-1]-(old[2][:,-1]+np.repeat(r['S0'],500))).max());assert error<1e-7;bench=[]
  for h in (.2,.1):
   rr=reference(29.,h,100.);mm=SpatialODE(h,14.,29.);out,sec=mm.evolve((rr['S0'][:,None]+rr['phi']).ravel(),20.);np.savez_compressed(RAW/f'null_benchmark_h{h}.npz',fields=out);bench.append(dict(h=h,T=20,seconds=sec,T2000_linear_seconds=sec*100,disclosed_preparation_only=True))
  for arm in ('A','C'):
   for eps in EPS:
    x=model.x;aa=14*(1+eps*np.sin(2*np.pi*x/100)) if arm=='A' else np.full(500,14.);bb=29*(1+eps*np.sin(2*np.pi*x/100)) if arm=='C' else np.full(500,29.);mm=SpatialODE(.2,aa,bb);state,residual,iterations=mm.steady();G,X,Y=state.reshape(3,500);gradient=np.fft.ifft(1j*2*np.pi*np.fft.fftfreq(500,d=.2)*np.fft.fft(G)).real;proxy=[float(bb.min()+bc(14)-bc(float(aa.max()))),float(bb.max()+bc(14)-bc(float(aa.min())))];np.savez_compressed(RAW/f'background_{arm}_{eps}.npz',fields=state.reshape(3,500),a=aa,b=bb,G_gradient=gradient)
    rows.append(dict(arm=arm,eps=eps,residual_inf=residual,newton_iterations=iterations,G_half_range=float(np.ptp(G)/2),X_half_range=float(np.ptp(X)/2),maximum_G_gradient=float(abs(gradient).max()),b_equiv_linear_proxy_extrema=proxy,proxy_in_historical_window=proxy[0]>=28.93 and proxy[1]<=29.40,nonlinear_window_not_certified=True));assert time.perf_counter()-start<180
  atomic_json(PUBLIC/'preflight.json',dict(backgrounds=rows,null_benchmarks=bench,constant_solver_max_error=error,jacobian_finite_difference_passed=True,bc14=bc(14),cached_reference_G_integral=float(r['phi'][0].sum()*.2),no_production_gradient_runs=True))
  assert binding()==plan['source'];atomic_json(PUBLIC/'attempt.json',dict(status='COMPLETE',seconds=time.perf_counter()-start,no_production_gradient_runs=True))
 except BaseException:atomic_json(PUBLIC/'attempt.json',dict(status='FAILED',seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':register() if sys.argv[1]=='register' else run()
