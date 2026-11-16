"""MG191 deterministic fields, core measurements and numerical stability checks."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse.linalg import eigs
from gravity_mg190 import SpatialODE
from runner import reference
from measure import fourier_shift
EPS=(.00025,.0005,.001);FINE=(.0005,.001)

def wrap(x):return (x+50)%100-50
def phase(j):return float(wrap(100*(j/8+1/32)))

def model_for(job):
 h=job['h'];x=(np.arange(round(100/h))-round(100/h)//2)*h;arm=job['arm'];epsilon=job.get('eps',0.)
 a=np.full(len(x),14.);b=np.full(len(x),29.)
 if arm in ('A','R','H'):a*=1+epsilon*np.sin(2*np.pi*x/100)
 elif arm=='C':b*=1+epsilon*np.sin(2*np.pi*x/100)
 elif arm in ('L','qa'):a*=1+epsilon
 elif arm=='qb':b*=1+epsilon
 elif arm!='N':raise ValueError(arm)
 return SpatialODE(h,a,b)

def core(field,background,h,reference_depth):
 v=np.asarray(field)-np.asarray(background);n=len(v);candidates=np.flatnonzero((v<np.roll(v,1))&(v<=np.roll(v,-1))&(v<-.5*reference_depth));count=len(candidates)
 if count!=1:return np.array([np.nan,np.nan,np.nan,float(count)])
 i=int(candidates[0]);left=i;right=i
 while left>i-n//2 and v[left%n]<=0:left-=1
 while right<i+n//2 and v[right%n]<=0:right+=1
 if v[left%n]<=0 or v[right%n]<=0:return np.array([np.nan,np.nan,np.nan,float(count)])
 l=(left-n//2-v[left%n]/(v[(left+1)%n]-v[left%n]))*h
 rr=(right-1-n//2-v[(right-1)%n]/(v[right%n]-v[(right-1)%n]))*h
 return np.array([wrap((l+rr)/2),-float(v[i]),rr-l,float(count)])

def reference_data(h):
 r=reference(29.,h,100.);c=core(r['phi'][1],np.zeros(len(r['x'])),h,abs(r['A']));assert np.isfinite(c[:3]).all()
 return r,dict(depth=float(c[1]),width=float(c[2]),wavelength=float(r['lam']))

def initialize(job,background):
 r,_=reference_data(job['h']);position=job.get('position',0.);return background.reshape(3,-1)+fourier_shift(r['phi'],job['h'],-position)

def evolve_chunk(model,initial,T,rtol=1e-8,atol=1e-10):
 samples=np.arange(round(T*2)+1)*.5;solution=solve_ivp(model.rhs,(0,T),np.asarray(initial).ravel(),method='BDF',jac=model.jac,rtol=rtol,atol=atol,max_step=.5,t_eval=samples)
 if not solution.success or solution.y.shape!=(3*model.n,len(samples)) or not np.isfinite(solution.y).all() or solution.y.min()<=0:raise RuntimeError('BDF failed/nonpositive/incomplete: '+solution.message)
 return solution.y.T.reshape(len(samples),3,model.n)

def stability(model,background):
 J=model.jac(0,background);v0=np.sin(np.arange(3*model.n)*np.sqrt(2))+np.cos(np.arange(3*model.n)*.17);start=time.perf_counter()
 values,vectors=eigs(J,k=6,which='LR',v0=v0,ncv=60,tol=1e-8,maxiter=15000);near,_=eigs(J,k=12,sigma=0,which='LM',v0=v0,ncv=60,tol=1e-9,maxiter=15000)
 residuals=[float(np.linalg.norm(J@vectors[:,i]-values[i]*vectors[:,i])/np.linalg.norm(vectors[:,i])) for i in range(6)]
 if max(residuals)>1e-5:raise RuntimeError('Eigenpair residual')
 maximum=max(float(values.real.max()),float(near.real.max()));local=[]
 for aa,bb in ((float(model.a.min()),float(model.b.min())),(float(model.a.max()),float(model.b.max()))):
  X=aa/.9;reaction=np.array([[-1,.1,0],[1,bb-1,X*X],[0,-bb,-X*X]]);modes=4/model.h**2*np.sin(np.pi*np.arange(model.n)/model.n)**2;disp=np.linalg.eigvals(reaction[None]-modes[:,None,None]*np.diag([1.,1.,12.])[None]).real.max();local.append(float(disp))
 return dict(passed=maximum < -1e-8 and max(local)<-1e-8,maximum_real=maximum,rightmost=[[float(z.real),float(z.imag)] for z in values],near_zero_maximum_real=float(near.real.max()),maximum_eigenpair_residual=max(residuals),local_extrema_maxima=local,seconds=time.perf_counter()-start,numerical_not_rigorous_spectral_certificate=True)

def uniform_stability(model):
 assert np.ptp(model.a)==0 and np.ptp(model.b)==0
 X=model.a[0]/.9;b=model.b[0];reaction=np.array([[-1,.1,0],[1,b-1,X*X],[0,-b,-X*X]]);modes=4/model.h**2*np.sin(np.pi*np.arange(model.n)/model.n)**2;maximum=float(np.linalg.eigvals(reaction[None]-modes[:,None,None]*np.diag([1.,1.,12.])[None]).real.max())
 return dict(passed=maximum < -1e-8,maximum_real=maximum,method='all discrete homogeneous Fourier modes')

def measure_series(fields,background,h,ref,previous=None,lost=False):
 values=np.array([core(v[1],background.reshape(3,-1)[1],h,ref['depth']) for v in fields]);valid=np.zeros(len(values),bool)
 for i,v in enumerate(values):
  if not np.isfinite(v[:3]).all() or v[3]!=1:lost=True
  elif previous is not None and abs(wrap(v[0]-previous))>ref['wavelength']/2:lost=True
  if not lost:valid[i]=True;previous=float(v[0])
 return values,valid,previous,lost

def velocity_blocks(observations,times,valid,ref):
 late=times>=500;shape=valid&(observations[:,1]>=.6*ref['depth'])&(observations[:,2]>=.75*ref['width'])&(observations[:,2]<=1.25*ref['width']);centers=np.full(len(times),np.nan)
 if valid.all():centers=np.unwrap(observations[:,0]*2*np.pi/100)*100/(2*np.pi)
 blocks=[]
 for left in range(500,2000,50):
  i=left*2;j=(left+50)*2;ok=bool(shape[i:j+1].all());blocks.append(dict(left=left,right=left+50,valid=ok,velocity=float((centers[j]-centers[i])/50) if ok else None))
 return dict(valid_late=bool(shape[late].all()),blocks=blocks,center_unwrapped=centers,shape_valid=shape)
