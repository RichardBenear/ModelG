"""Stationary Gaussian AR1 joint-profile variance-consistency diagnostic."""
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import chi2
D=998;N=64;M=16*D

def fit_reference(qspace,A,B,C,qfirst):
 def qtime(r):return qfirst+(A-2*r*B+r*r*C)/(1-r*r)
 def loglike(r,joint):
  qt=qtime(r);df=N*D+M if joint else N*D;ss=(qt+qspace)/df if joint else qt/df
  return -.5*(df*np.log(ss)+63*D*np.log1p(-r*r)+df)
 fits=[]
 for joint in (False,True):
  opt=minimize_scalar(lambda r:-loglike(r,joint),bounds=(-.99,.99),method='bounded',options={'xatol':1e-11});r=float(opt.x);ll=loglike(r,joint)
  if not joint:ll+=-.5*(M*np.log(qspace/M)+M)
  fits.append(dict(rho=r,scale=float((qtime(r)+(qspace if joint else 0))/(N*D+(M if joint else 0))),loglike=float(ll),success=bool(opt.success),interior=bool(abs(r)<.9899)))
 lr=max(0.,2*(fits[0]['loglike']-fits[1]['loglike']));return dict(LR=lr,p=float(chi2.sf(lr,1)),unrestricted=fits[0],shared=fits[1],valid=all(x['success'] and x['interior'] for x in fits),asymptotic_chi_square_df=1)
