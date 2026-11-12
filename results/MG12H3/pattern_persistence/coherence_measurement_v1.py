"""MG160 prepared-pattern coherence measurements; no particle-identity claim."""
import numpy as np
from measure import smooth_periodic
def homogenize(pop):
    """Integer mass-matched homogeneous control: each site differs by at most one."""
    pop=np.asarray(pop)
    if pop.ndim!=2 or pop.shape[0]!=7 or pop.dtype.kind!='i' or np.any(pop<0):raise ValueError('Invalid counts')
    n=pop.shape[1];out=np.zeros_like(pop)
    for s in range(7):
        total=sum(map(int,pop[s]));base,rem=divmod(total,n);out[s]=base+(np.arange(1,n+1)*rem//n-np.arange(n)*rem//n);out[s]=np.roll(out[s],(63*s)%n)
    assert np.array_equal(out.sum(1),pop.sum(1));return out
class Coherence:
    def __init__(self,reference,initial,omega,h=.2):
        reference=np.asarray(reference,float);initial=np.asarray(initial)
        if reference.shape!=initial.shape or reference.shape[0]!=7 or reference.ndim!=2 or not np.isfinite(reference).all() or omega<=0:raise ValueError('Invalid template')
        self.omega=omega;self.h=h;self.n=initial.shape[1];self.template=smooth_periodic(reference[1]-reference[1].mean(),h);self.norm=float(self.template@self.template)
        if self.norm<=0:raise ValueError('Flat reference')
        raw=self._measure(initial);self.initial_amplitude=raw['amplitude']
        if self.initial_amplitude<=0:raise ValueError('No initial feature')
    def _measure(self,pop):
        a=np.asarray(pop)
        if a.shape!=(7,self.n) or not np.isfinite(a).all() or np.any(a<0):raise ValueError('Invalid field')
        x=a[1]/self.omega;field=smooth_periodic(x-x.mean(),self.h);energy=float(field@field)
        if energy<=1e-28:return dict(amplitude=0.,correlation=0.,shift_sites=None,phase_ambiguous=False,participation_length=0.,residual_fraction=0.)
        # Explicit translations keep the lag convention auditable (500 sites only).
        dots=np.array([field@np.roll(self.template,k) for k in range(self.n)]);best=int(dots.argmax());peak=float(dots[best]);amplitude=peak/self.norm;correlation=peak/np.sqrt(energy*self.norm);tied=np.flatnonzero(abs(dots-peak)<=1e-10*max(abs(peak),1e-20));dist=np.minimum((tied-best)%self.n,(best-tied)%self.n);ambiguous=bool(np.any(dist>1));residual=field-amplitude*np.roll(self.template,best)
        return dict(amplitude=amplitude,correlation=float(correlation),shift_sites=int((best+self.n//2)%self.n-self.n//2),phase_ambiguous=ambiguous,participation_length=float(self.h*energy**2/np.sum(field**4)),residual_fraction=float(np.linalg.norm(residual)/np.sqrt(energy)))
    def measure(self,pop):
        out=self._measure(pop);out['relative_amplitude']=out['amplitude']/self.initial_amplitude;out['coherent']=bool(out['relative_amplitude']>=.5 and out['correlation']>=.8);return out
def continuity(times,metrics):
    t=np.asarray(times,float)
    if len(t)!=len(metrics) or len(t)<2 or t[0]!=0 or not np.isfinite(t).all() or np.any(np.diff(t)<=0):raise ValueError('Incomplete observation grid')
    failed=[i for i,m in enumerate(metrics) if not m['coherent']]
    if not metrics[0]['coherent']:return dict(status='initial feature not resolved',first_loss_interval=None)
    if not failed:return dict(status='sampled coherence retained',right_censored_at=float(t[-1]))
    i=failed[0];return dict(status='loss of initial coherence',first_loss_interval=[float(t[i-1]),float(t[i])])
