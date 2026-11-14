"""X/Y-only shared-geometry candidate with separate component and residual second-core checks."""
import numpy as np
from mg161_local_shape import LocalShape
class XYCore:
    def __init__(self,reference,initial,omega,cov):
        self.x=LocalShape(reference[1],initial[1],omega);m=self.x;self.cov=cov;self.inv=np.linalg.inv(cov);self.n=m.n
        raw=np.roll(reference[1:3],-m.anchor,axis=1);raw-=np.median(raw[:,abs(m.coords)>=25],axis=1)[:,None];self.raw=raw
        self.bank=np.array([self.signal(s,0) for s in m.scales]);self.g=self.bank-(self.bank@m.window/m.ws)[:,:,None];self.k=np.einsum('ij,sjn->sin',self.inv,self.g)*m.window;self.norm=np.sum(self.k*self.g,axis=(1,2));self.kfft=np.conj(np.fft.rfft(self.k,axis=2));self.ownnorm=np.sum(self.g*self.g*m.window,axis=2)
        initialfit=self.fit(initial[1:3]/omega);self.base=initialfit['amplitude'];self.initial_components=np.array(initialfit['individual_amplitudes'])
    def signal(self,s,k):
        m=self.x;raw=np.array([np.interp(m.coords/s,m.coords[m.order],r[m.order],left=0,right=0) for r in self.raw]);return np.roll(np.fft.irfft(np.fft.rfft(raw,axis=1)*m.smoother,n=m.n,axis=1),k,axis=1)
    def fit(self,field,already_smooth=False):
        m=self.x;f=field-field.mean(axis=1)[:,None]
        if not already_smooth:f=np.fft.irfft(np.fft.rfft(f,axis=1)*m.smoother,n=m.n,axis=1)
        dots=np.fft.irfft((np.fft.rfft(f,axis=1)[None,:,:]*self.kfft).sum(axis=1),n=m.n,axis=1);gain=np.maximum(dots,0)**2/self.norm[:,None];si,k=np.unravel_index(gain.argmax(),gain.shape);amp=dots[si,k]/self.norm[si];w=np.roll(m.window,k);fc=f-(f@w/m.ws)[:,None];gc=np.roll(self.g[si],k,axis=1);err=fc-amp*gc;res=np.einsum('in,ij,jn,n->',err,self.inv,err,w);individual=(fc*gc*w).sum(1)/self.ownnorm[si];deviation=abs(individual-amp)*np.sqrt(self.ownnorm[si]/np.diag(self.cov));per=(err*err)@w
        remainder=f-amp*np.roll(self.bank[si],k,axis=1);rdots=np.fft.irfft((np.fft.rfft(remainder,axis=1)[None,:,:]*self.kfft).sum(axis=1),n=m.n,axis=1);rgain=np.maximum(rdots,0)**2/self.norm[:,None];distance=np.minimum((np.arange(m.n)-k)%m.n,(k-np.arange(m.n))%m.n)*m.h;secondary=float(rgain[:,distance>m.background_exclusion].max())
        return dict(secondary_score=secondary,amplitude=float(amp),scale=float(m.scales[si]),center=int(k),score=float(gain[si,k]),residual=float(res),component_stat=float(deviation.max()),individual_amplitudes=individual.tolist(),component_deviations=deviation.tolist(),component_residuals=per.tolist(),bank_edge=bool(si in(0,len(m.scales)-1)))
    def label(self,f,c):
        if f['score']<=c['score'] or f['amplitude']/self.base<.5:return 'weak'
        if f['bank_edge'] or not .5<=f['scale']<=2:return 'scale'
        if np.any(np.array(f['individual_amplitudes'])/self.initial_components<.5):return 'component_weak'
        if f['component_stat']>c['component']:return 'component_mismatch'
        if f['secondary_score']>c['secondary']:return 'competing_core'
        if f['residual']>c['residual']:return 'shape_mismatch'
        return 'candidate_compatible'
