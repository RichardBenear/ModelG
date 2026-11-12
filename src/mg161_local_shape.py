"""MG161 fixed local template-bank fit; calibration accounts for search and noise."""
import numpy as np
class LocalShape:
    def __init__(self,reference,initial_x,omega,h=.2):
        reference=np.asarray(reference,float);self.n=len(reference);self.h=h;self.omega=omega;self.scales=np.arange(.25,3.0001,.125);self.anchor=int(np.argmax(abs(reference-reference.mean())));self.coords=((np.arange(self.n)+self.n//2)%self.n-self.n//2)*h;order=np.argsort(self.coords);far=abs(self.coords)>=self.n*h/4;raw=np.roll(reference,-self.anchor);raw=raw-np.median(raw[far]);self.raw=raw;self.order=order;freq=2*np.pi*np.fft.rfftfreq(self.n,d=h);self.smoother=np.exp(-.5*(freq*.4)**2)
        smooth=np.fft.irfft(np.fft.rfft(raw)*self.smoother,n=self.n);energy=smooth*smooth;radii=np.sort(np.unique(abs(self.coords)));self.radius95=float(next(r for r in radii if energy[abs(self.coords)<=r].sum()>=.95*energy.sum()));self.radius95=max(self.radius95,h);self.window=np.exp(-.5*(self.coords/(2*self.radius95))**2);self.ws=float(self.window.sum());self.window_fft=np.conj(np.fft.rfft(self.window));self.bank=np.array([np.fft.irfft(np.fft.rfft(self.signal(s,0))*self.smoother,n=self.n) for s in self.scales]);local_mean=(self.bank*self.window).sum(1)/self.ws;kernels=self.window[None,:]*(self.bank-local_mean[:,None]);self.norm=np.sum(kernels*self.bank,axis=1);assert np.all(self.norm>0);self.kernels_fft=np.conj(np.fft.rfft(kernels,axis=1));self.background_exclusion=max(2*self.radius95,2*h);self.initial_amplitude=self.fit(np.asarray(initial_x,float)/omega)['amplitude']
    def signal(self,scale,shift):
        core=np.interp(self.coords/scale,self.coords[self.order],self.raw[self.order],left=0.,right=0.);return np.roll(core,int(shift))
    def fit(self,field,noise_scale=1.):
        field=np.asarray(field,float)
        if field.shape!=(self.n,) or not np.isfinite(field).all() or noise_scale<=0:raise ValueError('Invalid field/noise scale')
        f=np.fft.irfft(np.fft.rfft(field-field.mean())*self.smoother,n=self.n);ft=np.fft.rfft(f);dot=np.fft.irfft(ft[None,:]*self.kernels_fft,n=self.n,axis=1);gain=np.maximum(dot,0.)**2/self.norm[:,None];s,k=np.unravel_index(int(gain.argmax()),gain.shape);amplitude=float(dot[s,k]/self.norm[s]);weighted_sum=np.fft.irfft(ft*self.window_fft,n=self.n);weighted_sq=np.fft.irfft(np.fft.rfft(f*f)*self.window_fft,n=self.n);local_energy=float(weighted_sq[k]-weighted_sum[k]**2/self.ws);residual=local_energy-float(gain[s,k]);assert residual>=-1e-8*max(1,local_energy);positions=np.arange(self.n);distance=np.minimum((positions-k)%self.n,(k-positions)%self.n)*self.h;outside=distance>self.background_exclusion;second=float(gain[:,outside].max()) if np.any(outside) else 0.
        return dict(amplitude=amplitude,scale=float(self.scales[s]),center_site=int(k),shift_sites=int((k-self.anchor+self.n//2)%self.n-self.n//2),score=float(gain[s,k]/noise_scale**2),residual=float(max(0,residual)/noise_scale**2),secondary_score=second/noise_scale**2,secondary_ratio=second/max(float(gain[s,k]),1e-30),bank_edge=bool(s in (0,len(self.scales)-1)))
    def classify(self,result,detection_threshold,residual_threshold):
        relative=result['amplitude']/self.initial_amplitude
        if result['score']<=detection_threshold or relative<.5:return 'absent_or_weak'
        if result['secondary_score']>detection_threshold and result['secondary_ratio']>=.5:return 'ambiguous_multiple_cores'
        if result['bank_edge'] or not .5<=result['scale']<=2.:return 'deformed_scale'
        if result['residual']>residual_threshold:return 'shape_mismatch'
        return 'core_compatible'
