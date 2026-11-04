"""Independent check of the post-MG030 centre diagnostic; no production imports or integrations.
Usage: python centre_audit.py CHECKOUT
"""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks

root=Path(sys.argv[1])/'results/MG9c'
reported=json.loads((root/'diagnostics/centre.json').read_text())
z=np.load(root/'snapshots.npz'); n=500; h=.2; lam=100/34
x=(np.arange(n)-n//2)*h
kernel=np.zeros(n); half=lam/(2*h)
for j in range(-10,11):
    kernel[j%n]=max(0,min(j+.5,half)-max(j-.5,-half))
kernel/=kernel.sum(); multiplier=np.fft.fft(kernel)**2
smooth=lambda f:np.fft.ifft(np.fft.fft(f)*multiplier)
env=lambda f:2*abs(smooth((f-smooth(f).real)*np.exp(-2j*np.pi*x/lam)))
k=2*np.pi*np.fft.fftfreq(n,d=h)
shift=lambda f,s:np.fft.ifft(np.fft.fft(f)*np.exp(-1j*k*s)).real

def minima(f,E):
    pk,_=find_peaks(-np.tile(f,3),prominence=E/2,distance=int(np.ceil(lam/(2*h))))
    pk=pk[(pk>=n)&(pk<2*n)]-n
    a=f[(pk-1)%n]; b=f[pk]; c=f[(pk+1)%n]
    return x[pk]+h*(a-c)/(2*(a-2*b+c))

out={'reviewed_revision':'2417833460d0209c6a6f256d3487ef877226524f','runs':{}}
errors=[]
for tag,rr in reported['runs'].items():
    rec=json.loads((root/'runs'/f'front_D02_{tag}.json').read_text())
    rows=rr['rows']; computed={}
    for row in rows:
        t=row['t']; f=z[f'front_D02_{tag}_t{t}'][n:2*n]; E=env(f)
        cs={}
        for name,w in [('amp',E**2),('phi2',f**2)]:
            C=float(x@w/sum(w)); var=float(((x-C)**2)@w/sum(w))
            vals={'C':C,'total':float(sum(w)*h),'var':var,'edge_weight_frac':float(sum(w[abs(x)>45])/sum(w))}
            for key,v in vals.items():
                errors.append(abs(v-row[name][key]));np.testing.assert_allclose(v,row[name][key],atol=1e-10,rtol=1e-12)
            cs[name]=C
        computed[t]=cs
    f=z[f'front_D02_{tag}_t15000'][n:2*n]; g=z[f'front_D02_{tag}_t30000'][n:2*n]
    # Coarse search checks uniqueness before refining the best alignment locally.
    obj=lambda s:float(np.sum((g-shift(f,s))**2))
    grid=np.linspace(-lam/4,lam/4,1001); vals=np.array([obj(s) for s in grid]); i=int(vals.argmin())
    assert 0<i<len(grid)-1
    fit=minimize_scalar(obj,bounds=(grid[i-1],grid[i+1]),method='bounded',options={'xatol':1e-12})
    s=float(fit.x); ratio=float(np.sqrt(fit.fun)/np.linalg.norm(g-f))
    np.testing.assert_allclose(s,rr['late']['alignment'],atol=1e-8)
    np.testing.assert_allclose(ratio,rr['late']['alignment_rem_of_change'],atol=1e-8)
    ma,mb=minima(f,rec['E_pat']),minima(g,rec['E_pat'])
    assert len(ma)==len(mb)
    dd=(mb-ma)[1:-1]
    m0=minima(z[f'front_D02_{tag}_t1800'][n:2*n],rec['E_pat'])[1:-1]
    ia=np.array([np.argmin(abs(ma-p)) for p in m0]); ib=np.array([np.argmin(abs(mb-p)) for p in m0])
    assert len(set(ia))==len(ia) and len(set(ib))==len(ib)
    anchored=mb[ib]-ma[ia]
    np.testing.assert_allclose(anchored.mean(),rr['late']['interior_minima'],atol=1e-12)
    # Additional interpolation sensitivity: locate minima of Fourier-interpolated profiles.
    def spectral_minima(field,seeds):
        coeff=np.fft.fft(field)/n
        def value(pos):return float(np.real(np.sum(coeff*np.exp(1j*k*(pos+50)))))
        return np.array([minimize_scalar(value,bounds=(p-h,p+h),method='bounded',options={'xatol':1e-11}).x for p in seeds])
    spectral_dd=(spectral_minima(g,mb)-spectral_minima(f,ma))[1:-1]
    shifts={'edge':rr['late']['C_edge'],'amp':computed[30000]['amp']-computed[15000]['amp'],
            'phi2':computed[30000]['phi2']-computed[15000]['phi2'],'minima':float(anchored.mean()),'alignment':s}
    out['runs'][tag]={'late_shifts':shifts,'unexplained_norm_fraction':ratio,
        'minima_sd':float(dd.std()),'max_minimum_deviation_from_mean':float(max(abs(dd-dd.mean()))),
        'late_order_subset_mean':float(dd.mean()),'original_anchored_subset_sd':float(anchored.std()),
        'original_anchored_subset_max_deviation':float(max(abs(anchored-anchored.mean()))),
        'spectral_minima_mean_shift':float(spectral_dd.mean()),
        'max_parabolic_vs_spectral_shift_difference':float(max(abs(dd-spectral_dd))),
        'max_relative_difference_from_alignment':max(abs(v-s)/abs(s) for v in shifts.values())}
out['max_center_moment_absolute_error']=max(errors)
Path(__file__).with_name('centre_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
