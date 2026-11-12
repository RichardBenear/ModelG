"""MG160 saved-only FFT measurement recomputation, separate from lag-loop code."""
import json
from pathlib import Path
import numpy as np
from mg10e2c1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG12H3/pattern_persistence'
def check():
    result=json.loads((BASE/'preparation_results.json').read_text())
    with np.load(BASE/'prepared_inputs.npz') as z:initial=z['initial'];null=z['homogeneous'];ref=z['reference'];omega=float(z['omega'])
    k=2*np.pi*np.fft.fftfreq(500,d=.2);kernel=np.exp(-.5*(k*.4)**2)
    def centered(v):return np.fft.ifft(np.fft.fft(v-v.mean())*kernel).real
    template=centered(ref[1]);norm=template@template
    for name,pop in (('initial_metrics',initial),('homogeneous_metrics',null)):
        f=centered(pop[1]/omega);corr=np.fft.ifft(np.fft.fft(f)*np.conj(np.fft.fft(template))).real;peak=float(corr.max());expected=result[name];amplitude=peak/norm;rho=peak/np.sqrt((f@f)*norm);extent=.2*(f@f)**2/np.sum(f**4);np.testing.assert_allclose([amplitude,rho,extent],[expected['amplitude'],expected['correlation'],expected['participation_length']],rtol=1e-10,atol=1e-12)
    np.testing.assert_array_equal(initial.sum(1),null.sum(1));atomic_json(BASE/'separate_measurement_recomputation.json',dict(passed=True,method='Explicit Fourier Gaussian smoothing and FFT cross-correlation versus translated-template loop; exact species totals. Same author/shared saved inputs, not external replication.'))
if __name__=='__main__':check()
