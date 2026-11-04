"""Synthetic MG027 audit only. Argument: reviewed ModelG checkout."""
import importlib.util
import json
import sys
from pathlib import Path
import numpy as np

repo = Path(sys.argv[1])
sys.path.insert(0, str(repo / "src"))
import run_mg9 as run

lam = 100/34

def average(f, h):
    # Independent physical cell-overlap quadrature, implemented as circular convolution.
    offsets = np.arange(-int(np.ceil(lam/h)), int(np.ceil(lam/h))+1)
    weights = np.maximum(0, np.minimum((offsets+.5)*h,lam/2) -
                         np.maximum((offsets-.5)*h,-lam/2))
    kernel = np.zeros(len(f))
    for j,w in zip(offsets,weights):
        kernel[j % len(f)] += w/lam
    return np.fft.ifft(np.fft.fft(f)*np.fft.fft(kernel))

def envelope(f,h):
    x=(np.arange(len(f))-len(f)//2)*h
    smooth=lambda a:average(average(a,h),h)
    return 2*abs(smooth((f-smooth(f).real)*np.exp(-2j*np.pi*x/lam)))

rows={}
for N,h in [(500,.2),(1000,.1),(1000,.2)]:
    x=(np.arange(N)-N//2)*h
    worst=0.
    cosine_lo,cosine_hi=2.,0.
    agreement=0.
    for phase in np.linspace(0,2*np.pi,40,endpoint=False):
        wave=np.cos(2*np.pi*x/lam+phase)
        e=envelope(wave,h)
        cosine_lo=min(cosine_lo,float(e.min()))
        cosine_hi=max(cosine_hi,float(e.max()))
        for a in (25.,17.3):
            for ell in (lam/2,lam/8,2*lam):
                f=wave*.5*(np.tanh((x+a)/ell)-np.tanh((x-a)/ell))
                e=envelope(f,h)
                agreement=max(agreement,float(np.max(abs(e-run.envelope(f,h)))))
                left=np.flatnonzero((e[:-1]<.5)&(e[1:]>=.5))
                right=np.flatnonzero((e[:-1]>=.5)&(e[1:]<.5))
                def crossing(indices,target):
                    points=x[indices]+h*(.5-e[indices])/(e[indices+1]-e[indices])
                    return points[np.argmin(abs(points-target))]
                worst=max(worst,abs(crossing(left,-a)+a)/lam,abs(crossing(right,a)-a)/lam)
    assert .98<=cosine_lo<=cosine_hi<=1.02 and worst<=.1 and agreement<1e-12
    harmonic=np.exp(4j*np.pi*x/lam)
    transfer=float(np.max(abs(average(harmonic,h))))
    rows[f"N{N}_h{h}"]=dict(cosine_range=[cosine_lo,cosine_hi],
                            worst_edge_error_in_wavelengths=worst,
                            independent_vs_production_max_error=agreement,
                            boxcar_2k_gain=transfer)
spec=importlib.util.spec_from_file_location("tests_mg9",repo/"tests/test_mg9.py")
tests=importlib.util.module_from_spec(spec)
spec.loader.exec_module(tests)
selected=["test_registered_rms_envelope_fails_its_validation",
          "test_envelope_validation_amendment_1","test_tracker_on_moving_patch",
          "test_classification_procedure","test_kernel_weights_sum_to_lam_over_h"]
for name in selected:
    getattr(tests,name)()
print(json.dumps(dict(domains=rows,repository_synthetic_tests_passed=selected),indent=2))
