"""Canonical MG049 spectral measurement, deliberately shared by analysis paths."""
import numpy as np

NAMES=("power_G","power_X","power_Y","cross_GX","cross_GY","cross_XY","real_X","imag_X")


def measure(populations,mean,omega,time,modes):
    n=populations.shape[-1]
    z=np.fft.rfft(populations-mean,axis=-1)[...,modes]/np.sqrt(n*omega*time)
    # Shape (...,species,mode). rFFT fixes the Nyquist imaginary component.
    g,x,y=z[...,0,:],z[...,1,:],z[...,2,:]
    values=np.stack([abs(g)**2,abs(x)**2,abs(y)**2,(g*x.conj()).real,
                     (g*y.conj()).real,(x*y.conj()).real,x.real,x.imag],axis=-1)
    return np.round(values,8),z
