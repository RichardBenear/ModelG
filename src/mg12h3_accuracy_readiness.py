"""MG158 saved-only paired accuracy planning; no stochastic trajectories."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import json,hashlib,math
import numpy as np
from scipy.stats import t,norm,chi2
from mg12h3_finer_precision import setup,power
from mg10e2c1_pool import atomic_json
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG12H3/tick_accuracy_readiness';PILOT=ROOT/'results/MG12H3/coupling_long_pilot'
def intervals(values,margins):
    a=np.asarray(values,float);m=np.asarray(margins,float)
    if a.ndim!=3 or a.shape[1:]!=(2,6) or len(a)<2 or m.shape!=(6,) or not np.isfinite(a).all() or not np.isfinite(m).all() or np.any(m<=0):raise ValueError('Invalid complete paired sample')
    mean=a.mean(0);sd=a.std(0,ddof=1);critical=float(t.ppf(1-.05/(2*24),len(a)-1));hw=critical*sd/math.sqrt(len(a));low=mean-hw;high=mean+hw
    labels=np.full((2,6),'unresolved',dtype='<U12');labels[(low>=-m)&(high<=m)&(sd>0)]='equivalent';labels[(low>m)|(high<-m)]='incompatible';labels[sd==0]='unresolved'
    overall='equivalent' if np.all(labels=='equivalent') else ('incompatible' if np.any(labels=='incompatible') else 'unresolved')
    return dict(mean=mean.tolist(),SD=sd.tolist(),critical=critical,lower=low.tolist(),upper=high.tolist(),labels=labels.tolist(),overall=overall)
def prepare():
    initial,contrast=setup();margins=np.r_[.05*abs(contrast.direct(initial)),[.05]*3];prior=json.loads((PILOT/'analysis_results.json').read_text());rows=[];inputs={};samples=[]
    # Deterministic classification and rejected-input checks; no RNG.
    x=np.tile(np.array([-1.,1.])[:,None,None],(8,2,6))*.000001
    assert intervals(x,margins)['overall']=='equivalent';bad=x.copy();bad[:,:,2]+=.1;assert intervals(bad,margins)['overall']=='incompatible';near=x.copy();near[:,:,2]+=margins[2];assert intervals(near,margins)['overall']=='unresolved';assert intervals(x*0,margins)['overall']=='unresolved';guards=4
    for invalid in (np.full_like(x,np.nan),np.full_like(x,np.inf),x[:1],x[:,:,:5]):
        try:intervals(invalid,margins)
        except ValueError:guards+=1
        else:raise AssertionError('Missing incomplete/invalid sample guard')
    for comp in range(2):
        data=[];fields=[]
        for r in range(8):
            path=PILOT/'production'/f'c{comp}_paired_r{r:02d}/observations.npz';inputs[path.relative_to(ROOT).as_posix()]=hashlib.sha256(path.read_bytes()).hexdigest()
            with np.load(path) as z:counts=z['counts']
            data.append([[*(contrast.direct(counts[f,0])-contrast.direct(counts[f,1])),*(power(counts[f,0],initial)-power(counts[f,1],initial))] for f in (5,10)])
            # Separate Parseval calculation of high-band power for all relevant fields.
            for f in (5,10):
                for side in (0,1):
                    delta=(counts[f,side,:3]-initial[:3])/10.;fft=np.fft.fft(delta,axis=1);freq=np.fft.fftfreq(500)*500;selected=np.abs(freq)>=125;other=10/np.array([140/9,140/9,29/(140/9)])*np.sum(abs(fft[:,selected])**2,axis=1)/500**2;np.testing.assert_allclose(other,power(counts[f,side],initial),rtol=1e-12,atol=1e-12)
        data=np.asarray(data);samples.append(data);mean=data.mean(0);sd=data.std(0,ddof=1);np.testing.assert_allclose(mean[:,:3],np.asarray(prior['comparisons'][comp]['mean_paired_difference'])[[5,10]],atol=1e-12,rtol=1e-12)
        grid=[]
        for n in (32,64,96,128,192,256):
            critical=float(t.ppf(1-.05/48,n-1));hw=critical*sd/math.sqrt(n);ratio=(abs(mean)+hw)/margins;gap=(margins-abs(mean))*math.sqrt(n)/sd;per_cell=norm.cdf(gap-critical);grid.append(dict(N=n,critical=critical,halfwidth=hw.tolist(),projected_edge_fraction=ratio.tolist(),worst_edge_fraction=float(ratio.max()),minimum_normal_approx_cell_pass_probability=float(per_cell.min())))
        # With eight pilot observations, SD and true mean are very uncertain.
        sd_multiplier=math.sqrt(7/chi2.ppf(.05,7));chosen=256 if comp==0 else 64;nominal=next(g for g in grid if g['N']==chosen);sensitivity=(abs(mean)+np.asarray(nominal['halfwidth'])*sd_multiplier)/margins
        rows.append(dict(comparison=comp,pilot_N=8,mean=mean.tolist(),SD=sd.tolist(),grid=grid,recommended_N=chosen,SD_upper_95_multiplier=sd_multiplier,sensitivity_worst_edge_fraction=float(sensitivity.max()),measured_pair_seconds=prior['comparisons'][comp]['mean_pair_cost_seconds']['paired']))
    np.savez_compressed(BASE/'saved_pilot_observables.npz',values=samples,margins=margins,times=[.005,.01])
    measured=sum(r['recommended_N']*r['measured_pair_seconds'] for r in rows);forecast=1.15*measured;atomic_json(BASE/'planning.json',dict(stage='MG158',scope='saved-only readiness, no accuracy outcome or new trajectories',margins=margins.tolist(),comparisons=rows,synthetic_checks=guards,measured_worker_seconds=measured,forecast_worker_seconds=forecast,forecast_wall_hours_at8=forecast/8/3600,proposed_worker_ceiling=300000,sessions=3,workers=8,job_cap_seconds=1200,input_hashes=inputs,limitations='Plug-in planning with eight pilot pairs; per-cell normal power is not joint power or a guarantee. Sensitivity holds pilot means fixed; true means and SD are uncertain. No pilot outcomes are scored as qualification.'))
    print(json.dumps(dict(comparisons=[dict(c=r['comparison'],N=r['recommended_N'],nominal_worst_edge=next(g for g in r['grid'] if g['N']==r['recommended_N'])['worst_edge_fraction'],sensitivity=r['sensitivity_worst_edge_fraction']) for r in rows],forecast_hours=forecast/8/3600,checks=guards)))
if __name__=='__main__':BASE.mkdir(parents=True,exist_ok=True);prepare()
