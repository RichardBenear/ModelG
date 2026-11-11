"""MG145 saved-input precision and archive-layout planning; no RNG or sampling."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import numpy as np,json,math,time,hashlib
from scipy.stats import t,norm
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG12H3/finer_tick_design'
S=np.array([140/9,140/9,29/(140/9)])
def run():
 BASE.mkdir(parents=True,exist_ok=True);clock=time.perf_counter();n=500;omega=10
 source=ROOT/'results/MG12H3/ring_feasibility/o0_a0_Omega10.npz'
 with np.load(source) as z:initial=z['counts']
 freq=np.fft.fftfreq(n)*n;mask=np.abs(freq)>=125;H=np.fft.ifft(np.fft.fft(np.eye(n),axis=0)*mask[:,None],axis=0).real;np.testing.assert_allclose(H@H,H,rtol=1e-10,atol=1e-10);assert int(mask.sum())==251
 rows=[];inputs={source.relative_to(ROOT).as_posix():hashlib.sha256(source.read_bytes()).hexdigest()}
 for species,name in enumerate(('G','X','Y')):
  p=ROOT/f'results/MG12H3/tick_resolution/o0_a0_Omega1_{name}.npz';inputs[p.relative_to(ROOT).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
  with np.load(p) as z:mu=z['eigenvalues'];U=z['eigenvectors'];exit_rate=float(z['exits'][0])
  template=initial[species]/omega/S[species];template=template-template.mean();template/=np.sqrt(np.mean(template**2));weight=template/(n*omega*S[species]);bytime={}
  for horizon in (.005,.01):
   variances={}
   for delta in (.0001,.00005,.00002):
    factor=-math.expm1(-exit_rate*delta)/exit_rate;steps=round(horizon/delta);assert abs(steps*delta-horizon)<1e-14;mult=(1+factor*mu)**steps;K=(U*mult)@U.T;assert K.min()>-1e-10
    counts=initial[species].astype(float);variance=float(counts@(K.T@(weight**2)-(K.T@weight)**2));C=np.diag(K@counts)-(K*counts[None,:])@K.T;HC=H@C@H;mean_delta=H@(K@counts-counts);power_variance=float((2*np.sum(HC**2)+4*mean_delta@HC@mean_delta)/(n*omega*S[species])**2);assert variance>=0 and power_variance>=0;variances[delta]=(variance,power_variance)
   for delta in (.0001,.00005):
    for reps in (32,48,64):
     for observable,k,margin in [('contrast',0,.01),('high_power',1,.05)]:
      v1,v2=variances[delta][k],variances[.00002][k];df=(v1+v2)**2/(v1*v1/(reps-1)+v2*v2/(reps-1));critical=float(t.ppf(1-.05/(2*24),df));se=math.sqrt((v1+v2)/reps);half=critical*se;power=max(0,float(norm.cdf(margin/se-critical)-norm.cdf(-margin/se+critical))) if margin>half else 0.;rows.append(dict(species=name,horizon=horizon,delta=delta,replicates=reps,observable=observable,margin=margin,transport_baseline_standard_error=se,nominal_halfwidth=half,approximate_zero_difference_equivalence_power=power))
 # Saved-data storage prototype, not a change to the old archive.
 directory=ROOT/'results/MG12H3/multitick_pilot/production/o0_a0_Omega10_d1';ids=[];outputs=[];available=[];total_old=0
 for file in sorted(directory.glob('tick_*.npz')):
  with np.load(file) as z:ids.append(np.bincount(z['channel_ids'],minlength=9000).astype(np.uint32));outputs.append(z['final']);available.append(z['available'])
  total_old+=file.stat().st_size
 compact=BASE/'saved_pilot_compact_layout.npz';np.savez_compressed(compact,channel_histograms=np.asarray(ids),final=np.asarray(outputs),available=np.asarray(available));storage=dict(ticks=len(ids),old_tick_trace_bytes=total_old,compact_histogram_state_bytes=compact.stat().st_size,excludes='RNG states, checkpoints and full sentinel traces; saved-layout comparison only, not exact new-writer forecast')
 pilot=json.loads((directory.parent/'o0_a0_Omega10_d1.json').read_text());bench=json.loads((ROOT/'results/MG12H3/sparse_development/benchmark.json').read_text());setup=max(x['setup_seconds'] for x in bench['records']);times={str(delta):pilot['elapsed_seconds']+max(0,round(.01/delta)-20)*setup for delta in (.0001,.00005,.00002)}
 budget=[dict(replicates=r,total_paths=3*r,planning_body_seconds=sum(times.values())*r,twice_body_core_hours=2*sum(times.values())*r/3600) for r in (32,48,64)]
 assert time.perf_counter()-clock<60;result=dict(complete=True,elapsed_seconds=time.perf_counter()-clock,precision=rows,storage=storage,budget=budget,source_hashes=inputs,random_draws=0,limitations='Exact transport-only linear contrast variance; Gaussian approximation for quadratic-power variance and rough Welch critical values. Not coupled-model power or a concrete run budget. Old pilot data unchanged.')
 (BASE/'planning.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(elapsed=result['elapsed_seconds'],storage=storage,budget=budget)))
if __name__=='__main__':run()
