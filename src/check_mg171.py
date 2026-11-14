"""Independent cache aggregation and componentwise spectral audit verification."""
import time,json,gzip,hashlib,math
from pathlib import Path
import numpy as np
from audit_mg171 import OUT,RAW,CLASSES,read,atomic_json,binding

def ls(values):
    values=np.asarray(values);m=float(values.max());return m+math.log(math.fsum(math.exp(float(v)-m) for v in values))
def weights(logs):
    a=np.exp(np.array(logs)-max(logs));return a/a.sum()
def summary(rows):
    if not rows:return dict(n=0)
    return dict(n=len(rows),median_retained_fraction=np.median([r['retained_fraction'] for r in rows],axis=0).tolist(),median_max_collapse_fraction=np.median([r['max_collapse_fraction'] for r in rows],axis=0).tolist(),median_topK_fraction=np.median([r['topK_fraction'] for r in rows],axis=0).tolist(),median_multiple_single_odds_shift=float(np.median([r['multiple_single_odds_shift'] for r in rows])),median_absent_single_odds_shift=float(np.median([r['absent_single_odds_shift'] for r in rows])),pruning_favors_multiple_frames=sum(r['multiple_single_odds_shift']>1e-10 for r in rows),winner_changed=sum(r['all_winner']!=r['retained_winner'] for r in rows),median_multiple_weight_all=float(np.median([r['all_class_weights'][4] for r in rows])),median_multiple_weight_retained=float(np.median([r['retained_class_weights'][4] for r in rows])),multiple_weight_ge95_all=sum(r['all_class_weights'][4]>=.95 for r in rows),multiple_weight_ge95_retained=sum(r['retained_class_weights'][4]>=.95 for r in rows))
def probe_summary(rows):
    if not rows:return dict(n=0)
    return dict(n=len(rows),median_selected_logBF=float(np.median([r['selected_logBF'] for r in rows])),median_full_logBF=float(np.median([r['full_uniform_logBF'] for r in rows])),median_log_inflation=float(np.median([r['selected_minus_full'] for r in rows])),min_log_inflation=min(r['selected_minus_full'] for r in rows),max_log_inflation=max(r['selected_minus_full'] for r in rows),selected_favors_double=sum(r['selected_logBF']>0 for r in rows),full_favors_double=sum(r['full_uniform_logBF']>0 for r in rows),selected_positive_full_negative=sum(r['selected_logBF']>0 and r['full_uniform_logBF']<0 for r in rows),median_two_term_coverage=float(np.median([np.exp(r['represented_uniform_logBF']-r['full_uniform_logBF']) for r in rows])),eligible_grid_sizes=sorted(set(r['eligible'] for r in rows)))
def run():
    start=time.perf_counter();attempt=OUT/'check_attempt.json';assert not attempt.exists();atomic_json(attempt,dict(status='RUNNING'));count=0;probecount=0;direct=0;allrows=[];allprobes=[]
    try:
        for p in sorted(OUT.glob('r*_rho*.json')):
            d=read(p);lookup={(r['replicate'],r['frame']):r for r in d['masses']};probes={(r['replicate'],r['frame']):r for r in d['probes']}
            for replicate in range(4):
                with np.load(RAW/f'{p.stem}_p{replicate}.npz') as z:
                    offsets=z['offsets'];logs=z['logs'];priors=z['priors'];kinds=z['kinds'];groups=z['groups'];retained=z['retained']
                    for t in range(24):
                        a,b=map(int,offsets[t:t+2]);v=logs[a:b];k=kinds[a:b];g=groups[a:b];keep=retained[t];assert abs(ls(priors[a:b]))<1e-9;maxima=np.full(int(g.max())+1,-np.inf);np.maximum.at(maxima,g,v);gkind=np.zeros(len(maxima),np.uint8);gkind[g]=k;assert np.all(gkind[g]==k);A=[];B=[];C=[]
                        for ci in range(5):
                            A.append(ls(v[k==ci]));B.append(ls(maxima[gkind==ci]));ids=keep[k[keep]==ci];C.append(ls(v[ids]));np.testing.assert_allclose(v[ids],maxima[g[ids]],atol=1e-12);top=np.sort(maxima[gkind==ci])[-(1 if ci==0 else 2):];np.testing.assert_allclose(np.sort(v[ids]),top,atol=1e-12)
                        row=lookup[(replicate,t)];np.testing.assert_allclose([A,B,C],[row['all_history_logmass'],row['max_history_logmass'],row['retained_logmass']],rtol=1e-11,atol=1e-9);np.testing.assert_allclose([weights(A),weights(B),weights(C)],[row['all_class_weights'],row['collapsed_class_weights'],row['retained_class_weights']],rtol=1e-10,atol=1e-10);np.testing.assert_allclose(np.array(C)-A,(np.array(B)-A)+(np.array(C)-B),atol=1e-12);np.testing.assert_allclose((C[4]-C[1])-(A[4]-A[1]),row['multiple_single_odds_shift'],atol=1e-9);count+=1;allrows.append(row)
                    for j,t in enumerate((0,10,15,23)):
                        row=probes[(replicate,t)];grid=z['grids'][j];center=row['primary_center'];distance=np.minimum((np.arange(500)-center)%500,(center-np.arange(500))%500)*.2;eligible=np.broadcast_to(distance[None,:]>6.4,grid.shape);N=int(eligible.sum());assert N==row['eligible'];selected=[grid[si,position] for position,si in row['secondary_indices']];full=ls(grid[eligible])-math.log(N);chosen=ls(selected)-math.log(2);represented=ls(selected)-math.log(N);np.testing.assert_allclose([chosen,full,represented],[row['selected_logBF'],row['full_uniform_logBF'],row['represented_uniform_logBF']],atol=1e-9);np.testing.assert_allclose(chosen-full,math.log(N/2)+represented-full,atol=1e-12);assert represented<=full+1e-10
                        e=z['residuals'][j];P=z['precision'];freq=z['freqweights'];var=row['variance']
                        for position,si in row['secondary_indices']+[(0,0)]:
                            mu=z['template_FT'][si]*np.exp(-2j*np.pi*np.arange(251)*position/500);ex,ey=e;mx,my=mu
                            linear=np.real(np.conj(ex)*(P[:,0,0]*mx+P[:,0,1]*my)+np.conj(ey)*(P[:,1,0]*mx+P[:,1,1]*my));quad=abs(mx)**2*P[:,0,0].real+abs(my)**2*P[:,1,1].real+2*np.real(np.conj(mx)*P[:,0,1]*my);value=float(math.fsum(float(f*(l-.5*q)) for f,l,q in zip(freq,linear,quad))/500/var);np.testing.assert_allclose(value,grid[si,position],rtol=1e-7,atol=2e-4);direct+=1
                        probecount+=1;allprobes.append(row)
            q=p.with_suffix('.json.gz');q.write_bytes(gzip.compress(p.read_bytes(),mtime=0));assert gzip.decompress(q.read_bytes())==p.read_bytes()
        assert binding()==read(OUT/'plan.json')['source'];strata=[]
        for when in(.1,1.):
            for rho in(0.,.95):
                rs=[r for r in allrows if r['regime']==when and r['rho']==rho and r['case'] in('static','moving','deforming')];ps=[r for r in allprobes if r['regime']==when and r['rho']==rho and r['all_scene_single']];accurate=[r for r in ps if r['primary_accurate']];split=[r for r in allprobes if r['regime']==when and r['rho']==rho and r['scene_split']];strata.append(dict(regime=when,rho=rho,intact_frame_masses=summary(rs),single_scene_probes=probe_summary(ps),accurate_primary_single_scene=probe_summary(accurate),true_split_probes=probe_summary(split)))
        atomic_json(OUT/'audit_summary.json',dict(classes=CLASSES,all_mass_summary=summary(allrows),all_probe_summary=probe_summary(allprobes),strata=strata,scope='One-step generated candidate trees conditioned on original pruned histories; conditional fixed-primary Bayes factors, not a replacement tracker or calibrated probabilities.'))
        atomic_json(OUT/'separate_check.json',dict(passed=True,prior_mass_checks=count,class_history_aggregation_checks=count,probe_integrals_checked=probecount,componentwise_secondary_logBF_checks=direct,all_source_bindings_unchanged=True,elapsed_seconds=time.perf_counter()-start,limits='Shared cache endpoint IDs/covariance/templates and same author; independent aggregation/formulas, not external replication.'))
        manifest={p.name:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in RAW.iterdir() if p.is_file()};atomic_json(OUT/'local_audit_manifest.json',dict(directory=str(RAW),files=manifest,total_bytes=sum(v['bytes'] for v in manifest.values()),full_cache_published=False));atomic_json(attempt,dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
    except BaseException:atomic_json(attempt,dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':run()
