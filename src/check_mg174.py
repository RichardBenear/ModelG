"""Independent structural/motion mass replay and explicit Gaussian history sums."""
import time,math,json,gzip,hashlib,collections
import numpy as np
from run_mg174 import OUT,RAW,read,atomic_json,binding,POSITIVES,NEGATIVES,SENS,geometry
from mg174_modes import PAIRS,MOTIONS
from analyze_mg170 import CLASSES

def ls(a):
    a=np.asarray(a);v=float(a.max());return v+math.log(math.fsum(math.exp(float(x)-v) for x in a))
def inner(a,b,P,w):return float(np.sum(w*np.real(np.conj(a[0])*(P[:,0,0]*b[0]+P[:,0,1]*b[1])+np.conj(a[1])*(P[:,1,0]*b[0]+P[:,1,1]*b[1])))/500)
def mu(kind,k,si,second,ft):
    f=ft[si]*np.exp(-2j*np.pi*np.arange(251)*k/500)
    if kind=='absent':return np.zeros_like(f)
    if kind=='Xonly':f=f.copy();f[1]=0
    if kind=='Yonly':f=f.copy();f[0]=0
    if kind=='multiple':f=f+ft[second[1]]*np.exp(-2j*np.pi*np.arange(251)*second[0]/500)
    return f

def check():
    start=time.perf_counter();attempt=OUT/'check_attempt.json';assert not attempt.exists();atomic_json(attempt,dict(status='RUNNING'));pairs_checked=0;frames=0;paths_checked=0;grids=0;scalars=0;qchecked=0;lookup={};diagnostics=[];modes=[]
    try:
        analysis=read(OUT/'analysis.json')
        for p in sorted(OUT.glob('r*_rho*.json')):
            d=read(p);case=d['case'];T=d['frames'];positive=case in POSITIVES or not d['primary'];onset=0 if case=='absent' else 10;row=next(r for r in analysis['rows'] if r['regime']==d['regime'] and r['rho']==d['rho'] and r['case']==case);success={k:[] for k in('modes','MG172')};false={k:[] for k in success}
            for path in d['paths']:
                for name,records in path['records'].items():
                    singles=0;accurate=0;previous=None;counts=collections.Counter();wronghigh=0;modeled=0;high=0
                    for t,r in enumerate(records):
                        truth=path['truth'][t]
                        # Independently reconstruct injected stop-start and sensitivity geometry.
                        basecase='static' if case in('weak_static','noise_scaled') else 'moving' if case=='weak_moving' else case;k=(path['seed']+(0 if basecase=='static' else t))%500;s=1.
                        if case=='deforming':k=(path['seed']+t//2)%500;s=[1.,1.125,1.25,1.125][t//10]
                        if case=='jump' and t>=10:k=(k+20)%500
                        if case=='stop_start':k=(path['seed']+min(max(t-7,0),8))%500
                        assert truth==[k,s];cw=r['class_weights'];assert abs(sum(cw.values())-1)<1e-9;winner=max(cw,key=cw.get);assert winner==r['winner'];label=winner if cw[winner]>=.95 and r['innovation_q']<=r['qcut'] else 'uncertain'
                        if label=='single' and(r['single_center_spread']>4 or r['single_width_spread']>.25):label='uncertain'
                        assert label==r['label'];counts[label]+=1
                        if name=='modes':
                            full=np.array(r['full_pair_logmass']);Z=ls(full);expected=np.exp(full-Z);np.testing.assert_allclose([r['pair_weights'][k+'/'+m] for k,m in PAIRS],expected,rtol=1e-9,atol=1e-12)
                            for ci,kind in enumerate(CLASSES):np.testing.assert_allclose(cw[kind],expected[2*ci:2*ci+2].sum(),rtol=1e-9,atol=1e-12)
                            singleZ=ls(full[2:4]);np.testing.assert_allclose([r['single_motion_weights'][m] for m in MOTIONS],np.exp(full[2:4]-singleZ),rtol=1e-9,atol=1e-12)
                            for i,(kind,motion) in enumerate(PAIRS):
                                states=[b for b in r['beam'] if b['kind']==kind and b['motion']==motion];selected=ls([b['merged_log'] for b in states]);np.testing.assert_allclose(np.exp(selected-full[i]),r['retained_pair_geometry_fraction'][kind+'/'+motion],rtol=1e-9,atol=1e-12);np.testing.assert_allclose(sum(np.exp(b['log_weight']) for b in states),expected[i],rtol=1e-9,atol=1e-12)
                                for b in states:np.testing.assert_allclose(b['log_weight'],b['merged_log']+full[i]-selected-Z,rtol=1e-10,atol=1e-9)
                                pairs_checked+=1
                            best=max([b for b in r['beam'] if b['kind']==winner],key=lambda b:b['merged_log']);assert best['center']==r['center'] and best['width']==r['width'] and best['motion']==r['representative_motion'];lookup[(p.stem,path['replicate'],t)]=r
                        if previous is not None:
                            for b in r['beam']:
                                if name=='modes' and(b['kind']=='absent' or b['motion']=='S'):assert any(b['center']==a['center'] and b['width']==a['width'] for a in previous)
                                else:assert any(min((b['center']-a['center'])%500,(a['center']-b['center'])%500)<=2 and abs(b['width']-a['width'])<=.12500001 for a in previous)
                        previous=r['beam'];singles+=label=='single';accurate+=label=='single' and min((r['center']-k)%500,(k-r['center'])%500)<=2 and abs(r['width']-s)<=.12500001;frames+=1
                        actual='single'
                        if case=='absent' or(case=='vanish' and t>=10):actual='absent'
                        elif case=='missingX' and t>=10:actual='Yonly'
                        elif case=='missingY' and t>=10:actual='Xonly'
                        elif case=='split' and t>=10:actual='multiple'
                        elif case in('displacedY','jump') and t>=10:actual=None
                        if actual is not None:modeled+=1;high+=cw[winner]>=.95;wronghigh+=cw[winner]>=.95 and winner!=actual
                    o=path['outcomes'][name];assert singles==o['single_frames'] and accurate==o['accurate_frames'];f=records[onset+5]['label']=='single' if not positive else None;ok=singles>=math.ceil(.9*T) and accurate>=math.ceil(.9*T) if positive else not f;assert bool(ok)==o['success'];success[name].append(ok);false[name].append(bool(f));paths_checked+=1;diagnostics.append(dict(regime=d['regime'],rho=d['rho'],case=case,primary=d['primary'],replicate=path['replicate'],method=name,label_counts=dict(counts),modeled_frames=modeled,high_weight_frames=int(high),wrong_high_weight_winner=int(wronghigh)))
                if d['primary'] and case in POSITIVES:
                    rec=path['records']['modes'];stats=[]
                    for a,b in((0,8),(8,16),(16,24)):
                        selected=[x['single_motion_weights']['S'] for x in rec[a:b] if x['label']=='single'];stats.append(dict(single_frames=len(selected),mean_stationary_given_single=float(np.mean(selected)) if selected else None))
                    modes.append(dict(regime=d['regime'],rho=d['rho'],case=case,replicate=path['replicate'],segments=stats))
            for name in success:
                assert sum(success[name])==row[name]['successful_paths'] and row[name]['screen_passed']==(sum(success[name])==d['n'])
                if not positive:assert sum(false[name])==row[name]['false_single_deadline']
            q=p.with_suffix('.json.gz');q.write_bytes(gzip.compress(p.read_bytes(),mtime=0));assert gzip.decompress(q.read_bytes())==p.read_bytes()
        for p in sorted(RAW.glob('*_snap*.npz')):
            tag,tt=p.stem.rsplit('_snap',1);t=int(tt);r=lookup[(tag,0,t)]
            with np.load(p) as z:
                ft=z['FT'];P=z['precision'];w=z['weights'];innovation=z['innovation'];rho=float(z['rho']);var=float(z['variance']);old=list(zip(z['old_kind'],z['old_motion'],z['old_center'],z['old_si'],z['old_mean'],z['old_log']));norms=np.array([inner(f,f,P,w) for f in ft]);K=np.einsum('fab,kbf->kaf',P,ft);pairlogs=[]
                for kind,motion in PAIRS:
                    stem=kind+'_'+motion;values=[]
                    for key,stored in zip(z[stem+'_keys'],z[stem+'_logs']):
                        k,si=map(int,key)
                        if kind!='multiple':
                            current=mu(kind,k,si,None,ft);v=[]
                            for ci,mi,pk,ps,pm,pl in old:
                                allowed=(k==pk and si==ps) if kind=='absent' or motion=='S' else min((k-pk)%500,(pk-k)%500)<=2 and abs(si-ps)<=1
                                if not allowed:continue
                                err=innovation-(current-rho*pm);q=inner(err,err,P,w)/var;v.append(pl+math.log(.98 if CLASSES[int(ci)]==kind else .005)+math.log(.98 if MOTIONS[int(mi)]==motion else .02)-.5*q-(math.log(5*(min(13,int(ps)+2)-max(0,int(ps)-1))) if kind!='absent' and motion=='M' else 0.))
                            value=ls(v);np.testing.assert_allclose(value,stored,rtol=1e-9,atol=2e-6);values.append(value);scalars+=1
                        else:
                            primary=mu('single',k,si,None,ft);distance=np.minimum((np.arange(500)-k)%500,(k-np.arange(500))%500)*.2;allowed_grid=np.broadcast_to(distance[None,:]>6.4,(13,500));N=int(allowed_grid.sum());combined=np.full((13,500),-np.inf)
                            for ci,mi,pk,ps,pm,pl in old:
                                allowed=(k==pk and si==ps) if motion=='S' else min((k-pk)%500,(pk-k)%500)<=2 and abs(si-ps)<=1
                                if not allowed:continue
                                e=innovation+rho*pm-primary;base=-.5*inner(e,e,P,w)/var;C=np.fft.irfft((e[None,:,:]*np.conj(K)).sum(axis=1),n=500,axis=1);nloc=1 if motion=='S' else 5*(min(13,int(ps)+2)-max(0,int(ps)-1));score=pl+math.log(.98 if CLASSES[int(ci)]=='multiple' else .005)+math.log(.98 if MOTIONS[int(mi)]==motion else .02)-math.log(nloc)-math.log(N)+base+(C-.5*norms[:,None])/var;combined=np.logaddexp(combined,np.where(allowed_grid,score,-np.inf))
                            np.testing.assert_allclose(combined[allowed_grid],stored[allowed_grid],rtol=1e-9,atol=2e-6);assert np.isneginf(stored[~allowed_grid]).all();values.append(ls(combined[allowed_grid]));grids+=1
                    pairlogs.append(ls(values))
                np.testing.assert_allclose(pairlogs,z['pairlogs'],rtol=1e-9,atol=2e-6);np.testing.assert_allclose(pairlogs,r['full_pair_logmass'],rtol=1e-9,atol=2e-6)
                best=max([b for b in r['beam'] if b['kind']==r['winner']],key=lambda b:b['merged_log']);si=int(round((best['width']-.5)/.125));current=mu(best['kind'],best['center'],si,best['second'],ft);incoming=[]
                for ci,mi,pk,ps,pm,pl in old:
                    allowed=(best['center']==pk and si==ps) if best['kind']=='absent' or best['motion']=='S' else min((best['center']-pk)%500,(pk-best['center'])%500)<=2 and abs(si-ps)<=1
                    if not allowed:continue
                    e=innovation-(current-rho*pm);q=inner(e,e,P,w)/var;score=pl+math.log(.98 if CLASSES[int(ci)]==best['kind'] else .005)+math.log(.98 if MOTIONS[int(mi)]==best['motion'] else .02)-.5*q
                    if best['kind']!='absent' and best['motion']=='M':score-=math.log(5*(min(13,int(ps)+2)-max(0,int(ps)-1)))
                    if best['kind']=='multiple':score-=math.log(5655)
                    incoming.append((score,q))
                np.testing.assert_allclose(max(incoming,key=lambda v:v[0])[1],r['innovation_q'],rtol=1e-8,atol=2e-5);qchecked+=1
        assert binding()==read(OUT/'plan.json')['source'];atomic_json(OUT/'separate_check.json',dict(passed=True,class_mode_mass_checks=pairs_checked,frame_labels_checked=frames,path_outcomes_checked=paths_checked,explicit_secondary_grids_checked=grids,explicit_scalar_history_sums=scalars,representative_quadratics_checked=qchecked,calibration_diagnostics=diagnostics,single_mode_segments=modes,elapsed_seconds=time.perf_counter()-start,limits='Independent residual-FFT/explicit history sums vs cached cross-correlations; shared stored precision/templates, same author, not external replication.'))
        manifest={p.name:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in RAW.iterdir() if p.is_file()};atomic_json(OUT/'local_raw_manifest.json',dict(directory=str(RAW),files=manifest,total_bytes=sum(v['bytes'] for v in manifest.values()),full_fields_cache_published=False));atomic_json(attempt,dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
    except BaseException:atomic_json(attempt,dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':check()
