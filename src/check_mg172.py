"""Independent full-grid/history aggregation and decision/weight replay."""
import time,json,gzip,hashlib,math,collections
import numpy as np
from run_mg172 import OUT,RAW,read,atomic_json,binding
from analyze_mg170 import CLASSES

def ls(values):
    values=np.asarray(values);m=float(values.max());return m+math.log(math.fsum(math.exp(float(v)-m) for v in values))
def inner(a,b,P,w):return float(np.sum(w*np.real(np.conj(a[0])*(P[:,0,0]*b[0]+P[:,0,1]*b[1])+np.conj(a[1])*(P[:,1,0]*b[0]+P[:,1,1]*b[1])))/500)
def reconstruct_mu(kind,k,si,second,ft):
    mu=ft[si]*np.exp(-2j*np.pi*np.arange(251)*k/500)
    if kind=='absent':return np.zeros_like(mu)
    if kind=='Xonly':mu=mu.copy();mu[1]=0
    if kind=='Yonly':mu=mu.copy();mu[0]=0
    if kind=='multiple':mu=mu+ft[second[1]]*np.exp(-2j*np.pi*np.arange(251)*second[0]/500)
    return mu

def check():
    start=time.perf_counter();attempt=OUT/'check_attempt.json';assert not attempt.exists();atomic_json(attempt,dict(status='RUNNING'));classchecks=0;framechecks=0;pathchecks=0;gridchecks=0;scalarchecks=0;qchecks=0;calibration=[]
    try:
        analysis=read(OUT/'analysis.json');lookup={}
        for p in sorted(OUT.glob('r*_rho*.json')):
            d=read(p);when=d['regime'];rho=d['rho'];case=d['case'];positive=case in('static','moving','deforming');onset=0 if case=='absent' else 10;summary=next(r for r in analysis['rows'] if r['regime']==when and r['rho']==rho and r['case']==case);ok={k:[] for k in('corrected','original')};false={k:[] for k in ok}
            for path in d['paths']:
                seed=path['seed'];previous={name:None for name in ok}
                for name,records in path['records'].items():
                    single=0;accurate=0;counts=collections.Counter();wrong_high=0;high=0;wrong_emitted=0;modeled=0
                    for t,rec in enumerate(records):
                        tr=path['truth'][t];k=(seed+(0 if case=='static' else t))%500;s=1.
                        if case=='deforming':k=(seed+t//2)%500;s=[1.,1.125,1.25,1.125][t//10]
                        if case=='jump' and t>=10:k=(k+20)%500
                        assert tr==[k,s];weights=rec['class_weights'];assert abs(sum(weights.values())-1)<1e-9;winner=max(weights,key=weights.get);assert winner==rec['winner'];label=winner if weights[winner]>=.95 and rec['innovation_q']<=rec['qcut'] else 'uncertain'
                        if label=='single' and(rec['single_center_spread']>4 or rec['single_width_spread']>.25):label='uncertain'
                        assert label==rec['label'];counts[label]+=1
                        if name=='corrected':
                            full=np.array(rec['full_class_logmass']);norm=ls(full);expected=np.exp(full-norm)
                            np.testing.assert_allclose([weights[c] for c in CLASSES],expected,rtol=1e-9,atol=1e-12)
                            for ci,kind in enumerate(CLASSES):
                                states=[b for b in rec['beam'] if b['kind']==kind];retained=ls([b['merged_log'] for b in states]);np.testing.assert_allclose(np.exp(retained-full[ci]),rec['retained_geometry_fraction'][kind],rtol=1e-9,atol=1e-12);np.testing.assert_allclose(sum(np.exp(b['log_weight']) for b in states),weights[kind],rtol=1e-9,atol=1e-12)
                                for state in states:np.testing.assert_allclose(state['log_weight'],state['merged_log']+full[ci]-retained-norm,rtol=1e-10,atol=1e-9)
                                classchecks+=1
                            best=max([b for b in rec['beam'] if b['kind']==winner],key=lambda b:b['merged_log']);assert best['center']==rec['center'] and best['width']==rec['width'];lookup[(p.stem,path['replicate'],t)]=rec
                        if previous[name] is not None:
                            for state in rec['beam']:assert any(min((state['center']-q['center'])%500,(q['center']-state['center'])%500)<=2 and abs(state['width']-q['width'])<=.12500001 for q in previous[name])
                        previous[name]=rec['beam'];single+=label=='single';accurate+=label=='single' and min((rec['center']-k)%500,(k-rec['center'])%500)<=2 and abs(rec['width']-s)<=.12500001;framechecks+=1
                        actual='single'
                        if case=='absent' or(case=='vanish' and t>=10):actual='absent'
                        elif case=='missingX' and t>=10:actual='Yonly'
                        elif case=='missingY' and t>=10:actual='Xonly'
                        elif case=='split' and t>=10:actual='multiple'
                        elif case in('displacedY','jump') and t>=10:actual=None
                        if actual is not None:
                            modeled+=1;high+=weights[winner]>=.95;wrong_high+=weights[winner]>=.95 and winner!=actual;wrong_emitted+=label not in('uncertain',actual)
                    outcome=path['outcomes'][name];assert single==outcome['single_frames'] and accurate==outcome['accurate_frames'];f=records[onset+5]['label']=='single' if not positive else None;passed=single>=22 and accurate>=22 if positive else not f;assert bool(passed)==outcome['success'];ok[name].append(passed);false[name].append(bool(f));pathchecks+=1;calibration.append(dict(regime=when,rho=rho,case=case,replicate=path['replicate'],method=name,label_counts=dict(counts),modeled_frames=modeled,high_weight_frames=int(high),high_weight_wrong_winner=int(wrong_high),wrong_emitted_class=int(wrong_emitted)))
            for name in ok:
                assert sum(ok[name])==summary[name]['successful_paths'] and summary[name]['screen_passed']==(sum(ok[name])==4)
                if not positive:assert sum(false[name])==summary[name]['false_single_deadline']
            q=p.with_suffix('.json.gz');q.write_bytes(gzip.compress(p.read_bytes(),mtime=0));assert gzip.decompress(q.read_bytes())==p.read_bytes()
        for p in sorted(RAW.glob('*_snap*.npz')):
            tag,t=p.stem.rsplit('_snap',1);t=int(t);rec=lookup[(tag,0,t)]
            with np.load(p) as z:
                ft=z['FT'];P=z['precision'];w=z['weights'];innovation=z['innovation'];rho=float(z['rho']);var=float(z['variance']);oldkind=z['old_kind'];oldcenter=z['old_center'];oldsi=z['old_si'];oldmean=z['old_mean'];oldlog=z['old_log'];Mkeys=z['multiple_keys'];Mlogs=z['multiple_logs'];norms=np.array([inner(a,a,P,w) for a in ft]);K=np.einsum('fab,kbf->kaf',P,ft);classlogs=[]
                for kind in CLASSES[:-1]:
                    values=[]
                    for (k,si),stored in zip(z[kind+'_keys'],z[kind+'_logs']):
                        mu=reconstruct_mu(kind,k,si,None,ft);v=[]
                        for ci,pk,ps,pm,pl in zip(oldkind,oldcenter,oldsi,oldmean,oldlog):
                            allowed=(k==pk and si==ps) if kind=='absent' else min((k-pk)%500,(pk-k)%500)<=2 and abs(si-ps)<=1
                            if not allowed:continue
                            err=innovation-(mu-rho*pm);q=inner(err,err,P,w)/var;score=pl+math.log(.98 if CLASSES[int(ci)]==kind else .005)-.5*q
                            if kind!='absent':score-=math.log(5*(min(13,int(ps)+2)-max(0,int(ps)-1)))
                            v.append(score)
                        value=ls(v);np.testing.assert_allclose(value,stored,rtol=1e-9,atol=2e-6);values.append(value);scalarchecks+=1
                    classlogs.append(ls(values))
                masses=[]
                for (k,si),stored in zip(Mkeys,Mlogs):
                    primary=ft[si]*np.exp(-2j*np.pi*np.arange(251)*k/500);distance=np.minimum((np.arange(500)-k)%500,(k-np.arange(500))%500)*.2;allowed=np.broadcast_to(distance[None,:]>6.4,(13,500));N=int(allowed.sum());combined=np.full((13,500),-np.inf)
                    for ci,pk,ps,pm,pl in zip(oldkind,oldcenter,oldsi,oldmean,oldlog):
                        if min((k-pk)%500,(pk-k)%500)>2 or abs(si-ps)>1:continue
                        err=innovation+rho*pm-primary;base=-.5*inner(err,err,P,w)/var;correlation=np.fft.irfft((err[None,:,:]*np.conj(K)).sum(axis=1),n=500,axis=1);ll=base+(correlation-.5*norms[:,None])/var;nloc=5*(min(13,int(ps)+2)-max(0,int(ps)-1));v=pl+math.log(.98 if CLASSES[int(ci)]=='multiple' else .005)-math.log(nloc)-math.log(N)+ll;combined=np.logaddexp(combined,np.where(allowed,v,-np.inf))
                    np.testing.assert_allclose(combined[allowed],stored[allowed],rtol=1e-9,atol=2e-6);assert np.isneginf(stored[~allowed]).all();masses.append(ls(combined[allowed]));gridchecks+=1
                classlogs.append(ls(masses));np.testing.assert_allclose(classlogs,z['classlogs'],rtol=1e-9,atol=2e-6);np.testing.assert_allclose(classlogs,rec['full_class_logmass'],rtol=1e-9,atol=2e-6)
                best=max([b for b in rec['beam'] if b['kind']==rec['winner']],key=lambda b:b['merged_log']);si=int(round((best['width']-.5)/.125));mu=reconstruct_mu(best['kind'],best['center'],si,best['second'],ft);incoming=[]
                for ci,pk,ps,pm,pl in zip(oldkind,oldcenter,oldsi,oldmean,oldlog):
                    allowed=(best['center']==pk and si==ps) if best['kind']=='absent' else min((best['center']-pk)%500,(pk-best['center'])%500)<=2 and abs(si-ps)<=1
                    if not allowed:continue
                    err=innovation-(mu-rho*pm);q=inner(err,err,P,w)/var;v=pl+math.log(.98 if CLASSES[int(ci)]==best['kind'] else .005)-.5*q
                    if best['kind']!='absent':v-=math.log(5*(min(13,int(ps)+2)-max(0,int(ps)-1)))
                    if best['kind']=='multiple':v-=math.log(5655)
                    incoming.append((v,q))
                np.testing.assert_allclose(max(incoming,key=lambda v:v[0])[1],rec['innovation_q'],rtol=1e-8,atol=2e-5);qchecks+=1
        assert binding()==read(OUT/'plan.json')['source'];atomic_json(OUT/'separate_check.json',dict(passed=True,class_mass_preservation_checks=classchecks,frame_labels_checked=framechecks,path_outcomes_checked=pathchecks,explicit_secondary_grids_checked=gridchecks,explicit_scalar_history_sums=scalarchecks,representative_quadratics_checked=qchecks,calibration_diagnostics=calibration,elapsed_seconds=time.perf_counter()-start,limits='Independent residual-FFT integration vs cached template cross-correlations; shared stored noise precision/templates, same author, not external replication.'))
        manifest={p.name:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in RAW.iterdir() if p.is_file()};atomic_json(OUT/'local_raw_manifest.json',dict(directory=str(RAW),files=manifest,total_bytes=sum(v['bytes'] for v in manifest.values()),full_fields_and_cache_published=False));atomic_json(attempt,dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
    except BaseException:atomic_json(attempt,dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':check()
