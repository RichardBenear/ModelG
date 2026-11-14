"""Independent weight/label/path replay and componentwise spectral score checks."""
import time,json,gzip,hashlib,collections
import numpy as np
from scipy.special import logsumexp
from analyze_mg170 import ROOT,OUT,INPUT,PRIOR,OLD,BENCH,Beam,binding,read,atomic_json
from mg167_xy import XYCore

def report():
    start=time.perf_counter();analysis=read(OUT/'analysis.json');bench=read(BENCH);metadata={};recordschecked=0;beamchecks=0;qchecks=0;pathschecked=0;labeltotals={};postdeadline=[]
    with np.load(OLD/'prepared_inputs.npz') as z:reference=z['reference'];initial=z['initial'];omega=float(z['omega'])
    samples=[]
    for r in(0,1):
        with np.load(OLD/'production'/f'homogeneous_r{r}'/'observations.npz') as z:samples.append(z['counts'][:,0,:3]/omega)
    def mean(model,state):
        if state['kind']=='absent':return np.zeros((2,251),complex)
        mu=model.signal(state['width'],state['center'])
        if state['kind']=='Xonly':mu[1]=0
        if state['kind']=='Yonly':mu[0]=0
        if state['kind']=='multiple' and state['second'] is not None:mu+=model.signal(.5+.125*state['second'][1],state['second'][0])
        return np.fft.rfft(mu,axis=1)
    for ri in(0,2):
        when=[.1,.5,1.][ri];model=XYCore(reference,initial,omega,np.array(bench['regime_inputs'][ri]['covariance'])[1:,1:]);data=np.array([s[round(when*1000)] for s in samples]);data-=data.mean(axis=2)[:,:,None];ft=np.fft.rfft(data,axis=2);S=np.einsum('rsf,rtf->fst',ft,np.conj(ft))/(2*500);S=np.stack([S[max(0,i-2):min(251,i+3)].mean(axis=0) for i in range(251)]);diag=np.zeros_like(S)
        for sp in range(3):diag[:,sp,sp]=S[:,sp,sp].real
        S=.9*S+.1*diag;ev,V=np.linalg.eigh(S);ev=np.maximum(ev,1e-8*ev.max(axis=1)[:,None]);L=V*np.sqrt(ev)[:,None,:];L[0]=0
        noise=Beam(model,L,0,100)
        def scalar_q(error,variance):
            ex,ey=error;P=noise.P;terms=abs(ex)**2*P[:,0,0].real+abs(ey)**2*P[:,1,1].real+2*np.real(np.conj(ex)*P[:,0,1]*ey);return float(sum(float(a*b) for a,b in zip(terms,noise.freqweights))/500/variance)
        for row in [r for r in analysis['rows'] if r['regime']==when]:
            rho=row['rho'];case=row['case'];tag=f'r{ri}_rho{rho:g}_{case}';d=read(OUT/(tag+'.json'));oldpath=PRIOR/(tag+'.json');published=PRIOR/(tag+'.json.gz');assert gzip.decompress(published.read_bytes())==oldpath.read_bytes();metadata[str(oldpath.relative_to(ROOT))]=dict(sha256=hashlib.sha256(oldpath.read_bytes()).hexdigest(),published_gzip_sha256=hashlib.sha256(published.read_bytes()).hexdigest(),matches_previously_published=True);olddata=read(oldpath)
            with np.load(INPUT/(tag+'.npz')) as z:fields=z['fields'][:4,:24]
            success=[];oldsuccess=[];false=[];allcounts=collections.Counter();positive=case in('static','moving','deforming');onset=0 if case=='absent' else 10;latefalse=0
            for path in d['paths']:
                r=path['replicate'];single=0;accurate=0;assert path['seed']==olddata['paths'][r]['initial_seed'];previousbeam=[dict(kind='single',center=path['seed'],width=1.,second=None)]
                for t,rec in enumerate(path['records']):
                    c=path['seed']+(0 if case=='static' else t);width=1.
                    if case=='deforming':c=path['seed']+t//2;width=[1.,1.125,1.25,1.125][t//10]
                    if case=='jump' and t>=10:c+=20
                    assert path['truth'][t]==[c%500,width]
                    weights=np.exp([x['log_weight'] for x in rec['beam']]);np.testing.assert_allclose(weights.sum(),1,rtol=1e-10,atol=1e-10);cw={k:float(sum(w for w,s in zip(weights,rec['beam']) if s['kind']==k)) for k in rec['class_weights']}
                    for k in cw:np.testing.assert_allclose(cw[k],rec['class_weights'][k],rtol=1e-10,atol=1e-10)
                    winner=max(cw,key=cw.get);assert winner==rec['winner'];best=max([b for b in rec['beam'] if b['kind']==winner],key=lambda b:b['log_weight']);assert best['center']==rec['center'] and best['width']==rec['width'];label=winner if cw[winner]>=.95 and rec['innovation_q']<=rec['qcut'] else 'uncertain';one=[s for s in rec['beam'] if s['kind']=='single'];spread=max(min((a['center']-b['center'])%500,(b['center']-a['center'])%500) for a in one for b in one);ws=max(abs(a['width']-b['width']) for a in one for b in one)
                    if label=='single' and(spread>4 or ws>.25):label='uncertain'
                    assert label==rec['label'];allcounts[label]+=1
                    for state in rec['beam']:
                        assert .5<=state['width']<=2 and any(min((state['center']-parent['center'])%500,(parent['center']-state['center'])%500)<=2 and abs(state['width']-parent['width'])<=.12500001 for parent in previousbeam)
                        if state['second'] is not None:assert min((state['second'][0]-state['center'])%500,(state['center']-state['second'][0])%500)*.2>model.x.background_exclusion
                        beamchecks+=1
                    if r==0 and t in(0,10,15,23):
                        field=fields[r,t];y=np.fft.rfft(field-field.mean(axis=1)[:,None],axis=1);mu=mean(model,best)
                        if t==0:poss=[scalar_q(y-mu,1)]
                        else:
                            priorfield=fields[r,t-1];yprev=np.fft.rfft(priorfield-priorfield.mean(axis=1)[:,None],axis=1);poss=[scalar_q(y-rho*yprev-(mu-rho*mean(model,parent)),1-rho*rho) for parent in previousbeam if (parent['center']==best['center'] and parent['width']==best['width']) or (best['kind']!='absent' and min((parent['center']-best['center'])%500,(best['center']-parent['center'])%500)<=2 and abs(parent['width']-best['width'])<=.12500001)]
                        assert min(abs(q-rec['innovation_q']) for q in poss)<2e-4;qchecks+=1
                    if label=='single':single+=1;accurate+=min((rec['center']-c)%500,(c-rec['center'])%500)<=2 and abs(rec['width']-width)<=.12500001
                    if not positive and t>=onset+5:latefalse+=label=='single'
                    previousbeam=rec['beam'];recordschecked+=1
                assert single==path['single_frames'] and accurate==path['accurate_frames'];f=path['records'][onset+5]['label']=='single' if not positive else None;ok=single>=22 and accurate>=22 if positive else not f;assert bool(ok)==path['success'];success.append(ok);false.append(bool(f))
                oldrec=olddata['paths'][r]['records']['local'][:24];n=sum(x['accepted'] for x in oldrec);acc=sum(x['accepted'] and x['center'] is not None and min((x['center']-tr[0])%500,(tr[0]-x['center'])%500)<=2 and abs(x['scale']-tr[1])<=.12500001 for x,tr in zip(oldrec,path['truth']));oldok=n>=22 and acc>=22 if positive else not oldrec[onset+5]['accepted'];assert bool(oldok)==path['old_success'];oldsuccess.append(oldok);pathschecked+=1
            assert sum(success)==row['successful_paths'] and sum(oldsuccess)==row['old_successful_paths'];assert row['screen_passed']==(sum(success)==4)
            if not positive:assert sum(false)==row['false_single_deadline']
            labeltotals[f'{when}_{rho}_{case}']=dict(allcounts);postdeadline.append(dict(regime=when,rho=rho,case=case,postdeadline_single_frames=latefalse if not positive else None))
    assert binding()==read(OUT/'plan.json')['source'];atomic_json(OUT/'separate_check.json',dict(passed=True,frame_labels_checked=recordschecked,beam_motion_and_weight_checks=beamchecks,independent_componentwise_spectral_q_checks=qchecks,path_outcomes_checked=pathschecked,label_totals=labeltotals,postdeadline=postdeadline,limits='Same author/shared reference and noise covariance construction; weights are not calibrated posteriors.'))
    import subprocess
    case_summary=PRIOR/'analysis.json';published_summary=subprocess.check_output(['git','show','3bbd570:'+str(case_summary.relative_to(ROOT)).replace('\\','/')],cwd=ROOT);assert published_summary.replace(b'\r\n',b'\n')==case_summary.read_bytes().replace(b'\r\n',b'\n');metadata[str(case_summary.relative_to(ROOT))]=dict(sha256=hashlib.sha256(case_summary.read_bytes()).hexdigest(),matches_pre_registration_Git_commit='3bbd570')
    atomic_json(OUT/'supplemental_metadata_manifest.json',dict(files=metadata,note='MG169 metadata seed/truth/comparator files were omitted from prospective input list. Their byte identity against previously published lossless gzip was verified post-analysis; independent truth/outcome replay passed. Original plan is not rewritten.'))
    for p in OUT.glob('r*_rho*.json'):
        q=p.with_suffix('.json.gz');q.write_bytes(gzip.compress(p.read_bytes(),mtime=0));assert gzip.decompress(q.read_bytes())==p.read_bytes()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(10,4.5),layout='constrained');conditions=[(.1,0.),(.1,.95),(1.,0.),(1.,.95)];bottom=np.zeros(4)
    for label,color in [('single','tab:blue'),('uncertain','tab:orange'),('multiple','tab:red'),('absent','gray'),('Xonly','tab:purple'),('Yonly','tab:green')]:
        values=[]
        for when,rho in conditions:values.append(sum(labeltotals[f'{when}_{rho}_{case}'].get(label,0) for case in('static','moving','deforming'))/288)
        ax.bar(np.arange(4),values,bottom=bottom,label=label,color=color);bottom+=values
    ax.set_xticks(np.arange(4),['Early / rho0','Early / rho0.95','Late / rho0','Late / rho0.95']);ax.set_ylabel('Label fraction on known intact single-core observations');ax.set_title('MG170: truth is single for every observation shown');ax.legend(ncol=3,fontsize=8);fig.savefig(OUT/'MG170_false_hypotheses.png',dpi=160);plt.close(fig)
    lines=['# MG170 - uncertainty and competing-hypothesis tracking diagnostic','',
    'Explicit alternative hypotheses did not produce a reliable tracker in this development test. The prototype handles early-noise examples well, but under late independent noise it wrongly labels nearly every intact single-core observation as multiple cores. Under late correlated noise it also mistakes some intact cores for absent ones. Its high normalized weights can therefore accompany a known wrong answer: adding an uncertainty label and competing models does not automatically calibrate uncertainty. The overall screen failed, and several late positive cases perform worse than the old local tracker on the same short paths. This is an inference-method failure on planted examples, not evidence of particle death or multiplication. G remains unresolved and MG163 is unchanged.','',
    '![Labels on known intact single cores](MG170_false_hypotheses.png)','','## Frozen prototype and scope','',
    'Five hypotheses are considered: no coherent core, one X/Y core, X-template-only, Y-template-only, and two X/Y cores. These concern template structure, not disappearance of all species tokens. The beam retains up to two histories/geometries per nonabsent class and one absent history, with bounded primary motion/width changes and no forced three-failure termination. A global secondary-core discovery step proposes competing structures. Secondary discovery does not establish an independently tracked second identity.','',
    'The procedure uses the full spatial/cross-species conditional Fourier covariance and the supplied synthetic ARcorrelation in innovation likelihoods. Ordinary amplitude1is known. Those are additional assumptions/information relative to MG169; any gain cannot be attributed solely to retaining hypotheses. Transition odds, model-weight threshold0.95, geometry-spread guard and conditional chi-square residual check were fixed before evaluation. Weights are normalized over data-selected, pruned proposals; they are not calibrated probabilities or physical confidence levels.','',
    'No new random fields or Model G trajectories were generated. The first four previously seen MG169paths percase and their first24observations are used, for early/late regimes0.1/1 andrho0/0.95:160sequences,3840frames. Initial geometry is supplied to every case, including no-core controls; the tracker receives no injected future position, width or case label. All comparisons are small development diagnostics, not independent real-noise validation.','',
    '## Intact-path results','',
    'Positive success requires>=22/24single labels and>=22/24estimates within2sites and0.125width of injected truth. Four of four paths are needed for the raw90%screen; this denominator cannot certify a population recall rate. The old local comparator is recomputed on the same four/24frame subset, not copied from its original16path/40frame verdict.','',
    '| Regime | rho | Static new/old successes | Moving new/old | Deforming new/old |','|---|---:|---:|---:|---:|']
    for when in(.1,1.):
        for rho in(0.,.95):
            cells=[]
            for case in('static','moving','deforming'):
                row=next(r for r in analysis['rows'] if r['regime']==when and r['rho']==rho and r['case']==case);cells.append(f"{row['successful_paths']}/4 / {row['old_successful_paths']}/4")
            lines.append(f'| {when:g} | {rho:g} | '+' | '.join(cells)+' |')
    lines+=['',
    'Early examples meet the intact and negative endpoint screens in this small subset. Under late independent noise, all12intact paths fail:276of288frames are wrongly labeled multiple, the other12uncertain, and none single. Under late correlated noise, the288intact frames receive175single,72absent,37uncertain and4multiple labels. Only5of12intact paths meet the whole-path criterion. The absent/multiple alternatives are demonstrably wrong on those planted intact frames, regardless of their numerical weight.','',
    '## Damage/no-core endpoint','',
    'A negative fails if the output is single at onset+5(frame5for initially absent,frame15for damaged cases). All four paths must avoid that false-single deadline. This endpoint is not identical to the original MG169 alive/stop deadline; the comparator here uses its accepted-frame label at the same observation. Full post-deadline histories are retained because a single deadline can hide later false continuation.','',
    '| Regime | rho | Absent | Vanish | X removed | Y removed | Y displaced | Split | Jump |','|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for when in(.1,1.):
        for rho in(0.,.95):
            cells=[next(r['false_single_deadline'] for r in analysis['rows'] if r['regime']==when and r['rho']==rho and r['case']==c) for c in('absent','vanish','missingX','missingY','displacedY','split','jump')];lines.append(f'| {when:g} | {rho:g} | '+' | '.join(f'{v}/4' for v in cells)+' |')
    lines+=['',
    'Most negative endpoints avoid a single label, but this is not evidence of a generally improved instrument when intact positives are also called absent/multiple. The late correlated split case retains one false-single deadline. No damaged/no-core result overrides the failed intact screen.','',
    '## What this failure diagnoses','',
    'The global secondary peaks are chosen using the same observation being scored. Assigning weight to a few best proposals is not the same as integrating over the full secondary-position/width prior. Beam pruning also changes retained class mass. These are plausible contributors to the spurious multiple-core preference and overconfident normalized weights; this study does not isolate their individual causal contributions. Model-relative weights can be high even when the entire comparison is poorly calibrated. More hypotheses alone are not a substitute for testing false assignments on known single cores.','',
    '## Independent checks, provenance and resources','',
    f"Independent replay reproduced {recordschecked}frame labels, {beamchecks}beam motion/weight checks and {pathschecked}path outcomes/comparator counts. Componentwise real/imaginary spectral sums reproduced {qchecks}selected innovation quadratics, with valid parent histories. Deterministic Parseval and saved pure-noise scale checks preceded classification. Same author/reference/noise inputs remain dependencies; this is not external replication.",'',
    'Prospective source and primary raw input hashes were published and upload confirmed before evaluation. A provenance omission was found during review: MG169seed/truth/comparator JSON metadata and case-denominator summary were not listed in the prospective manifest. Their byte-for-byte agreement with previously published lossless gzip exports was verified afterward, and true motion plus comparator outcomes were independently reconstructed. Supplemental metadata hashes are retained without rewriting the original plan. No outcomes were retuned.','',
    f"Evaluation took {analysis['elapsed_seconds']:.6g}seconds on one worker/library-thread1 within600active seconds. Full original fields remain in local_raw/MG169_tracking; no new full-field archive was created. Public compact beam/label/score/truth records are losslessly gzipped. G is UNRESOLVED; no persistence verdict or MG163rescore.",'',
    '## Recommended next step','',
    'Audit the single-versus-multiple model-selection calculation before building another tracker: compare proposal selection with integration over an explicit spatial/width prior or an independently calibrated search-volume penalty. Include known single-core and no-core cases specifically to measure false-multiple/false-absent assignments. Reassess beam-mass truncation and probability calibration separately from trajectory tracking, then test any correction on independent validation realizations. Do not loosen the0.95cutoff or reinterpret these high weights as physical evidence. No new tracker revision or kinetic run is authorized by this report.','']
    text='\n'.join(lines)
    for old,new in [('ARcorrelation','AR correlation'),('amplitude1is','amplitude 1 is'),('threshold0.95','threshold 0.95'),('MG169paths percase','MG169 paths per case'),('first24observations','first 24 observations'),('regimes0.1/1 andrho0/0.95:160sequences,3840frames','regimes 0.1/1 and rho0/0.95: 160 sequences, 3840 frames'),('requires>=22/24single','requires >=22/24 single'),('and>=22/24estimates within2sites and0.125width','and >=22/24 estimates within 2 sites and 0.125 width'),('raw90%screen','raw 90% screen'),('four/24frame','four-path/24-frame'),('original16path/40frame','original 16-path/40-frame'),('all12intact','all 12 intact'),('276of288frames','276 of 288 frames'),('other12uncertain','other 12 uncertain'),('the288intact','the 288 intact'),('175single,72absent,37uncertain and4multiple','175 single, 72 absent, 37 uncertain and 4 multiple'),('Only5of12intact','Only 5 of 12 intact'),('onset+5(frame5for','onset+5 (frame 5 for'),('absent,frame15for','absent, frame 15 for'),('}frame','} frame'),('}beam','} beam'),('}path','} path'),('}selected','} selected'),('MG169seed','MG169 seed'),('}seconds','} seconds'),('within600active','within 600 active'),('MG163rescore','MG163 rescore'),('the0.95cutoff','the 0.95 cutoff')]:text=text.replace(old,new)
    (OUT/'MG170_hypothesis_tracking_results.md').write_text(text,encoding='utf-8');atomic_json(OUT/'report_attempt_v2.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
if __name__=='__main__':report()
