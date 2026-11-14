"""Independent oracle closed forms and explicit sparse transition-node sums."""
import time,math,json,gzip,hashlib
import numpy as np
from audit_mg173 import OUT,RAW,INPUT,PRIOR,read,atomic_json,binding,MODES

def ls(values):
    values=np.asarray(values);m=float(values.max());return m+math.log(math.fsum(math.exp(float(v)-m) for v in values))
def inner(a,b,P,w):return float(math.fsum(float(z) for z in w*np.real(np.conj(a[0])*(P[:,0,0]*b[0]+P[:,0,1]*b[1])+np.conj(a[1])*(P[:,1,0]*b[0]+P[:,1,1]*b[1])))/500)
def mu(ft,si,k):return ft[si]*np.exp(-2j*np.pi*np.arange(251)*k/500)
def run():
    start=time.perf_counter();attempt=OUT/'check_attempt.json';assert not attempt.exists();atomic_json(attempt,dict(status='RUNNING'));oraclechecks=0;hmmchecks=0;nodes=0;snapshots=0;compressionchecks=0;summarychecks=0;closedform=[]
    try:
        data={};curve_lookup={};summ=read(OUT/'analysis.json')
        for p in sorted(OUT.glob('r*_rho*.json')):
            tag=p.stem;d=read(p);parts=tag.split('_');ri=int(parts[0][1:]);rho=d['rho'];case=d['case'];fields_path=INPUT/(tag+'.npz')
            with np.load(fields_path) as z:fields=z['fields']
            with np.load(INPUT/f'r{ri}_rho0_static_snap0.npz') as z:ft=z['FT'];P=z['precision'];w=z['weights']
            data[tag]=(fields,ft,P,w);old=json.loads(gzip.decompress((PRIOR/(tag+'.json.gz')).read_bytes()))
            for path in d['paths']:
                r=path['replicate'];seed=path['seed'];signal=mu(ft,4,seed);normmu=inner(signal,signal,P,w);ff=np.fft.rfft(fields[r]-fields[r].mean(axis=2)[:,:,None],axis=2);projection=np.array([inner(y,signal,P,w) for y in ff])
                for name,assumed in [('oracle',rho),('independence_assumed_oracle',0.)]:
                    logs=[]
                    for n in range(1,25):
                        effective=1+(n-1)*(1-assumed)/(1+assumed);weighted=projection[0] if n==1 else (projection[0]+projection[n-1]+(1-assumed)*projection[1:n-1].sum())/(1+assumed);logs.append(float(weighted-.5*normmu*effective))
                    np.testing.assert_allclose(logs,[x['logBF'] for x in path[name]['curve']],rtol=1e-9,atol=2e-6);oraclechecks+=24
                    D=normmu*(1+23*(1-assumed)/(1+assumed));np.testing.assert_allclose([D/2,np.sqrt(D)],[path[name]['expected_logBF_static'],path[name]['sd_logBF_static']],rtol=1e-10)
                # Independent 2-state forward recursion with direct componentwise quadratics.
                logits=np.log([.5,.5]);prob=[]
                for t in range(24):
                    rr=0 if t==0 else rho;v=1 if t==0 else 1-rr*rr;innovation=ff[t] if t==0 else ff[t]-rr*ff[t-1];out=[]
                    for current in(0,1):
                        terms=[]
                        for previous in(0,1):
                            e=innovation-(current-rr*previous)*signal;q=inner(e,e,P,w)/v;terms.append(logits[previous]+math.log((.98 if previous==current else .005)/.985)-.5*q)
                        out.append(ls(terms))
                    logits=np.array(out)-ls(out);prob.append(float(np.exp(logits[1])))
                np.testing.assert_allclose(prob,path['fixed_switching_hmm'],rtol=1e-8,atol=1e-10);np.testing.assert_allclose(prob,[x['present'] for x in path['curves']['fixed_exact']],rtol=1e-8,atol=1e-10);np.testing.assert_allclose(prob,[x['present'] for x in path['curves']['fixed_compressed']],rtol=1e-8,atol=1e-10);hmmchecks+=24
                for method,curves in path['curves'].items():
                    curve_lookup[(tag,method)]=d['paths'][0]['curves'][method]
                    for point in curves:assert abs(point['present']+point['absent']-1)<1e-10
                baseline=old['paths'][r]['records']['corrected'];bflog=np.array([v['full_class_logmass'][1]-v['full_class_logmass'][0] for v in baseline]);restricted=np.array([v['logodds'] for v in path['curves']['uniform_compressed']]);closedform.append(dict(regime=d['regime'],rho=rho,case=case,replicate=r,MG172_vs_restricted_max_logodds_difference=float(np.max(abs(bflog-restricted))),oracle_final_logBF=path['oracle']['curve'][-1]['logBF']))
            # Recompute all published final/threshold-count summaries.
            for row in [r for r in summ['rows'] if r['regime']==d['regime'] and r['rho']==rho and r['case']==case]:
                vals=[];counts=[]
                for path in d['paths']:
                    method=row['method']
                    if method in path['curves']:v=[x['present'] for x in path['curves'][method]]
                    elif method in('oracle','independence_assumed_oracle'):v=[x['present'] for x in path[method]['curve']]
                    elif method=='fixed_switching_hmm':v=path[method]
                    else:v=path['MG172_present_weights']
                    vals.append(v[-1]);counts.append(sum(x>=.95 for x in v))
                np.testing.assert_array_equal(vals,row['final_present_weights']);assert counts==row['frames_present_ge95'];summarychecks+=1
            q=p.with_suffix('.json.gz');q.write_bytes(gzip.compress(p.read_bytes(),mtime=0));assert gzip.decompress(q.read_bytes())==p.read_bytes()
        for p in sorted(RAW.glob('*.npz')):
            prefix,tt=p.stem.rsplit('_t',1);t=int(tt);parts=prefix.split('_');tag='_'.join(parts[:3]);method='_'.join(parts[3:]);fields,ft,P,w=data[tag]
            with np.load(p) as z:
                before=z['old'];new=z['new'];post=z['post'];alpha=float(z['alpha']);rr=float(z['rho']);var=float(z['variance']);seed=int(z['seed']);field=fields[0,t];y=np.fft.rfft(field-field.mean(axis=1)[:,None],axis=1)
                if t:
                    prevfield=fields[0,t-1];yprev=np.fft.rfft(prevfield-prevfield.mean(axis=1)[:,None],axis=1);innovation=y-rr*yprev
                else:innovation=y
                np.testing.assert_allclose(innovation,z['innovation'],rtol=1e-12,atol=1e-12)
                for kind in(0,1):
                    finite=np.flatnonzero(np.isfinite(new[kind].ravel()));chosen=list(finite[np.argsort(new[kind].ravel()[finite])[-2:]])
                    if np.isfinite(new[kind,4,seed]):chosen.append(4*500+seed)
                    for flat in sorted(set(chosen)):
                        si,k=np.unravel_index(int(flat),(13,500));current=np.zeros((2,251),complex) if kind==0 else mu(ft,si,k);terms=[]
                        for previous in(0,1):
                            if kind==0:parent_geometry=[(si,k,1.)]
                            else:
                                parent_geometry=[]
                                for ps in range(max(0,si-1),min(13,si+2)):
                                    nloc=5*(min(13,ps+2)-max(0,ps-1))
                                    for displacement in range(-2,3):
                                        probability=(1-alpha)/nloc+(alpha if ps==si and displacement==0 else 0.)
                                        if probability:parent_geometry.append((ps,(k-displacement)%500,probability))
                            for ps,pk,probability in parent_geometry:
                                prior=before[previous,ps,pk]
                                if not np.isfinite(prior):continue
                                past=np.zeros((2,251),complex) if previous==0 else mu(ft,ps,pk);err=innovation-(current-rr*past);quadratic=inner(err,err,P,w)/var;terms.append(prior+math.log((.98 if kind==previous else .005)/.985)+math.log(probability)-.5*quadratic)
                        value=ls(terms);np.testing.assert_allclose(value,new[kind,si,k],rtol=1e-9,atol=2e-6);nodes+=1
                masses=np.array([ls(new[k].ravel()) for k in(0,1)]);total=ls(masses);weights=np.exp(masses-total)
                for kind in(0,1):
                    np.testing.assert_allclose(np.exp(ls(post[kind].ravel())),weights[kind],rtol=1e-9,atol=1e-12)
                    if bool(z['compressed']):
                        live=post[kind][np.isfinite(post[kind])];K=min((1 if kind==0 else 2),np.isfinite(new[kind]).sum());assert len(live)==K;expected=np.sort(new[kind].ravel())[-K:];reweighted=expected+masses[kind]-ls(expected)-total;np.testing.assert_allclose(np.sort(live),reweighted,rtol=1e-10,atol=1e-8)
                    else:np.testing.assert_allclose(post[kind],new[kind]-total,rtol=1e-10,atol=1e-8)
                    compressionchecks+=1
                rec=curve_lookup[(tag,method)][t];np.testing.assert_allclose([weights[1],masses[1]-masses[0]],[rec['present'],rec['logodds']],rtol=1e-9,atol=2e-6);snapshots+=1
        assert binding()==read(OUT/'plan.json')['source'];atomic_json(OUT/'separate_check.json',dict(passed=True,oracle_prefix_closed_form_checks=oraclechecks,independent_two_state_HMM_checks=hmmchecks,explicit_transition_nodes_checked=nodes,grid_snapshots_checked=snapshots,mass_compression_checks=compressionchecks,summary_checks=summarychecks,MG172_restricted_comparisons=closedform,elapsed_seconds=time.perf_counter()-start,limits='Independent weighted-observation closed form and componentwise scalar Gaussian node sums; shared precision/templates, same author, not external replication.'))
        manifest={p.name:dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in RAW.iterdir() if p.is_file()};atomic_json(OUT/'local_audit_manifest.json',dict(directory=str(RAW),files=manifest,total_bytes=sum(v['bytes'] for v in manifest.values()),full_grids_published=False));atomic_json(attempt,dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start))
    except BaseException:atomic_json(attempt,dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':run()
