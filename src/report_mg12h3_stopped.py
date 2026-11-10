"""MG117 audit/report retained stopped run; no new integration or rescore."""
from prepare_mg12h3_population import *
from datetime import datetime,timezone
from mg11g2_pattern_v2_oracle import dense_projection
from mg11g2_pattern_v2 import reference as template
from mg10p1 import continuity
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    registration=read(BASE/'registration.json');assert registration['source_hashes']==hashes();assert hashlib.sha256((PREP/'results.json').read_bytes()).hexdigest()==registration['preparation_sha256'];ledger=read(BASE/'ledger.json');records={r['job']['id']:r for r in ledger['records']};details=[];fields={};times_by={};newvalues={};oldvalues={};other_values={};manifest=[];failed=None
    for job in registration['jobs']:
        p=BASE/'production'/(job['id']+'.npz')
        if not p.exists():details.append(dict(job=job,status='NOT STARTED'));continue
        case,ga,gb,k,mode,reference,parent=factory(job);n=ga['n']
        with np.load(p) as z:state=z['state'];times=z['times'];assert json.loads(str(z['job_json']))==job
        pop=state[:,:7*n].reshape(-1,7,n);imports=state[:,7*n:].reshape(-1,2,n);np.testing.assert_array_equal(pop[0],case['initial']);assert not np.any(imports[0]);expected=np.arange(int(job['duration']*4)+1)/4;complete=np.array_equal(times,expected);assert complete and np.all(np.isfinite(state));scale=np.maximum(1,case['initial'].max(axis=1))
        with np.load(reference) as z:old=z['state'][:len(times),:7*n].reshape(-1,7,n);np.testing.assert_array_equal(z['times'][:len(times)],times)
        error=abs(pop-old)/scale[None,:,None];maximum=float(error.max());location=np.unravel_index(np.argmax(error),error.shape);balance=0.
        for frame in range(len(times)):
            residual=sum(float(pop[frame,s].sum()) for s in range(7))-float(case['initial'].sum())-sum(float(imports[frame,s].sum()) for s in range(2));balance=max(balance,abs(residual)/float(case['initial'].sum()))
        minimum=float(state.min());assert minimum>=-1e-9 and balance<=1e-10
        if mode=='FINITE':assert not np.any(imports)
        status='VALIDATED COMPLETE' if job['id'] in records else 'REFERENCE CRITERION FAILED';assert (maximum<=1e-5)==(job['id'] in records)
        detail=dict(job=job,status=status,frames=len(times),max_reference_error=maximum,independent_balance=balance,minimum=minimum,max_location=dict(time=float(times[location[0]]),species=('G','X','Y','A','B','Z','W')[location[1]],site=int(location[2]),new=float(pop[location]),reference=float(old[location])),species_errors=error.max(axis=(0,2)).tolist())
        if status=='REFERENCE CRITERION FAILED':failed=detail;failed_errors=error;failed_times=times;failed_pop=pop;failed_old=old
        else:
            assert records[job['id']]['status']=='COMPLETE';np.testing.assert_allclose(maximum,records[job['id']]['check']['reference_error'],atol=1e-14)
        fields[job['id']]=pop;times_by[job['id']]=times;details.append(detail);manifest.append(dict(file=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
        if job['setting']==1:
            detail['reference_continuity']=continuity(pop[:,1],times+(10 if job['kind']=='FINITE' else 0))
            if job['kind']=='FINITE':
                primary=np.asarray(measure(job,state,times)['projection']);other=dense_projection(pop,template()[0]);reference_scores=dense_projection(old,template()[0]);np.testing.assert_allclose(primary,other,atol=1e-10,rtol=1e-10);keep=times>=10;key=job['origin'],job['start'],job['arm'];newvalues[key]=float(np.trapezoid(primary[keep,1],times[keep])/5);other_values[key]=float(np.sum((other[keep,1][1:]+other[keep,1][:-1])*.125)/5);oldvalues[key]=float(np.sum((reference_scores[keep,1][1:]+reference_scores[keep,1][:-1])*.125)/5)
    assert len(records)==22 and len(manifest)==23 and failed is not None;assert len(newvalues)==8
    refinement=[]
    for job in [j for j in jobs() if j['setting']==0]:
        other=job['id'][:-1]+'1'
        if job['id'] in fields and other in fields:
            scale=np.maximum(1,factory(job)[0]['initial'].max(axis=1));error=float(np.max(abs(fields[job['id']]-fields[other])/scale[None,:,None]));assert error<=1e-5;refinement.append(dict(base=job['id'],fine=other,error=error))
    def contrast(values):return sum((1 if s else -1)*(1 if a else -1)*values[o,s,a] for o in range(2) for s in range(2) for a in range(2))/2
    new=contrast(newvalues);old=contrast(oldvalues);independent=contrast(other_values);assert abs(new-independent)<1e-10 and abs(new-old)<1e-5
    first=float(failed_times[np.flatnonzero(failed_errors.max(axis=(1,2))>1e-5)[0]]);failed['first_failed_observation']=first;failed['relative_error_to_final_Z']=abs(failed['max_location']['new']-failed['max_location']['reference'])/abs(failed['max_location']['reference'])
    stop_at=datetime.fromtimestamp((BASE/'failure.json').stat().st_mtime,timezone.utc);before=read(BASE/'run_state.json');started=datetime.fromisoformat(before['started']);reconstructed=(stop_at-started).total_seconds()
    if not (BASE/'run_state_before_report.json').exists():atomic_json(BASE/'run_state_before_report.json',before)
    atomic_json(BASE/'run_state.json',dict(status='STOPPED ON REGISTERED REFERENCE CRITERION; no retry; audit/report ready',started=started.isoformat(),failure_record_mtime_utc=stop_at.isoformat(),reconstructed_elapsed_seconds=reconstructed,elapsed_note='Failure file timestamp minus recorded start; final failed job lacked a recorded performance-counter charge',validated_jobs=22,failed_jobs=1,not_started=1,failed_id=failed['job']['id']))
    result=dict(verdict='REGISTERED NUMERICAL REFERENCE CRITERION FAILED; STUDY INCOMPLETE',qualified=False,validated_complete=22,retained_full_trajectories=23,failed=1,not_started=1,details=details,available_refinement_pairs=refinement,finite_interaction_component=dict(new_fine=new,closed_reference=old,separate_projection=independent,error=abs(new-old),within_component_margin=True,scope='descriptive supported component; does not qualify incomplete global study'),failed_case=failed,all_retained_balance_and_nonnegative_checks_passed=True,source_reference_hashes_unchanged=True,known_complete_job_elapsed_seconds=ledger['elapsed_seconds'],reconstructed_elapsed_seconds=reconstructed,automatic_retry=False);atomic_json(BASE/'stopped_evaluation.json',result);atomic_json(BASE/'production_raw_manifest.json',dict(files=23,bytes=sum(x['bytes'] for x in manifest),records=manifest))
    fig,axes=plt.subplots(1,2,figsize=(11.2,4.5));zerror=failed_errors[:,5].max(axis=1);active=failed_errors[:,:3].max(axis=(1,2));axes[0].semilogy(failed_times,np.maximum(zerror,1e-16),label='Spent Z: largest scaled difference',color='#aa5748');axes[0].semilogy(failed_times,np.maximum(active,1e-16),label='Active G/X/Y: largest scaled difference',color='#187a77');axes[0].axhline(1e-5,color='#333333',linestyle='--',label='Registered limit');axes[0].set_xlabel('Time');axes[0].set_ylabel('Difference from archived fine reference');axes[0].set_title('Seeded isolated control, coarser setting');axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
    labels=['G','X','Y','A','B','Z','W'];axes[1].bar(labels,np.maximum(failed_errors.max(axis=(0,2)),1e-16),color=['#187a77']*5+['#aa5748','#8c8074']);axes[1].set_yscale('log');axes[1].axhline(1e-5,color='#333333',linestyle='--');axes[1].set_ylabel('Largest scaled difference over all frames/sites');axes[1].set_title('The registered failure is in accumulated Z');axes[1].grid(axis='y',alpha=.2);fig.tight_layout();fig.savefig(BASE/'MG12H3_stopped_summary.png',dpi=170);plt.close(fig)
    text=f'''# MG12-H3b deterministic reference bridge — stopped result (MG117)

The deterministic reference experiment recovered the earlier finite-reservoir trajectories and geometry-response measurement in the completed cases, but it did not pass the whole registered test. It stopped when one coarser isolated seeded control exceeded the accuracy limit for accumulated spent product Z. The active G/X/Y fields in that case remained close to the reference, and all retained trajectories passed positivity and population/import balance checks. This points to a numerical accuracy issue in an accumulated output rather than an observed loss of the pattern mechanism, but a tighter independent calculation is still needed to confirm the cause. The failed result is preserved: H3b is incomplete and is not qualified, and no automatic rerun has been made.

**Registered outcome: NUMERICAL REFERENCE CRITERION FAILED; STUDY INCOMPLETE. PI disposition pending.**

There are **22 validated jobs, one full trajectory that failed its reference check, and one job not started**. No source, threshold or numerical setting was changed after the failure.

![Stopped accuracy check and affected fields](MG12H3_stopped_summary.png)

## What stopped the run

The failed case was **ISOLATED_C1_s1_v0**: isolated finite reservoirs, seeded pattern, coarser numerical setting. Its trajectory reached T=25 and was retained. The next fine-setting seeded isolated job was not launched.

The registered criterion compares every species/frame/site against its archived fine reference, dividing by max(1, that species' maximum initial concentration). The limit is 1e-5. Z starts at zero, so its divisor is 1 and its criterion is effectively an absolute 1e-5 concentration difference, even after it accumulates substantially.

At T=25, site {failed['max_location']['site']}, Z was **{failed['max_location']['new']:.12f}**, versus reference **{failed['max_location']['reference']:.12f}**. The difference was **{failed['max_reference_error']:.12f}**, above 1e-5 by a factor of {failed['max_reference_error']/1e-5:.2f}. It is about **{failed['relative_error_to_final_Z']:.3g} of the local final Z value**. These are different normalizations: the latter small relative discrepancy does not override the frozen initial-scale criterion. The first saved observation exceeding the criterion was T={first:g}.

The largest scaled active-field difference in the failed case was **{max(failed['species_errors'][:3]):.12g}**, below 1e-5. W also stayed below the limit. Balance residual was **{failed['independent_balance']:.12g}**, and all fields were nonnegative. Accumulation of a small rate/integration difference is consistent with the increasing Z discrepancy; it has not been established by a new tighter-tolerance solve. No negative resources, lost tokens or physical sink interpretation is inferred.

## Supported completed components

All 16 finite prepared branches passed their reference checks, as did all four maintained-control jobs and both homogeneous isolated controls. Eleven available base/fine pairs passed the refinement criterion. The seeded isolated pair is unavailable because its fine job did not start.

The eight completed fine finite branches give the late-X interaction **{new:.12f}**, compared with recomputed closed G3 value **{old:.12f}**. Difference **{abs(new-old):.12g}** is within the registered 1e-5 component margin. A separate dense projection gives **{independent:.12f}**. This is a supported component of an incomplete deterministic reference study, not a global pass, new stochastic evidence, or a rescore of closed G3/G5.

The maintained seeded fine control retained sampled reference matching through T=25. Their A/B supplies are explicitly maintained through recorded external imports; this does not explain their physical source. Z/W remain outputs. Pattern projection or classifier loss is not physical particle death. All these trajectories are mathematical population-limit fields, not strict finite-tick token histories.

## Audit, execution and limitations

The independent reporting audit reread all **23 full raw trajectories**, checked canonical job/time metadata and starting states, compared all seven fields to frozen archives, reconstructed population-minus-import balance by a separate summation, verified positivity, paired the 11 available refinements and independently recomputed the finite response. Source/reference/preparation hashes match MG116 registration. The failed trajectory was included in the audit and retained rather than deleted or treated as absent.

The original driver recorded 22 completed jobs and **{ledger['elapsed_seconds']:.6f} seconds** through the last validated job. It did not write a performance-counter record for the failed case or replace its RUNNING status on the exception. The original failure traceback and run-state snapshot are preserved. Canonical run status is now explicitly stopped. From the failure file's UTC modification timestamp and recorded start, elapsed time is approximately **{reconstructed:.3f} seconds**; this is a reconstructed wall-time figure, not an exact final worker-charge ledger. The run was far below the approved 3,600-second cap. One local worker and one numerical-library thread were used; no timeout or budget stop occurred.

Registration/source/criteria were published at b42ab16 before new integrations. Richard approved the 24-job deterministic reference experiment, one-hour cap and full public raw/results. The new report/evaluation/figure are reporting-only; frozen execution code and scientific criteria remain unchanged. All raw/reference provenance, failed output, logs/metadata, source and numerical evidence are published and verified separately in production_retention_verified.json. Same-author separate arithmetic and archived references are not independent external replication.

## Recommendation

Preserve this run as failed/incomplete under its registered criteria. A separately registered numerical follow-up should first test tighter integration tolerances for the isolated seeded control against an independent high-accuracy solver, keeping the same output-error threshold. If that resolves the discrepancy, a continuation or replacement experiment must explicitly specify treatment of the 22 validated jobs, the failed trajectory and the unstarted fine job. Nothing is silently retried or rescored.

The evidence currently supports recovery of the completed deterministic finite-branch response and maintained pattern controls. It does not yet qualify the entire deterministic reference batch, the strict finite-tick stochastic particle process, a new higher-arity/3D particle or a vacuum supply mechanism. Formal mixed-law accuracy and efficient stochastic representation remain further H3 work. **No follow-up integration has been launched.**
'''
    (BASE/'MG12H3_stopped_results.md').write_text(text,encoding='utf-8')
if __name__=='__main__':main()
