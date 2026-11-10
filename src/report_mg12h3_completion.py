"""MG122 saved-only independent audit and corrected-matrix report."""
from mg12h3_completion import *
from mg11g2_pattern_v2_oracle import dense_projection
from mg11g2_pattern_v2 import reference as template
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def main():
    reg=read(BASE/'registration.json');e=read(BASE/'evaluation.json');run=read(BASE/'run_state.json');assert e['qualified'] and e['complete'] and e['completed']==24;assert completion_hashes()==reg['source_original_hashes'];assert matrix()==reg['entries'];verify_reuse(reg['entries']);manifest=[];audit=[];fields={};scores={};old_scores={};failed_original=OLD/'production/ISOLATED_C1_s1_v0.npz';assert failed_original.exists() and not (OLD/'production/ISOLATED_C1_s1_v1.npz').exists()
    for entry in reg['entries']:
        job=entry['job'];path=ROOT/entry['path'];case,ga,gb,k,mode,reference,parent=factory(job)
        with np.load(path) as z:state=z['state'];times=z['times'];assert json.loads(str(z['job_json']))==job
        expected=np.arange(int(job['duration']*4)+1)/4;np.testing.assert_array_equal(times,expected);pop=state[:,:3500].reshape(-1,7,500);imports=state[:,3500:].reshape(-1,2,500);scale=np.maximum(1,case['initial'].max(axis=1));np.testing.assert_array_equal(pop[0],case['initial']);assert not np.any(imports[0])
        with np.load(reference) as z:old=z['state'][:,:3500].reshape(-1,7,500);np.testing.assert_array_equal(z['times'],times)
        error=max(float(np.max(abs(pop[:,s]-old[:,s])))/scale[s] for s in range(7));balance=max(abs(sum(float(frame[s].sum()) for s in range(7))-float(case['initial'].sum())-sum(float(account[s].sum()) for s in range(2))) for frame,account in zip(pop,imports))/float(case['initial'].sum());assert error<=1e-5 and balance<=1e-10 and state.min()>=-1e-9 and np.all(np.isfinite(state))
        if mode=='FINITE':assert not np.any(imports)
        record=next(x for x in e['details'] if x['entry']['job']['id']==job['id']);np.testing.assert_allclose(error,record['check']['reference_error'],atol=1e-14)
        if entry['origin']=='fresh_corrected_replay':
            metadata=read(path.with_suffix('.json'));assert metadata['status']=='COMPLETE' and metadata['entry']==entry and metadata['sha256']==hashlib.sha256(path.read_bytes()).hexdigest()
            with np.load(ROOT/entry['diagnostic_reference']) as z:np.testing.assert_array_equal(state,z['state'])
        fields[job['id']]=pop;audit.append(dict(id=job['id'],reference_error=error,balance=balance,origin=entry['origin']));manifest.append(dict(file=entry['path'],job_id=job['id'],origin=entry['origin'],method=entry['method'],bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
        if job['kind']=='FINITE' and job['setting']==1:
            keep=times>=10;new=dense_projection(pop,template()[0]);old=dense_projection(old,template()[0]);key=job['origin'],job['start'],job['arm'];scores[key]=sum((new[keep,1][1:]+new[keep,1][:-1])*.125)/5;old_scores[key]=sum((old[keep,1][1:]+old[keep,1][:-1])*.125)/5
    for pair in e['pairs']:
        case=factory(next(j for j in jobs() if j['id']==pair['left']))[0];scale=np.maximum(1,case['initial'].max(axis=1));error=max(float(np.max(abs(fields[pair['left']][:,s]-fields[pair['right']][:,s])))/scale[s] for s in range(7));assert error<=1e-5;np.testing.assert_allclose(error,pair['error'],atol=1e-14)
    def contrast(values):return sum((1 if s else -1)*(1 if a else -1)*values[o,s,a] for o in range(2) for s in range(2) for a in range(2))/2
    new=float(contrast(scores));old=float(contrast(old_scores));np.testing.assert_allclose(new,e['interaction']['new'],atol=1e-10);assert abs(new-old)<=1e-5;assert completion_hashes()==reg['source_original_hashes'];atomic_json(BASE/'independent_check.json',dict(passed=True,records=audit,slots=24,reused=22,fresh=2,pairs=12,new_interaction=new,closed_reference=old,original_fine_job_still_unstarted=True,source_original_diagnostic_reference_hashes_unchanged=True,shared='same fields/fixtures and frozen dense projection, separate specieswise error and balance/pair summation; original/fresh metadata and byte-exact diagnostic replay checked; same author not external replication'));atomic_json(BASE/'raw_manifest.json',dict(files=24,new_files=2,reused_files=22,total_bytes=sum(r['bytes'] for r in manifest),new_bytes=sum(r['bytes'] for r in manifest if r['origin']=='fresh_corrected_replay'),records=manifest))
    categories=['Finite prepared','Maintained','Isolated H','Corrected isolated S'];maxima=[0.]*4;pairmax=[0.]*4
    def category(job):return 0 if job['kind']=='FINITE' else 1 if job['kind']=='MAINTAINED' else 2 if job['start']==0 else 3
    for detail in e['details']:i=category(detail['entry']['job']);maxima[i]=max(maxima[i],detail['check']['reference_error'])
    for pair in e['pairs']:job=next(j for j in jobs() if j['id']==pair['left']);i=category(job);pairmax[i]=max(pairmax[i],pair['error'])
    fig,axes=plt.subplots(1,2,figsize=(11.3,4.5))
    for ax,values,title in zip(axes,[maxima,pairmax],['Reference agreement: all 24 slots','Case-specific refinement: all 12 pairs']):
        ax.bar(np.arange(4),values,color=['#557f91']*3+['#187a77']);ax.set_yscale('log');ax.axhline(1e-5,color='#333333',linestyle='--',label='Unchanged limit');ax.set_xticks(np.arange(4),['Finite\nprepared','Maintained','Isolated\nH','Corrected\nisolated S']);ax.set_ylabel('Largest scaled field difference');ax.set_title(title);ax.grid(axis='y',alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Corrected deterministic completion: 22 validated reuse + 2 fresh replays');fig.tight_layout();fig.savefig(BASE/'MG12H3_completion_summary.png',dpi=170);plt.close(fig)
    maxref=max(r['reference_error'] for r in audit);maxpair=max(x['error'] for x in e['pairs']);maxbalance=max(x['balance'] for x in audit);maintained=next(x['reference_continuity']['matched_through'] for x in e['details'] if x['entry']['job']['id']=='MAINTAINED_s1_v1');fresh_records=read(BASE/'ledger.json')['records'];fresh_seconds=sum(x['elapsed_seconds'] for x in fresh_records)
    text=f'''# H3b corrected deterministic completion — MG122

The corrected completion study passed: the incidence-derived population-limit model recovers the tested earlier Model G trajectories and geometry-response measurement within the original accuracy requirements. We retained the 22 validated trajectories and generated two fresh isolated controls at the tighter settings supported by the independent-solver diagnostic. All 24 slots, all 12 refinement pairs and the population/import accounts passed. This completes the corrected deterministic reference bridge on the prescribed ring and resource graphs, while preserving the original stopped batch as failed and incomplete. It does not yet demonstrate particle behavior in the strict finite-tick stochastic engine, a particle on new higher-arity or 3D substrates, emergent space or a physical vacuum supply.

**Outcome: CORRECTED DETERMINISTIC REFERENCE MATRIX COMPLETE; REGISTERED CRITERIA SUPPORTED. PI scientific closure pending.**

![Reference and refinement checks](MG12H3_completion_summary.png)

## Evidence and provenance

This is the prospectively approved MG121 mixed-provenance matrix: **22 immutable validated original trajectories plus two fresh replays**, not 24 new independent runs. The original failed coarse isolated control is preserved and excluded from this new accepted matrix. The originally unstarted fine job remains unstarted in the original directory. No original result was overwritten or rescored.

The two new slots are isolated seeded controls to T25 at quarter-unit observations. Their BDF settings are rtol1e-11/atol1e-13/maxstep0.05 and rtol1e-12/atol1e-14/maxstep0.025. Reused slots retain their original numerical settings. Thus numerical levels are explicitly case-specific, rather than globally uniform. Physical equations, initial fields, reservoir conditions, archived references and acceptance thresholds remain unchanged.

Both fresh trajectories reproduced their corresponding fixed MG118 diagnostic arrays **exactly** on this platform. This is deterministic replay consistency, not independent physical replication. The separately assembled Radau result remains supporting diagnostic evidence, not an extra completion slot.

## Registered results

| Check | Observed result | Required limit |
|---|---:|---:|
| Largest scaled field/reference difference, all 24 slots | {maxref:.12g} | <=1e-5 |
| Largest paired field difference, all 12 context pairs | {maxpair:.12g} | <=1e-5 |
| Largest relative population/import residual | {maxbalance:.12g} | <=1e-10 |
| Fresh replay difference from diagnostic arrays | 0 for both | <=1e-10 |
| Fine late-X interaction difference from closed G3 | {abs(new-old):.12g} | <=1e-5 |

All fields were finite and nonnegative. Closed cases had zero external imports; maintained controls explicitly accounted for A/B replenishment. Z/W remain retained spent products, with no recycling mechanism added. Scaling remains max(1,initial_species_max), including divisor1 for initially absent outputs. No favorable final-Z scaling was introduced.

The completed finite branches gave interaction **{new:.12f}**, versus recomputed closed G3 **{old:.12f}**. A separate dense projection reproduced the value. This reproduces earlier deterministic evidence; it is not a new stochastic robustness result or an independent discovery of the same physical effect. G5 remains closed and unresolved.

The maintained seeded fine control retained sampled reference matching through **T={maintained:g}**. Its reservoirs are maintained through recorded external imports, not derived from a physical substrate. Finite prepared branches retain their earlier loss-of-reference-matching behavior; late projection is not physical lifetime, identity or death. The fixed ring is one-dimensional, and resource links are prescribed. No new geometry or quantum behavior is inferred.

## Independent saved-data audit

The complete-matrix evaluator rechecked all24slots with explicit origin/method provenance, archived references, balance/positivity, refinement, projection and held replay guards. A reporting-only second path reread all raw fields, used specieswise maxima and separate population/import summation, recomputed all12pair errors and the finite interaction, verified completion metadata/checksums and confirmed the fresh arrays match the pinned diagnostics exactly. Canonical plan/preparation/source/original/donor/diagnostic/reference hashes stayed unchanged.

Shared dependencies include raw fields, physical inputs, the legacy template/dense projection and some provenance helpers. These are same-author separate analysis paths, not external replication or a rigorous exact-solution error certificate. The archived fine references themselves have finite numerical accuracy, as documented by the preceding Radau diagnostic.

## Execution and retention

Richard approved the concrete reuse/two-fresh plan, numerical correction, unchanged criteria,300second cap and full public raw/provenance/results. Registration was published at **c72e397** before new integration. The source and all preexisting inputs were pinned; no new RNG or stochastic sample selection was involved.

One local worker, numerical-library thread1, completed the two fresh integrations and complete-matrix audit in **{run['elapsed_seconds']:.6f} seconds**, below300seconds and the conservative50.03second forecast. The two fresh solve/IO records sum to **{fresh_seconds:.6f} seconds**. Original and diagnostic costs remain separately recorded; reused evidence is not newly charged or described as fresh work. Finished **{run['finished']}** UTC. Reporting and publication are separate from the capped execution interval. No failure, retry, unplanned replacement/integration or threshold change occurred.

The complete24slot raw manifest references the22existing public trajectories and two newly published raw files without duplicating/altering the originals. It records each effective method, origin path, bytes and SHA256. Fresh fields/import accounts, prefix/completion metadata, exact plan/preparation/registration, evaluations, audit, report/figure and budget ledger are retained here. See raw_manifest.json and production_retention_verified.json for post-publication verification. OriginalMG117and diagnosticMG119remain separately visible.

## Scope and recommended next stage

This completes the corrected **deterministic population-limit reference experiment**, not all of H3. The mixed-law finite-tick engine still needs formal statistical qualification, and particle-scale execution needs an efficient exchangeable-token representation proven equivalent to the transparent identity-rich reference, including catalyst reservations. A later noisy-pattern study must specify adequate precision and its population/tick/horizon limits before production.

After PI acceptance of this corrected completion, the next useful H3 work is that mixed reaction/transport finite-tick qualification and efficiency bridge. New higher-arity or known-dimensional pattern substrates, spontaneously formed particles, evolving topology and substrate-generated supply require further separate tests. No next run is authorized by this report.
'''
    (BASE/'MG12H3_completion_results.md').write_text(text,encoding='utf-8')
if __name__=='__main__':main()
