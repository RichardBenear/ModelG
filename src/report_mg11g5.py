"""MG102 reporting-only raw/config/retention inventory and static scientific figure.

No simulation, seed selection, criterion changes or statistical rescoring.
"""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mg11g5 import ROOT,jobs,factory,OBS,TAUS,OMEGA,concentration_case
from mg11g5_supplement import hashes
from mg11g4_io import check_raw,digest
from mg11g3_oracle import accounting
from check_mg11g5 import audit
from mg10d1_pool import atomic_json
BASE=ROOT/'results/MG11G5'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def main():
    reg=read(BASE/'registration.json');assert hashes()==reg['all_source_hashes'];assert digest(BASE/'preparation_v2/results.json')==reg['preparation_sha256'];assert digest(BASE/'supplement/results.json')==reg['supplement_sha256']
    e=read(BASE/'evaluation.json');ledger=read(BASE/'production/ledger.json');assert e['complete'] and len(ledger['records'])==80 and all(r['status']=='COMPLETE' and r['returncode']==0 for r in ledger['records']);assert not ledger['stop_reason'] and not ledger['active_ids'] and not ledger['pending_ids']
    records=[];stats={v:dict(steps=0,adapted=0,max_probability=0.) for v in range(2)}
    for j in jobs():
        case,ga,gb,k=factory(j);p=BASE/'production'/(j['id']+'.npz');meta=read(p.with_suffix('.json'));assert meta['job']==j and meta['sha256']==digest(p)
        if j['kind']=='deterministic':
            with np.load(p) as z:
                assert json.loads(str(z['job_json']))==j;np.testing.assert_array_equal(z['times'],OBS);np.testing.assert_array_equal(z['state'][0,:3500],case['initial'].ravel()/OMEGA);assert accounting(z['state'],concentration_case(case),ga,gb)['passed']
            records.append(dict(file=p.relative_to(ROOT).as_posix(),sha256=digest(p),bytes=p.stat().st_size,kind='deterministic',accounting=True));continue
        for warm in (False,True):
            path=p.with_name(p.stem+'_warm.npz') if warm else p;check=check_raw(path)
            with np.load(path) as z:
                expected=dict(initial=case['initial'],a0=case['a0'],b0=case['b0'],diff=case['diff'],src=ga['src'],dst=ga['dst'],weight=ga['weight'],srcB=gb['src'],dstB=gb['dst'],weightB=gb['weight'],kappa=k,mask=np.ones(7),times=np.array([0.,.001]) if warm else OBS)
                for key,value in expected.items():np.testing.assert_array_equal(z[key],value)
                assert int(z['seed'])==j['warmup_seed' if warm else 'seed'] and float(z['omega'])==OMEGA and float(z['tau'])==TAUS[j['setting']] and int(z['engine'])==1 and int(z['requested_reps'])==int(z['reps_completed'])==1
                if not warm:
                    s=stats[j['setting']];s['steps']+=int(z['stats'][0,0]);s['adapted']+=int(z['stats'][0,1]);s['max_probability']=max(s['max_probability'],float(z['stats'][0,2]))
            records.append(dict(check,bytes=path.stat().st_size,kind='excluded_warmup' if warm else 'stochastic'))
    assert len(records)==144;atomic_json(BASE/'production_raw_manifest.json',dict(files=144,bytes=sum(x['bytes'] for x in records),records=records))
    checked=audit(BASE/'production',BASE/'evaluation.json',BASE/'report_independent_check.json');assert checked['passed']
    values={(r['job']['setting'],r['job']['block'],r['job']['origin'],r['job']['start'],r['job']['arm']):r['late_X'] for r in e['details'] if r['job']['kind']=='stochastic'};blocks=[];origins=[]
    for v in range(2):
        by_origin=[[values[v,b,o,1,1]-values[v,b,o,0,1]-values[v,b,o,1,0]+values[v,b,o,0,0] for b in range(4)] for o in range(2)];origins.append(by_origin);blocks.append(np.mean(by_origin,axis=0).tolist())
    charge=sum(r['charged_worker_seconds'] for r in ledger['records']);np.testing.assert_allclose(charge,ledger['charged_worker_seconds'],atol=1e-8)
    summary=dict(verdict='COMPLETE PILOT; STOCHASTIC RESPONSE AND STEP AGREEMENT UNRESOLVED',complete_jobs=80,blocks=blocks,per_origin_blocks=origins,settings=e['settings'],refinement=e['refinement'],adaptation=stats,charged_worker_seconds=charge,core_hours=charge/3600,ceiling_core_hours=24,finished=read(BASE/'run_state.json')['finished'],source_and_qualification_hashes_unchanged=True,raw_files=144,raw_bytes=sum(x['bytes'] for x in records),control_interactions=read(BASE/'secondary.json')['conditional_control_interaction']);atomic_json(BASE/'sampling_summary.json',summary)
    plt.rcParams.update({'font.size':10});fig,axes=plt.subplots(1,2,figsize=(11.6,4.6),gridspec_kw={'width_ratios':[1.25,1]})
    ax=axes[0]
    for v,row in enumerate(e['settings']):
        ax.scatter(blocks[v],np.full(4,v)+np.linspace(-.08,.08,4),color='#8798a5',s=28,label='Independent block' if v==0 else None,zorder=3)
        ax.errorbar(row['mean'],v,xerr=[[row['mean']-row['lower']],[row['upper']-row['mean']]],fmt='o',color='#146b8c',capsize=5,markersize=7,label='Mean and nominal 97.5% interval' if v==0 else None,zorder=4)
    control=float(summary['control_interactions']['0']);ax.axvline(0,color='#444444',linewidth=1);ax.axvline(control,color='#bc7930',linestyle='--',label='Rounded continuum control');ax.set_yticks([0,1],['Coarser max step\n1/24,000','Finer max step\n1/48,000']);ax.set_ylim(1.5,-.55);ax.set_xlabel('Signed X geometry-response interaction');ax.set_title('Positive means; direction unresolved');ax.legend(loc='lower left',fontsize=8);ax.grid(axis='x',alpha=.18)
    ax=axes[1];r=e['refinement'];ax.axvspan(-.02,.02,color='#d7e8dd',label='Registered agreement region');ax.axvline(0,color='#444444',linewidth=1);ax.errorbar(r['mean'],0,xerr=[[r['mean']-r['lower']],[r['upper']-r['mean']]],fmt='o',color='#146b8c',capsize=5);ax.set_ylim(-.6,.6);ax.set_yticks([]);ax.set_xlabel('Finer minus coarser mean response');ax.set_title('Step agreement unresolved\n90% Welch interval');ax.legend(loc='lower center',fontsize=8);ax.grid(axis='x',alpha=.18)
    fig.suptitle('MG11-G5: complete conditional stochastic pilot, four blocks per setting',fontsize=13);fig.tight_layout();fig.savefig(BASE/'MG11G5_summary.png',dpi=170);plt.close(fig)
if __name__=='__main__':main()
