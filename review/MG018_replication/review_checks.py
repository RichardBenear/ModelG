"""Independent MG5 record analysis. Does not import the simulator/evaluators.
Uses periodic KD-tree neighborhoods, explicit component traversal, and fresh
track/MSD/gate calculations. Does not re-integrate the ODE or verify smoothing.
"""
from pathlib import Path
import csv, json, platform
import numpy as np
import scipy
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'inputs/results/MG5'


def detect(pos,s,rows,theta,box):
    order=np.argsort(rows,kind='stable')
    pos,s,rows=pos[order],s[order],rows[order]
    if len(s)==0:return []
    tree=cKDTree(np.mod(pos,box),boxsize=box)
    candidates=[]
    for i in np.flatnonzero(s<theta):
        nb=np.asarray(tree.query_ball_point(pos[i],1.0))
        minimum=s[nb].min()
        tied=nb[s[nb]==minimum]
        if s[i]==minimum and rows[i]==rows[tied].min():candidates.append(i)
    if not candidates:return []
    ct=cKDTree(pos[candidates],boxsize=box)
    links=ct.query_ball_tree(ct,1.5)
    unseen=set(range(len(candidates))); answer=[]
    while unseen:
        todo=[min(unseen)]; group=[]; unseen.remove(todo[0])
        while todo:
            j=todo.pop();group.append(candidates[j])
            for k in links[j]:
                if k in unseen:unseen.remove(k);todo.append(k)
        i=min(group,key=lambda j:(s[j],rows[j]))
        nb=np.asarray(tree.query_ball_point(pos[i],.6))
        nb=nb[np.argsort(rows[nb])]
        d=pos[nb]-pos[i];d-=box*np.rint(d/box)
        w=np.clip(theta-s[nb],0,None)
        centre=np.mod(pos[i]+np.sum(w[:,None]*d,axis=0)/sum(w),box)
        answer.append(dict(centre=centre,depth=float(s[i]),row=int(rows[i])))
    return sorted(answer,key=lambda p:(p['depth'],p['row']))


def classify(parts,box):
    counts=[len(p) for p in parts[1:]]
    if len(parts[0])!=1:return dict(cls='unresolved')
    pairs=list(zip(counts[:-1],counts[1:]))
    if any(a>=2 and b>=2 for a,b in pairs):
        first=next(i for i,n in enumerate(counts) if n>=2)
        ref=parts[0][0]['centre']
        for p in parts[1:first+1]:
            if len(p)==1:ref=p[0]['centre']
        d=np.array([p['centre']-ref for p in parts[first+1]])
        d-=box*np.rint(d/box)
        second=np.sort(np.linalg.norm(d,axis=1))[1]
        qualifying=next(i for i,(a,b) in enumerate(pairs) if a>=2 and b>=2)
        return dict(cls='multiplied',sub='split' if second<=4 else 'remote',
                    first_multi_time=10*(first+1),confirmation_time=10*(qualifying+2))
    if any(a==0 and b==0 for a,b in pairs):return dict(cls='dissolved')
    if not all(n==1 for n in counts):return dict(cls='unresolved')
    centres=np.array([p[0]['centre'] for p in parts])
    steps=np.diff(centres,axis=0);steps-=box*np.rint(steps/box)
    if np.max(np.linalg.norm(steps,axis=1))>2:return dict(cls='unresolved')
    track=np.vstack([centres[0],centres[0]+np.cumsum(steps,axis=0)])
    md=float(np.max(np.linalg.norm(track-track[0],axis=1)))
    tau=np.arange(1,11)*10.
    msd=np.array([np.mean(np.sum((track[lag:]-track[:-lag])**2,axis=1)) for lag in range(1,11)])
    dp=float(np.dot(tau,msd)/(4*np.dot(tau,tau)))
    return dict(cls='single',mobility='mobile' if md>=1.5 else 'pinned',max_d=md,
                D_p=dp,track=track.tolist(),msd=msd.tolist(),max_step=float(np.max(np.linalg.norm(steps,axis=1))))


def read(path,theta):
    with np.load(path,allow_pickle=False) as z:d={k:z[k] for k in z.files}
    assert d['off'][0]==0 and d['off'][-1]==len(d['s'])
    assert np.all(np.isfinite(d['s'])) and np.all(d['s']< -1)
    box=np.array([float(d['Lx']),float(d['Ly'])]);parts=[]
    for i in range(len(d['t'])):
        a,b=map(int,d['off'][i:i+2]);r=d['rows'][a:b]
        assert np.all(np.diff(r)>0) and np.all((r>=0)&(r<d['N'][i]))
        assert len(np.unique(d['ids'][a:b]))==b-a
        if b>a:
            assert d['s'][a:b].min()==d['smin'][i]
        parts.append(detect(d['pos'][a:b],d['s'][a:b],r,theta,box))
    return d,parts,box


def main():
    evaluation=json.loads((DATA/'evaluation.json').read_text())
    out=dict(commit=json.loads((ROOT/'download_manifest.json').read_text())['commit'],
             environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),runs={})
    depths={};groups={};max_track_error=0.;max_dp_error=0.;max_depth_error=0.;audits=[];jumps=[]
    for path in sorted((DATA/'runs').glob('*.json')):
        rec=json.loads(path.read_text());key=rec['key']
        assert rec['status']=='COMPLETE' and not rec['flags'] and rec['t_reached']==rec['t_end']
        assert rec['dt']==.05 and rec['delta']==(1. if rec['nu']<=.04 else .2)
        assert np.allclose(rec['phi_shift'],[0,0,0 if rec['spec']['part']=='A' else -.9/14],rtol=0,atol=1e-15)
        expected_theta=-4.457581675528346 if rec['spec']['part']=='A' else -5.117354421620525
        assert rec['theta']==expected_theta
        d,parts,box=read(DATA/'checks'/f'{key}.npz',rec['theta'])
        assert np.array_equal(d['t'],np.arange(0,rec['t_end']+1,10))
        assert len(d['t'])==(51 if rec['spec']['part']=='A' else 31)
        mine=classify(parts,box);e=evaluation['runs'][key]['cls'];driver=rec['summary']
        assert mine['cls']==e['cls']==driver['cls']['cls']
        assert [len(p) for p in parts[1:]]==driver['counts']
        for p,sd,sc in zip(parts[1:],driver['depths'],driver['centres']):
            assert len(p)==len(sd)
            if p:
                derr=float(np.max(np.abs(np.array([v['depth'] for v in p])-sd)))
                delta=np.array([v['centre'] for v in p])-sc;delta-=box*np.rint(delta/box)
                assert np.max(np.abs(delta))<1e-9
                assert derr<1e-12;max_depth_error=max(max_depth_error,derr)
        if mine['cls']=='single':
            assert mine['mobility']==e['mobility']==driver['cls']['mobility']
            err=float(np.max(np.abs(np.array(mine['track'])-e['track'])))
            max_track_error=max(max_track_error,err);assert err<1e-9
            assert abs(mine['max_d']-e['max_d'])<1e-9
            dp_err=abs(mine['D_p']-e['D_p']);max_dp_error=max(max_dp_error,dp_err);assert dp_err<1e-12
        else:assert mine['sub']==e['sub']==driver['cls']['sub']
        K=round(rec['t_end']/rec['delta']);assert rec['n_factor']==K
        expected_ticks=np.arange(1,K)*rec['delta'] if rec['nu'] else np.empty(0)
        assert np.array_equal(d['tick_t'],expected_ticks)
        if len(expected_ticks):
            assert np.all(d['tick_n_prop_births']==d['tick_n_eff_births']+d['tick_n_ephemeral'])
            assert np.all(d['tick_n_prop_deaths']==d['tick_n_eff_deaths']+d['tick_n_ephemeral'])
            change=d['tick_n_eff_births']-d['tick_n_eff_deaths']
            assert np.array_equal(d['tick_N'],18901+np.cumsum(change))
            for i,t in enumerate(d['t']):
                j=int(np.searchsorted(d['tick_t'],t,side='left'))-1
                assert d['N'][i]==(18901 if j<0 else d['tick_N'][j])
            residual=np.max(np.abs(d['tick_I_after']-(d['tick_I_before']-d['tick_R']+d['tick_B']))/
                (1+np.abs(d['tick_I_before'])+np.abs(d['tick_R'])+np.abs(d['tick_B'])),axis=1)
            assert np.max(np.abs(residual-d['tick_jump_resid']))<1e-20
            assert np.all(residual<=1e-10) and np.all(d['tick_audit_max']<=1e-9)
            assert rec['max_jump_resid']==float(residual.max())
            assert rec['max_audit']==float(d['tick_audit_max'].max())
            jumps.extend(residual);audits.extend(d['tick_audit_max'])
        assert d['N'][-1]==rec['N_final']
        slim={k:v for k,v in mine.items() if k not in ('track','msd')}
        slim.update(counts=[len(p) for p in parts],smin_range=[float(d['smin'].min()),float(d['smin'].max())],
                    ticks=len(d['tick_t']),flags=len(rec['flags']))
        out['runs'][key]=slim
        depths[key]=np.array([p[0]['depth'] if len(p)==1 else np.nan for p in parts])
        groups.setdefault((rec['spec']['part'],rec['nu']),[]).append((key,mine))
        print(key,slim['cls'],slim.get('mobility'),slim.get('max_d'),flush=True)
    assert len(out['runs'])==24
    controls=groups['A',0]
    assert len(controls)==4 and all(r['cls']=='single' and r['max_d']<.1 for _,r in controls)
    rates={}
    for nu in [.002,.02,.2]:
        records=[r for _,r in groups['A',nu]]
        single=sum(r['cls']=='single' for r in records)
        mob=sum(r.get('mobility')=='mobile' for r in records)
        rates[str(nu)]=dict(single=single,pinned=single-mob,mobile=mob,
            mobility='MOBILE' if mob>=3 else 'PINNED' if single-mob>=3 else 'MIXED',
            max_d=[r['max_d'] for r in records],D_p=[r['D_p'] for r in records])
        assert rates[str(nu)]['mobility']==evaluation['partA']['rates'][str(nu)]['mobility']
    survive=[r['single']>=3 for r in rates.values()]
    gate='PASS' if all(survive) else 'PARTIAL' if any(survive) else 'FAIL'
    frozen=sum(r['cls']=='multiplied' for _,r in groups['B',0])
    fast_single=sum(r['cls']=='single' for _,r in groups['B',.2])
    fast_mult=sum(r['cls']=='multiplied' for _,r in groups['B',.2])
    b='A-void' if frozen<3 else 'A1' if fast_single>=3 else 'A2' if fast_mult>=3 else 'A3'
    assert gate==evaluation['gate'] and b==evaluation['partB']['outcome']
    out.update(gate=gate,partA=rates,partB=dict(outcome=b,frozen_multiplied=frozen,fast_single=fast_single),
        max_track_coordinate_error=max_track_error,max_depth_error=max_depth_error,max_Dp_error=max_dp_error,
        max_jump_residual=float(max(jumps)),max_interpolation_audit=float(max(audits)),edited_ticks=len(jumps))
    shape={}
    for key in depths:
        if key.startswith('A_nu') and not key.startswith('A_nu0_'):
            ctrl='A_nu0_'+key.split('_',2)[2]
            ratio=np.abs(depths[key]/depths[ctrl])
            shape[key]=dict(final_depth_ratio=float(ratio[-1]),min_depth_ratio=float(ratio.min()),max_depth_ratio=float(ratio.max()))
    out['shape_against_matched_frozen']=shape
    ledger=list(csv.DictReader((DATA/'budget_ledger.csv').open()))
    assert len(ledger)==24 and len({r['key'] for r in ledger})==24
    budget={}
    for row in ledger:
        rec=json.loads((DATA/'runs'/f"{row['key']}.json").read_text())
        assert abs(float(row['wall_s'])-rec['wall_s'])<.006
        group=f"{rec['spec']['part']}_nu{rec['nu']:g}"
        budget[group]=budget.get(group,0)+float(row['wall_s'])/3600
    hours=sum(budget.values());assert abs(hours-evaluation['core_h'])<1e-10
    out['ledger']=dict(core_h=hours,groups=budget,rows=len(ledger),lost_launch_core_h_reported_approx=65,
                       finished_plus_lost_core_h_approx=hours+65)
    p0={}
    for tag in ['step_dt0.05','step_dt0.025','tick_D0.2','tick_D0.1']:
        d,ps,box=read(DATA/'part0'/f'check_{tag}.npz',-4.457581675528346)
        assert all(len(p)==1 for p in ps) and np.array_equal(d['t'],np.arange(0,101,10))
        p0[tag]=ps
    out['part0']={}
    for label,a,b in [('step','step_dt0.05','step_dt0.025'),('tick','tick_D0.2','tick_D0.1')]:
        da=np.array([p[0]['centre'] for p in p0[a]]);db=np.array([p[0]['centre'] for p in p0[b]])
        d=da-db;d-=box*np.rint(d/box)
        depth_a=np.array([p[0]['depth'] for p in p0[a]]);depth_b=np.array([p[0]['depth'] for p in p0[b]])
        dc=float(np.max(np.linalg.norm(d,axis=1)));dd=float(np.max(np.abs(depth_a-depth_b)/np.abs(depth_a)))
        assert dc<=.1 and dd<=.02
        out['part0'][label]=dict(max_centre_difference=dc,max_relative_depth_difference=dd,counts_all_one=True)
    (ROOT/'review_checks.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k not in ('runs','shape_against_matched_frozen')},indent=2))


if __name__=='__main__':main()
