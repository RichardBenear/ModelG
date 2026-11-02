"""MG021 independent static review of committed MG6 Part 0 records.

Rebuilds candidates/node histories and density counts without simulator imports.
Re-detects on sparse smoothed X and diagnoses quadratic-centre sensitivity.
No ODE trajectories are integrated; alternative estimators never rescore Part 0.
"""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parent
P0=ROOT/'inputs/results/MG6/part0'


def seed(*key):
    return int(np.random.SeedSequence(20261027,spawn_key=tuple(key)+(0,)).generate_state(1,dtype=np.uint32)[0])


def wrap(d,box):return d-box*np.rint(d/box)


def detect(pos,s,rows,theta,box):
    tree=cKDTree(pos,boxsize=box);cand=[]
    for i in np.flatnonzero(s<theta):
        nb=np.asarray(tree.query_ball_point(pos[i],1.))
        mn=s[nb].min()
        if s[i]==mn and rows[i]==rows[nb[s[nb]==mn]].min():cand.append(i)
    if not cand:return []
    ct=cKDTree(pos[cand],boxsize=box);adj=ct.query_ball_tree(ct,1.5)
    unseen=set(range(len(cand)));out=[]
    while unseen:
        first=min(unseen);unseen.remove(first);todo=[first];group=[]
        while todo:
            a=todo.pop();group.append(cand[a])
            for b in adj[a]:
                if b in unseen:unseen.remove(b);todo.append(b)
        i=min(group,key=lambda j:(s[j],rows[j]))
        nb=np.asarray(tree.query_ball_point(pos[i],.6));nb=nb[np.argsort(rows[nb])]
        w=np.maximum(theta-s[nb],0)
        centre=np.mod(pos[i]+(w[:,None]*wrap(pos[nb]-pos[i],box)).sum(0)/w.sum(),box)
        out.append(dict(centre=centre,depth=float(s[i])))
    return sorted(out,key=lambda p:p['depth'])


def fit(pos,s,centre,box,radius):
    d=wrap(pos-centre,box);take=np.linalg.norm(d,axis=1)<radius
    x,y=d[take].T;v=s[take]
    design=np.column_stack([np.ones(len(x)),x,y,x*x,x*y,y*y])
    co,_,rank,_=np.linalg.lstsq(design,v,rcond=None)
    H=np.array([[2*co[3],co[4]],[co[4],2*co[5]]])
    offset=np.linalg.solve(H,-co[1:3]);eig=np.linalg.eigvalsh(H)
    assert rank==6 and eig.min()>0 and np.linalg.norm(offset)<radius
    return np.mod(centre+offset,box),dict(nodes=len(x),min_H_eigenvalue=float(eig.min()),
        condition_H=float(np.linalg.cond(H)),offset=float(np.linalg.norm(offset)),
        rmse=float(np.sqrt(np.mean((design@co-v)**2))))


def main():
    dec=json.loads((P0/'decisions.json').read_text());assert dec['proceed'] is False
    box=np.array([dec['Lx'],dec['Ly']]);A=float(np.prod(box));rho=18901/A
    rng=np.random.default_rng(seed(6,9,0,0,0))
    initial=np.column_stack([rng.uniform(0,box[0],18901),rng.uniform(0,box[1],18901)])
    rng=np.random.default_rng(seed(6,9,3,1,1));life0=rng.exponential(5.,18901)
    n=int(rng.poisson(.2*18901*np.exp(.4)*100))
    tb=np.sort(rng.random(n)*100)
    x,y,u,ell=rng.random(n)*box[0],rng.random(n)*box[1],rng.random(n),rng.exponential(5.,n)
    h=hashlib.sha256(np.array([18901],'<i8').tobytes())
    for a in (life0,tb,x,y,u,ell):h.update(np.asarray(a,dtype='<f8').tobytes())
    coords=np.vstack([initial,np.column_stack([x,y])]);death=np.concatenate([life0,tb+ell])
    birth=np.concatenate([np.zeros(18901),tb]);radii=np.array([1,1.5,2,3,5,8])
    ann=np.array([0,.5,1,1.5,2,2.5,3,4,5,6,8])
    specs={'d0.2_dt0.05':.2,'d0.2_dt0.025_forced':.2,'d0.1_dt0.05':.1,
           'd0.1_dt0.025_forced':.1,'d0.05_dt0.05':.05}
    runs={};out=dict(review_commit='ba765876dfde029388d5769a0cd690d25b38aeee',
        original_proceed=False,candidate_count=n,candidate_sha256=h.hexdigest(),runs={})
    for tag,delta in specs.items():
        with np.load(P0/f'check_{tag}.npz') as z:d={k:z[k] for k in z.files}
        with np.load(P0/f'accepted_{tag}.npz') as z:bits=np.unpackbits(z['bits'])[:int(z['n'])].astype(bool)
        assert n==int(d['n_candidates'])==len(bits)
        assert np.array_equal(bits,np.unpackbits(d['accepted_bits'])[:n])
        assert not np.any(bits & (tb+ell <= np.ceil(tb/delta)*delta))
        admitted=np.concatenate([np.ones(18901,bool),bits]);full=[];fit_centres={1.:[],1.5:[]};fit_details={1.:[],1.5:[]}
        centres=[];depth=[];d0fit={1.:[],1.5:[]}
        for i,t in enumerate(d['t']):
            cutoff=0. if t==0 else (round(t/delta)-1)*delta
            alive=admitted&(birth<=cutoff)&(death>cutoff)
            pos=coords[alive];full.append(pos)
            assert len(pos)==int(d['m_N'][i])
            a,b=map(int,d['off'][i:i+2]);sp,ss,rows,ids=d['pos'][a:b],d['s'][a:b],d['rows'][a:b],d['ids'][a:b]
            assert np.all(alive[ids]) and np.array_equal(sp,coords[ids])
            particles=detect(sp,ss,rows,-4.457581675528346,box)
            assert len(particles)==int(d['n_parts'][i])==1
            c=particles[0]['centre'];centres.append(c);depth.append(particles[0]['depth'])
            assert np.linalg.norm(wrap(c-d['centre'][i],box))<1e-12
            rr=np.linalg.norm(wrap(pos-c,box),axis=1)
            counts=np.array([(rr<R).sum() for R in radii])
            assert np.array_equal(counts,d['m_N_R'][i])
            assert np.allclose(counts-rho*np.pi*radii*radii,d['m_D0'][i],rtol=0,atol=1e-10)
            ac=np.array([((rr>=lo)&(rr<hi)).sum() for lo,hi in zip(ann[:-1],ann[1:])])
            assert np.array_equal(ac,d['m_ann_count'][i])
            for radius in fit_centres:
                fc,details=fit(sp,ss,c,box,radius);fit_centres[radius].append(fc);fit_details[radius].append(details)
                count=np.sum(np.linalg.norm(wrap(pos-fc,box),axis=1)<1.5)
                d0fit[radius].append(float(count-rho*np.pi*1.5**2))
        window=(d['t']>=50)&(d['t']<=100)
        runs[tag]=dict(d=d,full=full,centre=np.array(centres),depth=np.array(depth),
            fit_centres={r:np.array(v) for r,v in fit_centres.items()})
        out['runs'][tag]=dict(delta=delta,counts=d['n_parts'].tolist(),
            registered_D0_bar=float(d['m_D0'][window,1].mean()),P0_bar=float(d['m_P0'][window,1].mean()),
            max_recorded_phiG=float(d['tick_max_phiG'].max()),
            quadratic_fits={str(r):dict(details=fit_details[r],centres=fit_centres[r],
                D0_series=d0fit[r],D0_window_mean=float(np.array(d0fit[r])[window].mean()),
                max_absolute_count_change=float(np.max(np.abs(np.array(d0fit[r])-d['m_D0'][:,1])))) for r in fit_centres})
        print(tag,'counts / positions / sparse detection agree',flush=True)
    out['comparisons']={}
    for label,a,b in [('step_d0.2','d0.2_dt0.05','d0.2_dt0.025_forced'),
                      ('step_d0.1','d0.1_dt0.05','d0.1_dt0.025_forced'),
                      ('tick_0.2_0.1','d0.2_dt0.05','d0.1_dt0.05'),
                      ('tick_0.1_0.05','d0.1_dt0.05','d0.05_dt0.05')]:
        ra,rb=runs[a],runs[b]
        dis=np.linalg.norm(wrap(ra['centre']-rb['centre'],box),axis=1)
        dd=np.abs(ra['depth']-rb['depth'])/np.abs(ra['depth'])
        ref=out['runs'][b]['P0_bar'];raw_a=out['runs'][a]['registered_D0_bar'];raw_b=out['runs'][b]['registered_D0_bar']
        report=dict(centre_distances=dis.tolist(),max_centre=float(dis.max()),max_relative_depth=float(dd.max()),
            D0_window_means=[raw_a,raw_b],reference_P0=ref,
            relative_D0_difference=float(abs(raw_a-raw_b)/abs(ref)),part_b_pass=bool(abs(raw_a-raw_b)<=.25*abs(ref)),
            fitted_centres={})
        for radius in [1.,1.5]:
            fd=np.linalg.norm(wrap(ra['fit_centres'][radius]-rb['fit_centres'][radius],box),axis=1)
            fa=out['runs'][a]['quadratic_fits'][str(radius)]['D0_window_mean']
            fb=out['runs'][b]['quadratic_fits'][str(radius)]['D0_window_mean']
            report['fitted_centres'][str(radius)]=dict(distances=fd.tolist(),max_distance=float(fd.max()),
                D0_window_means=[fa,fb],difference=float(fa-fb))
        if label.startswith('tick'):
            ca=(ra['centre']+.5*wrap(rb['centre']-ra['centre'],box))%box
            common=[]
            for i,c in enumerate(ca):
                common.append([float(np.sum(np.linalg.norm(wrap(r['full'][i]-c,box),axis=1)<1.5)-rho*np.pi*1.5**2)
                               for r in [ra,rb]])
            common=np.array(common);w=(ra['d']['t']>=50)&(ra['d']['t']<=100)
            report['common_registered_midpoint_D0_window_means']=common[w].mean(0).tolist()
        out['comparisons'][label]=report
    base=runs['d0.1_dt0.05']['d'];ser=dict(zip(base['t'],base['m_D0'][:,1]))
    b1=float(np.mean([ser[60],ser[70]]));b2=float(np.mean([ser[90],ser[100]]))
    ref=out['runs']['d0.05_dt0.05']['P0_bar']
    out['amended_existing_checks']=dict(max_recorded_phiG=max(v['max_recorded_phiG'] for v in out['runs'].values()),
        margin_limit=.07,B1=b1,B2=b2,absolute_drift=abs(b2-b1),drift_limit=.25*abs(ref),
        chosen_window=[100,200] if abs(b2-b1)<=.25*abs(ref) else [150,250])
    for a,b,key in [('step_d0.2',0,'step_check'),('step_d0.1',1,'step_check'),('tick_0.2_0.1',0,'tick_check'),('tick_0.1_0.05',1,'tick_check')]:
        old=dec[key][b] if key=='step_check' else dec[key][b]['part_a']
        assert abs(out['comparisons'][a]['max_centre']-old['max_centre_diff'])<1e-12
        assert abs(out['comparisons'][a]['max_relative_depth']-old['max_depth_rel'])<1e-12
    def serial(v):
        if isinstance(v,np.ndarray):return v.tolist()
        raise TypeError(str(type(v)))
    (ROOT/'record_checks.json').write_text(json.dumps(out,indent=2,default=serial)+'\n')
    print(json.dumps(dict(comparisons=out['comparisons'],existing_checks=out['amended_existing_checks']),indent=2))


if __name__=='__main__':main()
