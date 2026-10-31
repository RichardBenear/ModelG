"""MG4c independent saved-field replication; no trajectory integration.
Registered builders reconstruct graph positions/edges. Measurements, FEM stiffness,
RHS, balances and refinement reconstruction below are independently written.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import sys,json,csv,gc
import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree,Delaunay
ROOT=Path(__file__).resolve().parent
IN=ROOT/'inputs'
sys.path.insert(0,str(IN/'src'))
import run_mg4c as registered
import graphs as original_graphs
EV=json.loads((IN/'results/MG4c/evaluation.json').read_text())
DEC=json.loads((IN/'results/MG4/part0/decisions.json').read_text())
TH29=-5.117354421620525

def delta(d,box):return d-box*np.rint(d/box)

def unweighted(g,h):
    i,j=g['I'],g['J'];n=g['N'];degree=np.bincount(np.r_[i,j],minlength=n)
    return sp.csr_matrix((np.r_[np.ones(2*len(i)),-degree],
        (np.r_[i,j,np.arange(n)],np.r_[j,i,np.arange(n)])),shape=(n,n))/h**2

def fem(g):
    # Unlike the implementation's centroid-in-tile selection, keep triangles
    # touching a central vertex and deduplicate by their physical node triples.
    n=g['N'];box=np.array([g['Lx'],g['Ly']]);pos=g['pos']
    shifts=np.array([(i,j) for i in [-1,0,1] for j in [-1,0,1]])*box
    q=np.concatenate([pos+s for s in shifts]);tri=Delaunay(q).simplices
    tri=tri[np.any(tri//n==4,axis=1)]
    _,ix=np.unique(np.sort(tri%n,axis=1),axis=0,return_index=True)
    tri=tri[ix];faces=tri%n;p=q[tri]
    u,v=p[:,1]-p[:,0],p[:,2]-p[:,0]
    det=u[:,0]*v[:,1]-u[:,1]*v[:,0];area=np.abs(det)/2
    grads=np.stack([np.column_stack([p[:,1,1]-p[:,2,1],p[:,2,0]-p[:,1,0]]),
        np.column_stack([p[:,2,1]-p[:,0,1],p[:,0,0]-p[:,2,0]]),
        np.column_stack([p[:,0,1]-p[:,1,1],p[:,1,0]-p[:,0,0]])],axis=1)/det[:,None,None]
    elem=-area[:,None,None]*np.einsum('fai,fbi->fab',grads,grads)
    rows=np.repeat(faces,3,axis=1).ravel();cols=np.tile(faces,(1,3)).ravel()
    K=sp.csr_matrix((elem.ravel(),(rows,cols)),shape=(n,n))
    mass=np.bincount(faces.ravel(),weights=np.repeat(area/3,3),minlength=n)
    L=sp.diags(1/mass)@K
    edge=np.sort(np.concatenate([faces[:,[0,1]],faces[:,[0,2]],faces[:,[1,2]]]),axis=1)
    unique,count=np.unique(edge,axis=0,return_counts=True)
    graph_edges=np.column_stack([g['I'],g['J']])
    checks=dict(n_faces=len(faces),expected_faces=2*n,total_area=float(area.sum()),
        box_area=float(np.prod(box)),min_node_area=float(mass.min()),
        edges_match=bool(np.array_equal(unique,graph_edges)),two_faces_each=bool(np.all(count==2)),
        weighted_symmetry_max=float(abs(K-K.T).max()),
        constant_rhs_max=float(abs(L@np.ones(n)).max()))
    source=original_graphs.cotan(g)
    diff=L-source['L']
    checks.update(operator_max_absolute_difference=float(abs(diff).max()),
        operator_max_entry=float(abs(source['L']).max()),
        mass_max_difference=float(abs(mass-source['A']).max()))
    return L.tocsr(),mass,checks

def particles(pos,box,tree,nd,near,s,theta):
    idx=np.arange(len(s));v=s[near]
    lower=(v<s[:,None])|((v==s[:,None])&(near<idx[:,None]))
    possible=np.flatnonzero((s<theta)&~np.any(lower&(nd<=1),axis=1));cand=[]
    for i in possible:
        nb=np.asarray(tree.query_ball_point(pos[i],1.0))
        if not np.any((s[nb]<s[i])|((s[nb]==s[i])&(nb<i))):cand.append(i)
    parent=list(range(len(cand)))
    def root(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for a in range(len(cand)):
        for b in range(a):
            if np.linalg.norm(delta(pos[cand[a]]-pos[cand[b]],box))<=1.5:parent[root(a)]=root(b)
    groups={}
    for a,i in enumerate(cand):groups.setdefault(root(a),[]).append(i)
    centres=[]
    for members in groups.values():
        i=min(members,key=lambda j:(s[j],j));nb=np.asarray(tree.query_ball_point(pos[i],.6))
        w=np.maximum(theta-s[nb],0)
        centres.append((pos[i]+(w[:,None]*delta(pos[nb]-pos[i],box)).sum(0)/w.sum())%box)
    return np.asarray(centres).reshape(-1,2)

def classify(pos,box,tree,nd,near,checks,theta):
    tracks=[particles(pos,box,tree,nd,near,s.astype(float),theta) for s in checks]
    counts=[len(p) for p in tracks];link=None
    if max(counts)>=2:cl='multiplied'
    elif counts[-1]==0:cl='decayed'
    elif min(counts)==1:
        link=float(max(np.linalg.norm(delta(tracks[k][0]-tracks[k-1][0],box)) for k in range(1,len(tracks))))
        cl='single' if link<=1 else 'unresolved'
    else:cl='unresolved'
    return dict(cls=cl,counts=counts,max_link=link),tracks[-1][0] if cl=='single' else None

def smooth(pos,box,tree,values):
    out=np.empty_like(values)
    for a in range(0,len(pos),256):
        b=min(a+256,len(pos));nb=tree.query_ball_point(pos[a:b],1.6)
        rows=np.repeat(np.arange(b-a),[len(x) for x in nb]);cols=np.concatenate(nb)
        d=delta(pos[cols]-pos[a:b][rows],box);w=np.exp(-np.sum(d*d,axis=1)/.32)
        w/=np.bincount(rows,weights=w,minlength=b-a)[rows]
        out[a:b]=sp.csr_matrix((w,(rows,cols)),shape=(b-a,len(pos)))@values
    return out

def rhs_and_balances(L,y,b,mass):
    G,X,Y=y;x0=14/(1-.1);y0=x0*b/x0**2
    rx=(x0+X)**2*(y0+Y)-x0**2*y0
    d=np.vstack([L@G-G+.1*X,L@X+G-(1+b)*X+rx,12*(L@Y)+b*X-rx])
    I=np.sum(mass*y,axis=1);di=np.sum(mass*d,axis=1)
    r1=di[0]+I[0]-.1*I[1];r2=sum(di)+.9*I[1]
    M1=abs(di[0])+abs(I[0])+.1*abs(I[1]);M2=sum(abs(di))+.9*abs(I[1])
    return dict(max_rhs=float(abs(d).max()),r1=float(r1),r2=float(r2),
        M1=float(M1),M2=float(M2),ok=bool(abs(r1)<=1e-10+1e-9*M1 and abs(r2)<=1e-10+1e-9*M2))

def shape(pos,box,centre,X,sX):
    radius=np.linalg.norm(delta(pos-centre,box),axis=1)
    bins=np.r_[0,np.arange(.5,8.000001,.25)];rc=(bins[:-1]+bins[1:])/2
    count=np.histogram(radius,bins)[0]
    def profile(v):
        val=np.histogram(radius,bins,weights=v)[0]
        return np.divide(val,count,out=np.full(len(val),np.nan),where=count>0)
    p,ps=profile(X),profile(sX);ok=np.isfinite(p);r,f=rc[ok],p[ok];z=[]
    for a in range(len(r)-1):
        if f[a]==0:z.append(r[a])
        elif f[a]*f[a+1]<0:z.append(r[a]-f[a]*(r[a+1]-r[a])/(f[a+1]-f[a]))
    z=np.asarray([x for x in z if 1<x<7]);lam=float(2*np.mean(np.diff(z))) if len(z)>=2 else None
    return dict(depth_s=float(ps[0]),depth_u=float(p[0]),lam=lam)

def main():
    recs=[json.loads(p.read_text()) for p in sorted((IN/'results/MG4c/runs').glob('*.json'))]
    out=dict(commit=json.loads((ROOT/'download_manifest.json').read_text())['commit'],runs={},tables={},fem={})
    groups={};thresholds={};shapes={};cls={}
    for r in recs:groups.setdefault((r['kind'],r['slot']),[]).append(r)
    # T is processed first to derive the thresholds independently of the records.
    for group in sorted(groups,key=lambda k:0 if k[0]=='T' else 1):
        rr=groups[group];g,h,_=registered.graph_for(rr[0],DEC);pos=g['pos'];box=np.array([g['Lx'],g['Ly']])
        tree=cKDTree(pos%box,boxsize=box);nd,near=tree.query(pos,k=16)
        arrays=[dict(np.load(IN/f"results/MG4c/fields/{r['key']}.npz")) for r in rr]
        sm=smooth(pos,box,tree,np.column_stack([f['final'][1] for f in arrays]))
        L=unweighted(g,h);wL,wMass,fc=fem(g) if group[0]=='T' or group[1]==0 else (None,None,None)
        if fc is not None:out['fem'][str(group)]=fc
        print('Graph',group,'runs',len(rr),flush=True)
        for k,(r,f) in enumerate(zip(rr,arrays)):
            minimum=float(sm[:,k].min());theta=TH29;eligible=None
            if group[0]=='T' and r['op'] is None:
                e,_=classify(pos,box,tree,nd,near,f['checks_sX'],-1)
                eligible=e['cls']=='single' and minimum<=-2
                if eligible:
                    candidate=.5*minimum;ct,_=classify(pos,box,tree,nd,near,f['checks_sX'],candidate)
                    eligible=ct['cls']=='single'
                    if eligible:theta=candidate
                thresholds[r['b']]=(theta,eligible,minimum)
            elif group[0]=='D' and r['op'] is None:theta=thresholds[r['b']][0]
            measured,centre=classify(pos,box,tree,nd,near,f['checks_sX'],theta)
            fixed,_=classify(pos,box,tree,nd,near,f['checks_sX'],TH29)
            bal=rhs_and_balances(wL if r['op'] else L,f['final'],r['b'],wMass if r['op'] else np.ones(g['N']))
            measured.update(theta=theta,smoothed_min=minimum,
                minimum_record_difference=minimum-r['min_smoothed_final'],
                final_smoothed_storage_max_error=float(abs(sm[:,k]-f['checks_sX'][-1]).max()),
                fixed_threshold_class=fixed['cls'],rhs_balance=bal,
                rhs_record_difference=bal['max_rhs']-r['max_dphidt'],stationary=bal['max_rhs']<1e-8)
            out['runs'][r['key']]=measured;cls[r['key']]=measured['cls']
            if centre is not None:shapes[r['key']]=shape(pos,box,centre,f['final'][1],sm[:,k])
        (ROOT/'replication.json').write_text(json.dumps(out,indent=2))
        del arrays,sm,tree,g;gc.collect()
    mism=[]
    all_single=[]
    for b,(theta,elig,tmin) in sorted(thresholds.items(),reverse=True):
        key=f'b{b:.6f}';names=[f'{key}_D{s}_p{p}' for s,p in [(0,0),(0,1),(1,0),(1,1)]]
        classes=[cls[k] for k in names];ev=EV['tables'][key]
        if classes!=list(ev['D'].values()) or elig!=ev['T_eligible'] or abs(theta-ev['theta'])>1e-10:mism.append(key)
        out['tables'][key]=dict(b=b,theta=theta,T_eligible=elig,T_min=tmin,D=classes)
        if all(c=='single' for c in classes):
            all_single.append(b)
            if elig:
                ts=shapes[f'{key}_T_p0'];ds=[shapes[k] for k in names]
                depth=[d['depth_s']/ts['depth_s'] for d in ds];lam=[d['lam']/ts['lam'] for d in ds]
                out['tables'][key]['shape_vs_T']=dict(depth_s_ratio=depth,lam_ratio=lam)
    # Reconstruct Case A decisions without importing driver plan/bisection functions.
    def all_s(b):return all(c=='single' for c in out['tables'][f'b{b:.6f}']['D'])
    grid=[29]+[28.75-.25*i for i in range(8)]
    yes=[i for i,b in enumerate(grid) if b!=29 and all_s(b)];top=min(yes);bottom=top
    while bottom+1<len(grid) and all_s(grid[bottom+1]):bottom+=1
    edges=[]
    for inside,outside in [(grid[top],grid[top-1]),(grid[bottom],grid[bottom+1])]:
        visits=[]
        for _ in range(6):
            mid=round((inside+outside)/2,9);truth=all_s(mid);visits.append([mid,truth])
            if truth:inside=mid
            else:outside=mid
        edges.append(dict(bracket=sorted([inside,outside]),visits=visits))
    weightedT=np.load(IN/'results/MG4c/fields/w_b29.000000_T_p0.npz')['final']
    ref=np.load(IN/'results/MG4/part0/steps/Tstate_122x142.npy').reshape(weightedT.shape)
    ledger=list(csv.DictReader((IN/'results/MG4c/budget_ledger.csv').open(encoding='utf-8')))
    out.update(class_threshold_mismatches=mism,all_single_b=all_single,independent_edges=edges,
        all_single_matches=all_single==EV['all_single_b'],
        edges_match=[x['bracket'] for x in edges]==[x['bracket'] for x in EV['edges']],
        T_control_maxdiff=float(abs(weightedT-ref).max()),
        ledger=dict(rows=len(ledger),unique_keys=len({r['key'] for r in ledger}),
                    core_h=sum(float(r['wall_s']) for r in ledger)/3600),
        unweighted_D_class_counts={c:sum(v['cls']==c for k,v in out['runs'].items() if '_D' in k and not k.startswith('w_')) for c in ['single','multiplied','decayed','unresolved']},
        max_smoothed_min_discrepancy=max(abs(v['minimum_record_difference']) for v in out['runs'].values()),
        max_rhs_discrepancy=max(abs(v['rhs_record_difference']) for v in out['runs'].values()),
        independent_balance_failures=[k for k,v in out['runs'].items() if not v['rhs_balance']['ok']],
        stationary_singles=sum(v['cls']=='single' and v['stationary'] for k,v in out['runs'].items() if '_D' in k and not k.startswith('w_')))
    shape_errors=[]
    for k,row in out['tables'].items():
        if 'shape_vs_T' in row:
            for field in ['depth_s_ratio','lam_ratio']:
                shape_errors.extend(abs(a-b) for a,b in zip(row['shape_vs_T'][field],EV['tables'][k]['shape_vs_T'][field]))
    out['max_shape_ratio_discrepancy']=max(shape_errors)
    out['window_status']='W1' if all_single else 'other'
    out['operator_status']='O1' if cls['w_b29.000000_D0_p0']=='single' and cls['w_b29.000000_T_p0']=='single' and out['T_control_maxdiff']<=1e-6 else 'other'
    (ROOT/'replication.json').write_text(json.dumps(out,indent=2))
    assert not mism and out['all_single_matches'] and out['edges_match']
    print(json.dumps({k:v for k,v in out.items() if k not in ['runs','tables']},indent=2))

if __name__=='__main__':main()
