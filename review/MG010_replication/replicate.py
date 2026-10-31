"""Independent measurements of hash-verified MG4 fields; no trajectory re-integration.
Uses registered graph builders only. Detector, links, smoothing, C3(b), Laplacian,
and high-precision balance arithmetic below are written separately.
Run: python review/MG010_replication/replicate.py
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import sys,json,math,gc
from decimal import Decimal,localcontext
import numpy as np
import scipy.sparse as sp
from scipy.spatial import cKDTree
from scipy.interpolate import CubicSpline
ROOT=Path(__file__).resolve().parent
IN=ROOT/'inputs'
sys.path.insert(0,str(IN/'src'))
import run_mg4_partA as builder
import graph_ode as original_ode
from lattice import steady_state,P17
DEC=json.loads((IN/'results/MG4/part0/decisions.json').read_text())
EV=json.loads((IN/'results/MG4/evaluation.json').read_text())
THETA=json.loads((IN/'results/MG4/partA/threshold.json').read_text())['theta']
CONT=dict(np.load(IN/'results/MG4/part0/continuum2d.npz'))

def delta(x,box):
    return x-box*np.rint(x/box)

def own_lap(g,h):
    i,j=g['I'],g['J']; n=g['N']; deg=np.bincount(np.r_[i,j],minlength=n)
    return sp.csr_matrix((np.r_[np.ones(2*len(i)),-deg],
                         (np.r_[i,j,np.arange(n)],np.r_[j,i,np.arange(n)])),shape=(n,n))/h**2

def detect(pos,box,tree,near,near_dist,s):
    ids=np.arange(len(s))
    # A subset of neighbours can only reject a genuine non-minimum. Check every
    # survivor against its full radius-1 neighbourhood, including index tie-break.
    lower=(s[near]<s[:,None])|((s[near]==s[:,None])&(near<ids[:,None]))
    cand=np.flatnonzero((s<THETA)&~np.any(lower&(near_dist<=1),axis=1))
    keep=[]
    for i in cand:
        nb=np.asarray(tree.query_ball_point(pos[i],1.0))
        if not np.any((s[nb]<s[i])|((s[nb]==s[i])&(nb<i))):keep.append(int(i))
    parent=list(range(len(keep)))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for a in range(len(keep)):
        for b in range(a):
            if np.linalg.norm(delta(pos[keep[a]]-pos[keep[b]],box))<=1.5:
                parent[find(a)]=find(b)
    groups={}
    for a,i in enumerate(keep):groups.setdefault(find(a),[]).append(i)
    centres=[]
    for group in groups.values():
        i=min(group,key=lambda j:(s[j],j))
        nb=np.asarray(tree.query_ball_point(pos[i],.6)); w=np.maximum(THETA-s[nb],0)
        centres.append((pos[i]+np.sum(w[:,None]*delta(pos[nb]-pos[i],box),axis=0)/sum(w))%box)
    return np.asarray(centres).reshape(-1,2)

def c3b(g,h,records,tree):
    pos=g['pos'];box=np.array([g['Lx'],g['Ly']]);L=own_lap(g,h)
    vals=[]; masks=[]
    for r in records:
        rn=np.linalg.norm(delta(pos-r['centre'],box),axis=1);masks.append(rn<6)
        for s in range(3):
            fv=CubicSpline(CONT['r'],CONT['phi'][s],bc_type=((1,0),'not-a-knot'))(rn)
            lc=CubicSpline(CONT['r'],CONT['lap'][s],bc_type=((1,0),'not-a-knot'))(rn)
            vals.extend([L@fv,lc])
    vals=np.column_stack(vals);smooth=np.empty_like(vals)
    # Independent chunked Gaussian smoother; avoids building a large full matrix.
    for a in range(0,g['N'],256):
        b=min(a+256,g['N']); nbs=tree.query_ball_point(pos[a:b],1.6)
        lens=np.array([len(x) for x in nbs]); rows=np.repeat(np.arange(b-a),lens)
        cols=np.concatenate(nbs);d=delta(pos[cols]-pos[a:b][rows],box)
        w=np.exp(-np.sum(d*d,axis=1)/.32)
        w/=np.bincount(rows,weights=w,minlength=b-a)[rows]
        S=sp.csr_matrix((w,(rows,cols)),shape=(b-a,g['N']))
        smooth[a:b]=S@vals
    out={}
    for k,r in enumerate(records):
        out[r['key']]={}
        for s,name in enumerate('GXY'):
            a,b=smooth[:,6*k+2*s:6*k+2*s+2].T;m=masks[k]
            out[r['key']][name]=float(np.linalg.norm((a-b)[m])/np.linalg.norm(b[m]))
    return out

def decimal_balance(g,h,y):
    # Windows longdouble has only 53 mantissa bits. Use actual 40-digit Decimal
    # edge fluxes and node reactions instead of assuming platform extended precision.
    with localcontext() as ctx:
        ctx.prec=40
        D=Decimal.from_float;zero=Decimal(0);one=Decimal(1);gg=D(.1)
        fields=[[D(float(x)) for x in row] for row in y]
        sums=[]; deriv=[];hh=D(float(h))**2
        x0=D(float(steady_state(P17)[1]));y0=D(float(steady_state(P17)[2]))
        laps=[]
        for f in fields:
            ll=[zero]*g['N']
            for i,j in zip(g['I'],g['J']):
                flux=(f[j]-f[i])/hh;ll[i]+=flux;ll[j]-=flux
            laps.append(ll);sums.append(sum(f,zero))
        dg=dx=dy=zero
        for i in range(g['N']):
            G,X,Y=(f[i] for f in fields)
            reaction=(x0+X)**2*(y0+Y)-x0*x0*y0
            dg+=laps[0][i]-G+gg*X
            dx+=laps[1][i]+G-Decimal(30)*X+reaction
            dy+=Decimal(12)*laps[2][i]+Decimal(29)*X-reaction
        return dict(r1=str(dg+sums[0]-gg*sums[1]),r2=str(dg+dx+dy+(one-gg)*sums[1]),precision_digits=40)

def main():
    records=[json.loads(p.read_text()) for p in sorted((IN/'results/MG4/partA/runs').glob('*.json'))]
    groups={}
    for r in records:groups.setdefault((r['fam'],r['slot'],r['torus']),[]).append(r)
    out=dict(commit=json.loads((ROOT/'download_manifest.json').read_text())['commit'],
             longdouble_mantissa_bits=int(np.finfo(np.longdouble).nmant+1),runs={},c3b={},balances={})
    for key,rr in groups.items():
        print('Reconstruct',key,flush=True)
        g,h=builder.get_graph(dict(fam=key[0],slot=key[1],torus=key[2]),DEC)
        box=np.array([g['Lx'],g['Ly']]);pos=g['pos']%box;tree=cKDTree(pos,boxsize=box)
        nd,near=tree.query(pos,k=16)
        assert g['N']==rr[0]['N'] and abs(h-rr[0]['h'])<1e-12
        for r in rr:
            f=np.load(IN/f"results/MG4/partA/fields_eval/{r['key']}.npz")
            checks=[detect(pos,box,tree,near,nd,s.astype(float)) for s in f['checks_sX']]
            counts=[len(c) for c in checks]
            link=None
            if max(counts)>=2:cl='multiplied'
            elif counts[-1]==0:cl='decayed'
            elif min(counts)==1:
                link=float(max(np.linalg.norm(delta(checks[k][0]-checks[k-1][0],box)) for k in range(1,len(checks))))
                cl='single' if link<=1 else 'unresolved'
            else:cl='unresolved'
            out['runs'][r['key']]=dict(cls=cl,counts=counts,max_link=link)
            if r['key'] in EV['balance_failures']:
                y=f['final'].astype(float);rhs,_,_=original_ode.make_system(g,h)
                raw=original_ode.balances(y.ravel(),rhs(300,y.ravel()),g['N'])
                precise=decimal_balance(g,h,y)
                out['balances'][r['key']]=dict(saved_float32_state_decimal=precise,float64_r1=float(raw['r1']),float64_r2=float(raw['r2']))
        if key[2]!='h03':out['c3b'].update(c3b(g,h,rr,tree))
        (ROOT/'replication.json').write_text(json.dumps(out,indent=2))
        print('Done',key,'counts',[(r['key'],out['runs'][r['key']]['counts'][-1]) for r in rr],flush=True)
        del g,tree;gc.collect()
    out['class_mismatches']=[k for k,v in out['runs'].items() if v['cls']!=EV['classes'][k]]
    out['c3b_medians']={fam:{s:float(np.median([v[s] for k,v in out['c3b'].items() if k.startswith(fam)])) for s in 'GXY'} for fam in 'DRT'}
    out['c3b_max_absolute_discrepancy']=max(abs(v[s]-(EV['T_C3b'] if k.startswith('T') else EV['families'][k[0]]['C3b'])['per_run'][k][s]) for k,v in out['c3b'].items() for s in 'GXY')
    # Analytic triangular dispersion and independent continuum Fourier-mode count.
    lx,ly=DEC['torus']['Lx'],DEC['torus']['Ly'];a=.2/math.sqrt(1.5)
    m,n=np.meshgrid(np.arange(-30,31),np.arange(-30,31));qx=2*np.pi*m/lx;qy=2*np.pi*n/ly
    nu=(6-2*np.cos(qx*a)-4*np.cos(qx*a/2)*np.cos(qy*a*math.sqrt(3)/2))/.04
    ratios=[int(np.sum(nu<=k*k*(1+1e-12)))/int(np.sum(qx*qx+qy*qy<=k*k*(1+1e-12))) for k in np.linspace(.9,1.1,21)*2*np.pi/3.08]
    out['analytic_T_C3a']=float(np.mean(ratios))
    out['spectral_panel_means_from_saved_ratios']={fam:float(np.mean([v['ratios'] for v in DEC['spectra'].values() if v['fam']==fam])) for fam in ['D','R','T','N3']}
    # Full precision Part 0 reference can support an independent stationarity check.
    g,h=builder.get_graph(dict(fam='T',slot=0),DEC);y=np.load(IN/'results/MG4/part0/steps/Tstate_122x142.npy').reshape(3,-1)
    L=own_lap(g,h);G,X,Y=y;x0,y0=steady_state(P17)[1:];rx=(x0+X)**2*(y0+Y)-x0*x0*y0
    dy=np.vstack([L@G-G+.1*X,L@X+G-30*X+rx,12*(L@Y)+29*X-rx])
    out['T_reference']=dict(depth_X=float(X.min()),max_rhs=float(abs(dy).max()),stored_dtype=str(y.dtype))
    (ROOT/'replication.json').write_text(json.dumps(out,indent=2))
    print(json.dumps({k:v for k,v in out.items() if k not in ['runs','c3b','balances']},indent=2),flush=True)

if __name__=='__main__':main()
