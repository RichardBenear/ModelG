"""Static MG6 Draft 2 checks: grid tiling and falloff-bin availability.
No MG6 integration, events, outcomes or evolving-mesh trajectories.
Uses the pinned MG014 graph builder and detector to reconstruct stored states.
"""
from pathlib import Path
import json,sys
import numpy as np
from scipy.spatial import Delaunay
ROOT=Path(__file__).resolve().parent
IN=ROOT/'inputs'
if not IN.exists():IN=ROOT.parent/'MG014_replication/inputs'
sys.path.insert(0,str(IN/'src'))
import run_mg4c as rm
import measure_graph as mg
dec=json.loads((IN/'results/MG4/part0/decisions.json').read_text())
out=dict(review_commit='92d8f2b8505526f1e6e559d450352682f574d6a3',states=[])
q=-np.expm1(-.2*.2)/(.2*.2)
out['ephemeral_candidates']=dict(nu=.2,delta=.2,
    mean_probability_survive_to_birth_tick=float(q),
    mean_probability_ephemeral=float(1-q),
    note='Uniform candidate time within a tick, independent Exp(nu) lifetime. '
         'Acceptance-rate checks must exclude ephemeral candidates or include this survival factor.')
edges=np.array([0,.5,1,1.5,2,2.5,3,4,5,6,8])
for slot,pl in [(0,0),(0,1),(1,0),(1,1)]:
    job=rm.job(28,'D',slot,pl);g,h,_=rm.graph_for(job,dec)
    F=np.load(IN/f"results/MG4c/fields/{job['key']}.npz")['final']
    box=np.array([g['Lx'],g['Ly']]);A=float(np.prod(box));rho=len(g['pos'])/A
    centre=mg.particles(g,mg.smoother(g)@F[1],-4.457581675528346)[0]['centre']
    shape=np.rint(box/.05).astype(int);spacing=box/shape;weight=A/np.prod(shape)
    xx=(np.arange(shape[0])+.5)*spacing[0];yy=(np.arange(shape[1])+.5)*spacing[1]
    grid=np.stack(np.meshgrid(xx,yy,indexing='ij'),axis=-1).reshape(-1,2)
    shifts=np.array([(x,y) for x in (-1,0,1) for y in (-1,0,1)])*box
    q=np.concatenate([g['pos']+s for s in shifts]);tri=Delaunay(q)
    si=tri.find_simplex(grid);assert np.all(si>=0)
    trans=tri.transform[si];b=np.einsum('ijk,ik->ij',trans[:,:2],grid-trans[:,2])
    w=np.column_stack([b,1-b.sum(1)])
    G=np.sum(w*F[0][tri.simplices[si]%len(g['pos'])],axis=1)
    d=grid-centre;d-=box*np.rint(d/box);r=np.linalg.norm(d,axis=1)
    means=[]
    for lo,hi in zip(edges[:-1],edges[1:]):
        use=(r>=lo)&(r<hi)
        means.append(float(G[use].mean()))
    # Locate the core's first sign change more finely; this only describes
    # the stored start, not the later coupled particle's eventual width.
    fine=np.arange(0,3.0001,.1)
    fm=np.array([G[(r>=lo)&(r<hi)].mean() for lo,hi in zip(fine[:-1],fine[1:])])
    crossing=next(i for i in range(len(fm)-1) if fm[i]<0 and fm[i+1]>=0)
    mids=(fine[:-1]+fine[1:])/2
    zero=float(mids[crossing]-fm[crossing]*(mids[crossing+1]-mids[crossing])/(fm[crossing+1]-fm[crossing]))
    bins=[[float(lo),float(hi)] for lo,hi in zip(edges[:-1],edges[1:]) if lo>=1 and hi<=zero]
    arithmetic=[]
    for name,n in [('floor',np.floor(box/.05).astype(int)),('round',np.rint(box/.05).astype(int)),('ceil',np.ceil(box/.05).astype(int))]:
        nominal=float(np.prod(n)*.0025)
        arithmetic.append(dict(rule=name,shape=n.tolist(),area_if_each_cell_00025=nominal,
                               area_error=nominal-A,equivalent_count_error=rho*(nominal-A)))
    out['states'].append(dict(key=job['key'],torus_lengths=box.tolist(),actual_area=A,
        length_over_005=(box/.05).tolist(),grid_examples=arithmetic,
        commensurate_midpoint_grid=dict(shape=shape.tolist(),spacing=spacing.tolist(),cell_area=float(weight),
            integrated_constant=float(len(grid)*weight)),annulus_edges=edges.tolist(),annulus_mean_G=means,
        first_G_zero_estimate=zero,eligible_falloff_annuli_from_r1=bins,
        possible_annuli_required_by_draft=3))
    print(job['key'],'first G zero',zero,'eligible full bins from r=1',len(bins),flush=True)
(ROOT/'draft2_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print('No MG6 dynamics performed.')
