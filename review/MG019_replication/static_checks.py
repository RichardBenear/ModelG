"""MG019: static checks on MG4c states; no MG6 integration or event stream.

Uses the inspected graph builder and MG4 detector for coordinate reconstruction.
Statistics and coupling-function arithmetic are calculated here independently.
"""
from pathlib import Path
import json,sys,hashlib
import numpy as np
from scipy.stats import t
ROOT=Path(__file__).resolve().parent
IN=ROOT/'inputs'
if not IN.exists():IN=ROOT.parent/'MG014_replication/inputs'
m=json.loads((ROOT/'download_manifest.json').read_text())
for e in m['files']:
    data=(IN/e['path']).read_bytes()
    assert hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==e['sha']
sys.path.insert(0,str(IN/'src'))
import run_mg4c as rm
import measure_graph as mg
dec=json.loads((IN/'results/MG4/part0/decisions.json').read_text())
out=dict(reference_commit=m['commit'],states=[],small_n_statistics=dict(
    independent_run_summaries=4,df=3,t_statistic=3,
    one_sided_p=float(t.sf(3,3)),two_sided_p=float(2*t.sf(3,3)),
    t_critical_two_sided_95=float(t.ppf(.975,3))))
for slot,pl in [(0,0),(0,1),(1,0),(1,1)]:
    job=rm.job(28,'D',slot,pl);graph,h,_=rm.graph_for(job,dec)
    F=np.load(IN/f"results/MG4c/fields/{job['key']}.npz")['final'];G=F[0]
    L=np.array([graph['Lx'],graph['Ly']]);A=float(np.prod(L));rho=len(G)/A
    s=mg.smoother(graph)@F[1]
    p=mg.particles(graph,s,-4.457581675528346)
    assert len(p)==1;centre=p[0]['centre']
    d=graph['pos']-centre;d-=L*np.rint(d/L);r=np.linalg.norm(d,axis=1)
    a=dict(key=job['key'],h=h,L=L.tolist(),area=A,rho=rho,G_min=float(G.min()),G_max=float(G.max()),
        G_node_sum=float(G.sum()),node_mean_area=A/len(G),h_squared=h*h,
        minimum_kappa_for_F2_clipping=-1/float(G.min()),
        F2_k4_negative_linear_nodes=int(np.sum(1+4*G<0)),
        F2_k4_raw_min=float(np.min(1+4*G)),
        F2_k4_clipped_total_excess=float(np.sum(np.maximum(0,1+4*G)-1)),
        well_integral_equal_area_start_estimate=float(np.maximum(-G,0).sum()/rho),
        halo_integral_equal_area_start_estimate=float(np.maximum(G,0).sum()/rho),F1={})
    for k in [2,4]:
        ex=np.expm1(k*G)
        a['F1'][str(k)]=dict(min_factor=float(np.exp(k*G.min())),max_factor=float(np.exp(k*G.max())),
            total_excess_equal_area_start_estimate=float(ex.sum()),radii={})
        for R in [1,1.5,2,3,5,8]:
            M=r<R;raw=float(ex[M].sum())
            # The proposal's initial equal-area quadrature is only a design
            # estimate on near-uniform starts; it is not valid on biased MG6 nodes.
            centred=raw-np.pi*R*R/A*float(ex.sum())
            a['F1'][str(k)]['radii'][str(R)]=dict(raw_target_count_excess=raw,
                excess_after_subtracting_target_global_density=centred,
                poisson_snapshot_sd=float(np.sqrt(rho*np.pi*R*R)))
    a['F1_k4_to_k2_R1p5_ratio']=a['F1']['4']['radii']['1.5']['raw_target_count_excess']/a['F1']['2']['radii']['1.5']['raw_target_count_excess']
    out['states'].append(a)
(ROOT/'static_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
