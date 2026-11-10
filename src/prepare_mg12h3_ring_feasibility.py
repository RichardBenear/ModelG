"""MG134 deterministic prepared-ring cost inventory; no RNG or integration."""
from mg12h3_population import ROOT,factory
from pathlib import Path
import numpy as np,json,time,hashlib
BASE=ROOT/'results/MG12H3/ring_feasibility'
def run():
    if (BASE/'inventory.json').exists():raise ValueError('Already inventoried')
    start=time.perf_counter();records=[];sources={}
    for origin in range(2):
      for arm in range(2):
        job=dict(id=f'f_o{origin}_s1_a{arm}_v1',kind='FINITE',origin=origin,start=1,arm=arm,setting=1,duration=15.)
        case,ga,gb,k,mode,reference,parent=factory(job);assert mode=='FINITE'
        sources[parent.relative_to(ROOT).as_posix()]=hashlib.sha256(parent.read_bytes()).hexdigest();n=case['initial'].shape[1]
        links=[('ActiveLink',i,j,1.) for i in range(n) for j in ((i-1)%n,(i+1)%n)]
        for kind,graph in [('ALink',ga),('BLink',gb)]:links.extend((kind,int(i),int(j),float(w)) for i,j,w in zip(graph['src'],graph['dst'],graph['weight']))
        for omega in (1,10):
          if time.perf_counter()-start>60:raise TimeoutError('Inventory deadline')
          exact=case['initial']*omega;assert np.all(exact>=0) and np.max(exact)<2**63;pop=np.floor(exact+.5).astype(np.int64);g,x,y,a,b,z,w=pop.astype(float);a0=case['a0']*omega;b0=case['b0']*omega;rates=np.array([14*omega*a/a0,g,.1*x,29*b*x/b0,x*(x-1)*y/omega**2,.9*x]);assert np.all(rates>=0);transport=np.zeros(5);channels=int(np.count_nonzero(rates));drift=np.zeros_like(pop);edge_rates=[]
          for kind,i,j,weight in links:
            species=range(3) if kind=='ActiveLink' else [3 if kind=='ALink' else 4]
            for s in species:
              coeff=case['diff'][s] if s<3 else k[s-3];rate=float(coeff*weight*pop[s,i]);transport[s]+=rate;edge_rates.append(rate);channels+=int(rate>0);drift[s,i]-=int(pop[s,i]);drift[s,j]+=int(pop[s,i])
          assert np.all(drift.sum(axis=1)==0);assert abs(sum(edge_rates)-transport.sum())<=1e-7*max(1,transport.sum())
          # Lower numeric payload only; Python tuples and boxed ints cost more.
          initial_rate=float(rates.sum()+transport.sum());row=dict(id=f'o{origin}_a{arm}_Omega{omega}',origin=origin,arm=arm,omega=omega,sites=n,counts_by_species=pop.sum(axis=1).tolist(),total_tokens=int(pop.sum()),quantization_error=float(abs(pop/omega-case['initial']).max()),channels=channels,reaction_rates=rates.sum(axis=1).tolist(),transport_rates=transport.tolist(),initial_rate=initial_rate,initial_activity_per_tick=initial_rate*.001,destination_payload_bytes=channels*7*n*8,count_channel_arrays_bytes=channels*2*7*n*8,transport_conservative=True)
          records.append(row);np.savez_compressed(BASE/(row['id']+'.npz'),counts=pop,concentrations=case['initial'],omega=omega,a0=a0,b0=b0,diff=case['diff'],coupling=k,links=json.dumps(links),job=json.dumps(job))
    files=['src/prepare_mg12h3_ring_feasibility.py','src/mg12h3_population.py','src/mg11g3.py','src/mg11g2_pattern_v2.py','src/mg12h3_exchangeable.py','src/mg12h1_reference.py','spec/MG12H3_ring_feasibility.md']
    sources.update({p:hashlib.sha256((ROOT/p).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for p in files})
    BASE.mkdir(parents=True,exist_ok=True);(BASE/'inventory.json').write_text(json.dumps(dict(complete=True,elapsed_seconds=time.perf_counter()-start,records=records,source_hashes=sources,random_draws=0,integrations=0),indent=2)+'\n')
if __name__=='__main__':BASE.mkdir(parents=True,exist_ok=True);run()
