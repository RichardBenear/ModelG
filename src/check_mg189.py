"""Separate saved-trace audit of MG189 causality, summaries and intervals."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import gzip,json,math,hashlib,time
from collections import Counter
import numpy as np
from scipy.stats import ttest_ind
from causal_mg189 import ROOT,PUBLIC,RAW,METHODS,METRICS,BINS,read,atomic_json,binding

def run():
 start=time.perf_counter();plan=read(PUBLIC/'plan.json');assert binding()==plan['source'];assert read(PUBLIC/'attempt.json')['status']=='COMPLETE';checked=0;events_checked=0;arrays=[]
 for mi,method in enumerate(METHODS):
  data=np.load(RAW/f'summaries_{mi}.npz');values=data['values'];assert values.shape==(8192,2,16);np.testing.assert_array_equal(data['seed'],plan['seeds'][method]);arrays.append(values)
  with gzip.open(RAW/f'traces_{mi}.jsonl.gz','rt',encoding='utf-8') as trace:
   for pi,line in enumerate(trace):
    row=json.loads(line);assert row['path']==pi and row['seed']==plan['seeds'][method][pi]
    init=[('G',0),('X',0),('X',0),('Y',0),('G',1),('X',1),('X',1),('Y',1),('G',2),('X',2),('X',3),('Y',3)];tokens={f'i{k}':(s,i,-1,-1) for k,(s,i) in enumerate(init)};ever=set(tokens);edge_lists=[[],[]];labels=[]
    for ei,event in enumerate(row['events']):
     rule=event['rule'];inputs=event['inputs'];labels.append(rule);parents=[Counter(),Counter()]
     expected_inputs={'R1':[],'R2':['G'],'R2r':['X'],'R3':['X'],'R4':['X','X','Y'],'R5':['X']}.get(rule,[rule[2:]])
     expected_outputs={'R1':['G'],'R2':['X'],'R2r':['G'],'R3':['Y','Z'],'R4':['X','X','X'],'R5':['W']}.get(rule,[rule[2:]])
     assert [v[1] for v in inputs]==expected_inputs and [v[1] for v in event['outputs']]==expected_outputs
     assert len({v[0] for v in inputs})==len(inputs) and all(v[2]==event['site'] for v in inputs) and all(v[2]==event['dest'] for v in event['outputs'])
     assert (event['dest'] in ((event['site']-1)%4,(event['site']+1)%4)) if rule.startswith('H_') else event['site']==event['dest']
     for rank,(uid,species,site) in enumerate(inputs):
      assert uid in tokens and tokens[uid][:2]==(species,site)
      for c in (0,1):
       producer=tokens[uid][2+c]
       if producer>=0 and not(c==1 and rule=='R4' and rank<2):parents[c][producer,species]+=1
     for c in (0,1):
      actual=sorted((p,s,n) for (p,s),n in parents[c].items());assert actual==[tuple(v) for v in event['edges'][c]];edge_lists[c].extend((p,ei) for p,s,n in actual)
     produced=[]
     for rank,(uid,species,site) in enumerate(event['outputs']):
      assert uid not in ever;ever.add(uid);shadow=tokens[inputs[rank][0]][3] if rule=='R4' and rank<2 else ei;produced.append((uid,(species,site,ei,shadow)))
     for uid,_,_ in inputs:del tokens[uid]
     tokens.update(produced);events_checked+=1
    for c in (0,1):
     n=len(labels);adj=np.zeros((n,n),bool)
     for p,q in edge_lists[c]:assert p<q;adj[p,q]=True
     levels=np.ones(n,int)
     for i in range(n):
      pp=np.flatnonzero(adj[:,i]);levels[i]=1+(int(levels[pp].max()) if len(pp) else 0)
     counts=Counter(labels);indegree=adj.sum(axis=0);outdegree=adj.sum(axis=1);metrics=[counts[r] for r in METRICS[:9]]+[n,int(levels.max()) if n else 0,max(Counter(levels).values(),default=0),float(indegree.mean()) if n else 0,int(indegree.max()) if n else 0,float(outdegree.mean()) if n else 0,int(outdegree.max()) if n else 0];np.testing.assert_allclose(metrics,values[pi,c],atol=1e-12,rtol=0)
    checked+=1
   assert pi+1==8192
  print('audited',method,flush=True)
 comparison=read(PUBLIC/'statistical_results.json')['comparisons'];number=128;eps=math.sqrt((13*math.log(2)+math.log(2*number/.05))/(2*8192))
 for row in comparison:
  mi=METHODS.index(row['method']);c=row['convention'];j=METRICS.index(row['metric']);a=arrays[0][:,c,j];b=arrays[mi][:,c,j]
  if a.var()+b.var():
   interval=ttest_ind(b,a,equal_var=False).confidence_interval(confidence_level=1-.05/number);np.testing.assert_allclose(row['CI95_family'],[interval.low,interval.high],atol=1e-10)
  bins_a=Counter(np.searchsorted(BINS,a,side='right')-1);bins_b=Counter(np.searchsorted(BINS,b,side='right')-1);tv=sum(abs(bins_a[k]-bins_b[k]) for k in range(13))/(2*8192);np.testing.assert_allclose(row['binned_TV'],tv);np.testing.assert_allclose(row['TV_bounds95_family'],[max(0,tv-2*eps),min(1,tv+2*eps)])
 atomic_json(PUBLIC/'verification.json',dict(passed=True,paths=checked,events=events_checked,stoichiometry_sites_and_unique_token_use=True,trace_edges_rederived=True,causal_metrics_recomputed_with_adjacency_matrices=True,mean_intervals_checked_with_Scipy_ttest=True,TV_recomputed_with_digitized_counter=True,all_comparisons=len(comparison),elapsed_seconds=time.perf_counter()-start,external_replication=False,shared_dependencies=['frozen fixture and summary definitions','NumPy/SciPy'],limits='Initial hazard/generator check and finite isomorphism witnesses are in original evaluator; independent saved-trace audit does not constitute external replication or continuum proof.'))
if __name__=='__main__':run()
