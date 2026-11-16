"""MG189 explicit-token small-system causal study; existing engines unchanged."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import sys,json,time,hashlib,itertools,math,gzip
from pathlib import Path
from collections import Counter,defaultdict
from dataclasses import dataclass
import numpy as np
from scipy.stats import t
from runner import seed_for
from run_mg176 import ROOT
from mg10e2c1_pool import atomic_json
PUBLIC=ROOT/'results/MG189_causal_invariance';RAW=ROOT/'local_raw/MG189_causal_invariance'
RULES=('R1','R2','R2r','R3','R4','R5','H_G','H_X','H_Y')
METHODS=('SSA','tick.02','tick.005','tick.001','sweep.02')
BINS=np.array([0,.1,.25,.5,1,2,3,4,6,10,20,40,80,np.inf])
METRICS=RULES+('events','depth','layer_width','mean_in','max_in','mean_out','max_out')
FILES=('src/causal_mg189.py','spec/MG189_causal_invariance.md','spec/causal_invariance_study_plan.md','src/mg184_open_bridge.py','src/mg12h3_sparse.py','src/lattice.py')

def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def binding():return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FILES}

@dataclass
class State:
 tokens:dict # id -> species,site,primary producer,sensitivity producer
 sources:tuple
 events:list

def fixture(name):
 definitions={'F1':[(0,'G'),(1,'G')],'F2':[(0,'G'),(0,'X'),(0,'X'),(0,'Y')],'F3':[(0,'G'),(0,'X'),(0,'X'),(1,'Y')],'F4':[(0,'G'),(0,'X'),(0,'X'),(0,'Y'),(1,'G'),(1,'X'),(1,'X'),(1,'Y'),(2,'G'),(2,'X'),(3,'X'),(3,'Y')]}
 n={'F1':2,'F2':1,'F3':2,'F4':4}[name];return State({f'i{k}':(s,i,-1,-1) for k,(i,s) in enumerate(definitions[name])},(0,)*n,[])

def matches(state):
 out=[];n=len(state.sources)
 for site in range(n):
  by={s:sorted(k for k,v in state.tokens.items() if v[0]==s and v[1]==site) for s in ('G','X','Y')};out.append(('R1',site,site,(),14.))
  for g in by['G']:out.append(('R2',site,site,(g,),1.))
  for x in by['X']:
   out.extend([('R2r',site,site,(x,),.1),('R3',site,site,(x,),29.),('R5',site,site,(x,),.9)])
  for x in itertools.combinations(by['X'],2):
   for y in by['Y']:out.append(('R4',site,site,x+(y,),2.))
  if n>1:
   for species,diff in (('G',1.),('X',1.),('Y',12.)):
    for uid in by[species]:
     for dest in sorted({(site-1)%n,(site+1)%n}):out.append(('H_'+species,site,dest,(uid,),diff))
 return out

def apply(state,match):
 rule,site,dest,inputs,_=match;assert all(uid in state.tokens for uid in inputs);tokens=state.tokens.copy();index=len(state.events);edges=[]
 for convention in (0,1):
  parents=Counter()
  for j,uid in enumerate(inputs):
   species,_,p,q=tokens[uid];producer=p if convention==0 else q
   if convention==1 and rule=='R4' and j<2:continue
   if producer>=0:parents[(producer,species)]+=1
  edges.append(tuple(sorted((p,s,c) for (p,s),c in parents.items())))
 output={'R1':('G',),'R2':('X',),'R2r':('G',),'R3':('Y','Z'),'R4':('X','X','X'),'R5':('W',)}.get(rule,(rule[2:],))
 source=list(state.sources);signature=(rule,site,dest,inputs,source[site] if rule=='R1' else 0);uidbase=hashlib.sha256(repr(signature).encode()).hexdigest();outputs=[]
 for rank,species in enumerate(output):
  uid=uidbase+':'+str(rank);shadow=state.tokens[inputs[rank]][3] if rule=='R4' and rank<2 else index;tokens[uid]=(species,dest,index,shadow);outputs.append((uid,species,dest))
 for uid in inputs:del tokens[uid]
 if rule=='R1':source[site]+=1
 event=dict(rule=rule,site=site,dest=dest,inputs=[(uid,state.tokens[uid][0],state.tokens[uid][1]) for uid in inputs],outputs=outputs,edges=edges)
 return State(tokens,tuple(source),state.events+[event])

def statekey(state,identity):
 if identity:return tuple(sorted((uid,v[0],v[1]) for uid,v in state.tokens.items()))
 return tuple(sorted(Counter((v[0],v[1]) for v in state.tokens.values()).items()))

def graphdata(state,convention):
 return tuple(e['rule'] for e in state.events),tuple((p,j,s,c) for j,e in enumerate(state.events) for p,s,c in e['edges'][convention])

def canonical(state,convention):
 labels,edges=graphdata(state,convention);n=len(labels);best=None
 for order in itertools.permutations(range(n)):
  inverse={old:new for new,old in enumerate(order)};candidate=(tuple(labels[i] for i in order),tuple(sorted((inverse[p],inverse[q],s,c) for p,q,s,c in edges)))
  if best is None or candidate<best:best=candidate
 return best

def statistics(state,convention):
 labels,edges=graphdata(state,convention);parents=[set() for _ in labels];children=[set() for _ in labels]
 for p,q,_,_ in edges:parents[q].add(p);children[p].add(q)
 depth=[]
 for p in parents:depth.append(1+max((depth[j] for j in p),default=0))
 count=Counter(labels);n=len(labels);return [count[r] for r in RULES]+[n,max(depth,default=0),max(Counter(depth).values(),default=0),np.mean([len(p) for p in parents]) if n else 0,max(map(len,parents),default=0),np.mean([len(p) for p in children]) if n else 0,max(map(len,children),default=0)]

def extensions(state,depth,deadline):
 allstates=[state];front=[state]
 for _ in range(depth):
  following=[]
  for current in front:
   if time.perf_counter()>deadline:raise TimeoutError('enumeration deadline')
   following.extend(apply(current,m) for m in matches(current))
  allstates.extend(following);front=following
  if len(allstates)>200000:raise RuntimeError('history cap')
 return allstates

def exact_tests(deadline):
 import networkx as nx
 results={};witnesses=[];conflicts=[];commute=0
 for name in ('F1','F4'):
  state=fixture(name);mm=matches(state)
  for a,b in itertools.combinations(mm,2):
   if set(a[3])&set(b[3]):continue
   first=apply(apply(state,a),b);second=apply(apply(state,b),a);assert statekey(first,True)==statekey(second,True)
   for c in (0,1):assert canonical(first,c)==canonical(second,c)
   commute+=1
 results['T2']=dict(passed=True,independent_pairs=commute,scope='all disjoint initial two-event match pairs F1/F4')
 for name in ('F2','F3'):
  initial=fixture(name);mm=matches(initial)
  cache={}
  for i,m in enumerate(mm):
   states=extensions(apply(initial,m),2,deadline);cache[i]=[set(statekey(v,c) for v in states) for c in (False,True)]
  for (i,a),(j,b) in itertools.combinations(enumerate(mm),2):
   shared=set(a[3])&set(b[3])
   if shared:conflicts.append(dict(fixture=name,rules=[a[0],b[0]],inputs=[a[3],b[3]],shared=sorted(shared),count_join_within2=bool(cache[i][0]&cache[j][0]),identity_join_within2=bool(cache[i][1]&cache[j][1])))
 results['T1']=dict(pairs=len(conflicts),count_not_joined_within2=sum(not r['count_join_within2'] for r in conflicts),identity_not_joined_within2=sum(not r['identity_join_within2'] for r in conflicts),no_global_unjoinability_claim=True)
 atomic_json(PUBLIC/'conflict_catalogue.json',conflicts)
 allresults=[]
 with gzip.open(RAW/'enumeration.jsonl.gz','wt',encoding='utf-8') as raw:
  for name,depth in (('F2',3),('F3',2),('F4',2)):
   states=extensions(fixture(name),depth,deadline);assert len({statekey(s,False) for s in states})+len({statekey(s,True) for s in states})<=1000000
   for state in states:raw.write(json.dumps(dict(fixture=name,state=list(statekey(state,True)),events=state.events))+'\n')
   for identity in (False,True):
    for convention in (0,1):
     for equal_depth in (False,True):
      groups=defaultdict(Counter);examples={}
      for state in states:
       key=(statekey(state,identity),len(state.events) if equal_depth else None);cg=canonical(state,convention);groups[key][cg]+=1;examples.setdefault((key,cg),state)
      reconverged=matched=0;witness=None
      for key,counts in groups.items():
       n=sum(counts.values());reconverged+=n*(n-1)//2;matched+=sum(v*(v-1)//2 for v in counts.values())
       if len(counts)>1 and witness is None:
        first,second=list(counts)[:2];one=examples[key,first];two=examples[key,second];witness=dict(fixture=name,identity=identity,convention=convention,equal_depth=equal_depth,state=repr(key[0]),branch1=one.events,branch2=two.events)
        # Independent isomorphism checker preserves species multiplicity labels.
        graphs=[]
        for ss in (one,two):
         labels,edges=graphdata(ss,convention);graph=nx.DiGraph();graph.add_nodes_from((i,dict(rule=r)) for i,r in enumerate(labels));ed=defaultdict(list)
         for p,q,s,c in edges:ed[p,q].append((s,c))
         for (p,q),v in ed.items():graph.add_edge(p,q,species=tuple(sorted(v)))
         graphs.append(graph)
        assert not nx.is_isomorphic(*graphs,node_match=lambda a,b:a['rule']==b['rule'],edge_match=lambda a,b:a['species']==b['species']);witnesses.append(witness)
      total=(sum(Counter(len(s.events) for s in states).values())*(len(states)-1)//2) if not equal_depth else sum(v*(v-1)//2 for v in Counter(len(s.events) for s in states).values())
      allresults.append(dict(fixture=name,depth=depth,histories=len(states),state_mode='identity' if identity else 'counts',convention=convention,equal_event_depth=equal_depth,unique_states=len(groups),all_branch_pairs=total,reconverging_pairs=reconverged,matching_graph_pairs=matched,reconvergence_fraction=reconverged/total if total else None,matching_fraction=matched/reconverged if reconverged else None,holds_in_tested_pairs=matched==reconverged))
 results['T3']=allresults;atomic_json(PUBLIC/'counterexamples.json',witnesses);return results

def proposals(state,delta,rng):
 proposed=[]
 for m in matches(state):
  if m[0]=='R1':
   for clock in rng.uniform(0,delta,size=rng.poisson(m[4]*delta)):proposed.append((float(clock),m))
  else:
   clock=float(rng.exponential(1/m[4]))
   if clock<delta:proposed.append((clock,m))
 return sorted(proposed,key=lambda v:v[0])

def path(method,seed,deadline):
 rng=np.random.default_rng(seed);state=fixture('F4');clock=0.;ticks=conflict_ticks=rejected=0
 if method=='SSA':
  while True:
   if time.perf_counter()>deadline:raise TimeoutError('sampling deadline')
   mm=matches(state);rates=np.array([m[4] for m in mm]);clock+=rng.exponential(1/rates.sum())
   if clock>=.1:break
   m=mm[min(int(np.searchsorted(np.cumsum(rates),rng.random()*rates.sum(),side='right')),len(mm)-1)];state=apply(state,m)
 else:
  delta=float(method[4:] if method.startswith('tick') else method[5:]);steps=round(.1/delta)
  for _ in range(steps):
   if time.perf_counter()>deadline:raise TimeoutError('sampling deadline')
   if method.startswith('tick'):
    old=state;used=set();accepted=[];bad=0
    for _,m in proposals(old,delta,rng):
     if used.intersection(m[3]):bad+=1;continue
     used.update(m[3]);accepted.append(m);state=apply(state,m)
    # Explicit event index permutation, avoiding raw event-label equality.
    reverse=old
    for m in reversed(accepted):reverse=apply(reverse,m)
    assert statekey(state,True)==statekey(reverse,True)
    base=len(old.events);mapping={i:i for i in range(base)};reverse_outputs={tuple(v[0] for v in e['outputs']):i for i,e in enumerate(reverse.events[base:],start=base)}
    mapping.update({i:reverse_outputs[tuple(v[0] for v in e['outputs'])] for i,e in enumerate(state.events[base:],start=base)})
    for uid,v in state.tokens.items():
     assert v[:2]==reverse.tokens[uid][:2]
     assert (-1 if v[2]<0 else mapping[v[2]])==reverse.tokens[uid][2]
    for c in (0,1):
     labels,edges=graphdata(state,c);rl,re=graphdata(reverse,c);assert all(labels[i]==rl[mapping[i]] for i in mapping);assert sorted((mapping[p],mapping[q],s,k) for p,q,s,k in edges)==sorted(re)
    ticks+=1;conflict_ticks+=bool(bad);rejected+=bad
   else:
    n=len(state.sources);channels=[(rule,site,site) for site in range(n) for rule in RULES[:6]]+[(rule,site,dest) for site in range(n) for dest in sorted({(site-1)%n,(site+1)%n}) for rule in RULES[6:]]
    for ci in rng.permutation(len(channels)):
     key=channels[ci];local=0.
     while True:
      mm=[m for m in matches(state) if m[:3]==key]
      if not mm:break
      rates=np.array([m[4] for m in mm]);local+=rng.exponential(1/rates.sum())
      if local>=delta:break
      m=mm[min(int(np.searchsorted(np.cumsum(rates),rng.random()*rates.sum(),side='right')),len(mm)-1)];state=apply(state,m)
   if len(state.events)>10000:raise RuntimeError('path event cap')
 return state,(ticks,conflict_ticks,rejected)

def register():
 assert not(PUBLIC/'plan.json').exists();PUBLIC.mkdir(exist_ok=True);used=set()
 def collect(value):
  if type(value)is int:used.add(value)
  elif isinstance(value,(dict,list)):
   for v in value.values() if isinstance(value,dict) else value:collect(v)
 def scan(value):
  if isinstance(value,dict):
   for k,v in value.items():collect(v) if 'seed' in k.lower() else scan(v)
  elif isinstance(value,list):
   for v in value:scan(v)
 for index,file in enumerate((ROOT/'results').rglob('*.json')):
  scan(read(file))
  if index%2000==0:print('seed audit files',index,flush=True)
 seeds={}
 for mi,method in enumerate(METHODS):
  values=[]
  for r in range(8192):
   nonce=0;s=seed_for(18900,mi,0,0,r)
   while s in used:nonce+=1;s=seed_for(18900,mi,0,nonce,r)
   used.add(s);values.append(s)
  seeds[method]=values
 atomic_json(PUBLIC/'plan.json',dict(source=binding(),seeds=seeds,paths_per_method=8192,workers=1,seconds_cap=1800,raw_bytes_cap=1000000000,approval_basis='PI requested reviewer study; D1-D4 adopted in reviewer outline',primary_predictions={'T2':'PASS','T3':'count-level exact invariance fails','T5':'agreement improves as tick shrinks; no forced precision pass'}))

def comparisons(arrays):
 n=arrays.shape[1];number=4*2*len(METRICS);out=[];eps=math.sqrt(((len(BINS)-1)*math.log(2)+math.log(2*number/.05))/(2*n))
 for m in range(1,5):
  for c in (0,1):
   for j,metric in enumerate(METRICS):
    a=arrays[0,:,c,j];b=arrays[m,:,c,j];va=a.var(ddof=1)/n;vb=b.var(ddof=1)/n;df=(va+vb)**2/(va*va/(n-1)+vb*vb/(n-1)) if va+vb else 2*n-2;half=float(t.ppf(1-.05/(2*number),df)*math.sqrt(va+vb));difference=float(b.mean()-a.mean());margin=.05*max(float(a.mean()),1);ci=[difference-half,difference+half];mean_status='EQUIVALENT' if ci[0]>=-margin and ci[1]<=margin else 'DIFFERENT' if ci[0]>margin or ci[1]<-margin else 'INCONCLUSIVE';ha=np.histogram(a,BINS)[0]/n;hb=np.histogram(b,BINS)[0]/n;tv=float(abs(ha-hb).sum()/2);lower=max(0,tv-2*eps);upper=min(1,tv+2*eps);status='EQUIVALENT' if upper<=.10 else 'DIFFERENT' if lower>.10 else 'INCONCLUSIVE'
    out.append(dict(method=METHODS[m],convention=c,metric=metric,SSA_mean=float(a.mean()),other_mean=float(b.mean()),difference=difference,CI95_family=ci,margin=margin,mean_status=mean_status,binned_TV=tv,TV_bounds95_family=[lower,upper],distribution_status=status))
 return out

def checks():
 from mg184_open_bridge import OpenBridge
 state=fixture('F4');population=np.zeros((7,4),np.int64);species={'G':0,'X':1,'Y':2,'Z':5,'W':6}
 for s,i,_,_ in state.tokens.values():population[species[s],i]+=1
 drift=np.zeros((3,4));hazard=0.
 for m in matches(state):
  after=apply(state,m);change=np.zeros((3,4))
  for s,i,_,_ in after.tokens.values():
   if s in ('G','X','Y'):change[species[s],i]+=1
  drift+=m[4]*(change-population[:3]);hazard+=m[4]
 model=OpenBridge(4,1,h=1);expected=np.zeros((7,4));total=0.
 for channel in model.channels:
  rate=channel.rate(population.ravel());total+=rate
  for pos,v in channel.change:expected.ravel()[pos]+=rate*v
 np.testing.assert_allclose(drift,expected[:3],atol=1e-12);np.testing.assert_allclose(hazard,total,atol=1e-12)
 return dict(open_bridge_hazard=float(total),maximum_drift_error=float(abs(drift-expected[:3]).max()),falling_factorial_pair_factor_checked=True)

def run():
 plan=read(PUBLIC/'plan.json');assert plan['source']==binding() and not RAW.exists();RAW.mkdir();start=time.perf_counter();deadline=start+1800;atomic_json(PUBLIC/'attempt.json',dict(status='RUNNING'))
 try:
  initialchecks=checks();exact=exact_tests(deadline);atomic_json(PUBLIC/'exact_results.json',exact);arrays=np.zeros((5,8192,2,len(METRICS)));tick_counts=np.zeros((5,3),np.int64)
  for mi,method in enumerate(METHODS):
   with gzip.open(RAW/f'traces_{mi}.jsonl.gz','wt',encoding='utf-8') as trace:
    for i,seed in enumerate(plan['seeds'][method]):
     state,ticks=path(method,seed,deadline);arrays[mi,i]=[statistics(state,c) for c in (0,1)];tick_counts[mi]+=ticks;trace.write(json.dumps(dict(path=i,seed=seed,events=state.events))+'\n')
     if i%256==0:
      atomic_json(PUBLIC/'progress.json',dict(method=method,completed=i+1,elapsed_seconds=time.perf_counter()-start));assert sum(p.stat().st_size for p in RAW.rglob('*') if p.is_file())<1000000000
   np.savez_compressed(RAW/f'summaries_{mi}.npz',values=arrays[mi],seed=np.array(plan['seeds'][method],dtype=np.uint32));print(method,'complete',time.perf_counter()-start,flush=True)
  comparisons_out=comparisons(arrays);verdict=[]
  for method in METHODS[1:]:
   rr=[r for r in comparisons_out if r['method']==method];statuses=[v for r in rr for v in (r['mean_status'],r['distribution_status'])];verdict.append(dict(method=method,classification='EQUIVALENT_IN_SCOPE' if all(s=='EQUIVALENT' for s in statuses) else 'DIFFERENT_IN_SOME_SUMMARIES' if 'DIFFERENT' in statuses else 'INCONCLUSIVE',mean_counts=dict(Counter(r['mean_status'] for r in rr)),distribution_counts=dict(Counter(r['distribution_status'] for r in rr))))
  atomic_json(PUBLIC/'statistical_results.json',dict(comparisons=comparisons_out,verdict=verdict,metrics=METRICS,paths_per_method=8192,tick_counts=[dict(method=m,ticks=int(v[0]),ticks_with_rejected_competition=int(v[1]),rejected_proposals=int(v[2]),conflict_fraction=None if not v[0] else float(v[1]/v[0])) for m,v in zip(METHODS,tick_counts)],initial_checks=initialchecks))
  atomic_json(PUBLIC/'attempt.json',dict(status='COMPLETE',elapsed_seconds=time.perf_counter()-start,raw_bytes=sum(p.stat().st_size for p in RAW.rglob('*') if p.is_file())))
 except BaseException:
  import traceback
  atomic_json(PUBLIC/'attempt.json',dict(status='FAILED',elapsed_seconds=time.perf_counter()-start,error=traceback.format_exc()));raise
if __name__=='__main__':register() if sys.argv[1]=='register' else run()
