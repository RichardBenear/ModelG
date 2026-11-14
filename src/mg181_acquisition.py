"""Fixed confirmation rules; observational acquisition metrics, not new labels."""
import math
from scipy.stats import beta

def episodes(records,consecutive):
 run=0;out=[]
 for t,r in enumerate(records):
  run=run+1 if r['emitted']=='single' else 0
  if run==consecutive:out.append(t)
 return out

def truth_state(case,t,onset):
 if case=='birth_moving':return t>=8
 if case in ['static','moving','stop_start','weak_static']:return True
 if case.startswith('damage_'):return t<onset
 return False

def metrics(records,truth,case,onset,consecutive):
 hits=episodes(records,consecutive);valid=[];false=[];overlap=[];late=[]
 for t in hits:
  states=[truth_state(case,j,onset) for j in range(t-consecutive+1,t+1)]
  if all(states):valid.append(t)
  elif any(states):overlap.append(t)
  else:
   false.append(t)
   if not case.startswith('damage_') or t-consecutive+1>=onset+5:late.append(t)
 first=hits[0] if hits else None;firstvalid=valid[0] if valid else None;acc=False
 if firstvalid is not None:
  r=records[firstvalid];g=truth[firstvalid];acc=min((r['center']-g[0])%500,(g[0]-r['center'])%500)<=2 and abs(r['width']-g[1])<=.12500001
 later=[t for t in range((firstvalid+1) if firstvalid is not None else len(records),len(records)) if truth_state(case,t,onset)];return dict(consecutive=consecutive,episode_times=hits,first_hit=first,first_valid_hit=firstvalid,valid_episodes=valid,false_episodes=false,overlap_episodes=overlap,response_allowed_false_episodes=late,acquired_by={str(t):first is not None and first<=t for t in [5,11,23,39]},first_valid_geometry_accurate=bool(acc),later_true_single_frames=len(later),later_single_labels=sum(records[t]['emitted']=='single' for t in later),birth_latency=(firstvalid-8 if firstvalid is not None and case=='birth_moving' else None))
def interval(k,n):return [0. if k==0 else float(beta.ppf(.025,k,n-k+1)),1. if k==n else float(beta.ppf(.975,k+1,n-k))]
