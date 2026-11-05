"""Static decision-rule counterexamples for MG9d Draft 2. No dynamics."""
import json
from pathlib import Path
import numpy as np

vf=5e-9
t=np.arange(0,10001,50,dtype=float)
# Equal, opposite resolved motions in the two scoring windows. Continuous positions.
d=np.where(t<=6000,4e-7*(t-2000),4e-7*4000-4e-7*(t-6000))
slope=lambda lo,hi:float(np.polyfit(t[(t>=lo)&(t<=hi)],d[(t>=lo)&(t<=hi)],1)[0])
a,b,v=slope(2000,6000),slope(6000,10000),slope(2000,10000)
assert abs(a)>=3*vf and abs(b)>=3*vf and a*b<0 and abs(v)<3*vf
out={'revision':'13997dca0a71869de4c5a8e8c82b9a03feb08e56',
 'reversal_cancellation':{'a':a,'b':b,'v':v,'C4':'REVERSAL',
    'overall_rule_3':'NO RESOLVED DRIFT despite resolved motion in both windows'},
 'invalid_reflection':{'R0':'VALID','R1':'VALID','R3':'INVALID','C1':'NUMERICS-AGREE',
    'issue':'C3 has no label for unmet prerequisites; overall rule 1 only checks ASYMMETRY-FLAG'},
 'unresolved_opposite_late_slope':{'a':-4e-7,'b':1e-9,'rho':-0.0025,
    'a_resolved':True,'b_resolved':False,'issue':'Reaches Otherwise in C4, contradicting rho >= 0'}}
Path(__file__).with_name('draft2_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
