"""MG149 illustrative precision projection; no equivalence verdict."""
from pathlib import Path
import numpy as np,json,math
from scipy.stats import t
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG12H3/finer_tick_precision_pilot'
def run():
 q=json.loads((BASE/'separate_recomputation.json').read_text());choice=json.loads((ROOT/'results/MG12H3/finer_tick_design/pattern_accuracy_choice.json').read_text());margins=np.array([x['margin_5percent'] for x in choice['records']]);plans=[]
 for level in (0,1):
  for time in range(2):
   for species in range(3):
    key='direct_variance' if species==2 and time==1 else 'integral_variance';v1=q['precision'][level][key][time][species];v2=q['precision'][2][key][time][species];n=48;df=(n-1)*(v1+v2)**2/(v1*v1+v2*v2);critical=float(t.ppf(1-.05/48,df));half=critical*math.sqrt((v1+v2)/n);plans.append(dict(level=level,time=(.005,.01)[time],species=('G','X','Y')[species],estimator='direct' if key=='direct_variance' else 'integral',assumed_replicates=48,halfwidth_from_pilot_variances=half,margin=float(margins[species]),halfwidth_to_margin=half/margins[species],scope='Planning from eight pilot paths, not a qualification interval or sample-count guarantee'))
 (BASE/'precision_planning.json').write_text(json.dumps(dict(planning_only=True,records=plans,selection='Illustrative fixed estimator by species/time: integral for G/X and early Y; direct for late Y. Not yet registered for qualification.',uncertainty='Pilot variances have only seven degrees of freedom percondition; point calculations can understate required sample counts.'),indent=2)+'\n')
if __name__=='__main__':run()
