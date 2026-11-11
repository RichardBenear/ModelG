"""MG154 primary/secondary fixed-blend validation figure."""
from pathlib import Path
import json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
BASE=Path(__file__).resolve().parents[1]/'results/MG12H3/Y_blend_validation'
def run():
 q=json.loads((BASE/'separate_recomputation.json').read_text());fig,axes=plt.subplots(1,2,figsize=(11,5.6),sharey=True)
 for ax,idx,title,color in [(axes[0],1,'Late Y, time 0.01: primary','#2377ae'),(axes[1],0,'Early Y, time 0.005: secondary','#23896b')]:
  row=q['statistics'][idx];ax.scatter(row['per_tick_blend_over_direct'],[0,1,2],s=55,color=color);v=row['blend_over_direct'];lo,hi=row['bootstrap_95_interval'];ax.errorbar(v,3,xerr=np.array([[v-lo],[hi-v]]),fmt='D',color=color,capsize=5,markersize=7);ax.axvline(1,color='#aa5155',ls='--');ax.set_xlim(0,1.06);ax.set_ylim(3.6,-.6);ax.set_title(title,loc='left',weight='bold');ax.set_xlabel('Blended variance / direct variance');ax.grid(axis='x',alpha=.2)
 axes[0].set_yticks([0,1,2,3],['Tick 0.0001','Tick 0.00005','Tick 0.00002','Pooled + 95% interval']);fig.text(.08,.035,'Ratios below 1 favor the blend. Only the pooled ratio has an interval shown.\nEight fresh paths per tick; paired bootstrap intervals have approximate small-sample coverage.',fontsize=10,color='#4c5c6b');fig.tight_layout(rect=(0,.10,1,1));fig.savefig(BASE/'Y_blend_variance.png',dpi=170,bbox_inches='tight');fig.savefig(BASE/'Y_blend_variance.svg',bbox_inches='tight')
if __name__=='__main__':run()
