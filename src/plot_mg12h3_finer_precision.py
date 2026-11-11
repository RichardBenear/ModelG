"""MG149 variance-ratio figure from pilot decision records."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np,json
BASE=Path(__file__).resolve().parents[1]/'results/MG12H3/finer_tick_precision_pilot'
def run():
 q=json.loads((BASE/'separate_recomputation.json').read_text());labels=[f'{s} / t={time}' for time in (.005,.01) for s in ('G','X','Y')];fig,ax=plt.subplots(figsize=(10,5.7));colors=['#2377ae','#219777','#bc7437']
 for level,delta in enumerate((.0001,.00005,.00002)):ax.scatter(np.asarray(q['precision'][level]['variance_ratio']).ravel(),np.arange(6)+(level-1)*.18,s=55,label=f'tick {delta:g}',color=colors[level],zorder=3)
 ax.axvline(1,color='#a34e59',ls='--',lw=1.4);ax.set_xscale('log');ax.set_yticks(np.arange(6),labels);ax.invert_yaxis();ax.set_xlabel('Direct variance ÷ integrated-estimator variance (log scale)');ax.set_title('G contrast gains precision; late Y does not',loc='left',weight='bold');ax.grid(axis='x',alpha=.2);ax.legend(loc='lower right');ax.text(.01,-.20,'Above 1: integrated estimator has lower sample variance. Below 1: direct estimate is less noisy.\nEight paths per condition; ratios are uncertain pilot estimates, not guaranteed gains.',transform=ax.transAxes,fontsize=10,color='#4e5b69');fig.tight_layout(rect=(0,.055,1,1));fig.savefig(BASE/'variance_ratio.png',dpi=170,bbox_inches='tight');fig.savefig(BASE/'variance_ratio.svg',bbox_inches='tight')
if __name__=='__main__':run()
