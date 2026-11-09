"""MG094 descriptive plotting from audited records; no new estimands/simulations."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11G3'
def main():
    e=json.loads((BASE/'evaluation.json').read_text(encoding='utf-8'));a=json.loads((BASE/'independent_check.json').read_text(encoding='utf-8'));assert a['passed'] and e['complete'];t=np.arange(101)/4
    def excess(o,arm):return np.array(e['records'][f'o{o}_s1_v1_a{arm}']['metrics']['projection'])[:,1]-np.array(e['records'][f'o{o}_s0_v1_a{arm}']['metrics']['projection'])[:,1]
    colors=['#007f8b','#d95f02','#7c3aed'];labels=['B uses local links','B uses cross-ring links','B transport off'];plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained');fig.suptitle('MG11-G3: ongoing B transport changes the conditional response',fontsize=15,fontweight='bold')
    for o,ax in enumerate(axes[0]):
        for arm,color,label in zip((0,1,2),colors,labels):ax.plot(t,excess(o,arm),color=color,label=label,lw=2)
        ax.set(title='Prepared with '+('local links' if o==0 else 'cross-ring links'),ylabel='X seed-excess projection');ax.legend(frameon=False,fontsize=9)
    ax=axes[1,0];d=[excess(o,1)-excess(o,0) for o in (0,1)]
    for o,color in zip((0,1),colors):ax.plot(t,d[o],color=color,label=('Local' if o==0 else 'Cross-ring')+' prepared')
    ax.plot(t,(d[0]+d[1])/2,color='#111827',lw=2,ls='--',label='Equal-origin average');ax.axhline(0,color='#999999',lw=.7);ax.set(title='B cross-ring minus local transport',ylabel='Difference in X projection');ax.text(.03,.96,f'Primary late average: {e["primary_GXY_by_setting"][1][1]:+.6f}\nRegistered threshold: 0.01',transform=ax.transAxes,va='top',fontsize=9);ax.legend(frameon=False,fontsize=9,loc='lower right')
    ax=axes[1,1]
    for o,color in zip((0,1),colors):ax.plot(t,excess(o,2)-excess(o,o),color=color,lw=2,label=('Local' if o==0 else 'Cross-ring')+' prepared')
    ax.axhline(0,color='#999999',lw=.7);ax.set(title='B-off minus unchanged transport (secondary)',ylabel='Difference in X projection');ax.legend(frameon=False,fontsize=9)
    for ax in axes.flat:ax.axvline(10,color='#6b7280',ls=':',lw=1.3);ax.axvspan(20,25,color='#d1d5db',alpha=.25);ax.set_xlabel('Absolute model time (switch at 10)');ax.grid(alpha=.15)
    fig.supxlabel('Conditional deterministic histories, not independent replicates. Higher projection is not better particle persistence.',fontsize=10);fig.savefig(BASE/'MG11G3_summary.png',dpi=170);plt.close(fig)
if __name__=='__main__':main()
