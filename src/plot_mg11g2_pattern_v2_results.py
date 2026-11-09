"""MG090 descriptive figure from audited results; no new scoring or simulations."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11G2/pattern_response_v2'

def main():
    e=json.loads((BASE/'evaluation.json').read_text(encoding='utf-8'));a=json.loads((BASE/'independent_check.json').read_text(encoding='utf-8'));assert a['passed'] and e['complete'];t=np.arange(101)/4
    colors=['#6b7280','#7c3aed','#007f8b','#d95f02'];labels=['Maintained','Isolated','Local links','Cross-ring links']
    def excess(c):return np.array(e['records'][f'c{c}_s1_v1']['metrics']['projection'])[:,1]-np.array(e['records'][f'c{c}_s0_v1']['metrics']['projection'])[:,1]
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained');fig.suptitle('MG11-G2: environmental connectivity changes pattern response',fontsize=16,fontweight='bold')
    for ax,conditions,title in ((axes[0,0],[0,1,2,3],'Lower resource capacity (C=1)'),(axes[0,1],[0,4,5,6],'Higher resource capacity (C=100)')):
        for c,color,label in zip(conditions,colors,labels):ax.plot(t,excess(c),color=color,label=label,lw=2)
        ax.axvspan(20,25,color='#d1d5db',alpha=.3);ax.set(title=title,xlabel='Model time',ylabel='X seed-excess projection');ax.legend(frameon=False,fontsize=9);ax.grid(alpha=.15)
    ax=axes[1,0];contrast=excess(3)-excess(2);base=[]
    for c in (3,2):base.append(np.array(e['records'][f'c{c}_s1_v0']['metrics']['projection'])[:,1]-np.array(e['records'][f'c{c}_s0_v0']['metrics']['projection'])[:,1])
    ax.plot(t,base[0]-base[1],label='Base numerical setting',color='#8b5cf6',lw=3,alpha=.65);ax.plot(t,contrast,label='Fine numerical setting',color='#111827',lw=1.4,ls='--');ax.axvspan(20,25,color='#d1d5db',alpha=.3);ax.axhline(0,color='#999999',lw=.8);ax.set(title='Primary contrast: cross-ring minus local, C=1',xlabel='Model time',ylabel='Difference in X projection');ax.text(.03,.94,f'Late-window mean: {e["primary_GXY_by_setting"][1][1]:+.6f}\nNumerically resolved; no statistical interval',transform=ax.transAxes,va='top',fontsize=9);ax.legend(frameon=False,fontsize=9,loc='lower right');ax.grid(alpha=.15)
    ax=axes[1,1];x=(np.arange(500)-250)*.2
    for c,color,label in zip((0,1,2,3),colors,labels):
        arrays=[]
        for s in (0,1):
            with np.load(BASE/f'production/c{c}_s{s}_v1.npz') as raw:arrays.append(raw['state'][-1,500:1000].copy())
        field=arrays[1]-arrays[0];field-=field.mean();ax.plot(x,field,color=color,label=label,lw=2)
    ax.set(title='Final centered X seed-excess profiles, C=1',xlabel='Fixed active-ring coordinate',ylabel='X concentration above matched H');ax.set_xlim(-10,10);ax.grid(alpha=.15)
    fig.supxlabel('Larger projection is not longer persistence: reference matching lasted local19, cross-ring16.75 time units at C=1.',fontsize=10)
    fig.savefig(BASE/'MG11G2_pattern_summary.png',dpi=170);plt.close(fig)
if __name__=='__main__':main()
