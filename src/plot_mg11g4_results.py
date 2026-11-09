"""MG098 descriptive qualification figure from audited frozen statistics."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'results/MG11G4'
def main():
    e=json.loads((BASE/'evaluation.json').read_text(encoding='utf-8'));a=json.loads((BASE/'independent_check.json').read_text(encoding='utf-8'));assert e['complete'] and a['passed']
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False});fig,axes=plt.subplots(1,2,figsize=(12,5.2),layout='constrained');fig.suptitle('MG11-G4: small-system stochastic coupling qualification',fontsize=16,fontweight='bold');x=np.arange(3);labels=['1 node\nChemistry','3 nodes\nScarce resources','3 nodes\nFast active transport']
    ax=axes[0];upper=[]
    for engine,ref,label,color,shift in ((3,0,'Finest leap vs exact','#007f8b',-.18),(0,4,'Exact vs independent exact','#8b5cf6',.18)):
        values=[max(row['upper'] for row in e['formal'] if row['cell']==c and row['engine']==engine and row['reference_engine']==ref) for c in range(3)];upper+=values;bars=ax.bar(x+shift,values,width=.34,label=label,color=color)
        ax.bar_label(bars,labels=[f'{v:.4f}' for v in values],padding=3,fontsize=9)
    ax.axhline(.08,color='#b91c1c',ls='--',label='Registered margin 0.08');ax.set_xticks(x,labels);ax.set_ylim(0,max(.10,max(upper)*1.3));ax.set_ylabel('Largest marginal CDF-distance upper bound');ax.set_title('All 196 registered bounds passed' if e['scientific_pass'] else 'Registered criterion outcome');ax.legend(frameon=True,facecolor='white',framealpha=.95,fontsize=8,loc='upper left');ax.grid(axis='y',alpha=.15)
    ax=axes[1]
    for c,color,label in zip(range(3),('#7c3aed','#007f8b','#d95f02'),('1-node chemistry','3-node scarce','3-node fast transport')):
        values=[]
        for engine in (1,2,3):values.append(max(row['distance'] for row in e['formal']+e['coarse_medium_descriptive'] if row['cell']==c and row['engine']==engine and row['reference_engine']==0))
        ax.plot(x,values,'o-',lw=2,color=color,label=label)
    ax.axvspan(1.8,2.2,color='#d1d5db',alpha=.3);ax.set_xticks(x,['0.001\nCoarse','0.0005\nMedium','0.00025\nPrespecified finest']);ax.set_xlim(-.15,2.25);ax.set_ylim(0,.026);ax.set_xlabel('Maximum numerical step (adaptive safety cap also applies)');ax.set_ylabel('Largest observed marginal CDF distance');ax.set_title('Coarse/medium comparisons are descriptive');ax.grid(alpha=.15);ax.legend(frameon=False,fontsize=9)
    fig.supxlabel('Marginal endpoint laws on three fixed small fixtures; not full temporal-law or large-lattice pattern validation.',fontsize=10);fig.savefig(BASE/'MG11G4_summary.png',dpi=170);plt.close(fig)
if __name__=='__main__':main()
