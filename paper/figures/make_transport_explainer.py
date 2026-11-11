from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle,FancyArrowPatch,FancyBboxPatch
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'svg.fonttype':'none'})
OUT=Path(__file__).resolve().parent;fig=plt.figure(figsize=(17,12),facecolor='#f7f9fc')
ink='#18304b';muted='#52657a';blue='#2874ae';green='#19866b';purple='#8053ad';orange='#db8b26'
fig.text(.045,.955,'How transport works in the current Model G runs',fontsize=24,weight='bold',color=ink)
fig.text(.045,.92,'500 fixed sites • 12 shown schematically • every connection permits movement in both directions',fontsize=13,color=muted)
angles=np.pi/2-np.arange(12)*2*np.pi/12;xy=np.c_[np.cos(angles),np.sin(angles)]
for index,(title,subtitle,kind,color) in enumerate([
 ('G / X / Y transport','Neighbors on the periodic ring','active',blue),
 ('A / B option: LOCAL','Ring + extra nearby connections','local',green),
 ('A / B option: CROSS_RING','Ring + opposite-site connections','cross',purple)]):
 ax=fig.add_axes([.035+index*.325,.545,.305,.345]);ax.set_aspect('equal');ax.axis('off');ax.set_xlim(-1.4,1.4);ax.set_ylim(-1.38,1.48)
 ax.text(0,1.39,title,ha='center',weight='bold',fontsize=15,color=ink);ax.text(0,1.17,subtitle,ha='center',fontsize=11,color=muted)
 edges=[(i,(i+1)%12) for i in range(12)]
 for i,j in edges:ax.plot(*np.array([xy[i],xy[j]]).T,color=color,lw=2,alpha=.75,zorder=1)
 extra=[]
 if kind=='local':extra=[(i,(i+3)%12) for i in range(0,12,2)]
 if kind=='cross':extra=[(i,i+6) for i in range(6)]
 for i,j in extra:ax.plot(*np.array([xy[i],xy[j]]).T,color=color,lw=1.8,alpha=.60,zorder=1)
 for i,(x,y) in enumerate(xy):
  ax.add_patch(Circle((x,y),.085,facecolor='white',edgecolor=color,lw=1.8,zorder=3));ax.text(x,y,str(i),ha='center',va='center',fontsize=9,color=ink,zorder=4)
 if kind=='active':
  ax.text(0,.23,'G and X: 25',ha='center',fontsize=13,color=blue,weight='bold');ax.text(0,-.04,'Y: 300',ha='center',fontsize=13,color=blue,weight='bold');ax.text(0,-.30,'per outgoing link',ha='center',fontsize=10,color=muted)
 else:ax.text(0,-1.25,'A and B: 1/3 per outgoing link',ha='center',fontsize=10.5,color=color)
fig.text(.045,.505,'All three drawings refer to the same sites. A and B can use different network options.',fontsize=12,color=ink)
fig.text(.045,.478,'These are fixed pair connections in the recent ring studies; line crossings do not create extra sites.',fontsize=11,color=muted)
# Site contents panel
ax=fig.add_axes([.045,.16,.33,.275]);ax.axis('off');ax.set_xlim(0,10);ax.set_ylim(0,7)
ax.text(0,6.65,'What occupies each site?',fontsize=16,weight='bold',color=ink)
box=FancyBboxPatch((.1,1.2),8.9,4.8,boxstyle='round,pad=.15,rounding_size=.3',facecolor='white',edgecolor='#b9c8d9',lw=1.5);ax.add_patch(box)
for label,x,y,col in [('G',1.6,4.7,blue),('X',3.6,4.7,blue),('Y',5.6,4.7,blue),('A',2.6,3.2,green),('B',4.6,3.2,green),('Z',6.6,2.0,'#8b96a4'),('W',8.0,2.0,'#8b96a4')]:
 ax.add_patch(Circle((x,y),.43,facecolor=col,edgecolor='none'));ax.text(x,y,label,ha='center',va='center',weight='bold',color='white',fontsize=12)
ax.text(.5,1.9,'Z / W stay at their site',fontsize=10,color=muted)
ax.text(.2,.56,'Tokens represent units of model species.',fontsize=11,color=ink)
ax.text(.2,.03,'Reactions happen locally; transport changes the site.',fontsize=10.5,color=muted)
# strict tick panel
ax=fig.add_axes([.405,.16,.55,.275]);ax.axis('off');ax.set_xlim(0,15);ax.set_ylim(0,7)
ax.text(0,6.65,'One participation per token per tick',fontsize=16,weight='bold',color=ink)
ax.text(0,6.05,'Example token A: a move is possible, not guaranteed.',fontsize=11,color=muted)
positions=[1.5,7,12.5]
for x,label in zip(positions,['Start of tick n','Accepted during tick n','Start of tick n + 1']):ax.text(x,5.25,label,ha='center',fontsize=11,weight='bold',color=ink)
for center in positions:
 for offset,label in [(-.85,'i'),(.85,'j')]:
  ax.add_patch(Circle((center+offset,3.9),.4,facecolor='white',edgecolor='#b9c8d9',lw=1.5));ax.text(center+offset,3.9,label,ha='center',va='center',color=muted,fontsize=10)
 ax.plot([center-.45,center+.45],[3.9,3.9],color='#c7d3df',lw=2,zorder=0)
for center,offset in [(positions[0],-.85),(positions[1],.85),(positions[2],.85)]:
 ax.add_patch(Circle((center+offset,3.9),.27,facecolor=green,edgecolor='white',zorder=4));ax.text(center+offset,3.9,'A',ha='center',va='center',color='white',weight='bold',fontsize=10,zorder=5)
for start,end in [(2.8,5.3),(8.3,10.8)]:ax.add_patch(FancyArrowPatch((start,3.9),(end,3.9),arrowstyle='-|>',mutation_scale=18,color=ink,lw=1.7))
for center,lines in zip(positions,[['Available in the','starting snapshot'],['Moved i → j','Reserved for this tick'],['Eligible again','in the next tick']]):
 ax.text(center,2.55,'\n'.join(lines),ha='center',va='top',fontsize=11,color=ink,linespacing=1.4)
ax.text(0,.65,'Arrivals, products and catalysts cannot act twice within a tick.',fontsize=11,color=muted)
ax.text(0,.04,'The next tick starts from the preceding tick’s output.',fontsize=11,color=muted)
fig.text(.045,.10,'Rates are per model-time unit. Tick duration in physical seconds has not been calibrated.',fontsize=12,color=ink)
fig.text(.045,.065,'The circular layout is a drawing of a one-dimensional periodic ring, not evidence of two-dimensional physical space.',fontsize=11,color=muted)
fig.savefig(OUT/'model_g_transport_explainer.png',dpi=170,facecolor=fig.get_facecolor(),bbox_inches='tight')
fig.savefig(OUT/'model_g_transport_explainer.svg',facecolor=fig.get_facecolor(),bbox_inches='tight')
print((OUT/'model_g_transport_explainer.png').resolve())
