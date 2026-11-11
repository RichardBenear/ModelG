from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch,Circle
import numpy as np
OUT=Path(__file__).resolve().parent;OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'svg.fonttype':'none'})
fig=plt.figure(figsize=(16,12),facecolor='#f7f9fc');ink='#17324d';muted='#51657a';blue='#2778af';green='#17866a';purple='#7953ac';gold='#c78928'
fig.text(.045,.95,'The spatial vocabulary of Model G',fontsize=27,weight='bold',color=ink)
fig.text(.045,.913,'A site is a location. Its state is the collection of token counts stored there.',fontsize=14,color=muted)
ax=fig.add_axes([.035,.565,.93,.30]);ax.set_xlim(-.5,17.3);ax.set_ylim(-1.45,3.0);ax.axis('off')
# Schematic row with labels that explicitly preserve number of sites.
positions=[0,1.4,2.8,4.2,5.6,7.0,12.4,13.8,15.2];labels=['Site 0','Site 1','Site 2','Site 3','Site 4','Site 5','Site 497','Site 498','Site 499']
for x,label in zip(positions,labels):
 ax.add_patch(FancyBboxPatch((x,.75),1.3,1.15,boxstyle='round,pad=.015,rounding_size=.07',facecolor='white',edgecolor=blue,lw=1.8));ax.text(x+.65,2.12,label,ha='center',fontsize=11,color=ink)
 if x<8:
  for xx,yy,col in [(x+.28,1.5,blue),(x+.73,1.5,green),(x+1.0,1.1,purple)]:ax.add_patch(Circle((xx,yy),.10,color=col))
 else:
  ax.text(x+.65,1.325,'token\ncounts',ha='center',va='center',fontsize=10,color=muted)
ax.text(10.2,1.3,'⋯',ha='center',va='center',fontsize=37,color=muted)
ax.text(10.2,.74,'Sites 6–496 omitted',ha='center',fontsize=10,color=muted)
# Integrated ruler: five intervals span the centers of six sites.
for i in range(5):
 x=.65+1.4*i
 ax.add_patch(FancyBboxPatch((x,.37),1.4,.20,boxstyle='square,pad=0',facecolor=green if i%2==0 else '#b7ded2',edgecolor='white',lw=1))
for i in range(6):
 x=.65+1.4*i
 ax.plot([x,x],[.30,.69],color=green,lw=1.3)
ax.annotate('',xy=(.65,.09),xytext=(2.05,.09),arrowprops=dict(arrowstyle='<->',color=blue,lw=1.4))
ax.text(1.35,-.16,'h = 0.2',ha='center',fontsize=11,weight='bold',color=blue)
ax.text(4.15,-.46,'1 model-length unit = 5 × 0.2',ha='center',fontsize=13,weight='bold',color=green)
# Wraparound route from last to first; label is not a spatial ruler.
ax.add_patch(FancyArrowPatch((15.85,.72),(.65,.72),connectionstyle='bar,fraction=-.11',arrowstyle='<->',mutation_scale=14,color=purple,lw=1.8))
ax.text(8.1,-.73,'PERIODIC: Site 499 connects back to Site 0',ha='center',fontsize=13,color=purple,weight='bold',bbox=dict(facecolor='#f7f9fc',edgecolor='none',pad=3))
ax.text(8.1,-1.31,'This is a connection, not an extra compartment or a diffusion cycle.',ha='center',fontsize=11,color=muted)
fig.text(.5,.535,'500 compartments × 0.2 model-length units = domain length L = 100',ha='center',fontsize=18,weight='bold',color=ink)
# Three vocabulary cards.
card_y=.265;card_h=.225
cards=[(.045,'SITE / COMPARTMENT',blue),(.365,'MODEL-LENGTH UNIT',green),(.685,'DOMAIN / PERIODIC',purple)]
for x,title,color in cards:
 rect=FancyBboxPatch((x,card_y),.285,card_h,transform=fig.transFigure,boxstyle='round,pad=.012,rounding_size=.012',facecolor='white',edgecolor='#ced9e5',lw=1.2,zorder=-1);fig.add_artist(rect);fig.text(x+.012,card_y+card_h-.033,title,color=color,fontsize=13,weight='bold')
fig.text(.057,.408,'A well-mixed location containing',fontsize=12,color=ink)
fig.text(.057,.381,'G, X, Y, A, B, Z and W.',fontsize=12,color=ink)
fig.text(.057,.346,'Site state = seven token counts:',fontsize=11.5,color=muted)
fig.text(.057,.315,'(nG, nX, nY, nA, nB, nZ, nW)',fontsize=12,color=blue,weight='bold')
fig.text(.057,.282,'Internal token positions are not modeled.',fontsize=10.5,color=muted)
# The ruler is now directly aligned with site centers above.
fig.text(.377,.408,'The green ruler spans five site spacings.',fontsize=11.5,color=ink)
fig.text(.377,.372,'Each spacing is 0.2 model-length units.',fontsize=11.5,color=ink)
fig.text(.377,.336,'Together: 1 model-length unit.',fontsize=12,color=green,weight='bold')
fig.text(.377,.288,'Physical size in meters is uncalibrated.',fontsize=10.5,color=muted)
fig.text(.697,.408,'Domain = the whole modeled region.',fontsize=12,color=ink)
fig.text(.697,.375,'500 sites; total length = 100 units.',fontsize=12,color=purple,weight='bold')
fig.text(.697,.337,'Periodic = the two ends connect.',fontsize=12,color=ink)
fig.text(.697,.300,'No endpoints. No required circulation.',fontsize=11,color=muted)
# Wavelength inset explicitly distinguish.
w=fig.add_axes([.045,.075,.53,.135]);w.set_xlim(0,10);w.set_ylim(-1.75,1.65);w.axis('off');x=np.linspace(0,10,600);lam=3.06;w.plot(x,.65*np.cos(2*np.pi*(x-1)/lam),color=gold,lw=2)
w.annotate('',xy=(1,1.02),xytext=(1+lam,1.02),arrowprops=dict(arrowstyle='<->',color=gold,lw=1.6));w.text(1+lam/2,1.32,'One pattern wavelength',ha='center',fontsize=12,weight='bold',color=gold)
w.text(0,-1.22,'Earlier characteristic wavelength ≈ 3.06 units ≈ 15.3 site spacings.',fontsize=11,color=ink)
w.text(0,-1.65,'Schematic concentration pattern; not the full 100-unit domain.',fontsize=10.5,color=muted)
fig.text(.63,.186,'WAVELENGTH ≠ DOMAIN LENGTH',fontsize=13,weight='bold',color=gold)
fig.text(.63,.15,'Wavelength: distance over which a pattern repeats.',fontsize=11.5,color=ink)
fig.text(.63,.117,'Diffusion: random spreading; it has no “full cycle.”',fontsize=11.5,color=ink)
fig.text(.045,.035,'Current fixed 1D lattice interpretation: compartments represent segments. Layout and marker sizes are schematic.',fontsize=11,color=muted)
fig.savefig(OUT/'model_g_spatial_terms_integrated.png',dpi=180,facecolor=fig.get_facecolor(),bbox_inches='tight');fig.savefig(OUT/'model_g_spatial_terms_integrated.svg',facecolor=fig.get_facecolor(),bbox_inches='tight')
print((OUT/'model_g_spatial_terms_integrated.png').resolve())
