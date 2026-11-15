"""Explicit open-boundary bridge, not the existing finite-reservoir model."""
import numpy as np
from mg12h3_sparse import Model,Channel
from lattice import P17
class OpenBridge(Model):
 def __init__(self,n,omega,h=.2,P=None):
  P=P17 if P is None else P;self.n=n;channels=[]
  def add(rule,reserve,change,k,factors):
   if k>0:channels.append(Channel(rule,tuple(reserve),tuple(change),float(k),tuple(factors)))
  for i in range(n):
   g,x,y=i,n+i,2*n+i
   add('R1_external',[],[(g,1)],omega*P['a'],[]);add('R1r_external',[(g,1)],[(g,-1)],P['q']-P['p'],[(g,1)]);add('R2',[(g,1)],[(g,-1),(x,1)],P['p'],[(g,1)]);add('R2r',[(x,1)],[(x,-1),(g,1)],P['g'],[(x,1)]);add('R3_external',[(x,1)],[(x,-1),(y,1)],P['b'],[(x,1)]);add('R4',[(x,2),(y,1)],[(x,1),(y,-1)],1/omega**2,[(x,2),(y,1)]);add('R5_external',[(x,1)],[(x,-1)],1-P['g'],[(x,1)])
   for s,d in enumerate((1.,P['dx'],P['dy'])):
    for j in ((i-1)%n,(i+1)%n):a=s*n+i;b=s*n+j;add('H_'+str(s),[(a,1)],[(a,-1),(b,1)],d/h**2,[(a,1)])
  self.channels=tuple(channels);deps=[[] for _ in range(7*n)]
  for j,c in enumerate(channels):
   for pos,_ in c.factors:deps[pos].append(j)
  self.dependencies=tuple(tuple(v) for v in deps)
def generator(pop,omega,h=.2):
 n=pop.shape[1];G,X,Y=pop[:3]/omega;P=P17;f4=X*np.maximum(0,X-1/omega)*Y;out=np.array([P['a']-P['q']*G+P['g']*X,P['p']*G-X-P['b']*X+f4,P['b']*X-f4])
 for s,d in enumerate((1.,P['dx'],P['dy'])):out[s]+=d/h**2*(np.roll(pop[s]/omega,1)+np.roll(pop[s]/omega,-1)-2*pop[s]/omega)
 return out
