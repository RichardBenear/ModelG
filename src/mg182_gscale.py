"""G-focused population scaling, direct observable and unrounded mean field."""
import json
import numpy as np
from scipy.sparse import coo_matrix
from mg12h3_sparse import Model
from run_mg176 import ROOT
PATH=ROOT/'results/MG12H3/ring_feasibility/o0_a0_Omega10.npz'

def setup(omega):
 with np.load(PATH) as z:
  ref=z['concentrations'].copy();old=float(z['omega']);a=z['a0']/old;b=z['b0']/old;links=json.loads(str(z['links']));diff=tuple(z['diff']);coupling=tuple(z['coupling'])
 pop=np.rint(ref*omega).astype(np.int64);model=Model(500,(float(omega),list(zip(a*omega,b*omega)),links,diff,coupling,(True,)*6));return dict(pop=pop,reference=ref,omega=float(omega),a=a,b=b,links=links,diff=diff,coupling=coupling,model=model)
def projection(field,reference):
 g=reference[0]-reference[0].mean();u=g/np.linalg.norm(g);f=field[0]-field[0].mean();q=float(f@u);return dict(G_projection=q,G_reference_norm=float(np.linalg.norm(g)),orthogonal_G_norm=float(np.linalg.norm(f-q*u)))
def rhs(flat,c,finite_omega=None):
 z=np.asarray(flat).reshape(7,500);G,X,Y,A,B=z[:5];f0=14*A/c['a'];f1=G;fr=.1*X;f3=29*B*X/c['b'];f4=X*X*Y if finite_omega is None else X*np.maximum(0,X-1/finite_omega)*Y;f5=.9*X;out=np.array([f0-f1+fr,f1-fr-f3+f4-f5,f3-f4,-f0,-f3,f3,f5])
 for s in range(3):out[s]+=c['diff'][s]*(np.roll(z[s],1)+np.roll(z[s],-1)-2*z[s])
 for kind,i,j,w in c['links']:
  if kind in ('ALink','BLink'):
   s=3 if kind=='ALink' else 4;flow=c['coupling'][s-3]*w*z[s,i];out[s,i]-=flow;out[s,j]+=flow
 return out.ravel()
def sparsity(model):
 rows=[];cols=[]
 for ch in model.channels:
  for dest,_ in ch.change:
   for source,_ in ch.factors:rows.append(dest);cols.append(source)
 return coo_matrix((np.ones(len(rows)),(rows,cols)),shape=(3500,3500)).tocsr()
