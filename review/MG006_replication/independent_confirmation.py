"""Independent Radau continuation and structure matching for stored MG3 fields (equations implemented separately)."""
import json,pickle,sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags,kron,coo_matrix
ROOT=Path(__file__).parent;INPUT=ROOT/'inputs'
class Limited(pickle.Unpickler):
    def find_class(self,m,n):
        d={('numpy._core.multiarray','_reconstruct'):np._core.multiarray._reconstruct,('numpy.core.multiarray','_reconstruct'):np._core.multiarray._reconstruct,('numpy','ndarray'):np.ndarray,('numpy','dtype'):np.dtype,('numpy._core.multiarray','scalar'):np._core.multiarray.scalar,('numpy.core.multiarray','scalar'):np._core.multiarray.scalar}
        if (m,n) not in d:raise ValueError((m,n))
        return d[m,n]
def confirm(phi,b,h=.2):
    N=phi.shape[1];x0=140/9;y0=b/x0;T=np.arange(200.,301.,10.)
    from scipy.sparse import lil_matrix
    L=diags([np.ones(N-1),-2*np.ones(N),np.ones(N-1)],[-1,0,1],format='lil')
    L[0,-1]=1;L[-1,0]=1;L=L.tocsc()/h**2
    D=kron(L,diags([1.,1.,12.]),format='csc')
    rows=np.repeat(np.arange(3*N),3);cols=np.repeat(np.arange(N),9)*3+np.tile(np.arange(3),3*N)
    def rhs(t,y):
        G,X,Y=y.reshape(N,3).T
        transfer=2*b*X+y0*X*X+(x0+X)**2*Y
        return D@y+np.column_stack([-G+.1*X,G-(1+b)*X+transfer,b*X-transfer]).ravel()
    def jac(t,y):
        f=y.reshape(N,3);X=x0+f[:,1];Y=y0+f[:,2]
        v=np.zeros((N,3,3));v[:,0,0]=-1;v[:,0,1]=.1;v[:,1,0]=1
        v[:,1,1]=-(1+b)+2*X*Y;v[:,1,2]=X*X;v[:,2,1]=b-2*X*Y;v[:,2,2]=-X*X
        return D+coo_matrix((v.ravel(),(rows,cols)),shape=D.shape).tocsc()
    sol=solve_ivp(rhs,[0,300],phi.T.ravel(),method='Radau',jac=jac,rtol=1e-9,atol=1e-11,max_step=.5,t_eval=T)
    if not sol.success:raise RuntimeError(sol.message)
    profiles=np.load(INPUT/'results/MG003/ref_profiles.npz');r=profiles[f'h{h}_b{b}'];A=r[2,N//2]
    xx=r[0];right=xx>=0;rf=xx[right];f=r[2,right];i=np.flatnonzero(f[:-1]*f[1:]<0)
    z=rf[i]-h*f[i]/(f[i+1]-f[i]);z=z[z<30];lam=2*np.mean(np.diff(z)[1:])
    minima=[]
    for k in range(len(T)):
        f=sol.y[:,k].reshape(N,3)[:,1]
        minima.append(np.flatnonzero((f<np.roll(f,1))&(f<=np.roll(f,-1))&(f<-.5*abs(A))))
    ends=[]
    for start in minima[0]:
        pos=int(start);valid=True
        for nxt in minima[1:]:
            if not len(nxt):valid=False;break
            dist=np.minimum(abs(nxt-pos),N-abs(nxt-pos))*h;j=int(np.argmin(dist))
            if dist[j]>lam:valid=False;break
            pos=int(nxt[j])
        if valid:ends.append(pos)
    last=sol.y[:,-1]
    cls='persistent' if ends else ('excursion' if max(abs(last))<1e-3 else 'unresolved')
    return dict(cls=cls,persistent_structures=len(set(ends)),deep_minima_end=len(minima[-1]),min_X=float(last.reshape(N,3)[:,1].min()),max_dt=float(max(abs(rhs(300,last)))))

def job(item):
    cell,rep,phi,b=item
    out=confirm(phi.astype(float),b)
    return dict(cell=cell,rep=rep,claimed='persistent',**out,agrees=out['cls']=='persistent')

if __name__=='__main__':
    jobs=[]
    for cell in ['B_150','B_120','B_100']:
        for p in sorted((INPUT/'results/MG3/runs'/cell).glob('*.pkl')):
            with p.open('rb') as f:r=Limited(f).load()
            if r['first'] is not None:jobs.append((cell,r['rep'],r['first']['field'],r['b']))
    out=[]
    with ProcessPoolExecutor(max_workers=2) as pool:
        for row in pool.map(job,jobs):
            out.append(row);print(row['cell'],row['rep'],row['cls'],flush=True)
    (ROOT/'independent_confirmation.json').write_text(json.dumps(out,indent=2))
    print('Agreement',sum(r['agrees'] for r in out),'/',len(out))
