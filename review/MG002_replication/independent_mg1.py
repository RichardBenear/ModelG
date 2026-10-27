import json, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags, kron, eye, block_diag

def run(h, dim=3):
    R=50.; N=round(R/h); r=np.arange(N)*h
    # Pointwise centered radial Laplacian, with the regular-origin limit n*f''(0).
    main=np.full(N,-2/h**2); upper=np.full(N-1,1/h**2); lower=np.full(N-1,1/h**2)
    main[0]=-2*dim/h**2; upper[0]=2*dim/h**2
    upper[1:] += (dim-1)/(2*h*r[1:-1])
    lower[:] -= (dim-1)/(2*h*r[1:])
    L=diags([lower,main,upper],[-1,0,1],format='csc')
    diffusion=kron(L,diags([1.,1.,12.]),format='csc')
    X0=140/9; Y0=261/140
    spatial=np.exp(-r*r/2)
    def rhs(t,y):
        f=y.reshape(N,3); G,X,Y=f.T
        # Expanded perturbation polynomial avoids subtracting background reactions.
        transfer=2*X0*Y0*X+Y0*X*X+(X0+X)**2*Y
        reaction=np.column_stack([-G+.1*X,G-30*X+transfer,-transfer+29*X])
        reaction[:,1]-=spatial*np.exp(-(t-10)**2/18)
        return diffusion@y+reaction.ravel()
    rows=np.repeat(np.arange(3*N),3); cols=(np.repeat(np.arange(N),9)*3+np.tile(np.arange(3),3*N))
    from scipy.sparse import coo_matrix
    def jac(t,y):
        f=y.reshape(N,3); X=X0+f[:,1]; Y=Y0+f[:,2]
        a=2*X*Y; b=X*X
        v=np.zeros((N,3,3));v[:,0,0]=-1;v[:,0,1]=.1
        v[:,1,0]=1;v[:,1,1]=-30+a;v[:,1,2]=b
        v[:,2,1]=29-a;v[:,2,2]=-b
        return diffusion+coo_matrix((v.ravel(),(rows,cols)),shape=(3*N,3*N)).tocsc()
    start=time.time()
    sol=solve_ivp(rhs,[0,100],np.zeros(3*N),method='Radau',jac=jac,rtol=2e-9,atol=2e-11,max_step=.4,t_eval=[100])
    if not sol.success:raise RuntimeError(sol.message)
    fields=sol.y[:,-1].reshape(N,3)
    wavelengths={}; pooled=[]; radii={}; zero_lists={}
    for i,name in enumerate(['G','X','Y']):
        f=fields[:,i]; ids=np.flatnonzero(f[:-1]*f[1:]<0)
        z=r[ids]-h*f[ids]/(f[ids+1]-f[ids]); radii[name]=float(z[0])
        z=z[z<30]; ds=np.diff(z)[1:];pool=ds.tolist();pooled+=pool
        wavelengths[name]=float(2*np.mean(ds));zero_lists[name]=z.tolist()
    result=dict(dim=dim,h=h,method='Radau',amplitudes=dict(zip(['G','X','Y'],fields[0].tolist())),radii=radii,
                wavelength=float(2*np.mean(pooled)),per_field=wavelengths,zeros=zero_lists,
                max_abs_dt=float(np.max(np.abs(rhs(100,sol.y[:,-1])))),seconds=time.time()-start)
    np.savez(f'review_work/independent_{dim}d_h{h}.npz',r=r,fields=fields)
    with open(f'review_work/independent_{dim}d_h{h}.json','w') as out:json.dump(result,out,indent=2)
    print(json.dumps(result),flush=True)
    return result

if __name__=='__main__':
    import sys
    for h in ([float(x) for x in sys.argv[1:]] or [.05,.025,.0125]):run(h)
