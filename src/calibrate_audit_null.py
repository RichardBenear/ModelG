"""MG005: large-sample null check of the within-cell audit (per-cell max studentised LOO, M = 12, 1404 tests)."""
import numpy as np
from scipy import stats
rng=np.random.default_rng(7); M=12; ntests=1404; nsim=200000
crit=stats.t.ppf(1-0.01/(2*ntests), M-2)
# per cell (12 runs): max studentized LOO; vectorised
x=rng.standard_normal((nsim,M))
S=x.sum(1,keepdims=True); SS=(x**2).sum(1,keepdims=True)
mr=(S-x)/(M-1); vr=((SS-x**2)-(M-1)*mr**2)/(M-2)
t=np.abs(x-mr)/np.sqrt(vr*(1+1/(M-1)))
p_cell=(t.max(1)>crit).mean()
print('crit',crit,'per-cell flag prob',p_cell,'-> gate (117 cells-projections)',1-(1-p_cell)**117)
print('single-test tail check', (t[:,0]>crit).mean(), 'expected', 0.01/(2*ntests)*2)
