"""Static MG9b Draft 1 checks. No ModelG integration.
Usage: python design_checks.py /path/to/ModelG
"""
import json,sys
from pathlib import Path
from unittest.mock import patch
import numpy as np
from scipy.optimize import minimize

root=Path(sys.argv[1]); sys.path.insert(0,str(root/'src'))
import confirm
captured={}
def fake_ode(P,N,h,t_end,**kw):
    captured.update(t_end=t_end,t_eval=kw['t_eval'].tolist())
    return np.arange(N)*h,kw['t_eval'],np.zeros((3*N,len(kw['t_eval']))),0.
with patch.object(confirm,'lattice_ode_1d',side_effect=fake_ode):
    confirm.confirm_field(np.zeros((3,20)),29.35,.2,-1.,3.,t_end=2000.)

# Deterministic illustrative count table, NOT ModelG data or a coverage study.
# Independent b intercepts; centered Omega gives the same slope-trend family.
bs=np.array([29.,29.2,29.35]); slopes=np.array([.064,.080,.020])
omega=np.tile([100.,115.,130.],3); bi=np.repeat(np.arange(3),3)
w=omega-115.; db=bs[bi]-29.175
X=np.column_stack([np.eye(3)[bi],-w,-w*db])
counts=np.rint(10000*np.exp(-slopes[bi]*w))
def fit(null=False):
    A=X[:,:4] if null else X
    def fun(beta):
        eta=A@beta; mu=np.exp(eta)
        return float(np.sum(mu-counts*eta)),A.T@(mu-counts)
    init=np.r_[np.repeat(np.log(10000),3),.06] if null else np.r_[np.repeat(np.log(10000),3),.06,0.]
    r=minimize(fun,init,jac=True,method='BFGS',options={'gtol':1e-7})
    return r
full,zero=fit(),fit(True)
out={'revision':'98d63b8659260879d7c674956610032707415142',
     'confirmation_t_end_only_probe':captured,
     'lack_of_fit_df':{'two_cells':0,'three_cells':1,'four_cells':2},
     'nonmonotone_example':{'b':bs.tolist(),'true_S':slopes.tolist(),'fitted_linear_S1':float(full.x[-1]),
         'likelihood_ratio_against_S1_zero':float(2*(zero.fun-full.fun)),
         'note':'Negative fitted trend despite an increase then decrease; intended only as a model-wording counterexample.'},
     'coverage_mc_se_at_0_93_R2000':float(np.sqrt(.93*.07/2000)),
     'duplicate_replacement_example':{'initial':[150,135,120],'both_lower_cells_fail_replacement_each':165,
         'distinct_cells_after_both_replacements':[150,165]}}
assert max(captured['t_eval'])==300 and captured['t_end']==2000
assert full.x[-1]<0 and out['nonmonotone_example']['likelihood_ratio_against_S1_zero']>3.8415
Path(__file__).with_name('design_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
