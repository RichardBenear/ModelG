"""Report the three prespecified lifetime fits; no additional model search."""
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import logsumexp
from lifetime_mg188 import OUT,read,atomic_json

def run():
 a=read(OUT/'fits.json');p=a['models']['exponential_in_Omega'];cohorts=a['cohorts'];rows=a['rows'];t=np.array([r['time'] for r in rows]);e=np.array([r['event'] for r in rows]);x=(np.array([r['omega'] for r in rows])-1000)/1000
 # Independent likelihood arithmetic from aggregated deaths and exposure.
 rate=np.array([np.exp(p['alpha']+p['beta']*(c['omega']-1000)/1000) for c in cohorts]);ll=sum(c['events']*np.log(rr)-rr*c['exposure'] for c,rr in zip(cohorts,rate));np.testing.assert_allclose(ll,p['loglik'],atol=1e-9);factor=np.exp(-p['beta']);np.testing.assert_allclose(factor,p['factor_per_1000']);zero=cohorts[-1];np.testing.assert_allclose(np.exp(-zero['exposure']/zero['zero_event_mean_lower95']),.05,atol=1e-12)
 atomic_json(OUT/'separate_arithmetic_check.json',dict(passed=True,aggregated_likelihood=float(ll),lifetime_factor=float(factor),zero_event_lower_bound_checked=True,external_replication=False))
 with (OUT/'observations.csv').open('w',newline='',encoding='utf-8') as file:
  w=csv.DictWriter(file,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 colors=plt.cm.viridis(np.linspace(.05,.9,5));fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained');grid=np.linspace(0,980,500)
 for color,c in zip(colors,cohorts):
  k=c['KM'];tt=[v['time'] for v in k];ss=[v['survival'] for v in k];axes[0].step(tt,ss,where='post',color=color,label=f"Ω={c['omega']} ({c['events']} losses)");axes[0].fill_between(tt,[v['lower'] for v in k],[v['upper'] for v in k],step='post',color=color,alpha=.08);rr=np.exp(p['alpha']+p['beta']*(c['omega']-1000)/1000);axes[0].plot(grid,np.exp(-rr*grid),color=color,ls='--',lw=.9)
 axes[0].set(xlabel='Model time',ylabel='Original-track survival probability',ylim=(-.03,1.03),xlim=(0,980),title='Kaplan–Meier curves and specified exponential fit');axes[0].legend(fontsize=8);axes[0].text(.38,.60,'Solid: observations with pointwise 95% bands\nDashed: fitted constant-hazard survival',transform=axes[0].transAxes,fontsize=8)
 om=np.linspace(100,3000,300);axes[1].plot(om,np.exp(p['alpha']+p['beta']*(om-1000)/1000),label='Fitted rate, Ω=100–3000');beyond=np.linspace(3000,10000,200);axes[1].plot(beyond,np.exp(p['alpha']+p['beta']*(beyond-1000)/1000),ls='--',color='tab:blue',label='Extrapolation beyond Ω=3000')
 for color,c in zip(colors,cohorts):
  if c['events']:axes[1].scatter(c['omega'],c['constant_hazard_rate'],color=color,s=35)
  else:
   upper=1/c['zero_event_mean_lower95'];axes[1].scatter(c['omega'],upper,marker='v',color=color,s=50);axes[1].annotate('0 losses: one-sided 95% rate bound',xy=(c['omega'],upper),xytext=(1800,.0015),arrowprops={'arrowstyle':'->'},fontsize=8)
 axes[1].set(yscale='log',xlabel='Population scale Ω',ylabel='Constant-hazard loss rate per model-time unit',title='Rate falls exponentially with Ω in the candidate law');axes[1].legend(fontsize=8,loc='lower left');fig.savefig(OUT/'MG188_lifetime.png',dpi=160);plt.close(fig)
 def f(v):return f'{v:.6g}'
 cif=p['factor_CI95'];secondary=a['secondary_solitary_exponential'];lines=['# MG188 lifetime law from saved trajectories','',
 'The saved trajectories support a much stronger statement than a survival count: a prespecified censored-lifetime fit estimates that mean original-track lifetime grows by about 9.4 times for every additional 1,000 in population scale Ω. Its model-conditional 95% interval is about 6.2–14.9 times. The exponential-in-Ω candidate passed the fixed diagnostic screens and substantially outperformed the prespecified power law; allowing a time-dependent Weibull hazard gave no decisive improvement. This is a useful empirical lifetime law for these prepared Phase 2 patterns, not a universal physical law or a measured lifetime at Ω=10,000. At that population no loss was observed, so its mean remains unmeasured.','',
 '## Quotable statement','',
 '> For the prepared Phase 2 Model G patterns studied here, a censored exponential-lifetime fit estimates a roughly 9.4-fold increase in mean original-track lifetime per additional 1,000 in Ω (model-conditional 95% interval: 6.2–14.9-fold).','',
 '“Original track” matters: this endpoint follows the original X-defined pattern until its recorded disappearance, even if replication occurred. It is neither a G-only shape-persistence test nor a guarantee that a single isolated particle remains intact. MG187 found that G and X usually fade together, which supports interpreting this as a coupled-pattern lifetime, but does not redefine the endpoint.','',
 '## Fixed analysis and censoring','',
 'All 150 prepared runs were included, 30 per population; controls have no prepared original track and were excluded. Fit choices and diagnostic rules were committed and published before computation, in [the specification](../../spec/MG188_lifetime_analysis.md). Prior survival counts were already known, so this remains an explicitly post hoc analysis. No new simulations or random resampling occurred.','',
 'The detector timestamps loss at the last matched observation and needs twenty later unmatched units to confirm it. Every run was therefore administratively censored at time 980, leaving the same confirmation window. Loss events at or before 980 count as events; other runs are right-censored. This avoids silently treating unconfirmable terminal disappearances as survived lifetimes. Observation spacing is 0.5 model-time units; there is no calibrated conversion to seconds.','',
 '## Empirical survival results','',
 '| Ω | Losses / 30 | Right-censored | KM median lifetime | Restricted mean through 980 | Losses / exposure rate |','|---:|---:|---:|---|---:|---:|']
 for c in cohorts:lines.append(f"| {c['omega']} | {c['events']}/30 | {c['censored']} | {'not reached' if c['KM_median'] is None else f(c['KM_median'])} | {f(c['restricted_mean_980'])} | {f(c['constant_hazard_rate'])} |")
 lines+=['','Restricted mean is the average time alive within the 980-unit observation window; it does not extrapolate beyond it. For Ω=100 and 300 all losses were observed, so it equals their sample mean lifetime. At Ω=3,000 and 10,000 the median is not reached; quoting 980 as their unrestricted mean would be wrong.','',
 '![Observed survival curves and population-dependent loss rate](MG188_lifetime.png)','',
 'The survival bands are pointwise Greenwood log-log intervals. Boundary bands at survival 1 with no events are degenerate and should not be read as certainty of perfect reliability; use the event-count intervals and censoring limits. The plotted zero-event rate bound is one-sided and assumes a constant time hazard.','',
 '## Fitted law','',
 f"The primary model is λ(Ω) = exp({f(p['alpha'])} − {f(-p['beta'])} × (Ω−1000)/1000), with S(t|Ω)=exp(−λt). Equivalently:",
 '',f"**Mean lifetime μ(Ω) ≈ {f(np.exp(-p['alpha']))} × exp({f(-p['beta'])} × (Ω−1000)/1000) model-time units.**",'',
 f"Each +1000 in Ω multiplies this fitted mean by {f(factor)}, with a 95% profile-likelihood interval [{f(cif[0])}, {f(cif[1])}]. This interval is conditional on independent trajectories and the chosen constant-hazard, log-linear-in-Ω model. It does not include model-selection uncertainty or prove extrapolation.",
 '', '| Prespecified candidate | Parameters | Log likelihood | AICc |','|---|---:|---:|---:|']
 for name,m in a['models'].items():lines.append(f"| {name} | {m['parameters']} | {f(m['loglik'])} | {f(m['AICc'])} |")
 weib=a['models']['Weibull_in_Omega'];delta=p['AICc']-weib['AICc'];lines+=['',f"The power-law candidate is worse by {f(a['models']['power_in_Omega']['AICc']-p['AICc'])} AICc units. Weibull improves AICc by only {f(delta)}, below the prespecified six-unit concern threshold. Its common shape is {f(weib['shape'])}, with approximate 95% Wald interval [{f(weib['shape_CI95'][0])}, {f(weib['shape_CI95'][1])}], which includes the constant-hazard value 1. Every primary-model predicted event fraction lies inside its population's observed exact binomial interval. These are compatibility screens, not proof of the law or calibrated global goodness-of-fit tests.",
 '', '| Ω | Observed losses | Predicted losses by 980, primary fit |','|---:|---:|---:|']
 for c in cohorts:lines.append(f"| {c['omega']} | {c['events']} | {f(c['predictions']['exponential_in_Omega']['expected_events'])} |")
 lines+=['','## What the high-population censoring can and cannot tell us','',
 f"At Ω=3,000, one observed loss gives a separate constant-hazard mean estimate of {f(cohorts[3]['constant_hazard_mean'])}; that extrapolated mean is very uncertain because it rests on just one event. At Ω=10,000 there is no finite separate maximum-likelihood estimate of mean lifetime. Under the constant-hazard assumption, its 30 censored runs provide only a one-sided 95% lower bound of {f(cohorts[4]['zero_event_mean_lower95'])} model-time units.",
 '', f"The population law extrapolates an Ω=10,000 mean of {f(cohorts[4]['predictions']['exponential_in_Omega']['extrapolated_mean'])} units—vastly beyond the 980-unit follow-up. **That number is a mathematical extrapolation, not an observed or validated lifetime, and should not be quoted as a persistence result.** The slope is principally informed by losses at Ω=100, 300, 1000 and the single loss at 3000. Zero losses at10000 are compatible with the trend but do not measure its continuation across the large population gap.",
 '', '## Original disappearance versus solitary-particle lifetime','',
 f"The prespecified secondary endpoint ends at the first loss, replication or merge. It gives an exponential-in-Ω factor of {f(secondary['factor_per_1000'])} per +1000, materially larger than the primary {f(factor)}. All30 Ω=100 runs had replication or merge events, compared with one at300 and none at1000 or above. Early replication shortens solitary-particle persistence even when the original track continues. The quotable primary statement must therefore retain its original-track endpoint; these two lifetime definitions cannot be interchanged.",
 '', '## Verification and conclusion','',
 'Event times were independently re-extracted from the immutable per-run summaries and matched to the MG185 exports. Every Kaplan–Meier restricted mean matched the direct mean of capped observation times. The analytic exponential profile likelihood matched a general numerical optimizer; a separate arithmetic path reconstructed its likelihood from population event/exposure totals and checked the zero-event lower bound. Censoring fixtures passed. Same author and shared legacy X detector/NumPy/SciPy remain dependencies; this is not independent external replication. Full precision, individual observation times, curves and all candidate fits are retained in [fits.json](fits.json), [observations.csv](observations.csv), [verification](verification.json) and [arithmetic check](separate_arithmetic_check.json).','',
 '**Conclusion:** “roughly ten times longer mean original-track lifetime per +1,000 Ω” is supported as a fitted empirical summary of this dataset, with the stated model dependence and uncertainty. It is stronger than raw survival counts because it uses when losses occurred and accounts for right-censoring. It remains a finite-data, endpoint-specific law for this implementation, not proof of indefinite persistence. No new simulations, further model search or historical rescoring follows this report. PI closure is pending.']
 (OUT/'MG188_results.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
if __name__=='__main__':run()
