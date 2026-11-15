# MG188 saved-data lifetime analysis — fixed fit specification

PI authorizes survival curves and loss-rate versus Omega using saved MG185 data, with fit specified before computation. This is post hoc model exploration: population counts and MG187 results are already known. No simulations, new tracker, resampling, threshold search or population exclusion. Freeze and publish this specification and analysis code before fitting. One deterministic analysis of 150 prepared trajectories; controls excluded because they have no prepared original particle.

## Endpoint and censoring

Primary lifetime: first recorded original-track `loss` event, in model-time units. Replication/merge are not original-track loss and do not end this primary endpoint. Report their incidence and a prespecified secondary endpoint, first loss/replication/merge, for solitary-particle persistence. Never pool attempts or replace failures. Use all30 original IDs per Omega100,300,1000,3000,10000.

The legacy loss requires20 unmatched time units and timestamps the last hit. To avoid unconfirmed terminal losses being mistaken for survival, set administrative censoring to980 for every trajectory, using only loss/event times at or before980. This discards the last20timeunits uniformly; it does not change MG185 verdicts. The detector samples every0.5unit. High-population no-loss runs are right-censored, not measured infinite or980-unit lifetimes.

## Survival and fits

Publish per-Omega Kaplan–Meier curves with log-log Greenwood pointwise95%bands, median (not reached where appropriate), restricted mean lifetime through980 (not extrapolated mean), total exposure and events. Report separate constant-hazard maximum-likelihood rates events/exposure. For the zero-event population, give the one-sided95% exponential-model lower mean-lifetime bound exposure/(-log .05); it is model-dependent.

Primary candidate: constant-in-time exponential survival with log loss rate alpha+beta*(Omega-1000)/1000. Mean lifetime is1/rate, so its multiplicative increase per+1000Omega is exp(-beta). Fit the full censored likelihood, not regression of survival percentages or measured lifetimes of failures alone. Report95%profile-likelihood slope/factor intervals using the usual chi-square1 cutoff; explicitly conditional on the parametric family.

Prespecified challengers: exponential survival with log rate alpha+beta*log(Omega/1000) (power law in population); Weibull survival with common shape k and log scale alpha+beta*(Omega-1000)/1000 (time-dependent hazard). Fit all three to the same150 observations, publish likelihood/AICc and predicted versus observed per-Omega event counts. Do not select a new population subset. AICc difference>=6 favoring a challenger means the primary law is not adequate to quote as an established lifetime law. If Weibull shape95%Wald interval excludes1, flag nonconstant time hazard. Also flag any fitted event count outside the corresponding exact95%binomial interval around observed event fraction (descriptive multiple-comparison screen). These are diagnostic screens, not proof that a surviving candidate is true. If primary fails a screen, report its factor only as a rejected/qualified model fit, not a measured universal rule. No fourth model or iterative instrument revision follows.

Repeat the primary candidate once for the secondary solitary-particle endpoint; disclose if its slope materially differs. Do not let replication-rich low populations silently redefine original-particle disappearance. Per-Omega rates, survival curves and restricted means remain the main empirical outputs if no simple law is adequate.

## Verification and limits

Bind inputs to published MG185 decision exports and immutable main/continuation summaries; independently re-extract event times from summaries and check rows. Verify restricted mean from KM integral against mean of capped observation times (common administrative censoring); analytic exponential profile versus general numerical likelihood optimizer; zero-event censoring fixture. Same author, shared legacy X detector; no independent replication or G-specific lifetime certification. No seconds conversion. Extrapolated mean at highOmega can be enormous and unsupported by follow-up; distinguish fit from observation. End with one report and stop; no new simulations authorized by this analysis.
