"""Logical counterexamples for MG9b Draft 2; no trajectories or coverage simulation."""
import json
from pathlib import Path

# Nine agreements and one incomplete continuation satisfy the written >=9 rule.
agreements=[True]*9+[False]
complete=[True]*9+[False]

# A finite mixture illustrates why passing pooled coverage does not validate
# every branch. These are assigned probabilities, NOT measured calibration.
weights=[.95,.05]
coverage=[.96,.50]
per_side_error=[.02,.25]
out={
    'revision':'9504aba071f8f7de88bc57786da35cd41bd602ff',
    'incomplete_control_case':{'agreements':sum(agreements),'completed':sum(complete),
        'passes_written_9_of_10_rule':sum(agreements)>=9,
        'passes_completion_required_rule':all(complete) and sum(agreements)>=9},
    'branch_mixture_example':{'branch_weights':weights,'branch_coverage':coverage,
        'branch_directional_error_each_side':per_side_error,
        'pooled_coverage':sum(w*c for w,c in zip(weights,coverage)),
        'pooled_directional_error_each_side':sum(w*c for w,c in zip(weights,per_side_error)),
        'expected_branch_counts_at_R2000':[2000*w for w in weights],
        'note':'Passes pooled 0.93 coverage and 0.035 per-side criteria despite poor rare-branch behavior.'}}
assert out['incomplete_control_case']['passes_written_9_of_10_rule']
assert out['branch_mixture_example']['pooled_coverage']>=.93
assert out['branch_mixture_example']['pooled_directional_error_each_side']<=.035
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
