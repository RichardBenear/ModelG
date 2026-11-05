# MG034 Draft 1 design checks

Reviewed revision: `98d63b8659260879d7c674956610032707415142`.

```text
python review/MG034_replication/design_checks.py /path/to/ModelG
```

Requires NumPy and SciPy. Writes `design_checks.json` beside the script. No trajectory integrations are performed. A mocked confirmation call records the sampling times used when t_end is set to 2000. An illustrative Poisson table demonstrates that a negative fitted linear slope trend can coexist with nonmonotone per-b slopes. The synthetic counts are not ModelG results or a calibration/power estimate. Additional outputs record degrees of freedom, Monte Carlo precision, and a duplicate replacement example.

Full review: `review/MG034_chatgpt.md`.
