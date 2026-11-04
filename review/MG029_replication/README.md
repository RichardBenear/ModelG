# MG029: MG9c Draft 1 design review

Reviewed revision `4466be3471056b6e2a9507fe6d98498b17de61ba`.

Run `python review/MG029_replication/design_checks.py CHECKOUT` with Python, NumPy and SciPy. CHECKOUT must contain the committed MG9 front records. Output is written beside the script as `design_checks.json`.

The script independently fits the existing tracks, checks their cell-count history, and constructs static counterexamples for the proposed bracket interpretation and opposite-sign motion rule. It imports no production modules and performs no new Model G integration. Hypothetical bracket outcomes are logical checks, not predictions or new experimental results.
