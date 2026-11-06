"""MG041: first MG10 qualification check, not scientific production.

At P17 q=p, A is decoupled: r'=-r/(C*T_ref). Check the specified
synchronous Euler updates against the exact ODE as required by MG040.
No active fields, particle outcomes, or recovery decisions are generated.
"""
from pathlib import Path
import hashlib
import json
import platform
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CAPACITIES = (1, 10, 100, 1000, 10000, 100000, 1000000)
STEPS_PER_UNIT = (1500, 3000, 6000)
T_REF = 300
TOLERANCE = 1e-8


def euler_checkpoints(capacities, steps_per_unit, horizon=300):
    """Literal float64 start-of-step Euler recurrence at integer times."""
    rates = 1.0 / (np.asarray(capacities, dtype=float) * T_REF)
    state = np.ones(len(capacities), dtype=float)
    saved = np.empty((horizon + 1, len(capacities)), dtype=float)
    saved[0] = state
    dt = 1.0 / steps_per_unit
    for second in range(1, horizon + 1):
        for _ in range(steps_per_unit):
            state -= dt * rates * state
        saved[second] = state
    return saved


def run(output):
    started = time.perf_counter()
    output.mkdir(parents=True, exist_ok=True)
    times = np.arange(T_REF + 1, dtype=float)
    capacities = np.array(CAPACITIES)
    exact = np.exp(-times[:, None] / (capacities[None, :] * T_REF))
    fields = {"times": times, "capacities": capacities, "exact": exact}
    rows = []
    for method, steps in enumerate(STEPS_PER_UNIT):
        actual = euler_checkpoints(capacities, steps)
        # A second oracle is the exact solution of the discrete recurrence.
        discrete = np.exp(times[:, None] * steps * np.log1p(
            -1.0 / (steps * capacities[None, :] * T_REF)))
        fields[f"E{method}"] = actual
        fields[f"E{method}_discrete_oracle"] = discrete
        for i, capacity in enumerate(capacities):
            error = np.abs(actual[:, i] - exact[:, i])
            k = int(np.argmax(error))
            rows.append({"method": f"E{method}", "capacity": int(capacity),
                         "steps_per_unit": steps,
                         "max_scaled_absolute_error": float(error[k]),
                         "worst_checkpoint": int(times[k]),
                         "discrete_oracle_max_error": float(np.max(
                             np.abs(actual[:, i] - discrete[:, i]))),
                         "pass": bool(np.max(error) <= TOLERANCE)})
    raw = output / "availability_checkpoints.npz"
    np.savez_compressed(raw, **fields)
    elapsed = time.perf_counter() - started
    report = {
        "stage": "MG041", "scope": "DEVELOPMENT_QUALIFICATION_ONLY",
        "verdict": "PASS" if all(r["pass"] for r in rows) else "QUALIFICATION FAILED",
        "criterion": "MG040 run plan section 3, mandatory check 4",
        "absolute_tolerance": TOLERANCE, "T_ref": T_REF,
        "production_dispatched": False, "external_review": "PI: No review needed. Go ahead.",
        "elapsed_worker_seconds": elapsed, "workers": 1,
        "qualification_ceiling_seconds": 900,
        "budget_exceeded": elapsed > 900,
        "platform": platform.platform(), "python": sys.version,
        "numpy": np.__version__, "raw_sha256": hashlib.sha256(raw.read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "run_plan_sha256": hashlib.sha256((ROOT / "spec/MG10_run_plan.md").read_bytes()).hexdigest(),
        "rows": rows,
        "unattempted": "Other qualification checks, full implementation, timing matrix and production.",
    }
    (output / "availability_qualification.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


if __name__ == "__main__":
    report = run(ROOT / "results/MG10/qualification")
    print(json.dumps({k: report[k] for k in (
        "verdict", "elapsed_worker_seconds", "production_dispatched")}, indent=2))
    for row in report["rows"]:
        if not row["pass"]:
            print(row)
    sys.exit(0 if report["verdict"] == "PASS" else 2)
