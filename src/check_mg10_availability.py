"""Separate MG041 recomputation: decimal arithmetic, no qualification imports."""
from decimal import Decimal, localcontext
from pathlib import Path
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def recompute(directory):
    record = json.loads((directory / "availability_qualification.json").read_text(encoding="utf-8"))
    raw = directory / "availability_checkpoints.npz"
    assert hashlib.sha256(raw.read_bytes()).hexdigest() == record["raw_sha256"]
    results = []
    with np.load(raw) as data, localcontext() as context:
        context.prec = 60
        assert np.array_equal(data["times"], np.arange(301))
        assert np.array_equal(data["capacities"], [1, 10, 100, 1000, 10000, 100000, 1000000])
        for row in record["rows"]:
            c = Decimal(row["capacity"])
            n = row["steps_per_unit"]
            column = list(data["capacities"]).index(int(c))
            observed = data[row["method"]][:, column]
            assert np.all(np.isfinite(observed)) and observed[0] == 1
            exact = np.array([float((-Decimal(t) / (c * 300)).exp()) for t in range(301)])
            # Integer-power discrete oracle, independently of log1p formulation.
            factor = Decimal(1) - Decimal(1) / (Decimal(n) * c * 300)
            discrete = np.array([float(factor ** (t * n)) for t in range(301)])
            error = float(np.max(np.abs(observed - exact)))
            passed = error <= 1e-8
            assert passed == row["pass"]
            assert abs(error - row["max_scaled_absolute_error"]) < 2e-15
            assert np.max(np.abs(observed - discrete)) < 1e-9
            results.append({"method": row["method"], "capacity": int(c),
                            "pass": passed, "decimal_exact_max_error": error,
                            "decimal_discrete_max_error": float(np.max(np.abs(observed - discrete)))})
    result = {"recomputation": "AGREES", "rows_checked": len(results),
              "failures": sum(not r["pass"] for r in results),
              "shared_dependencies": "Raw float64 checkpoints and NumPy array reductions; no qualification-code imports. Decimal oracles are separate.",
              "review_status": "Same-author recomputation, not external replication.", "rows": results}
    (directory / "availability_recomputed.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


if __name__ == "__main__":
    result = recompute(ROOT / "results/MG10/qualification")
    print({k: result[k] for k in ("recomputation", "rows_checked", "failures")})
