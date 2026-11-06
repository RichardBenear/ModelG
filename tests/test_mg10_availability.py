"""Meaningful exact-case tests for the MG041 early qualification check."""
from pathlib import Path
import sys
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from qualify_mg10_availability import euler_checkpoints


class AvailabilityTests(unittest.TestCase):
    def test_saved_times_and_discrete_recurrence(self):
        actual = euler_checkpoints([1, 10], 1500, horizon=2)
        self.assertEqual(actual.shape, (3, 2))
        for k in range(3):
            expected = np.exp(k * 1500 * np.log1p(-1 / (1500 * np.array([1, 10]) * 300)))
            np.testing.assert_allclose(actual[k], expected, rtol=0, atol=2e-13)

    def test_refinement_reduces_ode_error(self):
        errors = [abs(euler_checkpoints([1], n, horizon=2)[-1, 0] - np.exp(-2 / 300))
                  for n in (1500, 3000, 6000)]
        self.assertGreater(errors[0], errors[1])
        self.assertGreater(errors[1], errors[2])
        self.assertAlmostEqual(errors[0] / errors[1], 2, places=3)


if __name__ == "__main__":
    unittest.main()
