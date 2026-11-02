"""Saved-case and mathematical checks for independent mpmath comparisons."""

import csv
import gzip
import math
import unittest
import warnings
from pathlib import Path

import mpmath as mp

from sci_testing.case_csv import decode
from sci_testing.mpmath_general import MpmathGeneralComparisons, evaluate
from sci_testing.report import high_precision

ROOT = Path(__file__).resolve().parents[1] / "cases" / "generated"
GROUPS = (
    "complex", "constants", "linalg", "integrate", "differentiate", "stats",
    "optimize", "spatial", "signal", "filter", "cluster.vq", "fft",
)
INDICES = {0, 1, 2, 4, 20, 100, 500, 1000, 1500}
COMPLEX_APIS = {
    "abs", "acos", "acosh", "add", "arg", "asin", "asinh", "atan", "atanh",
    "cbrt", "conjugate", "cos", "cosh", "divide", "exp", "from_polar",
    "imag", "ln", "multiply", "negate", "norm_sqr", "powc", "powf",
    "powi", "powu", "real", "reciprocal", "sin", "sinh", "sqrt",
    "subtract", "tan", "tanh", "to_polar",
}


def saved_cases():
    for group in GROUPS:
        for path in sorted((ROOT / group).glob("*.csv*")):
            name = path.name.split(".csv")[0]
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "rt", newline="") as handle:
                for index, row in enumerate(csv.DictReader(handle)):
                    if index not in INDICES:
                        continue
                    yield {
                        "api": f"{group}.{name}", "id": row["id"],
                        "args": {key: decode(value) for key, value in row.items()
                                 if key not in ("id", "group")},
                    }


class MpmathGeneralTests(unittest.TestCase):
    def setUp(self):
        self.adapter = MpmathGeneralComparisons()

    def test_saved_case_dimensions_and_accuracy(self):
        """All published comparisons agree with full saved reference shapes."""
        supported = set()
        measured = 0
        for case in saved_cases():
            observed = self.adapter.evaluate(case)
            if not observed:
                continue
            with self.subTest(case=case["id"]):
                self.assertEqual(set(observed), {"mpmath"})
                self.assertEqual(set(observed["mpmath"]), {"values"})
                values = observed["mpmath"]["values"]
                self.assertTrue(values)
                self.assertTrue(all(isinstance(v, float) and math.isfinite(v)
                                    for v in values))
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    expected = high_precision(case)
                self.assertEqual(len(values), len(expected))
                for actual, reference in zip(values, expected):
                    self.assertTrue(math.isclose(actual, float(reference),
                                                  rel_tol=2e-8, abs_tol=2e-10),
                                    (actual, reference))
            supported.add(case["api"])
            measured += 1
        self.assertEqual({api.removeprefix("complex.") for api in supported
                          if api.startswith("complex.")}, COMPLEX_APIS)
        # All groups contribute independently evaluated operations.
        for group in GROUPS:
            self.assertTrue(any(api.startswith(group + ".") for api in supported),
                            group)
        self.assertGreaterEqual(len(supported), 156)
        self.assertGreaterEqual(measured, 1397)

    def test_shape_and_branch_regressions(self):
        with mp.workdps(80):
            self.assertEqual(evaluate({"api": "linalg.logm", "args": {
                "raw": {"api": "logm", "args": [[[2., 1.], [0., 3.]]]}}})[
                "mpmath"]["values"], [
                    float(mp.log(2)), float(mp.log(mp.mpf(3)/2)),
                    0., float(mp.log(3)), 0., 0., 0., 0.,
                ])
            self.assertEqual(evaluate({"api": "linalg.cholesky", "args": {
                "raw": {"api": "cholesky", "args": [
                    [[4., 2.], [2., 5.]], True]}}})[
                "mpmath"]["values"], [2., 0., 1., 2.])
            self.assertEqual(evaluate({"api": "linalg.svdvals", "args": {
                "raw": {"api": "svdvals", "args": [[[3., 0.], [0., 2.]]]}}})[
                "mpmath"]["values"], [3., 2.])
            below = evaluate({"api": "complex.sqrt",
                              "args": {"x": -4., "y": -0.}})["mpmath"]["values"]
            above = evaluate({"api": "complex.sqrt",
                              "args": {"x": -4., "y": 0.}})["mpmath"]["values"]
            self.assertEqual(below, [0., -2.])
            self.assertEqual(above, [0., 2.])
            self.assertEqual(evaluate({"api": "integrate.romb", "args": {
                "explicit": True, "values": [1., 2., 3., 4., 5., 6., 7., 8., 9.]}})[
                "mpmath"]["values"], [10.])
            self.assertAlmostEqual(evaluate({"api": "stats.ppf", "args": {
                "d": "Uniform", "p": [-2, 3], "x": [0.001]}})[
                "mpmath"]["values"][0], -1.995)
            self.assertAlmostEqual(evaluate({"api": "stats.ppf", "args": {
                "d": "Pareto", "p": [1, 3], "x": [0.001]}})[
                "mpmath"]["values"][0], float((1-mp.mpf('0.001'))**(-mp.mpf(1)/3)))
            self.assertAlmostEqual(evaluate({"api": "stats.cdf", "args": {
                "d": "Uniform", "p": [-2, 3], "x": [1.]}})[
                "mpmath"]["values"][0], .6)
            self.assertAlmostEqual(evaluate({"api": "stats.cdf", "args": {
                "d": "Pareto", "p": [1, 3], "x": [2.]}})[
                "mpmath"]["values"][0], .875)

    def test_unavailable_malformed_and_nonfinite_are_not_comparisons(self):
        unsupported = (
            {"api": "stats.unknown", "args": {"x": [1.]}},
            {"api": "linalg.det", "args": {"a": [[1., 2.], [3.]]}},
            {"api": "linalg.solve", "args": {"a": [[1., 2.], [2., 4.]],
                                            "b": [1., 2.]}},
            {"api": "linalg.inv", "args": {"raw": {"api": "inv",
                                                    "args": [[[1., 2.], [2., 4.]]]}}},
            {"api": "linalg.inv", "args": {"raw": {"api": "eigvalsh",
                                                    "args": [[[1.]]]}}},
            {"api": "stats.variance", "args": {"values": [math.inf, 1.]}},
            {"api": "complex.exp", "args": {"x": math.nan, "y": 1.}},
            {"api": "constants.lambda2nu", "args": {"x": 0.}},
            {"api": "stats.pdf", "args": {"d": "Unknown", "p": [], "x": [1.]}},
            {"api": "integrate.fixed_quad", "args": {"a": 0.}},
        )
        for case in unsupported:
            with self.subTest(case=case):
                self.assertEqual(evaluate(case), {})


if __name__ == "__main__":
    unittest.main()
