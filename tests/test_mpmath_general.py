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
from sci_testing.registered_cases import reference
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

    def test_five_saved_finite_endpoint_density_regressions(self):
        """Saved PDF inputs whose finite limits used to be incorrectly zero."""
        expected = {
            "stats.pdf/recipe-2484": ("Weibull", [1, 2], 0, 0.5),
            "stats.pdf/recipe-2447": ("Gamma", [1, 2], 0, 0.5),
            "stats.pdf/recipe-2374": ("Beta", [2, 1], 1, 2.),
            "stats.pdf/random-0689": ("Beta", [2, 1], 1, 2.),
            "stats.pdf/recipe-2336": ("Beta", [1, 3], 0, 3.),
        }
        seen = set()
        with gzip.open(ROOT / "stats" / "pdf.csv.gz", "rt", newline="") as handle:
            for row in csv.DictReader(handle):
                case_id = row["id"]
                if case_id not in expected:
                    continue
                with self.subTest(case=case_id):
                    args = {key: decode(value) for key, value in row.items()
                            if key not in ("id", "group")}
                    kind, parameters, x, target = expected[case_id]
                    self.assertEqual((args["d"], args["p"], args["x"]),
                                     (kind, parameters, [x]))
                    result = self.adapter.evaluate({"api": "stats.pdf", "args": args})
                    self.assertEqual(result, {"mpmath": {"values": [target]}})
                    self.assertEqual(reference({"api": "stats.pdf", "args": args}),
                                     [target])
                seen.add(case_id)
        self.assertEqual(seen, set(expected))

    def test_finite_support_density_limits(self):
        """Independent analytic limits including additional gamma-family cases."""
        finite = (
            ("Beta", [1, 3], 0, 3.),
            ("Beta", [2, 1], 1, 2.),
            ("Beta", [2, 3], 0, 0.),
            ("Beta", [2, 3], 1, 0.),
            ("Gamma", [1, 2], 0, .5),
            ("Gamma", [2, 2], 0, 0.),
            ("Erlang", [1, 2], 0, .5),
            ("ChiSquared", [2], 0, .5),
            ("Chi", [1], 0, math.sqrt(2/math.pi)),
            ("Chi", [3], 0, 0.),
            ("F", [2, 7], 0, 1.),
            ("F", [3, 7], 0, 0.),
            ("Weibull", [1, 2], 0, .5),
            ("Weibull", [1.5, 2], 0, 0.),
            ("Triangular", [-2, 3, -2], -2, .4),
            ("Triangular", [-2, 3, 3], 3, .4),
        )
        for kind, params, x, target in finite:
            with self.subTest(distribution=kind, params=params, x=x):
                case = {"api": "stats.pdf", "args": {"d": kind,
                        "p": params, "x": [x]}}
                actual = evaluate(case)
                self.assertIn("mpmath", actual)
                self.assertAlmostEqual(actual["mpmath"]["values"][0], target)

        for kind, params, x in (("Beta", [.5, 2], 0),
                                ("Beta", [2, .5], 1),
                                ("Gamma", [.5, 2], 0),
                                ("Weibull", [.5, 2], 0),
                                ("ChiSquared", [1], 0),
                                ("F", [1, 7], 0)):
            with self.subTest(infinite_density=kind, params=params):
                self.assertEqual(evaluate({"api": "stats.pdf", "args": {
                    "d": kind, "p": params, "x": [x]}}), {})

    def test_saved_distribution_quantile_endpoints(self):
        """Every finite saved q=0/1 endpoint is an independent comparison."""
        seen = finite = 0
        with gzip.open(ROOT / "stats" / "ppf.csv.gz", "rt", newline="") as handle:
            for row in csv.DictReader(handle):
                args = {key: decode(value) for key, value in row.items()
                        if key not in ("id", "group")}
                if args["x"] not in ([0], [1]):
                    continue
                seen += 1
                case = {"api": "stats.ppf", "args": args}
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    target = reference(case)[0]
                observed = evaluate(case)
                with self.subTest(case=row["id"]):
                    if not math.isfinite(float(target)):
                        self.assertEqual(observed, {})
                    else:
                        self.assertIn("mpmath", observed)
                        self.assertEqual(observed["mpmath"]["values"],
                                         [float(target)])
                        finite += 1
        self.assertEqual(seen, 51)
        self.assertEqual(finite, 51)

    def test_saved_inverse_quantile_precision_regressions(self):
        """Very small quantiles must not hit the absolute bisection floor."""
        ids = {
            "stats.ppf/recipe-5115",  # InverseGamma tiny scale, q=1e-100
            "stats.ppf/linspace-0299", "stats.ppf/linspace-0599",
            "stats.ppf/linspace-0899", "stats.ppf/random-0299",
            "stats.ppf/random-0599", "stats.ppf/random-0899",
            "stats.ppf/recipe-5065",  # Chi(1), q=1e-100
            "stats.ppf/recipe-0969",  # StudentT median is exactly zero
            "stats.ppf/recipe-5088",  # InverseGamma heavy upper tail
            "stats.ppf/recipe-5113", "stats.ppf/recipe-5114",
        }
        seen = set()
        with gzip.open(ROOT / "stats" / "ppf.csv.gz", "rt", newline="") as handle:
            for row in csv.DictReader(handle):
                if row["id"] not in ids:
                    continue
                args = {key: decode(value) for key, value in row.items()
                        if key not in ("id", "group")}
                case = {"api": "stats.ppf", "args": args}
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    target = float(reference(case)[0])
                with self.subTest(case=row["id"]):
                    result = evaluate(case)
                    self.assertIn("mpmath", result)
                    actual = result["mpmath"]["values"]
                    self.assertEqual(len(actual), 1)
                    if row["id"] == "stats.ppf/recipe-0969":
                        self.assertEqual(actual[0], 0.)
                    else:
                        # Deliberately no large absolute tolerance: a
                        # 1e-28 result is catastrophically wrong for 1e-203.
                        self.assertTrue(math.isclose(actual[0], target,
                                        rel_tol=2e-8, abs_tol=8*math.ulp(target)),
                                        (actual[0], target))
                    if args["d"] == "InverseGamma" and args["p"] == [1, 1e-200]:
                        with mp.workdps(90):
                            analytic = float(-mp.mpf(1e-200)/mp.log(
                                mp.mpf(args["x"][0])))
                        self.assertEqual(actual[0], analytic)
                    if args["d"] == "Chi" and args["p"] == [1]:
                        with mp.workdps(90):
                            analytic = float(mp.sqrt(2)*mp.erfinv(
                                mp.mpf(args["x"][0])))
                        self.assertEqual(actual[0], analytic)
                seen.add(row["id"])
        self.assertEqual(seen, ids)

    def test_every_saved_finite_distribution_quantile(self):
        """Audit all saved PPF inputs, including extreme tails, against Sci.

        Cached identical (distribution, parameters, probability) triples
        avoid recomputation but every original saved case ID is checked.
        """
        cache = {}
        total = compared = unavailable = 0
        with gzip.open(ROOT / "stats" / "ppf.csv.gz", "rt", newline="") as handle:
            for row in csv.DictReader(handle):
                args = {key: decode(value) for key, value in row.items()
                        if key not in ("id", "group")}
                key = (args["d"], tuple(args["p"]), tuple(args["x"]))
                if key not in cache:
                    case = {"api": "stats.ppf", "args": args}
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        expected = reference(case)
                    cache[key] = (evaluate(case), expected)
                observed, expected = cache[key]
                total += 1
                with self.subTest(case=row["id"]):
                    self.assertEqual(len(expected), 1)
                    target = float(expected[0])
                    self.assertTrue(math.isfinite(target))
                    if not observed:
                        # mpmath's regularized gamma does not reliably
                        # converge at Chi(100000) near its central mass.
                        self.assertEqual(args["d"], "Chi")
                        self.assertEqual(args["p"], [100000])
                        unavailable += 1
                        continue
                    actual = observed["mpmath"]["values"]
                    self.assertEqual(len(actual), 1)
                    self.assertTrue(math.isfinite(actual[0]))
                    if target == 0:
                        self.assertLessEqual(abs(actual[0]), 1e-12)
                    else:
                        self.assertTrue(math.isclose(
                            actual[0], target, rel_tol=2e-8,
                            abs_tol=8*math.ulp(target)), (actual[0], target))
                    compared += 1
        self.assertEqual(total, 2300)
        self.assertGreaterEqual(compared, 2286)
        self.assertLessEqual(unavailable, 14)

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
