"""SciPy-verified public Boost.Math/GSL summary statistics C ABI."""
import ctypes as ct
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
import warnings

import numpy as np
from scipy import stats

from sci_testing.stats_summary import StatsSummaryComparisons


ROOT = Path(__file__).resolve().parent.parent


class NativeStatsSummaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        temporary = tempfile.TemporaryDirectory()
        cls.addClassCleanup(temporary.cleanup)
        library_path = Path(temporary.name) / "stats-summary.so"
        compiler = os.environ.get("CXX", shutil.which("clang++") or "c++")
        environment = dict(os.environ)
        # Some AppImage launchers prepend a private LD_LIBRARY_PATH that
        # interferes with system compiler and virtualenv runtime components.
        environment.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            compiler, "-std=c++17", "-O2", "-fPIC", "-shared",
            str(ROOT / "native/stats_summary.cpp"),
            "-lgsl", "-lgslcblas", "-o", str(library_path),
        ], check=True, env=environment)
        cls.native = StatsSummaryComparisons(ct.CDLL(str(library_path)))

    @staticmethod
    def case(name, x, p=(), y=(), namespace="descriptive"):
        return {"api": f"stats.{namespace}.{name}",
                "args": {"x": list(x), "y": list(y), "p": list(p)}}

    def assert_values(self, name, case, expected, *, libraries=None, tolerance=5e-10):
        actual = self.native.evaluate(case)
        self.assertTrue(actual, f"{name}: no native competitors")
        if libraries is not None:
            self.assertEqual(set(actual), set(libraries))
        expected = np.atleast_1d(np.asarray(expected, dtype=float))
        for library, result in actual.items():
            with self.subTest(name=name, library=library):
                self.assertNotIn("error", result, result)
                observed = np.asarray(result["values"])
                self.assertEqual(len(observed), len(expected))
                self.assertTrue(
                    np.allclose(observed, expected, rtol=tolerance,
                                atol=tolerance, equal_nan=True),
                    f"{name} {library} got {observed} vs SciPy {expected}",
                )

    def test_univariate_bias_ddof_and_array_shapes(self):
        # A skewed, non-symmetric sample catches the different sample-SD
        # normalization of Boost/GSL skewness and GSL excess kurtosis.
        sample = [-3., -1., 0., 2., 4., 10., 12.]
        for correction in (0, 1, 3):
            self.assert_values("variance", self.case("variance", sample, [correction]),
                               np.var(sample, ddof=correction))
            self.assert_values("std", self.case("std", sample, [correction]),
                               np.std(sample, ddof=correction))
            self.assert_values("sem", self.case("sem", sample, [correction]),
                               stats.sem(sample, ddof=correction))
            self.assert_values("variation", self.case("variation", sample, [correction]),
                               stats.variation(sample, ddof=correction))
            self.assert_values("zscore", self.case("zscore", sample, [correction]),
                               stats.zscore(sample, ddof=correction))
        self.assert_values("mean", self.case("mean", sample), np.mean(sample))
        self.assert_values("median", self.case("median", sample), np.median(sample))
        for bias in (0, 1):
            self.assert_values("skew", self.case("skew", sample, [bias]),
                               stats.skew(sample, bias=bool(bias)))
            for fisher in (0, 1):
                self.assert_values("kurtosis",
                                   self.case("kurtosis", sample, [fisher, bias]),
                                   stats.kurtosis(sample, fisher=bool(fisher), bias=bool(bias)))
        for correction in (0, 1):
            for bias in (0, 1):
                description = stats.describe(sample, ddof=correction, bias=bool(bias))
                expected = [description.nobs, *description.minmax, description.mean,
                            description.variance, description.skewness, description.kurtosis]
                self.assert_values("describe", self.case("describe", sample, [correction, bias]),
                                   expected)

    def test_quantiles_and_median_abs_deviation(self):
        for sample in ([1., 2., 4., 8., 9.], [-3., -1., 0., 2., 4., 10.],
                       [4., 4., 4., 4.], [1.]):
            self.assert_values("median", self.case("median", sample), np.median(sample))
            for probability in (0., .1, .25, .5, .75, .9, 1.):
                self.assert_values("quantile", self.case("quantile", sample, [probability]),
                                   np.quantile(sample, probability), libraries=["GSL"])
                self.assert_values("scoreatpercentile",
                                   self.case("scoreatpercentile", sample, [probability * 100]),
                                   stats.scoreatpercentile(sample, probability * 100),
                                   libraries=["GSL"])
            self.assert_values("iqr", self.case("iqr", sample),
                               stats.iqr(sample), libraries=["GSL"])
            for scale in (1., 1.4826):
                self.assert_values("median_absolute_deviation",
                                   self.case("median_abs_deviation", sample, [scale]),
                                   stats.median_abs_deviation(sample, scale=scale))

    def test_correlations_include_correct_tail_probabilities(self):
        samples = [
            ([1., 2., 3., 4.], [2., 1., 5., 6.]),
            ([0., 0., 1., 1., 2., 3.], [5., 2., 1., 2., 4., 3.]),
            ([1., 2., 3., 4.], [4., 3., 2., 1.]),
        ]
        alternative = ("two-sided", "less", "greater")
        for x, y in samples:
            for code, tail in enumerate(alternative):
                test = stats.pearsonr(x, y, alternative=tail)
                self.assert_values("pearsonr",
                                   self.case("pearsonr", x, [code], y, "association"),
                                   [test.statistic, test.pvalue])
                test = stats.spearmanr(x, y, alternative=tail)
                self.assert_values("spearmanr",
                                   self.case("spearmanr", x, [code], y, "association"),
                                   [test.statistic, test.pvalue], libraries=["GSL"])
                # Sci spearmanrho uses Pearson's beta null tail after rank transform.
                test = stats.spearmanrho(x, y, alternative=tail)
                self.assert_values("spearmanrho",
                                   self.case("spearmanrho", x, [code], y, "association"),
                                   [test.statistic, test.pvalue], libraries=["GSL"])
        x = [0., 0., 1., 1.]
        y = [1., 2., 3., 4.]
        test = stats.pointbiserialr(x, y)
        self.assert_values("pointbiserialr",
                           self.case("pointbiserialr", x, [], y, "association"),
                           [test.statistic, test.pvalue])
        test = stats.spearmanrho([1., 2.], [2., 1.])
        self.assert_values("spearmanrho two-element", self.case(
            "spearmanrho", [1., 2.], [0], [2., 1.], "association"),
            [test.statistic, test.pvalue], libraries=["GSL"])

    def test_actual_saved_stats_recipes(self):
        # Exercise real, retained fixture case shapes (including aliases).
        path = ROOT / "sci/tests/scipy/statistics_cases.jsonl"
        selected = {"mean", "variance", "std", "sem", "variation", "skew",
                    "kurtosis", "median", "iqr", "quantile", "scoreatpercentile",
                    "median_abs_deviation", "pearsonr", "spearmanr",
                    "spearmanrho", "pointbiserialr", "describe", "zscore"}
        checked = 0
        for raw in map(json.loads, path.read_text().splitlines()):
            name = raw["api"]
            if name not in selected:
                continue
            args = {key: raw[key] for key in ("p", "x", "y")}
            case = {"api": "stats." + name, "args": args}
            # Sci 0.3 recipes may include invalid/constant samples; only
            # shared finite domains are compared numerically.
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    expected = self._scipy(name, args)
            except (ValueError, ZeroDivisionError):
                continue
            if not np.all(np.isfinite(np.asarray(expected, dtype=float))):
                continue
            results = self.native.evaluate(case)
            for library, result in results.items():
                with self.subTest(saved=name, library=library, index=checked):
                    self.assertNotIn("error", result, result)
                    np.testing.assert_allclose(result["values"], np.atleast_1d(expected),
                                               rtol=1e-8, atol=1e-10)
            checked += 1
        self.assertGreaterEqual(checked, 100)

    @staticmethod
    def _scipy(name, a):
        x, y, p = np.asarray(a["x"]), np.asarray(a["y"]), a["p"]
        if name == "mean": return np.mean(x)
        if name == "variance": return np.var(x, ddof=int(p[0]))
        if name == "std": return np.std(x, ddof=int(p[0]))
        if name == "sem": return stats.sem(x, ddof=int(p[0]))
        if name == "variation": return stats.variation(x, ddof=int(p[0]))
        if name == "zscore": return stats.zscore(x, ddof=int(p[0]))
        if name == "skew": return stats.skew(x, bias=bool(p[0]))
        if name == "kurtosis": return stats.kurtosis(x, fisher=bool(p[0]), bias=bool(p[1]))
        if name == "median": return np.median(x)
        if name == "quantile": return stats.quantile(x, p[0])
        if name == "scoreatpercentile": return stats.scoreatpercentile(x, p[0])
        if name == "iqr": return stats.iqr(x)
        if name == "median_abs_deviation": return stats.median_abs_deviation(x, scale=p[0])
        if name == "describe":
            r = stats.describe(x, ddof=int(p[0]), bias=bool(p[1]))
            return [r.nobs, *r.minmax, r.mean, r.variance, r.skewness, r.kurtosis]
        alt = ("two-sided", "less", "greater")
        kwargs = {} if name == "pointbiserialr" else {"alternative": alt[int(p[0])]}
        result = getattr(stats, name)(x, y, **kwargs)
        return [result.statistic, result.pvalue]

    def test_legacy_variance_and_unsupported_shapes(self):
        result = self.native.evaluate(
            {"api": "stats.variance", "args": {"values": [1., 2., 3.]}})
        self.assertEqual(set(result), {"Boost.Math", "GSL"})
        self.assert_values("historical variance",
                           {"api": "stats.variance", "args": {"values": [1., 2., 3.]}}, 1.)
        self.assertEqual(self.native.evaluate(self.case("bootstrap", [1, 2, 3])), {})
        self.assertEqual(self.native.evaluate(self.case("spearmanrho", [1, 2], [0], [2] )), {})
        self.assertEqual(self.native.evaluate(self.case("median", ["bad"])), {})
        self.assertEqual(self.native.evaluate(self.case("mean", [math.inf])), {})


if __name__ == "__main__":
    unittest.main()
