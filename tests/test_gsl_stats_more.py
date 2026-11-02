"""Exhaustive SciPy and GSL-derived statistical formula regressions."""

import csv
import ctypes as ct
import gzip
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
import warnings

import numpy as np
from scipy import stats

from sci_testing.case_csv import decode
from sci_testing.gsl_stats_more import StatsMoreGSLComparisons, SUPPORTED_NAMES
from sci_testing.registered_cases import reference


class StatsMoreGSLTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temp.cleanup)
        root = Path(__file__).resolve().parent.parent
        target = Path(cls.temp.name) / "gsl-stats-more.so"
        env = os.environ.copy()
        env.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2",
            "-shared", "-fPIC", str(root / "native/gsl_stats_more.cpp"),
            "-lgsl", "-lgslcblas", "-o", str(target),
        ], env=env, check=True)
        cls.native = StatsMoreGSLComparisons(ct.CDLL(str(target)))
        cls.generated = root / "cases/generated/stats"

    @staticmethod
    def case(name, *, x=None, y=None, p=None, g=None):
        return {"api": f"stats.{name}",
                "args": {"x": [] if x is None else x, "y": [] if y is None else y,
                         "p": [] if p is None else p, "g": [] if g is None else g,
                         "d": ""}}

    def result(self, case):
        returned = self.native.evaluate(case)
        self.assertEqual(set(returned), {"GSL"}, case)
        self.assertNotIn("error", returned["GSL"], (case, returned))
        return returned["GSL"]["values"]

    def test_all_saved_statistical_cases(self):
        counts = {}
        max_error = {}
        for name in sorted(SUPPORTED_NAMES):
            path = self.generated / (name + ".csv.gz")
            if not path.exists():
                path = self.generated / (name + ".csv")
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "rt", newline="") as handle:
                rows = list(csv.DictReader(handle))
            counts[name] = len(rows)
            max_error[name] = 0.0
            for record in rows:
                case = {"api": f"stats.{name}",
                        "args": {k: decode(v) for k, v in record.items()
                                 if k not in ("id", "group")}}
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    expected = reference(case)
                observed = self.result(case)
                with self.subTest(case=record["id"]):
                    self.assertEqual(len(observed), len(expected))
                    for actual, target in zip(observed, expected):
                        self.assertTrue(math.isfinite(float(target)), record["id"])
                        difference = abs(actual - target)
                        max_error[name] = max(
                            max_error[name], difference/max(abs(target), 1e-14))
                        self.assertTrue(math.isclose(actual, target,
                                                  rel_tol=3e-10, abs_tol=4e-12),
                                        (record["id"], actual, target, difference))
        expected_counts = {
            "gmean": 2003, "hmean": 2003, "circmean": 2003,
            "circstd": 2006, "circvar": 2003, "jarque_bera": 2002,
            "chisquare": 2006, "ttest_ind_from_stats": 2002,
            "f_oneway": 2003, "bartlett": 2003,
        }
        self.assertEqual(counts, expected_counts)
        self.assertEqual(sum(counts.values()), sum(expected_counts.values()))
        print("saved-case maximum relative errors:", max_error)

    def test_analytic_and_distribution_tails(self):
        self.assertAlmostEqual(self.result(self.case("gmean", x=[1, 4, 16]))[0], 4.)
        self.assertAlmostEqual(self.result(self.case("hmean", x=[1, 2, 4]))[0],
                               3/(1+0.5+0.25))
        self.assertAlmostEqual(self.result(self.case("circmean", x=[0, math.pi/2],
                                                p=[0, 2*math.pi]))[0], math.pi/4)
        for name in ("circstd", "circvar"):
            p = [0, 2*math.pi, 1] if name == "circstd" else [0, 2*math.pi]
            expected = getattr(stats, name)([0, .3, .7, 1.],
                                             high=2*math.pi, low=0,
                                             **({"normalize": True} if name == "circstd" else {}))
            self.assertAlmostEqual(self.result(self.case(name, x=[0, .3, .7, 1.], p=p))[0],
                                   expected, delta=2e-14)
        expected = stats.jarque_bera([-1, 0, 1, 2, 9, 10, 13])
        actual = self.result(self.case("jarque_bera", x=[-1, 0, 1, 2, 9, 10, 13]))
        np.testing.assert_allclose(actual, [expected.statistic, expected.pvalue],
                                   rtol=2e-14, atol=2e-14)
        chi = stats.chisquare([5, 15, 20, 30], [17.5]*4, ddof=1)
        np.testing.assert_allclose(self.result(self.case("chisquare", x=[5, 15, 20, 30],
                                                         y=[17.5]*4, p=[1])),
                                   [chi.statistic, chi.pvalue], rtol=2e-14)
        for equal in (False, True):
            for alt, word in enumerate(("two-sided", "less", "greater")):
                inputs = [1.5, 0.6, 8, 2.3, 1.1, 7, int(equal), alt]
                expected = stats.ttest_ind_from_stats(1.5, .6, 8, 2.3, 1.1, 7,
                                                     equal_var=equal, alternative=word)
                observed = self.result(self.case("ttest_ind_from_stats", p=inputs))
                np.testing.assert_allclose(observed[:2],
                                           [expected.statistic, expected.pvalue],
                                           rtol=2e-14, atol=2e-14)
        groups = [[1., 2., 3.], [4., 5., 6., 9.], [0., 1., 5., 8.]]
        for name in ("f_oneway", "bartlett"):
            expected = getattr(stats, name)(*groups)
            observed = self.result(self.case(name, g=groups))
            np.testing.assert_allclose(observed, [expected.statistic, expected.pvalue],
                                       rtol=1e-12, atol=1e-13)

    def test_invalid_input_and_abi_contract(self):
        library = self.native.library
        for name in SUPPORTED_NAMES:
            encoded = name.encode()
            self.assertEqual(library.gsl_stats_more_supported(0, encoded), 0)
            self.assertEqual(library.gsl_stats_more_supported(1, encoded),
                             3 if name == "ttest_ind_from_stats" else
                             2 if name in ("jarque_bera", "chisquare", "f_oneway", "bartlett")
                             else 1)
        self.assertEqual(self.native.evaluate(self.case("unsupported", x=[1, 2])), {})
        self.assertEqual(self.native.evaluate(self.case("gmean", x=[1, [2, 3]])), {})
        self.assertEqual(self.native.evaluate(self.case("gmean", x=[1, float("nan")])), {})
        self.assertEqual(self.native.evaluate(self.case("f_oneway", g=[[1, 2], [3, [4]]])), {})
        self.assertEqual(self.native.evaluate(self.case("circmean", x=[0, 1], p=[])), {})
        invalid = [
            self.case("gmean", x=[1, -1]),
            self.case("hmean", x=[1, -1]),
            self.case("circstd", x=[0, 1], p=[1, 1, 0]),
            self.case("circstd", x=[0, 1], p=[0, 2*math.pi, 5]),
            self.case("jarque_bera", x=[1, 1, 1]),
            self.case("chisquare", x=[1, 2, 3], y=[2, 2, 1], p=[0]),
            self.case("chisquare", x=[1, 2], y=[1, 2], p=[1]),
            self.case("ttest_ind_from_stats", p=[1, -2, 10, 2, 1, 10, 1, 0]),
            self.case("ttest_ind_from_stats", p=[1, 2, 10.5, 2, 1, 10, 1, 0]),
            self.case("f_oneway", g=[[1], [2, 3]]),
            self.case("bartlett", g=[[1, 1], [2, 3]]),
        ]
        for case in invalid:
            self.assertIn("error", self.native.evaluate(case)["GSL"], case)
        x = (ct.c_double * 3)(1., 2., 3.)
        output = (ct.c_double * 3)(91., 92., 93.)
        self.assertNotEqual(library.gsl_stats_more_evaluate(
            1, b"jarque_bera", x, 3, None, 0, None, 0, None, 0,
            None, 0, output, 1), 0)
        self.assertEqual(list(output), [91., 92., 93.])
        self.assertNotEqual(library.gsl_stats_more_evaluate(
            0, b"jarque_bera", x, 3, None, 0, None, 0, None, 0,
            None, 0, output, 3), 0)
        self.assertEqual(list(output), [91., 92., 93.])


if __name__ == "__main__":
    unittest.main()
