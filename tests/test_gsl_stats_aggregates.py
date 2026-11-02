"""GSL-derived aggregate checks against every saved SciPy comparison case."""

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
from sci_testing.gsl_stats_aggregates import (
    SUPPORTED_NAMES, StatsAggregatesGSLComparisons,
)
from sci_testing.registered_cases import reference


class StatsAggregatesGSLTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tmp.cleanup)
        root = Path(__file__).resolve().parent.parent
        target = Path(cls.tmp.name) / "gsl-stats-aggregates.so"
        environment = os.environ.copy()
        environment.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2", "-shared", "-fPIC",
            str(root / "native/gsl_stats_aggregates.cpp"), "-lgsl", "-lgslcblas",
            "-o", str(target),
        ], env=environment, check=True)
        cls.native = StatsAggregatesGSLComparisons(ct.CDLL(str(target)))
        cls.generated = root / "cases/generated/stats"

    @staticmethod
    def case(name, values, parameter):
        return {"api": "stats." + name,
                "args": {"x": list(values), "y": [], "p": [parameter],
                         "g": [], "d": ""}}

    def checked(self, case):
        observed = self.native.evaluate(case)
        self.assertEqual(set(observed), {"GSL"}, case)
        self.assertNotIn("error", observed["GSL"], (case, observed))
        return observed["GSL"]["values"]

    def test_every_saved_scipy_case(self):
        counts = {}
        maxima = {}
        for name in sorted(SUPPORTED_NAMES):
            path = self.generated / (name + ".csv.gz")
            with gzip.open(path, "rt", newline="") as handle:
                records = list(csv.DictReader(handle))
            counts[name] = len(records)
            maxima[name] = 0.0
            for row in records:
                case = {"api": "stats." + name,
                        "args": {k: decode(value) for k, value in row.items()
                                 if k not in ("id", "group")}}
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    expected = reference(case)
                observed = self.checked(case)
                with self.subTest(case=row["id"]):
                    self.assertEqual(len(observed), len(expected))
                    for actual, target in zip(observed, expected):
                        self.assertTrue(math.isfinite(float(target)),row["id"])
                        difference = abs(actual - target)
                        maxima[name] = max(maxima[name],difference/max(abs(target),1e-14))
                        self.assertTrue(math.isclose(actual,target,rel_tol=5e-11,
                                                  abs_tol=5e-12),
                                        (row["id"],actual,target,difference))
        self.assertEqual(counts, {
            "gstd": 2006, "gzscore": 2004, "pmean": 2018, "moment": 2024,
        })
        self.assertEqual(sum(counts.values()),8052)
        print("saved-case relative discrepancy:", maxima)

    def test_independent_power_and_moment_identities(self):
        x = [0.1, 0.5, 2., 7.]
        for power in (-8., -3., -1., -.5, 0., .25, .5, 1., 2., 3., 8., 16.):
            with self.subTest(power=power):
                actual = self.checked(self.case("pmean",x,power))[0]
                expected = float(stats.pmean(x,power))
                self.assertAlmostEqual(actual,expected,delta=4e-13*max(1.,expected))
        for power in (-3., -1., 0., 1., 2., 8.):
            with self.subTest(power=power, sample="includes zero"):
                value = self.checked(self.case("pmean",[0., 1., 4.],power))[0]
                self.assertAlmostEqual(value,float(stats.pmean([0., 1., 4.],power)),
                                       delta=5e-14)
        for x in ([1, 2, 4, 8, 9],[-3, 0, 2, 6, 12], [4, 4, 4, 4], [1]):
            for order in (0, 1, 2, 3, 4, 8, 12):
                with self.subTest(x=x,order=order):
                    actual = self.checked(self.case("moment",x,order))[0]
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", RuntimeWarning)
                        expected = float(stats.moment(x,order))
                    self.assertAlmostEqual(actual,expected,
                                           delta=1e-10*max(1.,abs(expected)))
        for ddof in (0, 1):
            actual = self.checked(self.case("gstd",[1,2,4,8],ddof))
            self.assertAlmostEqual(actual[0], stats.gstd([1,2,4,8],ddof=ddof))
            actual = self.checked(self.case("gzscore",[1,2,4,8],ddof))
            expected = stats.gzscore([1,2,4,8],ddof=ddof)
            np.testing.assert_allclose(actual,expected,rtol=3e-13,atol=3e-14)

    def test_rejected_domains_and_abi_shape(self):
        for name in SUPPORTED_NAMES:
            encoded = name.encode("ascii")
            self.assertEqual(self.native.library.gsl_stats_aggregates_supported(0,encoded),0)
            self.assertEqual(self.native.library.gsl_stats_aggregates_supported(1,encoded),1)
            self.assertEqual(self.native.library.gsl_stats_aggregates_output_count(
                encoded,4),4 if name=="gzscore" else 1)
        cases = (
            self.case("gstd",[1,-1,2],0), self.case("gstd",[1,2],2),
            self.case("gzscore",[2,2,2],0), self.case("gzscore",[1,2,3],4),
            self.case("moment",[1,2],-1), self.case("moment",[1,2],1.5),
            self.case("moment",[1,2],33),
            self.case("pmean",[-1,2,3],3), self.case("pmean",[1,2],40),
        )
        for case in cases:
            self.assertIn("error",self.native.evaluate(case)["GSL"],case)
        self.assertEqual(self.native.evaluate(self.case("unknown",[1,2],0)),{})
        self.assertEqual(self.native.evaluate(self.case("pmean",[1,math.nan],2)),{})
        self.assertEqual(self.native.evaluate(self.case("pmean",[1,[2]],2)),{})
        self.assertEqual(self.native.evaluate({"api":"stats.moment",
                                               "args":{"x":[1,2],"p":[]}}),{})
        self.assertEqual(self.native.evaluate({"api":"stats.moment",
                                               "args":{"x":[1,2],"p":[2],"g":[[1,2]]}}),{})
        # A failed call must leave every output component untouched.
        x = (ct.c_double * 4)(1,2,4,8)
        p = (ct.c_double * 1)(0)
        out = (ct.c_double * 4)(91,92,93,94)
        fn = self.native.library.gsl_stats_aggregates_evaluate
        self.assertNotEqual(fn(1,b"gzscore",x,4,None,0,p,1,None,0,out,3),0)
        self.assertEqual(list(out),[91,92,93,94])
        self.assertNotEqual(fn(0,b"gzscore",x,4,None,0,p,1,None,0,out,4),0)
        self.assertEqual(list(out),[91,92,93,94])


if __name__ == "__main__":
    unittest.main()
