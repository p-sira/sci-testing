"""Independent full-fixture regression tests for Boost-derived Bessel ratios."""

import ctypes as ct
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
import warnings

import mpmath as mp
from scipy import special

from sci_testing.bessel_ratio import BesselRatioComparisons
from sci_testing import special03_cases
from sci_testing.case_csv import read_case_file


ROOT = Path(__file__).resolve().parents[1]
NAMES = ("iv_ratio", "iv_ratio_c", "iv_ratioinv")


class BesselRatioTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        workspace = tempfile.TemporaryDirectory()
        cls.addClassCleanup(workspace.cleanup)
        binary = Path(workspace.name) / "libbessel-ratio.so"
        compiler = os.environ.get("CXX", shutil.which("clang++") or "c++")
        environment = dict(os.environ)
        environment.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            compiler, "-std=c++17", "-O2", "-Wall", "-Wextra", "-Werror",
            "-shared", "-fPIC", str(ROOT / "native/bessel_ratio.cpp"),
            "-o", str(binary),
        ], env=environment, check=True)
        cls.library = ct.CDLL(str(binary))
        cls.native = BesselRatioComparisons(cls.library)

    @staticmethod
    def case(name, v, argument):
        key = "r" if name == "iv_ratioinv" else "x"
        return {
            "api": f"special.bessel.{name}",
            "args": {"special03_case": True, "v": v, key: argument},
        }

    def get(self, name, v, arg):
        record = self.native.evaluate(self.case(name, v, arg))
        self.assertEqual(set(record), {"Boost.Math"})
        self.assertNotIn("error", record["Boost.Math"], record)
        values = record["Boost.Math"]["values"]
        self.assertEqual(len(values), 1)
        self.assertTrue(math.isfinite(values[0]))
        return values[0]

    def test_all_6024_actual_saved_cases_against_independent_mpmath(self):
        with mp.workdps(75):
            for name in NAMES:
                path = ROOT / f"cases/generated/special.bessel/{name}.csv"
                self.assertTrue(path.exists(), path)
                cases = read_case_file(path, f"special.bessel/{name}.csv")
                self.assertEqual(len(cases), 2008)
                for case in cases:
                    with self.subTest(fixture=case["id"]):
                        observed = self.native.evaluate(case)
                        self.assertEqual(set(observed), {"Boost.Math"})
                        self.assertNotIn("error", observed["Boost.Math"], observed)
                        values = observed["Boost.Math"]["values"]
                        self.assertEqual(len(values), 1)
                        actual = values[0]
                        expected = float(special03_cases.reference(case)[0])
                        self.assertTrue(math.isfinite(actual))
                        self.assertTrue(math.isclose(
                            actual, expected, rel_tol=4e-11, abs_tol=2e-14),
                            f"{case['id']}: {actual} != {expected}")

    def test_public_boost_bessel_i_and_special_half_order(self):
        # Boost can directly compute I_v and I_{v-1} at moderate x. Cross
        # check that this is a true public-Bessel comparison, not a Sci
        # implementation clone.
        for v in (0.5, 1., 1.5, 3., 8., 32., 128.):
            for x in (1e-8, 0.125, 1., 20., 200., 500.):
                actual = self.get("iv_ratio", v, x)
                if v == .5:
                    self.assertAlmostEqual(actual, math.tanh(x), delta=5e-15)
                else:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", RuntimeWarning)
                        scaled = special.ive(v, x) / special.ive(v - 1, x)
                    # Reference only when SciPy actually resolves both
                    # scaled Bessel values without numerical underflow.
                    if math.isfinite(scaled):
                        self.assertTrue(math.isclose(
                            actual, scaled, rel_tol=3e-11, abs_tol=2e-14))

    def test_complement_preserves_tiny_values(self):
        signed_zero = self.get("iv_ratio", 1., -0.0)
        self.assertEqual(signed_zero, 0.0)
        self.assertLess(math.copysign(1.0, signed_zero), 0.0)
        with mp.workdps(85):
            for v in (0.5, 1., 2., 8., 32., 128.):
                for x in (1e-200, 1., 1000., 4001., 1e5, 1e8, 1e12):
                    r = self.get("iv_ratio", v, x)
                    c = self.get("iv_ratio_c", v, x)
                    self.assertTrue(0 <= r <= 1)
                    self.assertTrue(0 <= c <= 1)
                    self.assertAlmostEqual(r + c, 1., delta=2e-15)
                    if v == .5:
                        expected = 2/(1+mp.exp(2*mp.mpf(x)))
                    else:
                        z = mp.mpf(x)
                        expected = 1 - mp.besseli(v, z)/mp.besseli(v-1, z)
                    self.assertTrue(math.isclose(
                        c, float(expected), rel_tol=3e-11, abs_tol=1e-300),
                        f"complement({v},{x}): {c} != {expected}")
        # Ordinary subtraction rounds to zero at x=1e8 but this adapter
        # exposes the true O(1/x) complement.
        self.assertGreater(self.get("iv_ratio_c", 1., 1e8), 0.)

    def test_inverse_roundtrip_and_unrepresentable_regions(self):
        with mp.workdps(85):
            for v in (.5, 1., 2., 8., 32., 128.):
                for r in (0., 1e-12, .001, .1, .5, .9, .97, .999, .99999):
                    x = self.get("iv_ratioinv", v, r)
                    expected = 0. if r == 0. else (
                        mp.tanh(x) if v == .5 else mp.besseli(v,x)/mp.besseli(v-1,x))
                    self.assertTrue(math.isclose(
                        float(expected), r, rel_tol=5e-11, abs_tol=1e-14))
                    if x > 0:
                        self.assertTrue(math.isclose(
                            self.get("iv_ratio", v, x), r,
                            rel_tol=5e-11, abs_tol=2e-14))
        for v, r in ((128., 1 - 1e-12), (8., 1 - 1e-14)):
            record = self.native.evaluate(self.case("iv_ratioinv", v, r))
            self.assertIn("error", record["Boost.Math"])

    def test_domain_rejection_no_fake_comparisons(self):
        for name in NAMES:
            for v in (math.nan, math.inf, -.5, 0.49, 129.):
                record = self.native.evaluate(self.case(name, v, 1.))
                self.assertIn("error", record["Boost.Math"], (name, v, record))
        for x in (-1., math.nan, math.inf, 1e13):
            for name in ("iv_ratio", "iv_ratio_c"):
                self.assertIn(
                    "error", self.native.evaluate(self.case(name, 2., x))["Boost.Math"])
        for r in (-.01, math.nan, math.inf, 1., 1.01):
            self.assertIn(
                "error", self.native.evaluate(
                    self.case("iv_ratioinv", 2., r))["Boost.Math"])
        for bad in (
            {"api": "special.bessel.iv_ratio", "args": {"v": 2.}},
            {"api": "special.bessel.iv_ratio", "args": {"v": 2., "x": [1., 2.]}},
            {"api": "special.bessel.unknown", "args": {"v": 2., "x": 3.}},
            {"api": "special.iv_ratio", "args": {"v": 2., "x": 3.}},
        ):
            self.assertEqual(self.native.evaluate(bad), {})

    def test_checked_c_abi_and_independent_link(self):
        for name in NAMES:
            self.assertEqual(
                self.library.bessel_ratio_supported(0, name.encode()), 2)
            self.assertEqual(
                self.library.bessel_ratio_supported(1, name.encode()), 0)
        self.assertEqual(self.library.bessel_ratio_supported(0, b"unknown"), 0)
        args = (ct.c_double * 2)(1., 1.)
        output = (ct.c_double * 2)(42., 43.)
        for library, name, count, capacity in (
            (1, b"iv_ratio", 2, 1),
            (0, b"iv_ratio", 1, 1),
            (0, b"iv_ratio", 2, 0),
            (0, b"nonexistent", 2, 1),
        ):
            status = self.library.bessel_ratio_evaluate(
                library, name, args, count, output, capacity)
            self.assertEqual(status, -1)
            self.assertEqual(list(output), [42., 43.])
            self.assertTrue(self.library.bessel_ratio_error())


if __name__ == "__main__":
    unittest.main()
