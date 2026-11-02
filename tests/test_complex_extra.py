"""Independent numeric checks of public GSL complex Gamma adapters."""

import ctypes as ct
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import mpmath as mp
from scipy import special

from sci_testing import special_complex03_cases as complex_cases
from sci_testing.complex_extra import ComplexExtraComparisons


ROOT = Path(__file__).resolve().parent.parent
NAMES = ("gamma", "loggamma", "rgamma")


class ComplexExtraTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        directory = tempfile.TemporaryDirectory()
        cls.addClassCleanup(directory.cleanup)
        library_path = Path(directory.name) / "complex-extra.so"
        compiler = os.environ.get("CXX", shutil.which("clang++") or "c++")
        environment = dict(os.environ)
        # Running from an AppImage may inject incompatible runtime libraries
        # into GCC/clang and the Python virtualenv.
        environment.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            compiler, "-std=c++17", "-O2", "-Wall", "-Wextra", "-Werror",
            "-fPIC", "-shared", str(ROOT / "native/complex_extra.cpp"),
            "-lgsl", "-lgslcblas", "-o", str(library_path),
        ], env=environment, check=True)
        cls.library = ct.CDLL(str(library_path))
        cls.native = ComplexExtraComparisons(cls.library)

    @staticmethod
    def case(name, real, imag):
        return {
            "api": f"special.complex.{name}",
            "args": {"special_complex03_case": True, "z": [real, imag]},
        }

    def get_value(self, name, real, imag):
        result = self.native.evaluate(self.case(name, real, imag))
        self.assertEqual(set(result), {"GSL"})
        self.assertNotIn("error", result["GSL"], (name, real, imag, result))
        observed = result["GSL"]["values"]
        self.assertEqual(len(observed), 2)
        return complex(*observed)

    def assert_complex_close(self, actual, reference, *, rtol=3e-11, atol=3e-12):
        expected = complex(reference)
        self.assertTrue(math.isclose(actual.real, expected.real, rel_tol=rtol, abs_tol=atol),
                        (actual, expected))
        self.assertTrue(math.isclose(actual.imag, expected.imag, rel_tol=rtol, abs_tol=atol),
                        (actual, expected))

    def test_independent_mpmath_and_scipy_reference(self):
        # Covers all quadrants, beyond +/-pi gamma phase, negative-real side
        # approaches, imaginary oscillation, and values with high log magnitude.
        points = [
            (0.25, 0.5), (0.5, -0.7), (1.0, 0.0), (3.5, 4.0),
            (-0.75, 0.25), (-0.75, -0.25), (-3.2, 1e-10),
            (-3.2, -1e-10), (1.0, 35.0), (-10.0, 11.0),
            (-20.0, 40.0), (170.0, 0.0), (1000.0, 0.5),
        ]
        with mp.workdps(95):
            for real, imag in points:
                z = mp.mpc(real, imag)
                for name in NAMES:
                    with self.subTest(name=name, z=(real, imag)):
                        if name == "gamma" and real == 1000.0:
                            self.assertIn(
                                "finite double range",
                                self.native.evaluate(
                                    self.case(name, real, imag))["GSL"]["error"])
                            continue
                        target = getattr(mp, name)(z)
                        observed = self.get_value(name, real, imag)
                        self.assert_complex_close(observed, target, rtol=1e-10)
                        scipy_result = getattr(special, name)(complex(real, imag))
                        if math.isfinite(scipy_result.real) and math.isfinite(scipy_result.imag):
                            self.assert_complex_close(observed, scipy_result, rtol=1e-10)

    def test_saved_complex_fixture_shapes_and_reference(self):
        # The official cases exercise real/imag pair keys, signed zero and
        # high-precision independent oracles rather than native Sci internals.
        for name in NAMES:
            selected = list(complex_cases.cases(
                12, 621, [f"special.complex.{name}"]))
            self.assertGreater(len(selected), 25)
            with mp.workdps(100):
                for case in selected:
                    with self.subTest(case=case["id"]):
                        z = case["args"]["z"]
                        if name in ("gamma", "loggamma") and z[0] == 0. and z[1] == 0.:
                            self.assertIn(
                                "pole", self.native.evaluate(case)["GSL"]["error"])
                            continue
                        expected = complex(*map(float, complex_cases.reference(case)))
                        actual = self.get_value(name, *z)
                        self.assert_complex_close(actual, expected, rtol=1e-10)

    def test_signed_negative_real_axis_loggamma(self):
        # Sci's analytic log-Gamma uses one-sided banks, not GSL's wrapped
        # principal Arg(Gamma). SciPy ignores some signed-zero variants, so
        # the independent mpmath-backed fixture oracle is authoritative.
        for real, expected_winding in ((-0.5, 1), (-1.5, 2), (-3.2, 4),
                                       (-7.5, 8)):
            for imag in (0.0, -0.0):
                with self.subTest(real=real, below=math.copysign(1., imag) < 0):
                    result = self.get_value("loggamma", real, imag)
                    side = 1.0 if math.copysign(1.0, imag) < 0 else -1.0
                    self.assertAlmostEqual(result.imag, side * math.pi * expected_winding,
                                           places=13)
                    self.assert_complex_close(
                        self.get_value("gamma", real, imag), mp.gamma(real),
                        rtol=2e-11)
                    self.assert_complex_close(
                        self.get_value("rgamma", real, imag), mp.rgamma(real),
                        rtol=2e-11)
            upper = self.get_value("loggamma", real, 1e-12)
            lower = self.get_value("loggamma", real, -1e-12)
            self.assertAlmostEqual(upper.imag, -lower.imag, places=10)

    def test_nonpositive_integer_poles_near_poles_and_reciprocal_zeros(self):
        for x in (0., -1., -2., -4., -15.):
            self.assertEqual(self.get_value("rgamma", x, 0.0), 0j)
            self.assertEqual(self.get_value("rgamma", x, -0.0), 0j)
            for name in ("gamma", "loggamma"):
                self.assertIn("pole", self.native.evaluate(
                    self.case(name, x, 0.))["GSL"]["error"])
        with mp.workdps(120):
            for imag in (1e-14, -1e-14, 1e-100, -1e-100, 1e-300):
                z = mp.mpc(-4, imag)
                for name in NAMES:
                    with self.subTest(name=name, imag=imag):
                        observed = self.get_value(name, -4., imag)
                        expected = complex(getattr(mp, name)(z))
                        # Near poles, compare the dominant imaginary component
                        # relatively, and the subdominant real component with a
                        # separate bound to detect exp(phase) cancellation.
                        self.assertTrue(math.isclose(
                            observed.imag, expected.imag, rel_tol=2e-11,
                            abs_tol=1e-290), (observed, expected))
                        if abs(expected.real) > 1e-300:
                            self.assertTrue(math.isclose(
                                observed.real, expected.real, rel_tol=2e-10,
                                abs_tol=1e-12), (observed, expected))

    def test_input_validation_unavailable_functions_and_abi_safety(self):
        self.assertEqual(self.library.complex_extra_supported(b"gamma"), 2)
        self.assertEqual(self.library.complex_extra_supported(b"loggamma"), 2)
        self.assertEqual(self.library.complex_extra_supported(b"rgamma"), 2)
        self.assertEqual(self.library.complex_extra_supported(b"digamma"), 0)
        self.assertEqual(self.native.evaluate(
            {"api": "special.complex.digamma", "args": {"z": [.5, .5]}}), {})
        self.assertEqual(self.native.evaluate(
            {"api": "special.gamma", "args": {"z": [.5, .5]}}), {})
        with self.assertRaisesRegex(ValueError, "expected a complex"):
            self.native.evaluate({"api": "special.complex.gamma", "args": {"z": [1.]}})
        for x, y in ((float("nan"), 1.), (1., float("inf"))):
            self.assertIn("nonfinite complex input",
                          self.native.evaluate(self.case("gamma", x, y))["GSL"]["error"])
        for name in NAMES:
            self.assertIn("comparison range",
                          self.native.evaluate(self.case(name, -3000., .5))["GSL"]["error"])

        input_values = (ct.c_double * 2)(.25, .5)
        out = (ct.c_double * 2)(17., 19.)
        for library, name, capacity in ((1, b"gamma", 1),
                                        (0, b"gamma", 2),
                                        (1, b"no-such-function", 2)):
            status = self.library.complex_extra_evaluate(
                library, name, input_values, 2, out, capacity)
            self.assertEqual(status, -1)
            self.assertEqual(list(out), [17., 19.])
            self.assertTrue(self.library.complex_extra_error())


if __name__ == "__main__":
    unittest.main()
