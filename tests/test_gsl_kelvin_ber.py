"""Exhaustive saved-case tests for GSL-derived Kelvin ber/bei derivatives.

Build the original GSL Kelvin integrator and the new lane-extraction wrapper
together. Compare all 8,004 saved scalar inputs against independent SciPy;
the underlying Kelvin result is used only for an exact tuple-layout check.
"""

import cmath
import ctypes as ct
import math
import os
from pathlib import Path
import subprocess
import tempfile
import time
import unittest

from scipy import special

from sci_testing.case_csv import read_case_file
from sci_testing.gsl_kelvin_ber import GSLKelvinBerComparisons, NAMES


ROOT = Path(__file__).resolve().parents[1]
LANES = {"ber": 0, "bei": 1, "berp": 4, "beip": 5}


class GSLKelvinBerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        temporary = tempfile.TemporaryDirectory()
        cls.addClassCleanup(temporary.cleanup)
        target = Path(temporary.name) / "gsl_kelvin_ber.so"
        environment = os.environ.copy()
        environment.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2",
            "-Wall", "-Wextra", "-Werror", "-Wpedantic", "-shared", "-fPIC",
            str(ROOT / "native/bessel_kelvin_k.cpp"),
            str(ROOT / "native/gsl_kelvin_ber.cpp"),
            "-lgsl", "-lgslcblas", "-o", str(target),
        ], env=environment, check=True)
        cls.library = ct.CDLL(str(target))
        cls.native = GSLKelvinBerComparisons(cls.library)
        cls.library.bessel_kelvin_k_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int)
        cls.library.bessel_kelvin_k_evaluate.restype = ct.c_int

    @staticmethod
    def case(name, x):
        return {"api": "special." + name, "args": {"p0": x}}

    def get(self, name, x):
        comparison = self.native.evaluate(self.case(name, x))
        self.assertEqual(set(comparison), {"GSL"}, (name, x, comparison))
        self.assertNotIn("error", comparison["GSL"], comparison)
        values = comparison["GSL"]["values"]
        self.assertEqual(len(values), 1)
        self.assertTrue(math.isfinite(values[0]))
        return values[0]

    def test_all_8004_actual_saved_cases_against_scipy(self):
        compared = 0
        started = time.perf_counter()
        for name in ("ber", "bei", "berp", "beip"):
            path = ROOT / "cases/generated/special" / (name + ".csv")
            cases = read_case_file(path, f"special/{name}.csv")
            self.assertEqual(len(cases), 2001, name)
            for case in cases:
                with self.subTest(case=case["id"]):
                    self.assertEqual(set(case["args"]), {"p0"})
                    x = case["args"]["p0"]
                    self.assertGreaterEqual(x, .1)
                    self.assertLessEqual(x, 8.)
                    actual = self.get(name, x)
                    independent = float(getattr(special, name)(x))
                    self.assertTrue(math.isclose(actual, independent,
                                                  rel_tol=2e-10, abs_tol=1e-12),
                                    f"{case['id']}: {actual} != {independent}")
                    compared += 1
        self.assertEqual(compared, 8004)
        # The original Kelvin integrator already meets a generous bounded
        # budget. An accidentally repeated full report regeneration would not.
        self.assertLess(time.perf_counter() - started, 30.)

    def test_exact_tuple_lane_mapping_and_amos_reference(self):
        # Check exact transport of the EXISTING ABI's lanes, rather than
        # copying its mathematical integration in this test or adapter.
        rotation = cmath.exp(1j * math.pi / 4)
        for x in (.1, .25, .4995, .5, .5005, 1., 2., 5., 8.):
            argument = (ct.c_double * 1)(x)
            tuple_out = (ct.c_double * 8)()
            status = self.library.bessel_kelvin_k_evaluate(
                1, b"kelvin", argument, 1, tuple_out, 8)
            self.assertEqual(status, 0)

            z = rotation * x
            complex_be = special.iv(0, z)
            complex_bep = rotation * special.iv(1, z)
            independent = (complex_be.real, complex_be.imag,
                           complex_bep.real, complex_bep.imag)
            for (name, lane), expected in zip(LANES.items(), independent):
                with self.subTest(api=name, x=x):
                    actual = self.get(name, x)
                    self.assertEqual(actual, tuple_out[lane])
                    self.assertTrue(math.isclose(actual, expected,
                                                  rel_tol=2e-10, abs_tol=3e-11))

    def test_full_validated_domain(self):
        for i in range(101):
            x = .1 + 7.9*i/100
            for name in NAMES:
                with self.subTest(api=name, x=x):
                    expected = float(getattr(special, name)(x))
                    actual = self.get(name, x)
                    self.assertTrue(math.isclose(actual, expected,
                                                  rel_tol=2e-10, abs_tol=3e-11))

    def test_explicit_checked_c_abi_and_error_recovery(self):
        library = self.library
        for name in LANES:
            encoded = name.encode("ascii")
            self.assertEqual(library.gsl_kelvin_ber_supported(1, encoded), 1)
            self.assertEqual(library.gsl_kelvin_ber_supported(0, encoded), 0)
            self.assertEqual(library.gsl_kelvin_ber_output_count(encoded), 1)
            input_x = (ct.c_double * 1)(.5)
            output = (ct.c_double * 2)(-999., -998.)
            for mode, args, count, cap in (
                (0, input_x, 1, 2), (1, input_x, 2, 2),
                (1, input_x, 1, 0), (1, None, 1, 1),
            ):
                status = library.gsl_kelvin_ber_evaluate(
                    mode, encoded, args, count, output, cap)
                self.assertNotEqual(status, 0)
                self.assertEqual(list(output), [-999., -998.])
                self.assertTrue(library.gsl_kelvin_ber_error())
            self.assertEqual(library.gsl_kelvin_ber_evaluate(
                1, encoded, input_x, 1, output, 2), 0)
            self.assertAlmostEqual(output[0], getattr(special, name)(.5),
                                   delta=2e-12)
            self.assertEqual(output[1], -998.)
            self.assertEqual(library.gsl_kelvin_ber_error(), b"")

        for name in (b"kelvin", b"ker", b"kei", b"keip", b"unknown"):
            self.assertEqual(library.gsl_kelvin_ber_supported(1, name), 0)
            self.assertEqual(library.gsl_kelvin_ber_output_count(name), 0)
        self.assertEqual(library.gsl_kelvin_ber_supported(1, None), 0)
        self.assertEqual(library.gsl_kelvin_ber_output_count(None), 0)

    def test_input_domains_and_unsupported_shapes(self):
        for name in LANES:
            for x in (-1., -.0, 0., .099999, 8.000001,
                      math.inf, -math.inf, math.nan):
                result = self.native.evaluate(self.case(name, x))
                with self.subTest(api=name, input=x):
                    self.assertEqual(set(result), {"GSL"})
                    self.assertIn("error", result["GSL"])
                    self.assertIn("0.1 <= x <= 8", result["GSL"]["error"])

        for case in (
            self.case("ber", complex(.5,.1)),
            self.case("bei", [.5, .6]),
            self.case("berp", {"x": .5}),
            {"api": "special.ber", "args": {"p0": .5, "p1": .6}},
            {"api": "special.bei", "args": {}},
            self.case("unknown", .5),
            {"api": "special.complex.ber", "args": {"p0": .5}},
        ):
            self.assertEqual(self.native.evaluate(case), {})


if __name__ == "__main__":
    unittest.main()
