"""Independent public-Boost.Math Bessel/Struve/Kelvin/Wright audit.

Compiles ONLY native/bessel_hypergeom.cpp and verifies every saved SciPy case
as well as closed forms/series that do not call Sci or its numerical backend.
"""

import cmath
import csv
import ctypes as ct
import math
import os
from pathlib import Path
import subprocess
import tempfile
import time
import unittest

import numpy as np
from scipy import special

from sci_testing.bessel_hypergeom import BesselHypergeomComparisons, NAMES


ROOT = Path(__file__).resolve().parents[1]
ARITIES = {
    "struve": 2, "modstruve": 2,
    "ber": 1, "bei": 1, "berp": 1, "beip": 1,
    "wright_bessel": 3, "log_wright_bessel": 3,
}


def _complex_j01(z):
    """Independent defining series for J0(z), J1(z), not a Boost wrapper."""
    z2 = z*z/4
    j0_term, j0 = 1.+0j, 1.+0j
    j1_term, j1 = z/2, z/2
    for k in range(1, 80):
        j0_term *= -z2/(k*k)
        j1_term *= -z2/(k*(k+1))
        j0 += j0_term
        j1 += j1_term
        if abs(j0_term) + abs(j1_term) < 1e-16:
            break
    return j0, j1


def _wright_half_series(b, x):
    """Log-gamma evaluated definition, independent of 0F2/pFq."""
    if x == 0:
        return 1/math.gamma(b)
    total = 0.
    for k in range(300):
        term = math.exp(k*math.log(x) - math.lgamma(k+1) - math.lgamma(k/2+b))
        total += term
        if k >= 10 and term < total*1e-16:
            break
    else:
        raise AssertionError("independent Wright series failed to converge")
    return total


class BesselHypergeomTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(temp.cleanup)
        shared = Path(temp.name) / "bessel-hypergeom.so"
        env = dict(os.environ)
        env.pop("LD_LIBRARY_PATH", None)  # CoS/AppImage may shadow host C++ runtimes.
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2", "-fPIC", "-shared",
            str(ROOT / "native/bessel_hypergeom.cpp"), "-o", str(shared),
        ], check=True, env=env)
        cls.comparison = BesselHypergeomComparisons(ct.CDLL(str(shared)))

    @staticmethod
    def case(name, *values):
        return {"api": "special." + name,
                "args": {f"p{i}": value for i, value in enumerate(values)}}

    def assert_scalar(self, name, args, reference, rtol=2e-12, atol=2e-14):
        actual = self.comparison.evaluate(self.case(name, *args))
        self.assertEqual(set(actual), {"Boost.Math"}, (name, args, actual))
        self.assertNotIn("error", actual["Boost.Math"], (name, args, actual))
        values = actual["Boost.Math"]["values"]
        self.assertEqual(len(values), 1)
        self.assertTrue(math.isfinite(values[0]), (name, args, values))
        self.assertTrue(math.isclose(values[0], reference,
                                     rel_tol=rtol, abs_tol=atol),
                        (name, args, values[0], reference))

    def test_exact_library_registration(self):
        self.assertEqual(set(ARITIES), NAMES)
        for name, expected in ARITIES.items():
            encoded = name.encode("ascii")
            self.assertEqual(self.comparison.library.bessel_hypergeom_supported(
                0, encoded), expected)
            self.assertEqual(self.comparison.library.bessel_hypergeom_supported(
                1, encoded), 0)
        self.assertEqual(self.comparison.evaluate(self.case("not_supported", .5)), {})
        self.assertEqual(self.comparison.evaluate(self.case("ber", 2.+1.j)), {})
        self.assertEqual(self.comparison.evaluate({"api": "stats.ber", "args": {"p0": 1.}}), {})

    def test_all_16008_saved_cases_against_scipy_and_performance(self):
        start = time.perf_counter()
        max_absolute, max_relative = 0., 0.
        total = 0
        for name, arity in ARITIES.items():
            path = ROOT / "cases/generated/special" / f"{name}.csv"
            with path.open(newline="") as source:
                rows = list(csv.DictReader(source))
            self.assertEqual(len(rows), 2001, name)
            expected_columns = {f"p{i}" for i in range(arity)}
            for row in rows:
                self.assertTrue(expected_columns.issubset(row), name)
            inputs = [[float(row[f"p{i}"].split(":", 1)[1])
                       for i in range(arity)] for row in rows]
            # SciPy is an INDEPENDENT validation source, evaluated as vectors.
            targets = np.asarray(getattr(special, name)(*np.asarray(inputs).T),
                                 dtype=float).ravel()
            self.assertEqual(len(targets), 2001)
            for args, target in zip(inputs, targets):
                observed = self.comparison.evaluate(self.case(name, *args))
                self.assertEqual(set(observed), {"Boost.Math"}, (name, args))
                result = observed["Boost.Math"]
                self.assertNotIn("error", result, (name, args, result))
                self.assertEqual(len(result["values"]), 1)
                value = result["values"][0]
                self.assertTrue(math.isfinite(value) and math.isfinite(float(target)),
                                (name, args, value, target))
                error = abs(value-target)
                max_absolute = max(max_absolute, error)
                max_relative = max(max_relative, error/max(abs(target), 1e-12))
                self.assertTrue(math.isclose(value, float(target), rel_tol=2e-12,
                                             abs_tol=2e-14),
                                (name, args, value, target))
            total += len(rows)
        elapsed = time.perf_counter() - start
        self.assertEqual(total, 8*2001)
        # Loose runaway-performance tripwire only, not a competitive benchmark:
        # the actual observed throughput is reported below for reproducibility.
        self.assertLess(elapsed, 60., f"16,008 native comparisons took {elapsed:.2f}s")
        print(f"\nBoost 0F2/0F3/pFq 16,008 saved evaluations: "
              f"{elapsed:.3f}s, {total/elapsed:.0f}/s, "
              f"max abs {max_absolute:.3g}, max scaled {max_relative:.3g}")

    def test_kelvin_against_complex_bessel_defining_series(self):
        angle = cmath.exp(3j*math.pi/4)
        for x in (-6., -2., -0.5, 0., .4995, .5, .5005, 2., 6.):
            j0, j1 = _complex_j01(x*angle)
            values = {
                "ber": j0.real,
                "bei": j0.imag,
                "berp": (-angle*j1).real,
                "beip": (-angle*j1).imag,
            }
            for name, target in values.items():
                with self.subTest(name=name, x=x):
                    self.assert_scalar(name, (x,), target, rtol=8e-12, atol=4e-12)

    def test_half_integer_struve_elementary_identities(self):
        for x in (.001, .2, 1.998, 2., 2.002, 4., 8.):
            prefactor = math.sqrt(2/(math.pi*x))
            self.assert_scalar("struve", (.5, x),
                               prefactor*(1-math.cos(x)), rtol=2e-10, atol=2e-14)
            self.assert_scalar("modstruve", (.5, x),
                               prefactor*(math.cosh(x)-1), rtol=2e-10, atol=2e-14)

    def test_wright_against_independent_positive_series(self):
        for b in (.5, 1., 1.5, 2.5, 4.):
            for x in (0., .01, .4995, .5, .5005, 1., 4., 8.):
                expected = _wright_half_series(b, x)
                self.assert_scalar("wright_bessel", (.5, b, x), expected)
                self.assert_scalar("log_wright_bessel", (.5, b, x),
                                   math.log(expected), atol=5e-14)

    def test_other_valid_struve_domains_against_scipy(self):
        for nu in (0., .5, 1., 1.75, 3., 4.):
            for x in (0., .25, 1., 3., 6., 8.):
                for name in ("struve", "modstruve"):
                    with self.subTest(nu=nu, x=x, name=name):
                        self.assert_scalar(name, (nu, x),
                                           float(getattr(special, name)(nu, x)),
                                           rtol=2e-10, atol=4e-13)

    def test_explicit_unsupported_domains_and_no_output_mutation(self):
        cases = [
            ("struve", (-1., 2.)), ("modstruve", (.5, -2.)),
            ("modstruve", (4.5, 2.)), ("struve", (.5, 8.1)),
            ("ber", (8.1,)), ("bei", (-8.1,)),
            ("berp", (float("nan"),)), ("beip", (float("inf"),)),
            ("wright_bessel", (0., 1.5, .5)),
            ("log_wright_bessel", (.5000000001, 1.5, .5)),
            ("wright_bessel", (.5, .4, .5)),
            ("log_wright_bessel", (.5, 1.5, -1.)),
            ("wright_bessel", (.5, 1.5, 9.)),
        ]
        for name, args in cases:
            with self.subTest(name=name, args=args):
                reported = self.comparison.evaluate(self.case(name, *args))
                self.assertIn("error", reported["Boost.Math"], reported)
                # Native failures must not expose a partially written value.
                raw_args = (ct.c_double * len(args))(*args)
                output = (ct.c_double * 1)(1234.5)
                status = self.comparison.library.bessel_hypergeom_evaluate(
                    0, name.encode(), raw_args, len(args), output, 1)
                self.assertNotEqual(status, 0)
                self.assertEqual(output[0], 1234.5)
        # ABI rejects GSL rather than masquerading as a second provider.
        args = (ct.c_double * 1)(.5)
        output = (ct.c_double * 1)(1234.5)
        self.assertNotEqual(self.comparison.library.bessel_hypergeom_evaluate(
            1, b"ber", args, 1, output, 1), 0)
        self.assertEqual(output[0], 1234.5)


if __name__ == "__main__":
    unittest.main()
