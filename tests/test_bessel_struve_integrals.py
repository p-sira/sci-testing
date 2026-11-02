"""Validate public Boost pFq-derived Struve integrals against SciPy."""

import csv
import ctypes as ct
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
from scipy import integrate, special

from sci_testing.bessel_struve_integrals import (
    BesselStruveIntegralsComparisons,
    SUPPORTED_NAMES,
)


class BesselStruveIntegralsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tmp.cleanup)
        root = Path(__file__).resolve().parent.parent
        cls.generated = root / "cases/generated/special"
        shared = Path(cls.tmp.name) / "struve-integrals.so"
        environment = os.environ.copy()
        environment.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2", "-shared", "-fPIC",
            str(root / "native/bessel_struve_integrals.cpp"),
            "-o", str(shared),
        ], env=environment, check=True)
        cls.native = BesselStruveIntegralsComparisons(ct.CDLL(str(shared)))

    def evaluate(self, name, x, *, namespace="special"):
        return self.native.evaluate({
            "api": f"{namespace}.{name}", "args": {"p0": x},
        })

    def assert_matches_scipy(self, name, x, rtol=2e-11, atol=3e-13):
        values = self.evaluate(name, x)
        self.assertEqual(set(values), {"Boost.Math"})
        result = values["Boost.Math"]
        self.assertNotIn("error", result, (name, x, result))
        self.assertEqual(len(result["values"]), 1)
        reference = getattr(special, name)(x)
        np.testing.assert_allclose(result["values"], [reference],
                                   rtol=rtol, atol=atol)
        return result["values"][0]

    def test_all_6003_saved_cases_against_independent_scipy(self):
        cases = 0
        for name in sorted(SUPPORTED_NAMES):
            with (self.generated / f"{name}.csv").open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2001, name)
            for row in rows:
                self.assertEqual(set(k for k in row if k.startswith("p")), {"p0"})
                kind, encoded = row["p0"].split(":", 1)
                self.assertEqual(kind, "s")
                x = float(encoded)
                self.assertTrue(1.998 <= x <= 2.002)
                with self.subTest(case=row["id"]):
                    self.assert_matches_scipy(name, x, rtol=3e-12, atol=3e-14)
                cases += 1
        self.assertEqual(cases, 6003)

    def test_integral_identities_parity_and_boundaries(self):
        for x in (-8., -5., -2., -0.5, -0.02, 0., 0.02, .5, 2., 5., 8.):
            for name in sorted(SUPPORTED_NAMES):
                with self.subTest(name=name, x=x):
                    value = self.assert_matches_scipy(name, x, rtol=5e-10,
                                                      atol=2e-11)
                    if name == "itstruve0":
                        integral = integrate.quad(lambda t: special.struve(0, t),
                                                  0., x, epsabs=1e-12)[0]
                    elif name == "itmodstruve0":
                        integral = integrate.quad(lambda t: special.modstruve(0, t),
                                                  0., x, epsabs=1e-12)[0]
                    else:
                        integral = math.pi/2 - integrate.quad(
                            lambda t: special.struve(0, t)/t if t else 2/math.pi,
                            0., x, epsabs=1e-12)[0]
                    self.assertAlmostEqual(value, integral, delta=2e-11)

        for x in (.1, 2., 6.):
            for name in ("itstruve0", "itmodstruve0"):
                self.assertEqual(self.evaluate(name, x), self.evaluate(name, -x))
            plus = self.evaluate("it2struve0", x)["Boost.Math"]["values"][0]
            minus = self.evaluate("it2struve0", -x)["Boost.Math"]["values"][0]
            self.assertAlmostEqual(plus+minus, math.pi, delta=5e-15)

    def test_input_errors_abi_bounds_and_unsupported_providers(self):
        native = self.native.library
        for name in SUPPORTED_NAMES:
            self.assertEqual(native.bessel_struve_integrals_supported(
                0, name.encode("ascii")), 1)
            self.assertEqual(native.bessel_struve_integrals_supported(
                1, name.encode("ascii")), 0)

        for x in (-9., 8.01, math.inf, -math.inf, math.nan):
            for name in SUPPORTED_NAMES:
                record = self.evaluate(name, x)
                self.assertEqual(set(record), {"Boost.Math"})
                self.assertIn("error", record["Boost.Math"])
        self.assertEqual(self.evaluate("itstruve0", [2., 0.]), {})
        self.assertEqual(self.evaluate("itstruve0", 2., namespace="special.complex"), {})
        self.assertEqual(self.evaluate("unknown", 2.), {})
        self.assertEqual(self.native.evaluate({"api": "special.itstruve0", "args": {}}), {})

        argument = (ct.c_double * 1)(2.)
        out = (ct.c_double * 1)(123.45)
        self.assertNotEqual(native.bessel_struve_integrals_evaluate(
            0, b"itstruve0", argument, 1, out, 0), 0)
        self.assertEqual(out[0], 123.45)
        self.assertNotEqual(native.bessel_struve_integrals_evaluate(
            1, b"itstruve0", argument, 1, out, 1), 0)
        self.assertEqual(out[0], 123.45)
        self.assertNotEqual(native.bessel_struve_integrals_evaluate(
            0, b"itstruve0", argument, 2, out, 1), 0)
        self.assertEqual(out[0], 123.45)


if __name__ == "__main__":
    unittest.main()
