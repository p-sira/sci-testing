"""Exhaustive saved-case and independent SciPy Hankel adapter checks."""

import csv
import ctypes as ct
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
from scipy import special

from sci_testing.bessel_hankel import BesselHankelComparisons, SUPPORTED_NAMES


class BesselHankelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tmp.cleanup)
        root = Path(__file__).resolve().parent.parent
        shared = Path(cls.tmp.name) / "bessel-hankel.so"
        environment = os.environ.copy()
        environment.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2", "-shared", "-fPIC",
            str(root / "native/bessel_hankel.cpp"), "-lgsl", "-lgslcblas",
            "-o", str(shared),
        ], env=environment, check=True)
        cls.native = BesselHankelComparisons(ct.CDLL(str(shared)))
        cls.generated = root / "cases/generated/special"

    @staticmethod
    def case(name, *arguments):
        return {"api": f"special.{name}",
                "args": {f"p{i}": x for i, x in enumerate(arguments)}}

    def assert_matches_scipy(self, name, *arguments, rtol=2e-10, atol=2e-13):
        results = self.native.evaluate(self.case(name, *arguments))
        self.assertEqual(set(results), {"Boost.Math", "GSL"})
        expected = complex(getattr(special, name)(*arguments))
        self.assertTrue(math.isfinite(expected.real) and math.isfinite(expected.imag))
        for library, result in results.items():
            with self.subTest(name=name, arguments=arguments, library=library):
                self.assertNotIn("error", result, result)
                self.assertEqual(len(result["values"]), 2)
                np.testing.assert_allclose(
                    result["values"], [expected.real, expected.imag],
                    rtol=rtol, atol=atol,
                )

    def test_all_12006_saved_special_cases(self):
        """Every saved real-scalar signature, including the 1000+1000 sweeps."""
        count = 0
        for name in sorted(SUPPORTED_NAMES):
            with (self.generated / f"{name}.csv").open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2001, name)
            self.assertEqual({row["p0"] for row in rows}, {"s:0.5"}, name)
            for row in rows:
                values = []
                for key in sorted(k for k in row if k.startswith("p")):
                    kind, encoded = row[key].split(":", 1)
                    self.assertIn(kind, ("s", "i"))
                    values.append(int(encoded) if kind == "i" else float(encoded))
                self.assertTrue(1.998 <= values[1] <= 2.002)
                self.assert_matches_scipy(name, *values, rtol=1e-11, atol=5e-14)
                count += 1
        self.assertEqual(count, 12006)

    def test_negative_orders_derivatives_scaling_and_real_input_edges(self):
        for order, x in ((0., 0.5), (0.5, 2.), (1., 3.), (2.3, 1.25),
                         (5., 20.), (-0.5, 0.25), (-1.5, 2.),
                         (-2.3, 6.), (-4., 3.), (8., 12.5)):
            for name in ("hankel1", "hankel2", "hankel1e", "hankel2e"):
                self.assert_matches_scipy(name, order, x, rtol=4e-10)
            for name in ("h1vp", "h2vp"):
                for n in (0, 1, 2, 3, 4):
                    self.assert_matches_scipy(name, order, x, n, rtol=4e-10)

    def test_invalid_domains_report_errors_without_output(self):
        for arguments in ((0.5, 0.), (0.5, -1.), (0.5, math.inf),
                          (257., 2.), (math.nan, 2.), (0.5, 1e7)):
            results = self.native.evaluate(self.case("hankel1", *arguments))
            self.assertEqual(set(results), {"Boost.Math", "GSL"})
            for value in results.values():
                self.assertIn("error", value)
        for n in (-1, 0.5, math.inf, 13):
            results = self.native.evaluate(self.case("h1vp", 0.5, 2., n))
            for value in results.values():
                self.assertIn("error", value)
        self.assertEqual(self.native.evaluate(self.case("hankel1", .5, [2., 0.])), {})
        self.assertEqual(self.native.evaluate(self.case("hankel1", [.5, 0.], 2.)), {})
        self.assertEqual(self.native.evaluate(self.case("hankel1", .5)), {})
        self.assertEqual(self.native.evaluate(self.case("jv", .5, 2.)), {})
        self.assertEqual(self.native.evaluate({"api": "special.complex.hankel1",
                                               "args": {"p0": .5, "p1": 2.}}), {})

        native = self.native.library
        self.assertEqual(native.bessel_hankel_supported(0, b"h1vp"), 3)
        self.assertEqual(native.bessel_hankel_supported(1, b"hankel1e"), 2)
        self.assertEqual(native.bessel_hankel_supported(2, b"hankel1"), 0)
        inputs = (ct.c_double * 2)(.5, 2.)
        output = (ct.c_double * 2)(73., 74.)
        self.assertNotEqual(native.bessel_hankel_evaluate(
            0, b"hankel1", inputs, 2, output, 1), 0)
        self.assertEqual(list(output), [73., 74.])


if __name__ == "__main__":
    unittest.main()
