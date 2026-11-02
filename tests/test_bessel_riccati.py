"""Riccati Bessel array ABI checks against independent SciPy reference."""

import ctypes as ct
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
from scipy import special

from sci_testing.bessel_riccati import BesselRiccatiComparisons
from sci_testing.case_csv import read_case_file


class BesselRiccatiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tmp.cleanup)
        library = Path(cls.tmp.name) / "bessel_riccati.so"
        source = Path(__file__).resolve().parent.parent / "native/bessel_riccati.cpp"
        environment = os.environ.copy()
        environment.pop("LD_LIBRARY_PATH", None)
        subprocess.run(
            [os.environ.get("CXX", "c++"), "-std=c++17", "-O2", "-shared",
             "-fPIC", str(source), "-lgsl", "-lgslcblas", "-o", str(library)],
            env=environment, check=True)
        cls.native = BesselRiccatiComparisons(ct.CDLL(str(library)))

    @staticmethod
    def case(name, n, x):
        return {"api": f"special.{name}", "args": {"p0": n, "p1": x}}

    def assert_matches(self, name, n, x, rtol=2e-10, atol=5e-13):
        expected = np.concatenate(getattr(special, name)(n, x))
        results = self.native.evaluate(self.case(name, n, x))
        self.assertEqual(set(results), {"Boost.Math", "GSL"})
        for library, result in results.items():
            with self.subTest(name=name, order=n, argument=x, library=library):
                self.assertNotIn("error", result, result)
                self.assertEqual(len(result["values"]), 2 * (n+1))
                np.testing.assert_allclose(result["values"], expected, rtol=rtol,
                                           atol=atol, err_msg=f"{name} {library}")

    def test_saved_signature_and_complete_fixture(self):
        root = Path(__file__).resolve().parent.parent / "cases/generated/special"
        for name in ("riccati_jn", "riccati_yn"):
            cases = read_case_file(root / f"{name}.csv", f"special/{name}.csv")
            self.assertEqual(len(cases), 2001)
            for case in cases:
                self.assertEqual(set(case["args"]), {"p0", "p1"})
                self.assert_matches(name, **{
                    "n": case["args"]["p0"], "x": case["args"]["p1"]})

    def test_orders_derivatives_and_negative_regular_argument(self):
        for name in ("riccati_jn", "riccati_yn"):
            for n in (0, 1, 2, 3, 5, 10):
                for x in (0.1, 0.5, 1.0, 3.0, 9.0):
                    self.assert_matches(name, n, x, rtol=2e-9, atol=1e-12)
        for n in (0, 1, 2, 3, 5):
            for x in (-0.1, -0.5, -3.0):
                self.assert_matches("riccati_jn", n, x, rtol=2e-9, atol=1e-12)

    def test_regular_origin_and_irregular_invalid_side(self):
        for n in (0, 1, 3, 10):
            self.assert_matches("riccati_jn", n, 0.)
            for x in (0., -0.5):
                results = self.native.evaluate(self.case("riccati_yn", n, x))
                self.assertEqual(set(results), {"Boost.Math", "GSL"})
                for result in results.values():
                    self.assertIn("error", result)
                    self.assertIn("sentinel", result["error"])

    def test_nonmatching_complex_and_invalid_shapes(self):
        for case in (
            self.case("riccati_jn", 1, 0.5+0.1j),
            self.case("riccati_yn", 3, [0.5, 1.]),
            self.case("riccati_jn", -1, 0.5),
            self.case("riccati_jn", 3.5, 0.5),
            self.case("riccati_jn", 513, 0.5),
            {"api": "special.riccati_jn", "args": {"p0": 2}},
            {"api": "special.riccati_jn", "args": {"p0": 2, "p1": .5, "p2": 0}},
            self.case("riccati_unknown", 3, 0.5),
        ):
            self.assertEqual(self.native.evaluate(case), {})

    def test_capabilities_and_guarded_capacity(self):
        lib = self.native.library
        for name in (b"riccati_jn", b"riccati_yn"):
            for library in (0, 1):
                self.assertEqual(lib.bessel_riccati_supported(library, name), 2)
            params = (ct.c_double * 2)(3., 0.5)
            self.assertEqual(lib.bessel_riccati_output_count(name, params, 2), 8)
            small = (ct.c_double * 7)()
            self.assertEqual(lib.bessel_riccati_evaluate(
                0, name, params, 2, small, 7), -1)
            self.assertIn(b"capacity", lib.bessel_riccati_error())
        self.assertEqual(lib.bessel_riccati_supported(3, b"riccati_jn"), 0)


if __name__ == "__main__":
    unittest.main()
