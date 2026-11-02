"""SciPy full-array reference and public Boost.Math Legendre ABI tests."""

import ctypes as ct
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
from scipy import special

from sci_testing.boost_special_arrays import BoostSpecialArraysComparisons
from sci_testing.case_csv import read_case_file


ROOT = Path(__file__).resolve().parents[1]


class BoostSpecialArraysTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(tmp.cleanup)
        library = Path(tmp.name) / "boost-special-arrays.so"
        environment = os.environ.copy()
        environment.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2", "-Wall", "-Wextra",
            "-Werror", "-shared", "-fPIC",
            str(ROOT / "native/boost_special_arrays.cpp"), "-o", str(library),
        ], env=environment, check=True)
        cls.native = BoostSpecialArraysComparisons(ct.CDLL(str(library)))

    @staticmethod
    def case(name, *args, prefix="special"):
        return {"api": f"{prefix}.{name}",
                "args": {f"p{i}": value for i, value in enumerate(args)}}

    def assert_matches(self, name, *args, atol=3e-13, rtol=4e-12):
        result = self.native.evaluate(self.case(name, *args))
        self.assertEqual(set(result), {"Boost.Math"})
        self.assertNotIn("error", result["Boost.Math"], result)
        expected = np.concatenate(getattr(special, name)(*args))
        actual = np.asarray(result["Boost.Math"]["values"])
        self.assertEqual(actual.shape, expected.shape)
        np.testing.assert_allclose(actual, expected, rtol=rtol, atol=atol)

    def test_all_saved_real_scipy_cases(self):
        for name in ("roots_legendre", "lqn"):
            cases = read_case_file(ROOT / f"cases/generated/special/{name}.csv",
                                   f"special/{name}.csv")
            self.assertEqual(len(cases), 2001)
            for case in cases:
                with self.subTest(case=case["id"]):
                    args = case["args"]
                    self.assertEqual(set(args), {"p0"} if name == "roots_legendre"
                                     else {"p0", "p1"})
                    values = (args["p0"],) if name == "roots_legendre" else (
                        args["p0"], args["p1"])
                    self.assert_matches(name, *values)

    def test_roots_symmetry_positive_weights_and_exact_sum(self):
        for n in (1, 2, 3, 4, 7, 16, 24, 48, 100):
            self.assert_matches("roots_legendre", n)
            output = self.native.evaluate(self.case("roots_legendre", n))["Boost.Math"]["values"]
            nodes, weights = output[:n], output[n:]
            self.assertTrue(all(a < b for a, b in zip(nodes, nodes[1:])))
            self.assertTrue(all(w > 0 for w in weights))
            self.assertTrue(math.isclose(sum(weights), 2., abs_tol=6e-15))
            for i in range(n):
                self.assertAlmostEqual(nodes[i], -nodes[n - 1 - i], delta=4e-15)
                self.assertAlmostEqual(weights[i], weights[n - 1 - i], delta=4e-15)

    def test_legendre_q_interior_exterior_and_derivatives(self):
        for n in (0, 1, 2, 3, 4, 8, 16, 32):
            for x in (-.99, -.5, 0., .3, .99):
                self.assert_matches("lqn", n, x, rtol=5e-11, atol=1e-12)
        for n in (0, 1, 2, 3, 4):
            for x in (-3., -2., -1.2, 1.2, 2., 3.):
                self.assert_matches("lqn", n, x, rtol=2e-8, atol=1e-12)
        # SciPy's negative exterior forward recurrence loses accuracy rapidly
        # at higher orders; compare to the mathematical real-Q continuation.
        for n in (8, 16, 32):
            for x in (1.2, 2., 3.):
                self.assert_matches("lqn", n, x, rtol=5e-10, atol=5e-14)
        import mpmath as mp
        with mp.workdps(75):
            for n in (8, 16, 32, 128):
                for x in (-1.2, -2., -3.) if n != 128 else (-1.0001, 1.0001):
                    results = self.native.evaluate(self.case("lqn", n, x))
                    self.assertNotIn("error", results["Boost.Math"])
                    actual = results["Boost.Math"]["values"]
                    z = mp.mpf(x)
                    reference = [mp.log((z + 1) / (z - 1)) / 2]
                    reference.append(z * reference[0] - 1)
                    for k in range(1, n):
                        reference.append(((2*k + 1) * z * reference[k]
                                          - k * reference[k-1]) / (k + 1))
                    for k, value in enumerate(reference[:n+1]):
                        expected = float(value)
                        self.assertTrue(math.isclose(actual[k], expected,
                                                      rel_tol=5e-10, abs_tol=5e-14))

    def test_rejects_unvalidated_shapes_and_domains(self):
        cases = (
            self.case("roots_legendre", 0),
            self.case("roots_legendre", 1.5),
            self.case("roots_legendre", 129),
            self.case("roots_legendre", [3]),
            self.case("roots_legendre", 3, 0.2),
            self.case("lqn", 2),
            self.case("lqn", -1, 0.3),
            self.case("lqn", 3, 1.0),
            self.case("lqn", 3, float("nan")),
            self.case("lqn", 9, 1.00000001),
            self.case("lqn", 2, 0.3 + 1j),
            self.case("root_legendre", 3),
        )
        for case in cases:
            with self.subTest(case=case):
                self.assertEqual(self.native.evaluate(case), {})

    def test_abi_capacity_and_library_selection(self):
        library = self.native.library
        self.assertEqual(library.boost_special_arrays_supported(0, b"roots_legendre"), 2)
        self.assertEqual(library.boost_special_arrays_supported(0, b"lqn"), 2)
        self.assertEqual(library.boost_special_arrays_supported(1, b"lqn"), 0)
        params = (ct.c_double * 1)(3.)
        self.assertEqual(library.boost_special_arrays_output_count(b"roots_legendre", params, 1), 6)
        short = (ct.c_double * 5)()
        self.assertEqual(library.boost_special_arrays_evaluate(
            0, b"roots_legendre", params, 1, short, 5), -1)
        self.assertIn(b"capacity", library.boost_special_arrays_error())
        self.assertEqual(library.boost_special_arrays_evaluate(
            1, b"roots_legendre", params, 1, short, 5), -1)


if __name__ == "__main__":
    unittest.main()
