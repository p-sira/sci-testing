"""Accuracy and capability checks for the independent special_extra C ABI."""

import ctypes as ct
import csv
import gzip
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
from scipy import special

from sci_testing.special_extra import SpecialExtraComparisons


class SpecialExtraTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tmp.cleanup)
        shared = Path(cls.tmp.name) / "special-extra.so"
        root = Path(__file__).resolve().parent.parent
        # CoS' AppImage library path can shadow the host GCC toolchain.
        environment = os.environ.copy()
        environment.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2", "-shared", "-fPIC",
            str(root / "native/special_extra.cpp"), "-lgsl", "-lgslcblas",
            "-o", str(shared)], env=environment, check=True)
        cls.native = SpecialExtraComparisons(ct.CDLL(str(shared)))

    @staticmethod
    def case(name, *arguments, path="special"):
        return {"api": f"{path}.{name}",
                "args": {f"p{i}": value for i, value in enumerate(arguments)}}

    def assert_matches_scipy(self, name, *arguments, libraries=None, path="special", rtol=2e-10, atol=5e-13):
        actual = self.native.evaluate(self.case(name, *arguments, path=path))
        self.assertTrue(actual, name)
        if libraries is not None:
            self.assertEqual(set(actual), set(libraries))
        expected = np.asarray(getattr(special, name)(*arguments), dtype=float).ravel()
        for library, result in actual.items():
            with self.subTest(name=name, arguments=arguments, library=library):
                self.assertNotIn("error", result, result)
                self.assertEqual(len(result["values"]), len(expected))
                np.testing.assert_allclose(result["values"], expected, rtol=rtol, atol=atol)

    def test_derivative_recurrences_and_negative_orders(self):
        for name in ("jvp", "yvp", "ivp", "kvp"):
            for n in (0, 1, 2, 3, 4):
                for order, x in ((0.5, 2.), (2., 1.25), (3.5, 5.)):
                    self.assert_matches_scipy(name, order, x, n,
                                              libraries=("Boost.Math", "GSL"))

    def test_bessel_root_vectors(self):
        for order, count in ((0, 1), (1, 3), (2, 5), (10, 5)):
            self.assert_matches_scipy("jn_zeros", order, count,
                                      libraries=("Boost.Math", "GSL"))
            self.assert_matches_scipy("yn_zeros", order, count,
                                      libraries=("Boost.Math",))

    def test_orthogonal_and_shifted_polynomials(self):
        for degree in range(5):
            for x in (-.7, .3, .8):
                self.assert_matches_scipy("eval_chebyc", degree, x,
                                          libraries=("Boost.Math",))
                self.assert_matches_scipy("eval_chebys", degree, x,
                                          libraries=("Boost.Math",))
            self.assert_matches_scipy("eval_sh_jacobi", degree, 2.5, 1.5, .3,
                                      libraries=("Boost.Math",))
        for degree, order, x in ((0, 0, .3), (3, 1, .3), (4, 0, -.5),
                                 (5, 3, .8), (3, -1, .3), (5, -3, .8), (3, 4, .3)):
            self.assert_matches_scipy("legendre_p", degree, x,
                                      libraries=("Boost.Math", "GSL"))
            self.assert_matches_scipy("assoc_legendre_p", degree, order, x,
                                      libraries=("Boost.Math", "GSL"))
            # Both Sci 0.3 standalone and potential nested orthogonal API.
            self.assert_matches_scipy("legendre_p", degree, x,
                                      path="special.orthogonal")

    def test_elementary_and_combinatorics(self):
        for x in (-4., -0.5, 0., 1e-12, .5, 4.):
            self.assert_matches_scipy("sinc", x, libraries=("Boost.Math", "GSL"))
            self.assert_matches_scipy("exprel", x, libraries=("GSL",))
        for n, k in ((5., 2), (5.5, 2), (8., 0), (10., 9),
                     (3.5, 4), (-1., 0), (0., 0)):
            self.assert_matches_scipy("perm", n, k, libraries=("Boost.Math",))
        self.assertEqual(self.native.evaluate(self.case("softplus", .5)), {})

    def test_modified_spherical_bessel_values_and_derivatives(self):
        for name in ("spherical_in", "spherical_kn"):
            for n in (0, 1, 2, 4):
                for x in (0.5, 1.5, 6.):
                    for derivative in (False, True):
                        self.assert_matches_scipy(name, n, x, derivative,
                                                  libraries=("Boost.Math", "GSL"),
                                                  rtol=2e-9)
        for n in (0, 1, 2, 3):
            for x in (-1.5, -.5, 0.):
                for derivative in (False, True):
                    self.assert_matches_scipy("spherical_in", n, x, derivative,
                                              libraries=("Boost.Math", "GSL"))

    def test_argument_validation_and_abi_capabilities(self):
        for name, inputs in (("jn_zeros", (2., 0.)),
                             ("yn_zeros", (-1., 2.)),
                             ("jvp", (.5, 2., -1.)),
                             ("spherical_kn", (2., -.1, 0.)),
                             ("assoc_legendre_p", (2., 1.5, .5))):
            result = self.native.evaluate(self.case(name, *inputs))
            self.assertTrue(result)
            for value in result.values():
                self.assertIn("error", value)
        for case in (self.case("jvp", 1.+2.j, 2., 1.),
                     self.case("jvp", [1., 2.], 2., 1.),
                     self.case("unknown", 1.),
                     {"api": "special.jvp", "args": {"p0": 1}}):
            self.assertEqual(self.native.evaluate(case), {})

    def test_sampled_generated_special_cases(self):
        """Exercise the CSV input format used in the actual report pipeline."""
        root = Path(__file__).resolve().parent.parent / "cases/generated/special"
        names = (
            "jvp", "yvp", "ivp", "kvp", "jn_zeros", "yn_zeros",
            "eval_chebyc", "eval_chebys", "eval_sh_jacobi", "legendre_p",
            "assoc_legendre_p", "sinc", "exprel", "perm",
            "spherical_in", "spherical_kn",
        )
        for name in names:
            plain = root / f"{name}.csv"
            compressed = root / f"{name}.csv.gz"
            self.assertTrue(plain.exists() or compressed.exists(), name)
            opener = open if plain.exists() else gzip.open
            target = plain if plain.exists() else compressed
            with opener(target, mode="rt", newline="") as file:
                rows = list(csv.DictReader(file))
            for record in (rows[0], rows[len(rows)//2], rows[-1]):
                arguments = []
                for key in sorted(k for k in record if k.startswith("p")):
                    kind, raw = record[key].split(":", 1)
                    arguments.append(int(raw) if kind in ("i", "b") else float(raw))
                self.assert_matches_scipy(name, *arguments, rtol=2e-9)


if __name__ == "__main__":
    unittest.main()
