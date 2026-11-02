"""Independent SciPy and saved-case checks of combinatorics_extra C ABI."""

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

from sci_testing.combinatorics_extra import CombinatoricsExtraComparisons


class CombinatoricsExtraTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temporary.cleanup)
        root = Path(__file__).resolve().parent.parent
        shared = Path(cls.temporary.name) / "combinatorics-extra.so"
        environment = os.environ.copy()
        environment.pop("LD_LIBRARY_PATH", None)  # CoS AppImage can shadow host GCC.
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2", "-fPIC", "-shared",
            str(root / "native/combinatorics_extra.cpp"), "-lgsl", "-lgslcblas",
            "-o", str(shared)], env=environment, check=True)
        cls.comparisons = CombinatoricsExtraComparisons(ct.CDLL(str(shared)))

    @staticmethod
    def case(name, *arguments):
        return {"api": f"special.{name}", "args": {
            f"p{i}": argument for i, argument in enumerate(arguments)}}

    @staticmethod
    def reference(name, *arguments):
        if name == "bernoulli":
            return np.asarray(special.bernoulli(int(arguments[0])), dtype=float).ravel()
        if name == "harmonic":
            n = int(arguments[0])
            return np.array([float(special.digamma(n + 1) + np.euler_gamma) if n else 0.])
        if name == "gen_harmonic":
            n, exponent = arguments
            return np.array([math.fsum(j ** (-exponent) for j in range(1, int(n) + 1))])
        raise AssertionError(name)

    def assert_comparison(self, name, *arguments, expected_libraries, atol=2e-12, rtol=8e-12):
        actual = self.comparisons.evaluate(self.case(name, *arguments))
        self.assertEqual(set(actual), set(expected_libraries), (name, arguments, actual))
        reference = self.reference(name, *arguments)
        for library, result in actual.items():
            with self.subTest(api=name, args=arguments, library=library):
                self.assertNotIn("error", result, result)
                np.testing.assert_allclose(result["values"], reference,
                                           atol=atol, rtol=rtol)

    def test_complete_bernoulli_vector_and_finite_range(self):
        for n in (0, 1, 2, 3, 4, 8, 16, 32, 100, 200, 258):
            self.assert_comparison("bernoulli", n, expected_libraries=("Boost.Math",))
        result = self.comparisons.evaluate(self.case("bernoulli", 260))
        self.assertIn("error", result["Boost.Math"])
        self.assertEqual(self.comparisons.library.combinatorics_extra_supported(
            1, b"bernoulli"), 0)

    def test_harmonic_integer_edges_up_to_one_million(self):
        for n in (0, 1, 2, 8, 16, 100, 1000, 1000000):
            self.assert_comparison("harmonic", n, expected_libraries=("Boost.Math", "GSL"))

    def test_generalized_harmonic_mathematically_supported_region(self):
        for n, exponent in (
            (0, -4.), (0, 2.), (1, 1.), (2, 1.), (2, 1.998),
            (5, 2.), (100, 2.), (1000, 2.), (1000000, 2.),
            (5, 1000.), (1000000, 1000.), (5, 1.1),
        ):
            self.assert_comparison("gen_harmonic", n, exponent,
                                   expected_libraries=("GSL",))
        for n, exponent in ((5, .5), (100, -1.5), (1000000, 0),
                            (1, 1.00000001)):
            result = self.comparisons.evaluate(self.case("gen_harmonic", n, exponent))
            self.assertIn("error", result["GSL"], result)

    def test_all_actual_generated_inputs(self):
        root = Path(__file__).resolve().parent.parent / "cases/generated/special"
        for name in ("bernoulli", "harmonic", "gen_harmonic"):
            with (root / f"{name}.csv").open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            self.assertGreaterEqual(len(rows), 2000)
            for row in rows:
                arguments = []
                for key in sorted(key for key in row if key.startswith("p")):
                    tag, raw = row[key].split(":", 1)
                    arguments.append(int(raw) if tag == "i" else float(raw))
                result = self.comparisons.evaluate(self.case(name, *arguments))
                expected_libraries = ("Boost.Math",) if name == "bernoulli" else (
                    ("Boost.Math", "GSL") if name == "harmonic" else ("GSL",))
                self.assertEqual(set(result), set(expected_libraries), row)
                if name == "gen_harmonic" and arguments[0] != 0 and arguments[1] < 1.:
                    self.assertIn("error", result["GSL"], row)
                    continue
                expected = self.reference(name, *arguments)
                for library, actual in result.items():
                    with self.subTest(case=row["id"], library=library):
                        self.assertNotIn("error", actual, (row, actual))
                        np.testing.assert_allclose(actual["values"], expected,
                                                   rtol=2e-11, atol=3e-12)

    def test_invalid_input_and_capacity_checks(self):
        for name, arguments in (
            ("bernoulli", (-1,)), ("bernoulli", (1.5,)),
            ("harmonic", (-1,)), ("harmonic", (1000001,)),
            ("gen_harmonic", (2.5, 2.)), ("gen_harmonic", (5, math.inf)),
        ):
            result = self.comparisons.evaluate(self.case(name, *arguments))
            self.assertTrue(result)
            for value in result.values():
                self.assertIn("error", value)
        for case in (self.case("bernoulli", [2, 4]),
                     self.case("harmonic", 2 + 3j),
                     {"api": "special.bernoulli", "args": {}},
                     self.case("factorial", 4)):
            self.assertEqual(self.comparisons.evaluate(case), {})
        native = self.comparisons.library
        parameters = (ct.c_double * 1)(8.)
        undersized = (ct.c_double * 1)()
        status = native.combinatorics_extra_evaluate(
            0, b"bernoulli", parameters, 1, undersized, 1)
        self.assertNotEqual(status, 0)
        self.assertIn(b"capacity", native.combinatorics_extra_error())


if __name__ == "__main__":
    unittest.main()
