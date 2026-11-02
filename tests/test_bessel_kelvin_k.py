"""Checks of the public-GSL-quadrature-derived Kelvin K-family ABI.

Tests compare against independent SciPy 1.18.1 special-function values,
not other derived adapters in the repository.
"""

import ctypes as ct
import cmath
import math
import os
from pathlib import Path
import subprocess
import tempfile
import time
import unittest

from scipy import special

from sci_testing.bessel_kelvin_k import BesselKelvinKComparisons
from sci_testing.case_csv import read_case_file


FUNCTIONS = ("ker", "kei", "kerp", "keip", "kelvin")


class BesselKelvinKTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temporary.cleanup)
        source = (Path(__file__).resolve().parent.parent /
                  "native/bessel_kelvin_k.cpp")
        shared = Path(cls.temporary.name) / "bessel_kelvin_k.so"
        environment = os.environ.copy()
        environment.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2",
            "-shared", "-fPIC", str(source), "-lgsl", "-lgslcblas",
            "-o", str(shared)], env=environment, check=True)
        cls.native = BesselKelvinKComparisons(ct.CDLL(str(shared)))

    @staticmethod
    def make_case(name, x):
        return {"api": f"special.{name}", "args": {"p0": x}}

    @staticmethod
    def reference(name, x):
        result = getattr(special, name)(x)
        if name == "kelvin":
            return [component for value in result for component in
                    (float(value.real), float(value.imag))]
        return [float(result)]

    def check_case(self, name, x, *, rtol=2e-10, atol=1e-12):
        actual = self.native.evaluate(self.make_case(name, x))
        self.assertEqual(set(actual), {"GSL"})
        self.assertNotIn("error", actual["GSL"], actual)
        values = actual["GSL"]["values"]
        reference = self.reference(name, x)
        self.assertEqual(len(values), len(reference))
        for index, (observed, expected) in enumerate(zip(values, reference)):
            self.assertTrue(
                math.isclose(observed, expected, rel_tol=rtol, abs_tol=atol),
                f"{name}(x={x}) index={index}: GSL derived {observed}, SciPy {expected}")

    def test_every_generated_saved_case(self):
        root = Path(__file__).resolve().parent.parent / "cases/generated/special"
        comparisons = 0
        for name in FUNCTIONS:
            cases = read_case_file(root / f"{name}.csv", f"special/{name}.csv")
            self.assertEqual(len(cases), 2001, name)
            for case in cases:
                self.assertEqual(set(case["args"]), {"p0"})
                self.check_case(name, case["args"]["p0"])
                comparisons += 8 if name == "kelvin" else 1
        self.assertEqual(comparisons, 8004 + 16008)

    def test_extended_validated_domain(self):
        # Probe the entire supported interval, not only the narrow saved sweep.
        for x in (.1, .15, .25, .49, .5, 1., 2., 3., 5., 8.):
            for name in FUNCTIONS:
                with self.subTest(api=name, x=x):
                    self.check_case(name, x, atol=3e-11)

    def test_dense_validated_domain(self):
        # Catch oscillatory sign changes and near-zero components away from
        # the dense, but narrow, generated sweep around x=0.5.
        for index in range(401):
            x = .1 + 7.9 * index / 400
            for name in FUNCTIONS:
                self.check_case(name, x, atol=3e-11)

    def test_independent_complex_amos_bessel_identities(self):
        # SciPy's complex-argument AMOS Bessels provide an independent
        # reference from the real Kelvin routines used in other tests.
        rotation = cmath.exp(1j*math.pi/4)
        for x in (.1, .25, .5, 1., 2., 5., 8.):
            z = rotation*x
            be = special.iv(0, z)
            ke = special.kv(0, z)
            bep = rotation*special.iv(1, z)
            kep = -rotation*special.kv(1, z)
            expected = [coordinate for item in (be,ke,bep,kep)
                        for coordinate in (item.real, item.imag)]
            result = self.native.evaluate(self.make_case("kelvin", x))
            self.assertIn("values", result["GSL"], result)
            for index, (actual, target) in enumerate(zip(result["GSL"]["values"], expected)):
                self.assertTrue(
                    math.isclose(actual, target, rel_tol=2e-10, abs_tol=3e-11),
                    f"AMOS discrepancy x={x} index={index}: {actual} != {target}")

    def test_unsupported_domain_and_nonfinite_inputs(self):
        for name in FUNCTIONS:
            for x in (-1., 0., .099, 8.001, float("inf"), float("nan")):
                with self.subTest(api=name, x=x):
                    result = self.native.evaluate(self.make_case(name, x))
                    self.assertEqual(set(result), {"GSL"})
                    self.assertIn("error", result["GSL"])
                    self.assertIn("0.1 <= x <= 8", result["GSL"]["error"])

    def test_complex_vector_and_nonmatching_shapes(self):
        for case in (
            self.make_case("ker", .5+.01j),
            self.make_case("kei", [.5, 1.]),
            self.make_case("kelvin", {"x": .5}),
            {"api": "special.ker", "args": {"p0": .5, "p1": 1.}},
            {"api": "special.ker", "args": {}},
            self.make_case("unknown", .5),
            {"api": "special.complex.ker", "args": {"p0": .5}},
        ):
            self.assertEqual(self.native.evaluate(case), {})

    def test_supported_libraries_and_exact_abi_capacity(self):
        library = self.native.library
        for name in FUNCTIONS:
            encoded = name.encode()
            count = 8 if name == "kelvin" else 1
            self.assertEqual(library.bessel_kelvin_k_supported(0, encoded), 0)
            self.assertEqual(library.bessel_kelvin_k_supported(1, encoded), count)
            self.assertEqual(library.bessel_kelvin_k_output_count(encoded), count)
            x = (ct.c_double * 1)(.5)
            output = (ct.c_double * 8)(*([-99.] * 8))
            self.assertEqual(library.bessel_kelvin_k_evaluate(
                1, encoded, x, 1, output, count-1), -1)
            self.assertIn(b"capacity", library.bessel_kelvin_k_error())
            self.assertEqual(list(output), [-99.] * 8)
            self.assertEqual(library.bessel_kelvin_k_evaluate(
                1, encoded, x, 2, output, 8), -1)
            self.assertEqual(library.bessel_kelvin_k_evaluate(
                1, encoded, None, 1, output, 8), -1)
            self.assertEqual(library.bessel_kelvin_k_evaluate(
                1, encoded, x, 1, None, 8), -1)
        self.assertEqual(library.bessel_kelvin_k_supported(1, b"ber"), 0)
        self.assertEqual(library.bessel_kelvin_k_output_count(b"ber"), 0)

    def test_bounded_2000_case_performance(self):
        start = time.perf_counter()
        for i in range(2000):
            result = self.native.evaluate(
                self.make_case("kelvin", .4995 + i*0.001/1999))
            self.assertIn("values", result["GSL"], result)
        elapsed = time.perf_counter()-start
        # Generous CI-safe ceiling; observed native runtime is ~0.1s.
        self.assertLess(elapsed, 15., f"2000 Kelvin vectors took {elapsed:.3f}s")


if __name__ == "__main__":
    unittest.main()
