"""Independent SciPy checks of GSL 2.7+ fixed quadrature roots adapters."""

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

from sci_testing.case_csv import read_case_file
from sci_testing.gsl_quadrature_roots import (
    GSLQuadratureRootsComparisons, ROOT_ARITIES,
)


class GSLQuadratureRootsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tmp.cleanup)
        library = Path(cls.tmp.name) / "gsl-quadrature-roots.so"
        source = Path(__file__).resolve().parent.parent / (
            "native/gsl_quadrature_roots.cpp")
        env = os.environ.copy()
        env.pop("LD_LIBRARY_PATH", None)
        subprocess.run(
            [os.environ.get("CXX", "c++"), "-std=c++17", "-O2",
             "-shared", "-fPIC", str(source), "-lgsl", "-lgslcblas",
             "-o", str(library)],
            env=env, check=True,
        )
        cls.adapter = GSLQuadratureRootsComparisons(ct.CDLL(str(library)))

    @staticmethod
    def case(name, *arguments):
        return {
            "api": f"special.{name}",
            "args": {f"p{i}": value for i, value in enumerate(arguments)},
        }

    def compare_reference(self, name, *arguments, rtol=2e-10, atol=2e-12):
        result = self.adapter.evaluate(self.case(name, *arguments))
        self.assertEqual(set(result), {"GSL"})
        self.assertNotIn("error", result["GSL"], (name, arguments, result))
        expected_nodes, expected_weights = getattr(special, name)(*arguments)
        actual = result["GSL"]["values"]
        n = int(arguments[0])
        self.assertEqual(len(actual), 2*n)
        np.testing.assert_allclose(
            actual[:n], expected_nodes, rtol=rtol, atol=atol,
            err_msg=f"{name}{arguments} nodes")
        np.testing.assert_allclose(
            actual[n:], expected_weights, rtol=rtol, atol=atol,
            err_msg=f"{name}{arguments} weights")
        self.assertTrue(all(w > 0. and math.isfinite(w) for w in actual[n:]))
        self.assertTrue(all(a < b for a, b in zip(actual[:n-1], actual[1:n])))
        # Independently confirm that weights have the measure's correct
        # normalization, without claiming optional mu=True is in saved API.
        _, _, mu = getattr(special, name)(*arguments, mu=True)
        self.assertAlmostEqual(sum(actual[n:])/float(mu), 1., delta=5e-10)

    def test_all_saved_cases_independent_scipy(self):
        root = Path(__file__).resolve().parent.parent / "cases/generated/special"
        total_cases = 0
        for name, arity in ROOT_ARITIES.items():
            cases = read_case_file(root / f"{name}.csv", f"special/{name}.csv")
            self.assertEqual(len(cases), 2001, name)
            for case in cases:
                self.assertEqual(set(case["args"]), {
                    f"p{i}" for i in range(arity)})
                parameters = [case["args"][f"p{i}"] for i in range(arity)]
                self.compare_reference(name, *parameters)
                total_cases += 1
        self.assertEqual(total_cases, 2001*15)

    def test_more_orders_and_parameter_domains(self):
        cases = {
            "roots_legendre": [(1,), (2,), (5,), (15,), (32,)],
            "roots_jacobi": [(3,.5,.3), (10,.2,.7), (12,1.5,2.3)],
            "roots_gegenbauer": [(3,.4995), (10,.25), (12,1.5)],
            "roots_laguerre": [(1,), (5,), (16,)],
            "roots_genlaguerre": [(2,.2), (8,1.5), (12,2.0)],
            "roots_hermite": [(1,), (5,), (16,)],
            "roots_hermitenorm": [(1,), (5,), (16,)],
            "roots_chebyt": [(1,), (5,), (16,)],
            "roots_chebyu": [(1,), (5,), (16,)],
            "roots_chebyc": [(1,), (5,), (16,)],
            "roots_chebys": [(1,), (5,), (16,)],
            "roots_sh_legendre": [(1,), (5,), (16,)],
            "roots_sh_jacobi": [(2,2.5,1.5), (10,1.7,.9)],
            "roots_sh_chebyt": [(1,), (5,), (16,)],
            "roots_sh_chebyu": [(1,), (5,), (16,)],
        }
        for name, families in cases.items():
            for args in families:
                with self.subTest(name=name, args=args):
                    self.compare_reference(name, *args, rtol=1e-9, atol=1e-11)

    def test_abi_supported_and_capacity_guards(self):
        native = self.adapter.library
        for name, arity in ROOT_ARITIES.items():
            encoded = name.encode("ascii")
            self.assertEqual(native.gsl_quadrature_roots_supported(0, encoded), 0)
            self.assertEqual(native.gsl_quadrature_roots_supported(1, encoded), arity)
            params = (ct.c_double * arity)(3., *([.5] * (arity-1)))
            count = native.gsl_quadrature_roots_output_count(encoded, params, arity)
            self.assertEqual(count, 6)
            output = (ct.c_double * 6)(*([-99.] * 6))
            self.assertEqual(native.gsl_quadrature_roots_evaluate(
                1, encoded, params, arity, output, 5), -1)
            self.assertIn(b"capacity", native.gsl_quadrature_roots_error())
            self.assertEqual(list(output), [-99.] * 6)
            self.assertEqual(native.gsl_quadrature_roots_evaluate(
                1, encoded, params, arity+1, output, 6), -1)
            self.assertEqual(native.gsl_quadrature_roots_evaluate(
                1, encoded, None, arity, output, 6), -1)
            self.assertEqual(native.gsl_quadrature_roots_evaluate(
                1, encoded, params, arity, None, 6), -1)
        self.assertEqual(native.gsl_quadrature_roots_supported(1, b"unknown"), 0)

    def test_invalid_degree_domains_and_input_shapes(self):
        for n in (-1., 0., 1.5, 129., float("inf"), float("nan")):
            result = self.adapter.evaluate(self.case("roots_legendre", n))
            self.assertEqual(set(result), {"GSL"})
            self.assertIn("error", result["GSL"])
        for name, args in (
            ("roots_jacobi", (3., -1., .5)),
            ("roots_jacobi", (3., .5, -1.)),
            ("roots_genlaguerre", (3., -1.)),
            ("roots_gegenbauer", (3., 0.)),
            ("roots_sh_jacobi", (3., 2.5, 0.)),
            ("roots_sh_jacobi", (3., .5, 2.)),
        ):
            result = self.adapter.evaluate(self.case(name, *args))
            self.assertIn("error", result["GSL"])
        for case in (
            self.case("roots_legendre", [3]),
            self.case("roots_legendre", 3+1j),
            self.case("roots_legendre", {"degree": 3}),
            self.case("roots_jacobi", 3, .5),
            {"api": "special.roots_legendre", "args": {"p0":3, "mu":True}},
            self.case("roots_unknown", 3),
        ):
            self.assertEqual(self.adapter.evaluate(case), {})

    def test_2000_case_performance(self):
        start = time.perf_counter()
        for index in range(2000):
            actual = self.adapter.evaluate(
                self.case("roots_jacobi", 3, .5, .2997 + index*.0006/1999))
            self.assertIn("values", actual["GSL"])
        self.assertLess(time.perf_counter()-start, 10.)


if __name__ == "__main__":
    unittest.main()
