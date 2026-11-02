"""GSL 2.7 public polynomial, Steffen, and interval-search interpolation ABI."""

import ctypes as ct
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
from scipy import interpolate

from sci_testing.case_csv import read_case_file
from sci_testing.gsl_interpolation import (
    GslInterpolationComparisons, LIBRARY, _interpn, _simple)


ROOT = Path(__file__).resolve().parents[1]
POLYNOMIAL = ("barycentric_interpolate", "krogh_interpolate")


def saved_cases(name):
    file = ROOT / f"cases/generated/interpolate/{name}.csv.gz"
    if not file.is_file():
        file = file.with_suffix("")
    return read_case_file(file, f"interpolate/{name}.csv")


def independent_scipy(case):
    args = case["args"]
    name = case["api"].split(".", 1)[-1]
    if name == "interpn":
        axes, values, query, method, bounds, fill = args["raw"]["args"]
        shape = tuple(len(ax) for ax in axes)
        return np.asarray(interpolate.interpn(
            axes, np.asarray(values).reshape(shape), query,
            method=method, bounds_error=bounds, fill_value=fill)).ravel()
    return np.asarray(getattr(interpolate, name)(
        args["nodes"], args["values"], args["query"])).ravel()


class GslInterpolationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        folder = tempfile.TemporaryDirectory()
        cls.addClassCleanup(folder.cleanup)
        path = Path(folder.name) / "gsl_interpolation.so"
        environ = os.environ.copy()
        environ.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2", "-Wall",
            "-Wextra", "-Werror", "-fPIC", "-shared",
            str(ROOT / "native/gsl_interpolation.cpp"), "-lgsl",
            "-lgslcblas", "-o", str(path)], check=True, env=environ)
        cls.adapter = GslInterpolationComparisons(ct.CDLL(str(path)))

    def assert_matching(self, case, target, rtol=2e-10, atol=5e-12):
        result = self.adapter.evaluate(case)
        self.assertEqual(set(result), {LIBRARY}, case["id"])
        self.assertNotIn("error", result[LIBRARY], (case["id"], result))
        got = np.asarray(result[LIBRARY]["values"], dtype=float)
        want = np.asarray(target, dtype=float).ravel()
        self.assertEqual(got.shape, want.shape, case["id"])
        np.testing.assert_allclose(
            got, want, rtol=rtol, atol=atol, err_msg=case["id"])

    def test_all_4100_polynomial_saved_cases(self):
        count = 0
        for name in POLYNOMIAL:
            records = saved_cases(name)
            self.assertEqual(len(records), 2050)
            for case in records:
                with self.subTest(case=case["id"]):
                    self.assert_matching(case, independent_scipy(case))
                    count += 1
        self.assertEqual(count, 4100)

    def test_all_2014_multidimensional_saved_cases(self):
        records = saved_cases("interpn")
        self.assertEqual(len(records), 2014)
        count = 0
        for case in records:
            with self.subTest(case=case["id"]):
                self.assert_matching(case, independent_scipy(case))
                count += 1
        self.assertEqual(count, 2014)

    def test_steffen_not_misrepresented_as_scipy_pchip(self):
        # Every saved varying-slope PCHIP case must be excluded rather than
        # silently replaced with a different GSL monotone interpolant.
        records = saved_cases("pchip_interpolate")
        self.assertEqual(len(records), 2000)
        for case in records:
            with self.subTest(case=case["id"]):
                self.assertEqual(self.adapter.evaluate(case), {})
        affine = {
            "id": "affine-steffen",
            "api": "interpolate.pchip_interpolate",
            "args": {"nodes": [0., 1., 2., 3.],
                     "values": [2., 4., 6., 8.],
                     "query": [0.25, 1.5, 2.75]},
        }
        self.assert_matching(affine, independent_scipy(affine))
        knots = {
            "id": "nonaffine-exact-knots",
            "api": "interpolate.pchip_interpolate",
            "args": {"nodes": [0., 1., 2., 3.],
                     "values": [1., 1.5, 1.2, 2.],
                     "query": [0., 1., 2., 3.]},
        }
        self.assert_matching(knots, independent_scipy(knots))

    def test_domain_guards(self):
        valid = saved_cases("barycentric_interpolate")[0]
        invalid_args = [
            {"nodes": [0., 0., 1.], "values": [1., 2., 3.], "query": [0.2]},
            {"nodes": [1., 0., 2.], "values": [1., 2., 3.], "query": [0.2]},
            {"nodes": [0., 1.], "values": [1.], "query": [0.2]},
            {"nodes": [0., 1.], "values": [1., 2.], "query": []},
            {"nodes": [0., 1.], "values": [1., 2.], "query": [float("inf")]},
            {"nodes": [0., 1.], "values": [1., 2.], "query": ["bad"]},
        ]
        for args in invalid_args:
            case = dict(valid, args=args)
            self.assertEqual(self.adapter.evaluate(case), {})
        wrong = dict(valid, api="interpolate.griddata")
        self.assertEqual(self.adapter.evaluate(wrong), {})
        raw = saved_cases("interpn")[0]["args"]
        malformed = dict(raw, raw=dict(raw["raw"], args=[[]]))
        self.assertEqual(self.adapter.evaluate(
            {"api": "interpolate.interpn", "args": malformed}), {})
        outside = {
            "fields": [], "suite": "interpolation", "suite_name": "interpolation",
            "raw": {"api": "interpn", "args": [
                [[0., 1.]], [0., 1.], [[2.]], "linear", True, None]},
        }
        self.assertEqual(self.adapter.evaluate(
            {"api": "interpolate.interpn", "args": outside}), {})

    def test_abi_capacity_error_and_library_select(self):
        lib = self.adapter.library
        for name in (*POLYNOMIAL, "pchip_interpolate", "interpn"):
            encoded = f"interpolate.{name}".encode()
            self.assertEqual(lib.gsl_interpolation_supported(1, encoded), 2)
            self.assertEqual(lib.gsl_interpolation_supported(0, encoded), 0)
        params = (ct.c_double * 10)(
            3., 2., 0., 1., 2., 1., 2., 3., 0.5, 1.5)
        expected_count = lib.gsl_interpolation_output_count(
            b"interpolate.barycentric_interpolate", params, 10)
        self.assertEqual(expected_count, 2)
        output = (ct.c_double * 2)()
        self.assertEqual(lib.gsl_interpolation_evaluate(
            1, b"interpolate.barycentric_interpolate",
            params, 10, output, 2), 0)
        np.testing.assert_allclose(list(output), [1.5, 2.5])
        self.assertEqual(lib.gsl_interpolation_evaluate(
            1, b"interpolate.barycentric_interpolate",
            params, 10, output, 1), -1)
        self.assertIn(b"capacity", lib.gsl_interpolation_error())
        self.assertEqual(lib.gsl_interpolation_evaluate(
            0, b"interpolate.barycentric_interpolate",
            params, 10, output, 2), -1)
        self.assertEqual(lib.gsl_interpolation_output_count(
            b"interpolate.barycentric_interpolate", params, 7), 0)


if __name__ == "__main__":
    unittest.main()
