"""End-to-end saved-case comparisons of public Boost.Numeric.uBLAS operations."""

import ctypes as ct
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
import mpmath as mp
from scipy import linalg

from sci_testing.boost_ublas import BoostUblasComparisons, LIBRARY
from sci_testing.case_csv import read_case_file
from sci_testing.sparse_cases import reference as sparse_reference


ROOT = Path(__file__).resolve().parents[1]
DENSE = ("det", "solve", "dot", "matvec", "matmul", "transpose")
SPARSE = ("add", "multiply", "matvec", "matmul", "transpose", "from_dense", "to_dense")


def cases(module, name):
    path = ROOT / f"cases/generated/{module}/{name}.csv.gz"
    if not path.is_file():
        path = path.with_suffix("")
    return read_case_file(path, f"{module}/{name}.csv")


def dense_reference(name, args):
    if name == "dot":
        return [np.dot(args["u"], args["v"])]
    if name == "matvec":
        return np.asarray(args["a"], dtype=float) @ args["u"]
    if name == "matmul":
        return np.asarray(args["a"], dtype=float) @ np.asarray(args["b"], dtype=float)
    if name == "transpose":
        return np.asarray(args["a"], dtype=float).T
    if name == "det":
        return [linalg.det(np.asarray(args["a"], dtype=float))]
    return linalg.solve(np.asarray(args["a"], dtype=float), np.asarray(args["b"], dtype=float))


class BoostUblasTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        temporary = tempfile.TemporaryDirectory()
        cls.addClassCleanup(temporary.cleanup)
        library = Path(temporary.name) / "boost_ublas.so"
        environment = os.environ.copy()
        environment.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2", "-Wall",
            "-Wextra", "-Werror", "-shared", "-fPIC",
            str(ROOT / "native/boost_ublas.cpp"), "-o", str(library),
        ], check=True, env=environment)
        cls.native = BoostUblasComparisons(ct.CDLL(str(library)))

    def check_reference(self, case, expected):
        result = self.native.evaluate(case)
        self.assertEqual(set(result), {LIBRARY}, case["id"])
        self.assertNotIn("error", result[LIBRARY], (case["id"], result))
        actual = np.asarray(result[LIBRARY]["values"], dtype=float)
        target = np.asarray(expected, dtype=float).ravel()
        self.assertEqual(actual.shape, target.shape, case["id"])
        np.testing.assert_allclose(
            actual, target, rtol=2e-10, atol=2e-12,
            err_msg=case["id"])

    def test_all_dense_saved_cases_against_scipy(self):
        for name in DENSE:
            records = cases("linalg", name)
            self.assertEqual(len(records), 2058 if name in ("det", "solve") else 2000)
            for case in records:
                with self.subTest(case=case["id"]):
                    expected = dense_reference(name, case["args"])
                    if name == "solve":
                        a, b = case["args"]["a"], case["args"]["b"]
                        condition = np.linalg.cond(np.asarray(a, dtype=float))
                        if condition > 1e9:
                            # SciPy/LAPACK's binary64 solve can lose ~1e-7
                            # on the saved 8x8 Hilbert case. A wider Boost
                            # uBLAS LU legitimately differs from that result.
                            # Confirm the reference with independent 100-digit
                            # LU rather than rewarding a less accurate answer.
                            with mp.workdps(100):
                                exact = mp.lu_solve(mp.matrix(a), mp.matrix(b))
                                target = [float(v) for v in exact]
                            np.testing.assert_allclose(
                                expected, target,
                                rtol=25*np.finfo(float).eps*condition,
                                atol=2e-12, err_msg=case["id"] + " SciPy conditioning")
                            expected = target
                    self.check_reference(case, expected)

    def test_all_sparse_saved_cases_and_exact_reference_tags(self):
        for name in SPARSE:
            records = cases("sparse", name)
            self.assertEqual(len(records), 2003 if name == "matmul" else 2002)
            for case in records:
                with self.subTest(case=case["id"]):
                    self.check_reference(case, sparse_reference(case))

    def test_dense_singular_and_shape_guards(self):
        for case in (
            {"api": "linalg.solve", "args": {"a": [[1., 2.], [2., 4.]], "b": [1., 2.]}},
            {"api": "linalg.det", "args": {"a": [[1., 2., 3.]]}},
            {"api": "linalg.solve", "args": {"a": [[1., 2.], [0., 1.]], "b": [1.]}},
            {"api": "linalg.dot", "args": {"a": [], "b": [], "u": [1.], "v": [1., 2.], "explicit": True}},
            {"api": "linalg.matvec", "args": {"a": [[1., 2.]], "b": [], "u": [1.], "v": [], "explicit": True}},
            {"api": "linalg.det", "args": {"a": [[float("inf")]]}},
            {"api": "linalg.det", "args": {"a": [[1+2j]]}},
            {"api": "linalg.inverse", "args": {"a": [[1.]]}},
            {"api": "sparse.linalg.cg", "args": {}},
        ):
            self.assertEqual(self.native.evaluate(case), {}, case)

    def test_sparse_coo_duplicates_rectangles_and_rejections(self):
        for name in SPARSE:
            for case in cases("sparse", name)[-1:]:
                self.check_reference(case, sparse_reference(case))
        for case in (
            {"api": "sparse.add", "args": {"raw": {"api": "add", "rows": 2, "cols": 3,
                                                   "a": [[0, 4, 1.]], "brows": 2,
                                                   "bcols": 3, "b": []},
                                             "fields": [], "suite": "sparse", "suite_name": "sparse"}},
            {"api": "sparse.matvec", "args": {"raw": {"api": "matvec", "rows": 2, "cols": 3,
                                                      "a": [[0, 0, 1.]], "x": [1.]},
                                                "fields": [], "suite": "sparse", "suite_name": "sparse"}},
        ):
            self.assertEqual(self.native.evaluate(case), {})

    def test_binary_abi_and_output_capacity(self):
        lib = self.native.library
        for name in (b"linalg.det", b"linalg.solve", b"linalg.dot",
                     b"linalg.matvec", b"linalg.matmul", b"linalg.transpose",
                     b"sparse.add", b"sparse.multiply", b"sparse.matvec",
                     b"sparse.matmul", b"sparse.transpose",
                     b"sparse.from_dense", b"sparse.to_dense"):
            self.assertEqual(lib.boost_ublas_supported(0, name), 2)
            self.assertEqual(lib.boost_ublas_supported(1, name), 0)
        input_values = (ct.c_double * 6)(2., 2., 1., 2., 3., 4.)
        self.assertEqual(lib.boost_ublas_output_count(
            b"linalg.transpose", input_values, 6), 4)
        short = (ct.c_double * 3)()
        self.assertEqual(lib.boost_ublas_evaluate(
            0, b"linalg.transpose", input_values, 6, short, 3), -1)
        self.assertIn(b"capacity", lib.boost_ublas_error())
        self.assertEqual(lib.boost_ublas_evaluate(
            1, b"linalg.transpose", input_values, 6, short, 3), -1)


if __name__ == "__main__":
    unittest.main()
