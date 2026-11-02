"""All saved SciPy sparse cases versus real public GSL sparse routines."""
import ctypes as ct
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np

from sci_testing.case_csv import read_case_file
from sci_testing.gsl_sparse import GslSparseComparisons, NAMES
from sci_testing.sparse_cases import reference as independent_scipy

ROOT = Path(__file__).resolve().parents[1]
COUNTS = {
    "matvec": 2002, "matrix_transpose": 2002,
    "transpose": 2002, "from_dense": 2002, "to_dense": 2002,
    "add": 2002, "matmul": 2003,
}


def cases(name):
    path = ROOT / "cases/generated/sparse" / (name + ".csv.gz")
    return read_case_file(path, "sparse/" + name + ".csv")


class GslSparseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        folder = tempfile.TemporaryDirectory()
        cls.addClassCleanup(folder.cleanup)
        output = Path(folder.name) / "gsl-sparse.so"
        env = os.environ.copy()
        env.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2",
            "-Wall", "-Wextra", "-Werror", "-shared", "-fPIC",
            str(ROOT / "native/gsl_sparse.cpp"), "-lgsl", "-lgslcblas",
            "-o", str(output),
        ], check=True, env=env)
        cls.adapter = GslSparseComparisons(ct.CDLL(str(output)))

    def test_all_saved_independent_scipy_cases(self):
        failures = []
        counts = {}
        for name in sorted(NAMES):
            seen = 0
            for case in cases(name):
                seen += 1
                result = self.adapter.evaluate(case)
                if "GSL" not in result or "error" in result["GSL"]:
                    failures.append((case["id"], result))
                    continue
                actual = np.asarray(result["GSL"]["values"], dtype=float)
                expected = np.asarray(independent_scipy(case), dtype=float)
                if actual.shape != expected.shape or not np.allclose(
                        actual, expected, rtol=2e-10, atol=2e-12):
                    failures.append((case["id"], actual.tolist(), expected.tolist()))
            counts[name] = seen
            self.assertEqual(seen, COUNTS[name], name)
        self.assertFalse(failures,
                         f"{len(failures)} cases fail: {failures[:12]}")
        self.assertEqual(sum(counts.values()), 14015)

    def test_abi_capacity_provider_and_malformed_input(self):
        lib = self.adapter.library
        for name in NAMES:
            self.assertEqual(lib.gsl_sparse_supported(
                1, ("sparse." + name).encode()), 1)
            self.assertEqual(lib.gsl_sparse_supported(
                0, ("sparse." + name).encode()), 0)
        self.assertEqual(lib.gsl_sparse_supported(1,b"sparse.multiply"), 0)
        self.assertEqual(self.adapter.evaluate({"api":"sparse.linalg.cg","args":{}}),{})
        sample = cases("transpose")[-1]
        self.assertIn("values", self.adapter.evaluate(sample)["GSL"])
        invalid_cases = [
            {"api":"sparse.matvec", "rows":2, "cols":3,
             "a":[[0,0,1.]],"x":[1.]},
            {"api":"sparse.add", "rows":2, "cols":3, "a":[],
             "brows":3, "bcols":3, "b":[]},
            {"api":"sparse.matmul", "rows":2, "cols":3,"a":[[0,4,1.]],
             "brows":3,"bcols":4,"b":[]},
            {"api":"sparse.transpose", "rows":2, "cols":3,
             "a":[[0,1,float("nan")]]},
            {"api":"sparse.transpose", "rows":2, "cols":3,
             "a":[[0,1,1.],[0,1,float("inf")]]},
        ]
        for raw in invalid_cases:
            name=raw["api"].removeprefix("sparse.")
            case = {"api":"sparse."+name,
                    "args":{"fields":[],"raw":{**raw,"api":name},
                            "suite":"sparse","suite_name":"sparse"}}
            self.assertEqual(self.adapter.evaluate(case),{},case)
        # Ensure native code rejects insufficient capacity without overruns.
        data=(ct.c_double * 6)(2.,2.,1.,0.,0.,2.)
        short=(ct.c_double * 2)()
        self.assertEqual(lib.gsl_sparse_evaluate(
            1,b"sparse.transpose",data,6,short,2),-1)
        self.assertIn(b"capacity",lib.gsl_sparse_error())

    def test_coo_duplicates_and_explicit_zeros(self):
        for name in ("transpose", "matrix_transpose", "from_dense", "to_dense"):
            case=cases(name)[-1]
            raw=case["args"]["raw"]
            self.assertTrue(any(triple[2] == 0. for triple in raw["a"]),name)
            self.assertNotEqual(len(raw["a"]),
                                len({(t[0],t[1]) for t in raw["a"]}),name)
            expected=independent_scipy(case)
            result=self.adapter.evaluate(case)["GSL"]
            self.assertNotIn("error",result,(name,result))
            np.testing.assert_allclose(result["values"],expected,
                                       rtol=2e-10,atol=2e-12)


if __name__=="__main__":
    unittest.main()
