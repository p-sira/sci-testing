"""Full independent SciPy validation of saved GSL linalg cases and native ABI."""
import csv
import ctypes as ct
import gzip
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

import mpmath as mp
import numpy as np
from scipy import linalg

from sci_testing.gsl_linalg import NAMES, GslLinalgComparisons

ROOT = Path(__file__).resolve().parents[1]
COUNTS = {
    "det": 2058, "solve": 2058, "inv": 2003, "dot": 2000,
    "matvec": 2000, "matmul": 2000, "cholesky": 2004,
    "cho_solve": 2004, "eigvalsh": 2006, "svdvals": 2005,
}


def _vec(value):
    assert value.startswith("v:"), value
    return [float(n) for n in value[2:].split()]


def _mat(value):
    assert value.startswith("m:"), value
    return [[float(n) for n in row.split()] for row in value[2:].split(";")]


def cases(name):
    path = ROOT / "cases" / "generated" / "linalg" / f"{name}.csv.gz"
    with gzip.open(path, mode="rt", newline="") as stream:
        for row in csv.DictReader(stream):
            if "raw" in row:
                args = {"raw": json.loads(row["raw"].removeprefix("j:"))}
            else:
                args = {"a": _mat(row["a"])}
                if "b" in row:
                    args["b"] = _mat(row["b"]) if row["b"].startswith("m:") else _vec(row["b"])
                if "u" in row:
                    args["u"] = _vec(row["u"])
                if "v" in row:
                    args["v"] = _vec(row["v"])
            yield row["id"], {"api": f"linalg.{name}", "args": args}


def independent(name, args):
    if "raw" in args:
        saved = args["raw"]["args"]
        a = np.asarray(saved[0], dtype=float)
        if name == "cholesky":
            return linalg.cholesky(a, lower=saved[-1]).ravel()
        if name == "cho_solve":
            factor = linalg.cho_factor(a, lower=saved[-1])
            return linalg.cho_solve(factor, saved[1]).ravel()
        if name == "eigvalsh":
            return linalg.eigvalsh(a, lower=saved[-1]).ravel()
        if name == "svdvals":
            return linalg.svdvals(a).ravel()
        return linalg.inv(a).ravel()
    a = np.asarray(args["a"], dtype=float)
    if name == "det":
        return np.atleast_1d(linalg.det(a))
    if name == "solve":
        return linalg.solve(a, args["b"]).ravel()
    if name == "dot":
        return np.atleast_1d(np.dot(args["u"], args["v"]))
    if name == "matvec":
        return (a @ np.asarray(args["u"])).ravel()
    if name == "matmul":
        return (a @ np.asarray(args["b"])).ravel()
    raise ValueError(name)


class GslLinalgTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(tmp.cleanup)
        target = Path(tmp.name) / "gsl-linalg.so"
        env = os.environ.copy()
        env.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2", "-fPIC",
            "-shared", str(ROOT / "native/gsl_linalg.cpp"),
            "-lgsl", "-lgslcblas", "-o", str(target),
        ], check=True, env=env)
        cls.adapter = GslLinalgComparisons(ct.CDLL(str(target)))

    def test_full_saved_scipy(self):
        failures = []
        measured = {}
        illconditioned = 0
        for name, count in COUNTS.items():
            seen = 0
            for identifier, case in cases(name):
                seen += 1
                result = self.adapter.evaluate(case)
                if "GSL" not in result or "error" in result["GSL"]:
                    failures.append((identifier, result))
                    continue
                actual = np.asarray(result["GSL"]["values"])
                expected = np.asarray(independent(name, case["args"]), dtype=float)
                if name == "solve":
                    a = np.asarray(case["args"]["a"], dtype=float)
                    if np.linalg.cond(a) > 1e8:
                        # The saved 8x8 Hilbert matrix amplifies binary64 LU
                        # rounding by ~1e10. Against exact input floats, GSL
                        # is more accurate than SciPy here: compare against an
                        # independent 110-decimal-digit LU solution instead
                        # of forcing agreement with SciPy's less accurate LU.
                        illconditioned += 1
                        with mp.workdps(110):
                            high_a = mp.matrix([[mp.mpf(float(v)) for v in row]
                                                for row in a])
                            high_b = mp.matrix([mp.mpf(float(v))
                                                for v in case["args"]["b"]])
                            reference = np.array([float(v) for v in
                                                  mp.lu_solve(high_a, high_b)])
                        # This Hilbert solve has condition number ~1.5e10.
                        # GSL's binary64 LU differs by ~1e-7 relative across
                        # platforms even with a ~1e-17 backward residual.
                        if actual.shape != reference.shape or not np.allclose(
                                actual, reference, atol=2e-11, rtol=2e-7):
                            failures.append((identifier, actual.tolist(),
                                             reference.tolist(), "high-precision LU"))
                        continue
                if actual.shape != expected.shape or not np.allclose(
                        actual, expected, atol=2e-11, rtol=3e-10):
                    failures.append((identifier, actual.tolist(), expected.tolist()))
            measured[name] = seen
            self.assertEqual(seen, count, name)
        self.assertFalse(failures,
                         f"{len(failures)} saved-case mismatches: {failures[:20]}")
        self.assertEqual(sum(measured.values()), 20138)
        self.assertEqual(illconditioned, 1, "unexpected high-condition cases")

    def test_abi_errors_and_unsupported(self):
        library = self.adapter.library
        for name in NAMES:
            self.assertEqual(library.gsl_linalg_supported(1,name.encode()), 1)
            self.assertEqual(library.gsl_linalg_supported(0,name.encode()), 0)
        self.assertEqual(library.gsl_linalg_supported(1,b"unsupported"), 0)
        self.assertEqual(self.adapter.evaluate({"api": "linalg.unsupported",
                                               "args": {"a": [[1.]]}}), {})
        self.assertEqual(self.adapter.evaluate({"api": "stats.det",
                                               "args": {"a": [[1.]]}}), {})
        # Native should reject singular matrices, rather than abort on GSL's
        # default error handler.
        invalid = self.adapter.evaluate({
            "api": "linalg.inv",
            "args": {"raw": {"args": [[[1.,2.],[2.,4.]]]}},
        })
        self.assertIn("error", invalid["GSL"], invalid)
        # The Python schema rejects non-rectangular and nonfinite input.
        self.assertEqual(self.adapter.evaluate({
            "api": "linalg.det", "args": {"a": [[1.], [2.,3.]]}}), {})
        self.assertEqual(self.adapter.evaluate({
            "api": "linalg.det", "args": {"a": [[float("nan")]]}}), {})


if __name__ == "__main__":
    unittest.main()
