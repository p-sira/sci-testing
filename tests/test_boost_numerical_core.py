"""Full independent verification of five public Boost.Math numerical oracles.

Every saved fixture is checked against SciPy and/or exact mathematical
derivatives, including the complete vector/matrix output shape. Native
compilation uses no GSL and changes no shared build/dispatch files.
"""

import csv
import ctypes as ct
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
from scipy import differentiate, integrate

from sci_testing.boost_numerical_core import (
    ARITIES, OUTPUT_COUNTS, BoostNumericalCoreComparisons,
)

ROOT = Path(__file__).resolve().parents[1]
COUNTS = {
    "integrate.fixed_quad": 2068,
    "differentiate.derivative": 2201,
    "differentiate.gradient": 2076,
    "differentiate.jacobian": 2000,
    "differentiate.hessian": 2000,
}


def cases(api):
    module, name = api.split(".")
    with (ROOT / "cases/generated" / module / (name + ".csv")).open(
            newline="") as source:
        raw = list(csv.DictReader(source))
    assert len(raw) == COUNTS[api], (api, len(raw))
    result = []
    for row in raw:
        if api == "integrate.fixed_quad":
            args = {"a": float(row["a"][2:]), "b": float(row["b"][2:])}
        elif api == "differentiate.derivative":
            args = {"x": float(row["x"][2:]), "step": float(row["step"][2:])}
        elif api == "differentiate.gradient":
            args = {"point": [float(x) for x in row["point"][2:].split()],
                    "step": float(row["step"][2:])}
        else:
            args = {"point": [float(x) for x in row["point"][2:].split()],
                    "explicit": row["explicit"] == "b:1"}
        result.append({"id": row["id"], "api": api, "args": args})
    return result


class BoostNumericalCoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(tmp.cleanup)
        native_path = Path(tmp.name) / "boost-numerical-core.so"
        env = dict(os.environ)
        env.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2",
            "-Wall", "-Wextra", "-Werror", "-shared", "-fPIC",
            str(ROOT / "native/boost_numerical_core.cpp"),
            "-o", str(native_path),
        ], check=True, env=env)
        cls.adapter = BoostNumericalCoreComparisons(ct.CDLL(str(native_path)))

    def evaluate(self, case):
        result = self.adapter.evaluate(case)
        self.assertEqual(set(result), {"Boost.Math"}, (case["id"], result))
        self.assertNotIn("error", result["Boost.Math"], (case["id"], result))
        values = np.asarray(result["Boost.Math"]["values"], dtype=float)
        self.assertEqual(len(values), OUTPUT_COUNTS[case["api"]])
        self.assertTrue(np.all(np.isfinite(values)))
        return values

    def test_all_saved_fixed_gauss_quadrature(self):
        fixtures = cases("integrate.fixed_quad")
        maximum_relative = maximum_absolute = 0.
        for case in fixtures:
            a, b = case["args"]["a"], case["args"]["b"]
            actual = self.evaluate(case)[0]
            scipy = integrate.fixed_quad(np.exp, a, b, n=16)[0]
            with self.subTest(id=case["id"]):
                absolute = abs(actual-scipy)
                maximum_absolute = max(absolute,maximum_absolute)
                maximum_relative = max(maximum_relative,
                                       absolute/max(abs(scipy),1e-300))
                self.assertTrue(math.isclose(actual,scipy,
                                              rel_tol=3e-13,abs_tol=2e-12))
                if a == b:
                    self.assertEqual(actual,0.)
                if a != b:
                    self.assertEqual(actual>0,b>a)
        self.assertLess(maximum_relative,3e-13)
        self.assertLess(maximum_absolute,2e-12)

    def test_all_saved_sine_derivatives(self):
        fixtures=cases("differentiate.derivative")
        observed=np.asarray([self.evaluate(case)[0] for case in fixtures])
        x=np.asarray([case["args"]["x"] for case in fixtures])
        expected=np.cos(x)
        np.testing.assert_allclose(observed,expected,rtol=8e-12,atol=7e-14)
        # SciPy's independent adaptive method accepts a vector of these x
        # values; verify numerical derivatives and status for every saved x.
        scipy=differentiate.derivative(
            np.sin,x,tolerances={"atol":1e-10,"rtol":1e-10})
        self.assertTrue(np.all(scipy.success))
        np.testing.assert_allclose(observed,scipy.df,rtol=8e-10,atol=3e-10)

    def test_all_saved_gradient_cases(self):
        fixtures=cases("differentiate.gradient")
        points=np.asarray([case["args"]["point"] for case in fixtures])
        observed=np.asarray([self.evaluate(case) for case in fixtures])
        x,y=points.T
        expected=np.column_stack((2*x+3*y,3*x+4*y))
        np.testing.assert_allclose(observed,expected,rtol=1e-14,atol=1e-14)
        scipy=differentiate.jacobian(
            lambda p: np.asarray([p[0]**2+3*p[0]*p[1]+2*p[1]**2]),
            points.T,tolerances={"atol":1e-10,"rtol":1e-10})
        self.assertTrue(np.all(scipy.success))
        np.testing.assert_allclose(observed,scipy.df[0].T,
                                   rtol=5e-9,atol=5e-9)

    def test_all_saved_jacobian_cases(self):
        fixtures=cases("differentiate.jacobian")
        points=np.asarray([case["args"]["point"] for case in fixtures])
        observed=np.asarray([self.evaluate(case) for case in fixtures])
        x,y=points.T
        expected=np.column_stack((y,x,np.ones_like(x),np.ones_like(x)))
        np.testing.assert_allclose(observed,expected,rtol=1e-14,atol=1e-14)
        scipy=differentiate.jacobian(
            lambda p: np.asarray([p[0]*p[1],p[0]+p[1]]),
            points.T,tolerances={"atol":1e-10,"rtol":1e-10})
        self.assertTrue(np.all(scipy.success))
        np.testing.assert_allclose(observed,scipy.df.transpose(2,0,1).reshape(-1,4),
                                   rtol=5e-9,atol=5e-9)

    def test_all_saved_hessian_cases(self):
        fixtures=cases("differentiate.hessian")
        points=np.asarray([case["args"]["point"] for case in fixtures])
        observed=np.asarray([self.evaluate(case) for case in fixtures])
        expected=np.tile([2.,3.,3.,4.],(len(fixtures),1))
        np.testing.assert_allclose(observed,expected,rtol=0,atol=0)
        scipy=differentiate.hessian(
            lambda p: p[0]**2+3*p[0]*p[1]+2*p[1]**2,
            points.T,tolerances={"atol":1e-10,"rtol":1e-10})
        self.assertTrue(np.all(scipy.success))
        np.testing.assert_allclose(observed,scipy.ddf.transpose(2,0,1).reshape(-1,4),
                                   rtol=3e-8,atol=3e-8)

    def test_provider_abi_and_guards(self):
        library=self.adapter.library
        for api, count in ARITIES.items():
            name=api.rsplit(".",1)[-1].encode()
            self.assertEqual(library.boost_numerical_core_supported(0,name),count)
            self.assertEqual(library.boost_numerical_core_supported(1,name),0)
            self.assertEqual(library.boost_numerical_core_supported(-1,name),0)
        self.assertEqual(self.adapter.evaluate(
            {"api":"differentiate.unmapped","args":{"x":1.}}),{})
        bad = [
            ("integrate.fixed_quad",{"a":float("nan"),"b":1.}),
            ("differentiate.derivative",{"x":0.,"step":0.}),
            ("differentiate.derivative",{"x":0.,"step":-1e-3}),
            ("differentiate.gradient",{"point":[1.,2.],"step":0.}),
            ("differentiate.gradient",{"point":[1.],"step":1e-5}),
            ("differentiate.jacobian",{"point":[1.,2.],"explicit":False}),
            ("differentiate.hessian",{"point":[1.,2.,3.],"explicit":True}),
            ("integrate.fixed_quad",{"a":0.,"b":float("inf")}),
        ]
        for api,args in bad:
            with self.subTest(api=api,args=args):
                r=self.adapter.evaluate({"api":api,"args":args})
                self.assertIn("error",r["Boost.Math"])
        args=(ct.c_double*2)(0.,1.)
        out=(ct.c_double*1)()
        self.assertEqual(
            library.boost_numerical_core_evaluate(
                0,b"fixed_quad",args,2,out,0),-1)
        self.assertIn(b"capacity",library.boost_numerical_core_error())
        self.assertEqual(
            library.boost_numerical_core_evaluate(
                1,b"fixed_quad",args,2,out,1),-1)
        self.assertEqual(
            library.boost_numerical_core_output_count(
                b"fixed_quad",args,1),0)


if __name__=="__main__":
    unittest.main()
