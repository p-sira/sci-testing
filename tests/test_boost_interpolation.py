"""Independent SciPy numerical tests for public Boost.Math interpolators."""
import ctypes as ct
import math
import os
from pathlib import Path
import subprocess
import tempfile
import time
import unittest

from scipy import interpolate

from sci_testing.boost_interpolation import (
    BoostInterpolationComparisons, SUPPORTED,
)
from sci_testing.case_csv import read_case_file

ROOT = Path(__file__).resolve().parent.parent


class BoostInterpolationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tempdir = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tempdir.cleanup)
        libpath = Path(cls.tempdir.name) / "boost_interpolation.so"
        env = {k: v for k, v in os.environ.items() if k != "LD_LIBRARY_PATH"}
        subprocess.run([
            os.environ.get("CXX", "/usr/bin/clang++"), "-std=c++17", "-O2",
            "-fPIC", "-shared", str(ROOT / "native/boost_interpolation.cpp"),
            "-o", str(libpath)], check=True, env=env)
        cls.adapter = BoostInterpolationComparisons(ct.CDLL(str(libpath)))

    @staticmethod
    def make_case(name, nodes, values, query):
        return {"api": "interpolate." + name,
                "args": {"nodes": nodes, "values": values, "query": query}}

    @staticmethod
    def scipy_values(case):
        name = case["api"].rsplit(".", 1)[-1]
        args = case["args"]
        return list(getattr(interpolate, name)(
            args["nodes"], args["values"], args["query"]))

    def assert_matches_scipy(self, case, *, rtol=2e-10, atol=3e-10):
        result = self.adapter.evaluate(case)
        self.assertEqual(set(result), {"Boost.Math"}, case["api"])
        self.assertNotIn("error", result["Boost.Math"], (case["api"], result))
        actual = result["Boost.Math"]["values"]
        expected = self.scipy_values(case)
        self.assertEqual(len(actual), len(expected))
        for i, (got, ref) in enumerate(zip(actual, expected)):
            self.assertTrue(math.isclose(got, float(ref),
                rel_tol=rtol, abs_tol=atol), (case["id"] if "id" in case else case["api"],
                                              i, got, float(ref)))

    def test_all_saved_cases(self):
        count = 0
        for name in sorted(SUPPORTED):
            cases = read_case_file(
                ROOT / "cases/generated/interpolate" / (name + ".csv.gz"),
                "interpolate/" + name + ".csv",
            )
            self.assertEqual(len(cases), 2000 if name == "pchip_interpolate" else 2050)
            for case in cases:
                self.assert_matches_scipy(case)
                count += 1
        self.assertEqual(count, 6100)

    def test_pchip_endpoint_shape_and_extrapolation_guard(self):
        x=[-2.,-1.,0.,1.,2.]
        y=[1.,1.9,3.,2.,1.]
        for points in ([x[0],x[-1]],[-1.7,-.75,.75,1.7],[-1.,0.,1.]):
            self.assert_matches_scipy(self.make_case(
                "pchip_interpolate",x,y,points),atol=2e-12)
        for query in ([-2.00001], [2.001], [-3.,0.]):
            response=self.adapter.evaluate(self.make_case(
                "pchip_interpolate",x,y,query))
            self.assertIn("error",response["Boost.Math"])
            self.assertIn("extrapolate",response["Boost.Math"]["error"])

    def test_pchip_zero_secants_and_clipped_boundaries(self):
        x=[-3.,-2.25,-.5,.25,2.,3.]
        datasets=[
            [0.,0.,3.,3.,4.,4.],
            [0.,2.,2.,2.,0.,-5.],
            [-4.,-2.,1.,3.,3.,3.],
            [0.,100.,-50.,-49.,-48.,30.],
            [5.,5.,5.,5.,5.,5.],
        ]
        for y in datasets:
            for q in ([x[0],x[1],x[-2],x[-1]],[-1.8,0.,1.,2.5]):
                with self.subTest(y=y,q=q):
                    self.assert_matches_scipy(
                        self.make_case("pchip_interpolate", x,y,q),
                        rtol=2e-12,atol=1e-12,
                    )

    def test_polynomial_and_distinct_krogh_extrapolation(self):
        x = [-3.,-1.5,-.25,.5,2.1,3.5]
        y = [v**4-2*v+3 for v in x]
        q = [-3.5,-3.,-2.125,-1.,.0,2.,3.5,4.1]
        for name in ("barycentric_interpolate","krogh_interpolate"):
            self.assert_matches_scipy(self.make_case(name,x,y,q))

    def test_guard_duplicate_nonascending_and_invalid_shapes(self):
        for name in sorted(SUPPORTED):
            for invalid_nodes in (
                [1.,2.,2.,4.,5.], [0.,2.,1.,3.,4.],
            ):
                response=self.adapter.evaluate(self.make_case(
                    name,invalid_nodes,[1.,2.,3.,4.,5.],[1.]))
                self.assertIn("error",response["Boost.Math"])
            for invalid in (
                {"api": "interpolate."+name,
                 "args": {"nodes":[0.,1.],"values":[2.],"query":[.5]}},
                self.make_case(name,[0.,1.,2.],[[1.],[2.],[3.]],[.5]),
                self.make_case(name,[0.,1.,2.],[1.,2.,3.],[.5+.3j]),
                self.make_case(name,[0.,1.,float("nan")],[1.,2.,3.],[.5]),
                self.make_case(name,[0.,1.,2.],[1.,2.,3.],[[.5]]),
            ):
                self.assertEqual(self.adapter.evaluate(invalid), {})
        self.assertEqual(self.adapter.evaluate(
            self.make_case("spline_interpolate",[0.,1.],[0.,1.],[.5])), {})

    def test_c_abi_capacity_and_errors(self):
        lib=self.adapter.library
        x=(ct.c_double*5)(-2.,-1.,0.,1.,2.)
        y=(ct.c_double*5)(1.,1.9,3.,2.,1.)
        q=(ct.c_double*2)(-1.,0.)
        for name in sorted(SUPPORTED):
            key=name.encode()
            self.assertEqual(lib.boost_interpolation_supported(key),1)
            out=(ct.c_double*2)(-99.,-99.)
            self.assertEqual(lib.boost_interpolation_evaluate(
                key,x,y,5,q,2,out,1),-1)
            self.assertIn(b"capacity",lib.boost_interpolation_error())
            self.assertEqual(list(out),[-99.,-99.])
            self.assertEqual(lib.boost_interpolation_evaluate(
                key,x,y,5,q,2,None,2),-1)
            self.assertEqual(lib.boost_interpolation_evaluate(
                key,x,y,5,None,2,out,2),-1)
            self.assertEqual(lib.boost_interpolation_evaluate(
                key,None,y,5,q,2,out,2),-1)
            self.assertEqual(list(out),[-99.,-99.])
        self.assertEqual(lib.boost_interpolation_supported(b"interp1d"),0)
        # Bounds are explicit; no accidental large-order polynomial claim.
        many=(ct.c_double*17)(*range(17))
        out=(ct.c_double*2)(-99.,-99.)
        self.assertEqual(lib.boost_interpolation_evaluate(
            b"krogh_interpolate",many,many,17,q,2,out,2),-1)
        self.assertEqual(list(out),[-99.,-99.])

    def test_2000_case_speed(self):
        cases=read_case_file(ROOT / "cases/generated/interpolate/pchip_interpolate.csv.gz",
                             "interpolate/pchip_interpolate.csv")
        start=time.perf_counter()
        for case in cases:
            actual=self.adapter.evaluate(case)
            self.assertIn("values",actual["Boost.Math"])
        self.assertLess(time.perf_counter()-start,15.)


if __name__=="__main__":
    unittest.main()
