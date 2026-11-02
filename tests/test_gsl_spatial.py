"""Independent saved SciPy fixtures and GSL-derived spatial BLAS ABI tests.

GSL has no matching *named* spatial.distance functions; all implementation
paths are explicitly grounded in public gsl_blas or gsl_stats primitives.
"""
import ctypes as ct
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
import warnings

import numpy as np
from scipy.spatial import distance

from sci_testing.case_csv import read_case_file
from sci_testing.gsl_spatial import GSLSpatialComparisons, NAMES

ROOT=Path(__file__).resolve().parents[1]
COUNTS={
    "euclidean":2064,
    "sqeuclidean":2004,
    "cityblock":2004,
    "chebyshev":2003,
    "cosine":2064,
    "correlation":2003,
}


def vectors(case):
    args=case["args"]
    if "raw" not in args:
        return np.asarray(args["u"],dtype=float),np.asarray(args["v"],dtype=float),None
    raw=args["raw"]
    convert=lambda s: np.asarray([float(x) for x in s.split(";")],dtype=float
                                 ) if s else np.asarray([],dtype=float)
    return convert(raw["u"]),convert(raw["v"]),convert(raw["w"]) if raw["w"] else None


class GSLSpatialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        workspace=tempfile.TemporaryDirectory()
        cls.addClassCleanup(workspace.cleanup)
        path=Path(workspace.name)/"gsl_spatial.so"
        env=os.environ.copy()
        env.pop("LD_LIBRARY_PATH",None)
        subprocess.run([
            os.environ.get("CXX","c++"),"-std=c++17","-O2",
            "-Wall","-Wextra","-Werror","-Wpedantic","-fPIC","-shared",
            str(ROOT/"native/gsl_spatial.cpp"),"-lgsl","-lgslcblas",
            "-o",str(path),
        ],env=env,check=True)
        cls.adapter=GSLSpatialComparisons(ct.CDLL(str(path)))

    @staticmethod
    def case(name,u,v,w=None):
        raw={"u":";".join(map(str,u)),"v":";".join(map(str,v)),
             "w":";".join(map(str,w)) if w is not None else ""}
        return {"api":"spatial."+name,"args":{"raw":raw}}

    def evaluated(self,case):
        answer=self.adapter.evaluate(case)
        self.assertEqual(set(answer),{"GSL"},case)
        self.assertNotIn("error",answer["GSL"],(case,answer))
        values=answer["GSL"]["values"]
        self.assertEqual(len(values),1)
        self.assertTrue(math.isfinite(values[0]))
        return values[0]

    def test_every_12142_saved_spatial_case(self):
        confirmed=0
        metrics={}
        for name in sorted(NAMES):
            path=ROOT/"cases/generated/spatial"/(name+".csv.gz")
            fixtures=read_case_file(path,"spatial/"+name+".csv")
            self.assertEqual(len(fixtures),COUNTS[name])
            maximum_absolute=maximum_relative=0.
            for case in fixtures:
                u,v,w=vectors(case)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore",RuntimeWarning)
                    reference=float(getattr(distance,name)(u,v,w=w))
                observed=self.evaluated(case)
                absolute=abs(observed-reference)
                relative=absolute/max(abs(reference),1e-300)
                maximum_absolute=max(maximum_absolute,absolute)
                maximum_relative=max(maximum_relative,relative)
                with self.subTest(case=case["id"]):
                    self.assertTrue(math.isclose(observed,reference,
                                                  rel_tol=3e-12,abs_tol=3e-13),
                                    (case["id"],observed,reference))
                confirmed+=1
            metrics[name]=(maximum_absolute,maximum_relative)
        self.assertEqual(confirmed,12142)
        for name,(_,relative) in metrics.items():
            self.assertLess(relative,5e-12,(name,metrics))

    def test_all_vector_sizes_and_custom_weighted_cases(self):
        for n in (0,1,2,3,5,8,32,127,256):
            u=np.linspace(-2.,2.,n).tolist()
            v=(np.linspace(3.,-1.,n)+.1).tolist()
            w=(np.arange(n,dtype=float)%4).tolist()
            for name in NAMES:
                if n==0 and name in ("chebyshev","cosine","correlation"):
                    continue
                if name=="correlation" and n<3:
                    continue
                if name=="cosine" and n==1 and u[0]==0:
                    continue
                for weights in (None,w):
                    if weights is not None and not any(weights) and name in (
                            "cosine","correlation"):
                        continue
                    with self.subTest(name=name,n=n,weighted=weights is not None):
                        case=self.case(name,u,v,weights)
                        actual=self.evaluated(case)
                        expected=float(getattr(distance,name)(u,v,w=weights))
                        self.assertTrue(math.isclose(actual,expected,
                                                     rel_tol=3e-11,abs_tol=3e-12),
                                        (name,n,actual,expected))

    def test_weight_rules_and_scipy_cosine_overflow_case(self):
        u=[-2.,1.,4.]
        v=[1.,-3.,2.]
        for weights in ([.5,2.,1.],[1.,0.,3.],[0.,0.,0.]):
            for name in NAMES:
                if name in ("cosine","correlation") and not any(weights):
                    continue
                with self.subTest(name=name,weights=weights):
                    actual=self.evaluated(self.case(name,u,v,weights))
                    expected=float(getattr(distance,name)(u,v,w=weights))
                    self.assertAlmostEqual(actual,expected,delta=2e-12)
        # SciPy's sqrt(dot(u,u)*dot(v,v)) overflows although dnrm2 stays
        # finite: reproduce the saved SciPy 1.0 convention, not 1.54554.
        u=[1e100,2e100,-3e100]
        v=[-1e100,1e100,2e100]
        self.assertEqual(self.evaluated(self.case("cosine",u,v)),1.)

    def test_unsupported_wrong_shapes_and_invalid_inputs(self):
        self.assertEqual(self.adapter.evaluate(
            {"api":"spatial.jaccard","args":{"u":[1],"v":[2]}}),{})
        for name in NAMES:
            self.assertEqual(self.adapter.library.gsl_spatial_supported(
                1,name.encode()),1)
            self.assertEqual(self.adapter.library.gsl_spatial_supported(
                0,name.encode()),0)
        cases=[
            ("euclidean",[1.,2.],[2.],None),
            ("sqeuclidean",[1.,2.],[2.,3.],[1.,2.,3.]),
            ("cityblock",[1.,2.],[3.,4.],[-1.,2.]),
            ("chebyshev",[],[],None),
            ("cosine",[0.,0.],[1.,2.],None),
            ("correlation",[1.,1.,1.],[2.,1.,0.],None),
            ("correlation",[1.,2.,3.],[2.,3.,4.],[0.,0.,0.]),
            ("euclidean",[float("nan"),2.],[2.,3.],None),
        ]
        for name,u,v,w in cases:
            with self.subTest(name=name,u=u,v=v,w=w):
                response=self.adapter.evaluate(self.case(name,u,v,w))
                self.assertIn("error",response["GSL"])

    def test_native_capacity_provider_and_error_integrity(self):
        lib=self.adapter.library
        u=(ct.c_double*3)(1.,2.,3.)
        v=(ct.c_double*3)(2.,4.,1.)
        out=(ct.c_double*2)(987.,654.)
        self.assertEqual(lib.gsl_spatial_output_count(b"euclidean",3),1)
        self.assertEqual(lib.gsl_spatial_output_count(b"euclidean",-1),0)
        self.assertEqual(lib.gsl_spatial_output_count(b"unknown",3),0)
        for library,capacity in ((1,0),(0,1)):
            self.assertEqual(lib.gsl_spatial_evaluate(
                library,b"euclidean",u,v,3,None,0,out,capacity),-1)
            self.assertEqual(list(out),[987.,654.])
            self.assertTrue(lib.gsl_spatial_error())
        self.assertEqual(lib.gsl_spatial_evaluate(
            1,b"euclidean",u,v,3,None,0,out,1),0)
        self.assertAlmostEqual(out[0],3.)
        self.assertEqual(out[1],654.)
        self.assertEqual(lib.gsl_spatial_error(),b"")


if __name__=="__main__":
    unittest.main()
