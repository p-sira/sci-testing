"""GSL 2.7/2.8 public direct and derived special-function adapters.

All saved 20,015 SciPy fixture inputs are checked, including the full four-
lane Airy zero collections, variable-length multinomial counts, signed
quantiles and independent mathematical/forward-function identities.
"""
import csv
import ctypes as ct
import math
import os
import subprocess
import tempfile
import unittest
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy import special

from sci_testing.gsl_special_derived import GslSpecialDerivedComparisons, SUPPORTED


ROOT = Path(__file__).resolve().parents[1]
NAMES = (
    "ai_zeros", "bi_zeros", "erfinv", "erfcinv", "erfi",
    "erfcx", "gdtria", "elliprg", "ln_multinomial", "log_ndtr",
)
ARITY = {
    "ai_zeros": 1, "bi_zeros": 1, "erfinv": 1, "erfcinv": 1,
    "erfi": 1, "erfcx": 1, "gdtria": 3, "elliprg": 3,
    "ln_multinomial": 2, "log_ndtr": 1,
}


@lru_cache(maxsize=512)
def reference(name, args):
    if name == "ln_multinomial":
        n, counts = args
        return np.array([special.gammaln(n + 1) -
                         sum(special.gammaln(k+1) for k in counts)])
    observed = getattr(special, name)(*args)
    if isinstance(observed, tuple):
        return np.concatenate([np.asarray(item).ravel() for item in observed])
    return np.asarray(observed).reshape(-1)


class GslSpecialDerivedTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        temporary = tempfile.TemporaryDirectory()
        cls.addClassCleanup(temporary.cleanup)
        path = Path(temporary.name) / "gsl-special-derived.so"
        env = dict(os.environ)
        env.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2",
            "-Wall", "-Wextra", "-Werror", "-fPIC", "-shared",
            str(ROOT / "native/gsl_special_derived.cpp"),
            "-lgsl", "-lgslcblas", "-o", str(path),
        ], check=True, env=env)
        cls.native = GslSpecialDerivedComparisons(ct.CDLL(str(path)))

    @staticmethod
    def case(name, args):
        return {"api": "special." + name,
                "args": {f"p{i}": v for i, v in enumerate(args)}}

    def evaluated(self, name, args):
        case = self.case(name, args)
        result = self.native.evaluate(case)
        self.assertEqual(set(result), {"GSL"}, (name,args,result))
        self.assertNotIn("error", result["GSL"], (name,args,result))
        return np.array(result["GSL"]["values"])

    def test_all_saved_20015_cases(self):
        for name in NAMES:
            with (ROOT / "cases/generated/special" / (name + ".csv")).open(
                    newline="") as source:
                rows = list(csv.DictReader(source))
            self.assertEqual(len(rows), 2006 if name == "ln_multinomial" else 2001)
            for row in rows:
                values = []
                for i in range(ARITY[name]):
                    raw = row[f"p{i}"]
                    kind, value = raw.split(":", 1)
                    if kind == "i": values.append(int(value))
                    elif kind == "s": values.append(float(value))
                    elif kind == "v": values.append(tuple(int(v) for v in value.split()))
                    else: self.fail(f"unknown fixture encoding {raw}")
                args = tuple(values)
                observed = self.evaluated(name, args)
                expected = reference(name, args)
                with self.subTest(row=row["id"]):
                    self.assertEqual(observed.size, expected.size)
                    tolerance = 3e-10 if name in ("ai_zeros","bi_zeros") else 2e-11
                    np.testing.assert_allclose(observed, expected,
                                               rtol=3e-10, atol=tolerance)
                    if name in ("ai_zeros","bi_zeros"):
                        self.assertEqual(observed.size, 4*args[0])

    def test_airy_roots_and_companion_identity(self):
        for name in ("ai_zeros", "bi_zeros"):
            first = name == "ai_zeros"
            parts = self.evaluated(name, (64,)).reshape(4,64)
            zeros, derivative_zeros, function_at_derivative, derivative_at_zeros = parts
            self.assertTrue(np.all(zeros<0))
            self.assertTrue(np.all(derivative_zeros<0))
            self.assertTrue(np.all(np.diff(zeros)<0))
            self.assertTrue(np.all(np.diff(derivative_zeros)<0))
            for i in range(64):
                with self.subTest(api=name,rank=i+1):
                    f = special.airy(zeros[i])
                    df = special.airy(derivative_zeros[i])
                    value_index = 0 if first else 2
                    derivative_index = 1 if first else 3
                    self.assertLess(abs(f[value_index]), 2e-7)
                    self.assertLess(abs(df[derivative_index]),2e-7)
                    self.assertAlmostEqual(function_at_derivative[i],
                                           df[value_index],delta=2e-9)
                    self.assertAlmostEqual(derivative_at_zeros[i],
                                           f[derivative_index],delta=2e-9)

    def test_airy_full_supported_rank_and_inverse_tails(self):
        for name in ("ai_zeros","bi_zeros"):
            with self.subTest(name=name):
                values=self.evaluated(name,(256,))
                np.testing.assert_allclose(values,reference(name,(256,)),
                                           rtol=2e-10,atol=4e-10)
        for name,number in (
            ("erfinv",-.99999999999999),
            ("erfinv",.99999999999999),
            ("erfcinv",1e-14),
            ("erfcinv",2-1e-14),
            ("log_ndtr",30.),
            ("erfcx",25.),
            ("erfi",-25.),
        ):
            with self.subTest(name=name,x=number):
                observed=self.evaluated(name,(number,))
                np.testing.assert_allclose(observed,reference(name,(number,)),
                                           rtol=2e-10,atol=1e-11)

    def test_derivation_identities_and_zero_inputs(self):
        for x in (-6.,-2.,-.5,0.,.3,3.,6.):
            with self.subTest(x=x):
                if -1<x<1:
                    inverse=self.evaluated("erfinv",(x,))[0]
                    self.assertAlmostEqual(special.erf(inverse),x,delta=1e-13)
                if 0<x<2:
                    inverse=self.evaluated("erfcinv",(x,))[0]
                    self.assertAlmostEqual(special.erfc(inverse),x,delta=1e-13)
                for name in ("erfi","erfcx","log_ndtr"):
                    observed=self.evaluated(name,(x,))[0]
                    expected=reference(name,(x,))[0]
                    np.testing.assert_allclose(observed,expected,
                                               rtol=2e-13,atol=3e-11)
        for p,b,x in ((.01,.5,.3),(.3,2.5,.5),(.99,8.0,20.)):
            output=self.evaluated("gdtria",(p,b,x))[0]
            self.assertAlmostEqual(float(special.gdtr(output,b,x)),p,delta=1e-12)
        for args in ((0.,0.,0.),(0.,0.,5.),(0.,2.,3.),
                     (2.5,1.5,.5),(1.,1.,1.)):
            with self.subTest(elliprg=args):
                np.testing.assert_allclose(self.evaluated("elliprg",args),
                                           reference("elliprg",args),atol=2e-12)
        for total,counts in ((0,()),(0,(0,0)),(4,(2,1,1)),
                             (20,(5,5,5,5)),(100,(50,50))):
            np.testing.assert_allclose(self.evaluated("ln_multinomial",(total,counts)),
                                       reference("ln_multinomial",(total,counts)),atol=1e-11)

    def test_errors_library_ownership_and_buffer_capacity(self):
        bridge=self.native.library
        for name in NAMES:
            arity=ARITY[name]
            self.assertEqual(bridge.gsl_special_derived_supported(1,name.encode()),arity)
            self.assertEqual(bridge.gsl_special_derived_supported(0,name.encode()),0)
        self.assertEqual(self.native.evaluate(self.case("does_not_exist",(3,))),{})
        for name,args in (
            ("ai_zeros",(0,)),("bi_zeros",(257,)),("ai_zeros",(2.5,)),
            ("erfinv",(-1.5,)),("erfinv",(1.,)),("erfcinv",(0.,)),
            ("erfcinv",(2.5,)),("erfi",(27.,)),("erfcx",(-27.,)),
            ("gdtria",(-.1,2.,.3)),("gdtria",(.4,2.,0.)),
            ("elliprg",(-1.,1.,1.)),("ln_multinomial",(4,(2,1))),
            ("ln_multinomial",(3,(1,-1,3))),("log_ndtr",(float("nan"),)),
        ):
            result=self.native.evaluate(self.case(name,args))
            with self.subTest(name=name,args=args):
                self.assertIn("error",result["GSL"])
        raw=(ct.c_double*1)(3.)
        too_small=(ct.c_double*1)()
        status=bridge.gsl_special_derived_evaluate(
            1,b"ai_zeros",raw,1,too_small,1)
        self.assertNotEqual(status,0)
        self.assertIn(b"capacity",bridge.gsl_special_derived_error())


if __name__ == "__main__":
    unittest.main()
