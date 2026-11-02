"""Audit the saved SciPy vector Bessel-zero conventions via public Boost."""

import csv
import ctypes as ct
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from functools import lru_cache

import numpy as np
from scipy import special

from sci_testing.bessel_zeros import BesselZerosComparisons, SUPPORTED


ROOT = Path(__file__).resolve().parents[1]
SPECS = {
    "jnp_zeros": (2, 1, 2004),
    "ynp_zeros": (2, 1, 2004),
    "y0_zeros": (1, 2, 2001),
    "y1_zeros": (1, 2, 2001),
    "y1p_zeros": (1, 2, 2001),
    "jnjnp_zeros": (1, 4, 2001),
    "jnyn_zeros": (2, 4, 2001),
    "ai_zeros": (1, 4, 2001),
    "bi_zeros": (1, 4, 2001),
}


@lru_cache(maxsize=1024)
def scipy_values(name, args):
    value = getattr(special, name)(*args)
    parts = value if isinstance(value, tuple) else (value,)
    # Only the real-axis Sci facade is being compared. SciPy y0/y1 returned
    # values are dtype complex even though their imaginary parts are zero.
    for part in parts:
        if np.iscomplexobj(part):
            if np.any(np.abs(np.imag(part)) > 1e-13):
                raise AssertionError(f"complex Bessel zeros in {name}")
    return np.concatenate([np.asarray(part).real.ravel() for part in parts])


class BesselZeroTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        temporary = tempfile.TemporaryDirectory()
        cls.addClassCleanup(temporary.cleanup)
        target = Path(temporary.name) / "bessel-zeros.so"
        env = dict(os.environ)
        env.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            "c++", "-std=c++17", "-O2", "-Wall", "-Wextra", "-Werror",
            "-shared", "-fPIC", str(ROOT / "native/bessel_zeros.cpp"),
            "-o", str(target),
        ], env=env, check=True)
        cls.native = BesselZerosComparisons(ct.CDLL(str(target)))

    @staticmethod
    def case(name, args):
        return {"api": f"special.{name}",
                "args": {f"p{i}": v for i, v in enumerate(args)}}

    def checked(self, name, args):
        value = self.native.evaluate(self.case(name, args))
        self.assertEqual(set(value), {"Boost.Math"}, (name,args,value))
        self.assertNotIn("error",value["Boost.Math"], (name,args,value))
        actual = np.array(value["Boost.Math"]["values"])
        self.assertEqual(actual.size, int(args[-1])*SPECS[name][1], (name,args))
        return actual

    def test_all_saved_case_vectors_and_metadata(self):
        for name, (arity, lanes, count) in SPECS.items():
            with (ROOT / "cases/generated/special" / f"{name}.csv").open(
                    newline="") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), count, name)
            for row in rows:
                args = tuple(int(row[f"p{i}"].split(":",1)[1]) for i in range(arity))
                expected = scipy_values(name, args)
                observed = self.checked(name,args)
                with self.subTest(case=row["id"]):
                    self.assertEqual(len(expected),len(observed))
                    # SciPy's Y1' root solver produces ~1.9e-8 error by
                    # the fourth zero; root residuals below certify Boost's
                    # refined answer rather than inheriting that error.
                    tolerance = 1e-7 if name == "y1p_zeros" else 4e-9
                    np.testing.assert_allclose(observed,expected,atol=tolerance,rtol=2e-8)
                    if name == "jnjnp_zeros":
                        nt=int(args[-1])
                        self.assertTrue(np.array_equal(observed[nt:],
                                                       np.rint(observed[nt:])))

    def test_first_positive_derivative_roots_and_root_residuals(self):
        for n in (0,1,2,3,8,20,50,100):
            for name, derivative in (
                ("jnp_zeros",lambda x:special.jvp(n,x,1)),
                ("ynp_zeros",lambda x:special.yvp(n,x,1)),
            ):
                zeros = self.checked(name,(n,6))
                self.assertTrue(np.all(np.diff(zeros)>0))
                self.assertTrue(np.all(zeros>0))
                self.assertLess(np.max(abs(derivative(zeros))),1e-9)
                np.testing.assert_allclose(zeros,scipy_values(name,(n,6)),
                                           atol=2e-9,rtol=2e-11)
                if n == 0:
                    known = special.jn_zeros(1,6) if name == "jnp_zeros" else special.yn_zeros(1,6)
                    np.testing.assert_allclose(zeros,known,atol=2e-11)

    def test_paired_values_and_joint_root_contract(self):
        for name, base, companion in (
            ("y0_zeros",lambda x:special.yv(0,x),lambda x:special.yv(1,x)),
            ("y1_zeros",lambda x:special.yv(1,x),lambda x:special.yvp(1,x,1)),
            ("y1p_zeros",lambda x:special.yvp(1,x,1),lambda x:special.yv(1,x)),
        ):
            values=self.checked(name,(9,))
            zeros, paired=np.split(values,2)
            self.assertTrue(np.all(np.diff(zeros)>0))
            self.assertLess(np.max(abs(base(zeros))),1e-9)
            np.testing.assert_allclose(paired,companion(zeros),atol=1e-11,rtol=1e-10)
        jnyn=self.checked("jnyn_zeros",(3,6)).reshape(4,6)
        for name,values in zip(("jn_zeros","jnp_zeros","yn_zeros","ynp_zeros"),jnyn):
            np.testing.assert_allclose(values,getattr(special,name)(3,6),atol=1e-10)

        merged=self.checked("jnjnp_zeros",(50,)).reshape(4,50)
        zeros,orders,ranks,flags=merged
        self.assertEqual((zeros[0],orders[0],ranks[0],flags[0]),(0.,0.,0.,1.))
        self.assertTrue(np.all(np.diff(zeros)>=0))
        self.assertTrue(np.all(orders==np.rint(orders)))
        self.assertTrue(np.all(ranks==np.rint(ranks)))
        self.assertTrue(np.all(np.isin(flags,(0.,1.))))
        for x,n,k,derivative in zip(zeros[1:],orders[1:],ranks[1:],flags[1:]):
            residual=special.jvp(n,x,1) if derivative else special.jv(n,x)
            self.assertLess(abs(residual),1e-9)
            reference=(special.jnp_zeros(int(n),int(k)) if derivative
                       else special.jn_zeros(int(n),int(k)))
            self.assertAlmostEqual(x,reference[-1],delta=2e-10)
        self.assertAlmostEqual(zeros[4],zeros[5],delta=2e-14)
        self.assertEqual((orders[4],orders[5],flags[4],flags[5]),(1.,0.,0.,1.))

    def test_airy_roots_derivative_roots_and_paired_values(self):
        for name,fn,fprime in (
            ("ai_zeros",lambda x:special.airy(x)[0],lambda x:special.airy(x)[1]),
            ("bi_zeros",lambda x:special.airy(x)[2],lambda x:special.airy(x)[3]),
        ):
            values=self.checked(name,(20,)).reshape(4,20)
            roots,prime_roots,fn_values,prime_values=values
            self.assertTrue(np.all(roots<0))
            self.assertTrue(np.all(prime_roots<0))
            self.assertTrue(np.all(np.diff(roots)<0))
            self.assertTrue(np.all(np.diff(prime_roots)<0))
            self.assertLess(np.max(abs(fn(roots))),1e-10)
            self.assertLess(np.max(abs(fprime(prime_roots))),1e-10)
            np.testing.assert_allclose(fn_values,fn(prime_roots),atol=1e-11)
            np.testing.assert_allclose(prime_values,fprime(roots),atol=1e-11)

    def test_high_rank_and_joint_boundary(self):
        # The merged sequence has duplicate J1 / J0' roots; sorting by
        # near-equal x must preserve all four metadata arrays up to cap 128.
        for name,args,atol in (
            ("jnjnp_zeros",(128,),2e-10),
            ("jnp_zeros",(64,15),2e-10),
            ("ynp_zeros",(64,15),2e-10),
            ("jnyn_zeros",(64,15),2e-10),
            ("ai_zeros",(128,),3e-10),
            ("bi_zeros",(128,),3e-10),
        ):
            with self.subTest(name=name,args=args):
                actual=self.checked(name,args)
                expected=scipy_values(name,args)
                np.testing.assert_allclose(actual,expected,atol=atol,rtol=2e-10)

    def test_native_abi_library_support_and_buffer_lengths(self):
        library=self.native.library
        for name,(arity,lanes,_) in SPECS.items():
            encoded=name.encode()
            self.assertEqual(library.bessel_zeros_supported(0,encoded),arity)
            self.assertEqual(library.bessel_zeros_supported(1,encoded),0)
            args=(ct.c_double*arity)(*(2,3) if arity==2 else (3,))
            self.assertEqual(library.bessel_zeros_output_count(encoded,args,arity),3*lanes)
            insufficient=(ct.c_double*1)()
            status=library.bessel_zeros_evaluate(
                0,encoded,args,arity,insufficient,1)
            if lanes*3>1:
                self.assertNotEqual(status,0)
                self.assertIn(b"capacity",library.bessel_zeros_error())

    def test_unsupported_and_invalid_requests_explicit(self):
        for name in ("jn_zeros","yn_zeros",
                     "kelvin_zeros","erf_zeros","fresnel_zeros"):
            self.assertEqual(self.native.evaluate(self.case(name,(2,3))),{})
        for name,args in (
            ("jnp_zeros",(-1,3)),("jnp_zeros",(1,0)),
            ("ynp_zeros",(257,3)),("ynp_zeros",(2.2,4)),
            ("jnyn_zeros",(2,257)),("jnjnp_zeros",(129,)),
            ("y0_zeros",(0,)),("y1_zeros",(-1,)),("y1p_zeros",(3.5,)),
        ):
            result=self.native.evaluate(self.case(name,args))
            self.assertIn("error",result["Boost.Math"],(name,args,result))
        self.assertIn("error",self.native.evaluate(self.case("jnp_zeros",(2,)))["Boost.Math"])


if __name__ == "__main__":
    unittest.main()
