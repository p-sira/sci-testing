"""Independent GSL-derived full saved-case Bessel/Neumann zero verification."""

import csv
import ctypes as ct
from functools import lru_cache
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
from scipy import special

from sci_testing.gsl_bessel_zeros import GSLBesselZerosComparisons


ROOT = Path(__file__).resolve().parents[1]
SPECS = {
    "jnp_zeros": (2, 1, 2004),
    "yn_zeros": (2, 1, 2004),
    "ynp_zeros": (2, 1, 2004),
    "y0_zeros": (1, 2, 2001),
    "y1_zeros": (1, 2, 2001),
    "y1p_zeros": (1, 2, 2001),
    "jnjnp_zeros": (1, 4, 2001),
    "jnyn_zeros": (2, 4, 2001),
}


@lru_cache(maxsize=1024)
def reference(name, args):
    """Use independent SciPy zero functions, preserving complete tuple layout."""
    value = getattr(special, name)(*args)
    parts = value if isinstance(value, tuple) else (value,)
    for part in parts:
        assert np.max(np.abs(np.asarray(part).imag)) < 1e-12, name
    return np.concatenate([np.asarray(part).real.ravel() for part in parts])


class GSLBesselZerosTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        workspace = tempfile.TemporaryDirectory()
        cls.addClassCleanup(workspace.cleanup)
        target = Path(workspace.name) / "gsl-bessel-zeros.so"
        env = dict(os.environ)
        env.pop("LD_LIBRARY_PATH", None)  # CoS AppImage can shadow compiler.
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2",
            "-Wall", "-Wextra", "-Werror", "-Wpedantic", "-shared", "-fPIC",
            str(ROOT / "native/gsl_bessel_zeros.cpp"), "-lgsl",
            "-lgslcblas", "-o", str(target),
        ], env=env, check=True)
        cls.native = GSLBesselZerosComparisons(ct.CDLL(str(target)))

    @staticmethod
    def case(name, args):
        return {
            "api": "special." + name,
            "args": {f"p{i}": value for i, value in enumerate(args)},
        }

    def result(self, name, args):
        observed = self.native.evaluate(self.case(name, args))
        self.assertEqual(set(observed), {"GSL"}, (name, args, observed))
        self.assertNotIn("error", observed["GSL"], (name,args,observed))
        values = np.asarray(observed["GSL"]["values"], dtype=float)
        self.assertEqual(len(values), int(args[-1])*SPECS[name][1])
        self.assertTrue(np.all(np.isfinite(values)))
        return values

    def test_all_16017_generated_saved_cases(self):
        total = 0
        for name, (arity, _, expected_rows) in SPECS.items():
            with (ROOT / "cases/generated/special" / (name + ".csv")).open(
                    newline="") as source:
                rows = list(csv.DictReader(source))
            self.assertEqual(len(rows), expected_rows, name)
            for row in rows:
                args = tuple(int(row[f"p{i}"].split(":", 1)[1])
                             for i in range(arity))
                actual = self.result(name, args)
                expected = reference(name, args)
                with self.subTest(case=row["id"]):
                    # scipy.special's Y1' zero root cache loses 1e-8 at
                    # some higher ranks; independently assert residuals below.
                    atol = 1e-7 if name == "y1p_zeros" else 4e-9
                    np.testing.assert_allclose(actual, expected,
                                               rtol=2e-8, atol=atol)
                    if name == "jnjnp_zeros":
                        n = args[0]
                        self.assertTrue(np.array_equal(actual[n:],
                                                       np.rint(actual[n:])))
                total += 1
        self.assertEqual(total, 16017)

    def test_values_are_actual_zeros_and_companion_contract(self):
        for n in (0, 1, 2, 3, 8, 20, 50, 64):
            for name, derivative in (
                ("jnp_zeros", special.jvp),
                ("yn_zeros", special.yv),
                ("ynp_zeros", special.yvp),
            ):
                zeros = self.result(name, (n, 6))
                self.assertTrue(np.all(zeros > 0))
                self.assertTrue(np.all(np.diff(zeros) > 0))
                residual = derivative(n, zeros)
                if name != "yn_zeros":
                    residual = derivative(n, zeros, 1)
                self.assertLess(max(abs(residual)), 1e-9, (name,n,residual))
                np.testing.assert_allclose(zeros, reference(name, (n, 6)),
                                           rtol=2e-10, atol=2e-9)
        for name, family, other in (
            ("y0_zeros", lambda x:special.yv(0,x),lambda x:special.yv(1,x)),
            ("y1_zeros", lambda x:special.yv(1,x),lambda x:special.yvp(1,x,1)),
            ("y1p_zeros",lambda x:special.yvp(1,x,1),lambda x:special.yv(1,x)),
        ):
            values = self.result(name,(20,)).reshape(2,20)
            zeros, companions=values
            self.assertLess(max(abs(family(zeros))),1e-9)
            np.testing.assert_allclose(companions,other(zeros),
                                       rtol=2e-10,atol=2e-11)

    def test_global_joint_order_and_four_arrays(self):
        for nt in (1, 16, 50, 128):
            roots, orders, ranks, derivative = self.result(
                "jnjnp_zeros", (nt,)).reshape(4,nt)
            self.assertEqual((roots[0],orders[0],ranks[0],derivative[0]),
                             (0.,0.,0.,1.))
            self.assertTrue(np.all(np.diff(roots) >= 0.))
            for metadata in (orders,ranks,derivative):
                np.testing.assert_array_equal(metadata,np.rint(metadata))
            self.assertTrue(np.all(np.isin(derivative,(0.,1.))))
            for root,n,rank,flag in zip(roots[1:],orders[1:],
                                        ranks[1:],derivative[1:]):
                residual = special.jvp(n,root,1) if flag else special.jv(n,root)
                self.assertLess(abs(residual),1e-9)
                expected = (special.jnp_zeros(int(n),int(rank))[-1] if flag
                            else special.jn_zeros(int(n),int(rank))[-1])
                self.assertAlmostEqual(root,expected,delta=1e-10)
            np.testing.assert_allclose(
                np.ravel([roots,orders,ranks,derivative]),
                reference("jnjnp_zeros",(nt,)),rtol=2e-11,atol=3e-10)
        self.assertEqual(self.result("jnjnp_zeros",(6,)).reshape(4,6)[3,4:6].tolist(),
                         [0.,1.])

    def test_joint_jn_yn_layout(self):
        n,nt=10,15
        values=self.result("jnyn_zeros",(n,nt)).reshape(4,nt)
        for part,name in zip(values,("jn_zeros","jnp_zeros","yn_zeros","ynp_zeros")):
            np.testing.assert_allclose(
                part,getattr(special,name)(n,nt),rtol=2e-11,atol=1e-10)

    def test_abi_guards_and_no_false_library_attribution(self):
        lib = self.native.library
        for name,(arity,lanes,_) in SPECS.items():
            encoded=name.encode()
            self.assertEqual(lib.gsl_bessel_zeros_supported(1,encoded),arity)
            self.assertEqual(lib.gsl_bessel_zeros_supported(0,encoded),0)
            inputs=(ct.c_double*arity)(*(2,3) if arity==2 else (3,))
            self.assertEqual(lib.gsl_bessel_zeros_output_count(
                encoded,inputs,arity),3*lanes)
            output=(ct.c_double*1)(123.)
            status=lib.gsl_bessel_zeros_evaluate(
                1,encoded,inputs,arity,output,1)
            self.assertNotEqual(status,0)
            self.assertEqual(output[0],123.)
            self.assertIn(b"capacity",lib.gsl_bessel_zeros_error())
        self.assertEqual(lib.gsl_bessel_zeros_supported(1,b"kelvin_zeros"),0)
        for name,args in (
            ("jnp_zeros",(-1,3)),("jnp_zeros",(2,0)),
            ("ynp_zeros",(65,3)),("yn_zeros",(2.5,3)),
            ("jnyn_zeros",(2,129)),("jnjnp_zeros",(129,)),
            ("y0_zeros",(0,)),("y1_zeros",(-1,)),("y1p_zeros",(3.5,)),
            ("jnp_zeros",(math.nan,3)),("yn_zeros",(2,math.inf)),
        ):
            observed=self.native.evaluate(self.case(name,args))
            self.assertIn("error",observed["GSL"],(name,args,observed))
        for unsupported in ("jn_zeros","bei_zeros","erf_zeros"):
            self.assertEqual(self.native.evaluate(self.case(unsupported,(3,))),{})
        self.assertIn("error",self.native.evaluate(
            {"api":"special.yn_zeros","args":{"p0":1}})["GSL"])


if __name__ == "__main__":
    unittest.main()
