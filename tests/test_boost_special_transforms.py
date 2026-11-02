"""Independent tests for Boost-derived real scalar transforms.

Checks every saved SciPy case (11 functions x 2,001 cases), separate closed
identities, guard behavior and the standalone native C ABI. Does not link Sci
implementations or regenerate reference reports.
"""

import csv
import ctypes as ct
import math
import os
from pathlib import Path
import subprocess
import tempfile
import time
import unittest

import numpy as np
from scipy import special

from sci_testing.boost_special_transforms import BoostSpecialTransformsComparisons, NAMES


ROOT = Path(__file__).resolve().parents[1]
ARITIES = {
    "multigammaln": 2, "boxcox": 2, "boxcox1p": 2,
    "softplus": 1, "log_expit": 1, "nrdtrimn": 3, "nrdtrisd": 3,
    "ndtri_exp": 1, "log_ndtr": 1, "factorialk": 2, "agm": 2,
}


class BoostSpecialTransformsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(temp.cleanup)
        library_path = Path(temp.name) / "boost-special-transforms.so"
        env = os.environ.copy()
        env.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2", "-Wall", "-Wextra",
            "-Werror", "-shared", "-fPIC",
            str(ROOT / "native/boost_special_transforms.cpp"), "-o", str(library_path),
        ], check=True, env=env)
        cls.adapter = BoostSpecialTransformsComparisons(ct.CDLL(str(library_path)))

    @staticmethod
    def case(name, *args):
        return {"api": "special."+name,
                "args": {f"p{i}": float(arg) for i, arg in enumerate(args)}}

    def value(self, name, *args):
        result = self.adapter.evaluate(self.case(name, *args))
        self.assertEqual(set(result), {"Boost.Math"}, (name, args, result))
        entry = result["Boost.Math"]
        self.assertNotIn("error", entry, (name, args, entry))
        self.assertEqual(len(entry["values"]), 1)
        self.assertTrue(math.isfinite(entry["values"][0]))
        return entry["values"][0]

    def test_registration_and_wrapper_contract(self):
        self.assertEqual(set(ARITIES), NAMES)
        for name, arity in ARITIES.items():
            self.assertEqual(self.adapter.library.boost_special_transforms_supported(
                0, name.encode()), arity)
            self.assertEqual(self.adapter.library.boost_special_transforms_supported(
                1, name.encode()), 0)
        self.assertEqual(self.adapter.evaluate(self.case("not_an_api", 2.)), {})
        self.assertEqual(self.adapter.evaluate({"api": "stats.softplus",
                                               "args": {"p0": 2.}}), {})
        self.assertEqual(self.adapter.evaluate({"api": "special.softplus",
                                               "args": {}}), {})
        self.assertEqual(self.adapter.evaluate({"api": "special.softplus",
                                               "args": {"p0": "oops"}}), {})

    def test_all_22011_saved_SciPy_cases(self):
        count = 0
        greatest_scaled, greatest_absolute = 0., 0.
        started = time.perf_counter()
        for name, arity in ARITIES.items():
            path = ROOT / "cases/generated/special" / (name + ".csv")
            with path.open(newline="") as stream:
                cases = list(csv.DictReader(stream))
            self.assertEqual(len(cases), 2001, name)
            parameters = np.asarray([
                [float(case[f"p{i}"].split(":", 1)[1]) for i in range(arity)]
                for case in cases], dtype=float)
            if name == "factorialk":
                expected = np.asarray([
                    special.factorialk(int(values[0]), int(values[1]))
                    for values in parameters], dtype=float)
            elif name == "multigammaln":
                expected = np.asarray([
                    special.multigammaln(values[0], int(values[1]))
                    for values in parameters], dtype=float)
            else:
                expected = np.asarray(getattr(special, name)(*parameters.T), dtype=float)
            self.assertEqual(expected.shape, (2001,))
            for args, target in zip(parameters, expected):
                with self.subTest(name=name, args=tuple(args)):
                    value = self.value(name, *args)
                    self.assertTrue(math.isfinite(float(target)))
                    absolute = abs(value - float(target))
                    scaled = absolute / max(abs(float(target)), 1e-14)
                    greatest_absolute = max(greatest_absolute, absolute)
                    greatest_scaled = max(greatest_scaled, scaled)
                    self.assertTrue(math.isclose(value, float(target), rel_tol=2e-11,
                                                 abs_tol=3e-13),
                                    (name, tuple(args), value, target))
            count += len(cases)
        elapsed = time.perf_counter() - started
        self.assertEqual(count, 22011)
        self.assertLess(elapsed, 60., "Boost special transform case evaluation stalled")
        self.assertLess(greatest_scaled, 2e-11)
        print(f"\nBoost-derived transforms: {count} saved SciPy cases / {elapsed:.3f}s; "
              f"max absolute {greatest_absolute:.3g}, max scaled {greatest_scaled:.3g}")

    def test_independent_transform_identities(self):
        for x in (1e-9, .2, .5, 2., 3., 10.):
            for lam in (-2., -.5, 0., .3, 1., 2.):
                target = math.log(x) if lam == 0. else math.expm1(lam*math.log(x))/lam
                self.assertAlmostEqual(self.value("boxcox", x, lam), target, delta=3e-13)
                if x < 3.:
                    target1p = math.log1p(x) if lam == 0. else math.expm1(lam*math.log1p(x))/lam
                    self.assertAlmostEqual(self.value("boxcox1p", x, lam), target1p,
                                           delta=3e-13)
        for x in (-100., -5., -.5, 0., .5, 5., 100.):
            soft = max(x, 0.) + math.log1p(math.exp(-abs(x)))
            self.assertAlmostEqual(self.value("softplus", x), soft, delta=2e-14)
            self.assertAlmostEqual(self.value("log_expit", x),
                                   -(max(-x,0.) + math.log1p(math.exp(-abs(x)))), delta=2e-14)

    def test_multivariate_gamma_and_multifactorial(self):
        for d in range(1, 8):
            for a in (3.5, 5.5, 10.5):
                if a <= (d-1)/2: continue
                target = d*(d-1)/4 * math.log(math.pi) + sum(
                    math.lgamma(a + (1-j)/2) for j in range(1,d+1))
                self.assertAlmostEqual(self.value("multigammaln", a, d), target,
                                       delta=3e-12)
        for k in (1,2,3,4,5,8):
            for n in range(0,41):
                exact = math.prod(range(n,0,-k))
                self.assertTrue(math.isclose(self.value("factorialk",n,k), exact,
                                             rel_tol=5e-13, abs_tol=3e-13))

    def test_normal_cdf_inverse_and_parameter_roundtrips(self):
        for p in (1e-7, .01, .1, .3, .8, .99, 1-1e-7):
            mean, std, x = .2, 2.3, 1.5
            shifted_x = self.value("nrdtrimn",p,std,x)
            self.assertAlmostEqual(.5*math.erfc(-(x-shifted_x)/(std*math.sqrt(2))),p,
                                   delta=2e-14)
            # Choose x below mean for p<0.5 and above for p>0.5.
            if p != .5:
                position = mean + (1 if p>.5 else -1)*1.7
                fitted_std = self.value("nrdtrisd",mean,p,position)
                self.assertGreater(fitted_std,0.)
                cdf=.5*math.erfc(-(position-mean)/(fitted_std*math.sqrt(2)))
                self.assertAlmostEqual(cdf,p,delta=2e-14)
        for y in (-699., -100., -20., -2., -.5, -.1, -1e-8):
            x=self.value("ndtri_exp",y)
            # Compare against SciPy's separate numerically stable log-CDF
            # and check that inversion is consistent even in the far tail.
            self.assertAlmostEqual(float(special.log_ndtr(x)),y,delta=1e-12)
        for x in (-35., -30., -10., -2., 0., 2., 10., 35.):
            self.assertAlmostEqual(self.value("log_ndtr",x),
                                   float(special.log_ndtr(x)),delta=3e-12)

    def test_agm_against_independent_mean_iteration(self):
        for a,b in ((1.,2.), (1.,1.), (1e-6,1e-4), (.125,128.), (3.,5.), (100.,101.)):
            arithmetic, geometric=a,b
            for _ in range(100):
                new_a=(arithmetic+geometric)/2
                new_g=math.sqrt(arithmetic*geometric)
                if abs(new_a-new_g)<=2e-15*max(new_a,new_g):
                    arithmetic, geometric=new_a,new_g
                    break
                arithmetic, geometric=new_a,new_g
            self.assertTrue(math.isclose(self.value("agm",a,b),arithmetic,
                                         rel_tol=2e-12,abs_tol=1e-14))

    def test_guarded_domains_and_untouched_abi_output_on_failure(self):
        invalid = [
            ("multigammaln", (1.,3.)), ("multigammaln", (2.,0.)),
            ("multigammaln", (2.,2.5)), ("boxcox",(0.,.3)),
            ("boxcox",(-1.,0.)), ("boxcox1p",(-2.,.5)),
            ("nrdtrimn",(.3,0.,1.)), ("nrdtrimn",(.0,2.,1.)),
            ("nrdtrisd",(1.,.5,2.)), ("nrdtrisd",(1.,.3,2.)),
            ("ndtri_exp",(.01,)), ("ndtri_exp",(-701.,)),
            ("log_ndtr",(-36.,)), ("factorialk",(4.5,2.)),
            ("factorialk",(5.,0.)), ("agm",(-1.,2.)),
            ("agm",(1.,0.)), ("agm",(1e-8,1e8)),
            ("softplus",(float('nan'),)),
            ("log_expit",(float('inf'),)),
        ]
        for name,args in invalid:
            with self.subTest(name=name,args=args):
                got=self.adapter.evaluate(self.case(name,*args))
                self.assertIn("error",got.get("Boost.Math",{}),got)
                parameters=(ct.c_double*len(args))(*args)
                out=(ct.c_double*1)(13579.5)
                status=self.adapter.library.boost_special_transforms_evaluate(
                    0,name.encode(),parameters,len(args),out,1)
                self.assertNotEqual(status,0)
                self.assertEqual(out[0],13579.5)
        params=(ct.c_double*1)(1.)
        out=(ct.c_double*1)(13579.5)
        self.assertNotEqual(self.adapter.library.boost_special_transforms_evaluate(
            1,b"softplus",params,1,out,1),0)
        self.assertEqual(out[0],13579.5)


if __name__ == "__main__":
    unittest.main()
