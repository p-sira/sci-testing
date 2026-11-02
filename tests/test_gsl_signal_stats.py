"""Independent saved-case SciPy tests of GSL-derived signal statistics."""

import ctypes as ct
import math
import os
from pathlib import Path
import subprocess
import tempfile
import time
import unittest
import warnings

import numpy as np
from scipy import signal

from sci_testing.case_csv import read_case_file
from sci_testing.gsl_signal_stats import GSLSignalStatsComparisons, NAMES


ROOT = Path(__file__).resolve().parents[1]
COUNTS = {"detrend": 2006, "medfilt": 2003, "medfilt2d": 2003, "wiener": 2005}


def independent(case):
    name = case["api"].rsplit(".", 1)[-1]
    a = case["args"]["arguments"]
    x = np.asarray(a[0], dtype=np.float64)
    if name == "detrend":
        result = signal.detrend(x, type="linear" if a[1] else "constant")
    elif name == "medfilt":
        result = signal.medfilt(x, kernel_size=a[1])
    elif name == "medfilt2d":
        result = signal.medfilt2d(x, kernel_size=a[1])
    else:
        result = signal.wiener(x, mysize=a[1], noise=a[2])
    return np.asarray(result, dtype=float).ravel()


class GSLSignalStatsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        temporary = tempfile.TemporaryDirectory()
        cls.addClassCleanup(temporary.cleanup)
        target = Path(temporary.name) / "gsl-signal-stats.so"
        environment = os.environ.copy()
        environment.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2",
            "-Wall", "-Wextra", "-Werror", "-Wpedantic", "-shared", "-fPIC",
            str(ROOT / "native/gsl_signal_stats.cpp"),
            "-lgsl", "-lgslcblas", "-o", str(target),
        ], env=environment, check=True)
        cls.native = GSLSignalStatsComparisons(ct.CDLL(str(target)))

    @staticmethod
    def case(name, *args):
        return {"api": "signal." + name, "args": {"arguments": list(args)}}

    def checked(self, case):
        response = self.native.evaluate(case)
        self.assertEqual(set(response), {"GSL"}, (case, response))
        entry = response["GSL"]
        self.assertNotIn("error", entry, (case, response))
        values = np.asarray(entry["values"], dtype=float)
        self.assertTrue(np.all(np.isfinite(values)))
        return values

    def test_every_8017_saved_scipy_case(self):
        checked = 0
        maxima = {}
        started = time.perf_counter()
        for name, count in COUNTS.items():
            source = ROOT / "cases/generated/signal" / (name + ".csv.gz")
            if not source.exists():
                source = source.with_suffix("")
            cases = read_case_file(source, "signal/" + name + ".csv")
            self.assertEqual(len(cases), count, name)
            maxima[name] = 0.
            for case in cases:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    expected = independent(case)
                with self.subTest(saved_id=case["id"]):
                    actual = self.checked(case)
                    self.assertEqual(actual.shape, expected.shape)
                    np.testing.assert_allclose(actual, expected,
                                               rtol=2e-10, atol=3e-12)
                maxima[name] = max(
                    maxima[name], float(np.max(np.abs(actual-expected))))
                checked += 1
        self.assertEqual(checked, 8017)
        self.assertLess(time.perf_counter()-started, 45.)
        print(f"GSL signal saved cases: {checked}, max absolute discrepancy {maxima}")

    def test_detrend_linear_and_constant_nontrivial_profiles(self):
        for values in (
            [1.], [2.,3.], [1.,3.,5.,7.,9.],
            [3.,-1.,4.,1.5,-2.,8.,0.],
            [1e9,1e9+1,1e9+3,1e9+6],
        ):
            for linear in (False, True):
                with self.subTest(values=values, linear=linear):
                    case = self.case("detrend", values, linear)
                    expected = independent(case)
                    actual = self.checked(case)
                    np.testing.assert_allclose(actual, expected,
                                               rtol=2e-10, atol=2e-6 if max(values)>1e8 else 1e-12)
                    if linear and len(values)>1:
                        self.assertLess(abs(np.sum(actual)), 1e-5)

    def test_zero_padded_median_1d_and_2d(self):
        x = [-9.,3.,0.,8.,-3.,2.,5.,-1.]
        for kernel in (1,3,5,7,9):
            actual = self.checked(self.case("medfilt",x,kernel))
            expected = signal.medfilt(x,kernel_size=kernel)
            np.testing.assert_array_equal(actual,expected)
        arr = [[1.,-2.,0.,3.],[-1.,2.,4.,5.],[3.,-1.,9.,0.]]
        for kernel in ([1,1],[3,3],[3,5],[5,7]):
            actual = self.checked(self.case("medfilt2d",arr,kernel))
            expected = signal.medfilt2d(np.asarray(arr),kernel_size=kernel).ravel()
            np.testing.assert_array_equal(actual,expected)

    def test_wiener_local_mean_variance_noise_estimation(self):
        x = [1.,2.,-1.,0.,3.,2.,4.,-2.,1.,0.,2.,3.]
        for size in (1,3,5,7):
            for noise in (0.,.2,1.):
                with self.subTest(size=size,noise=noise):
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore",RuntimeWarning)
                        expected = signal.wiener(np.asarray(x),mysize=size,noise=noise)
                    if not np.all(np.isfinite(expected)):
                        result=self.native.evaluate(self.case("wiener",x,size,noise))
                        self.assertIn("error",result["GSL"],result)
                        continue  # SciPy's 0/0 is not a finite numerical oracle.
                    actual = self.checked(self.case("wiener",x,size,noise))
                    np.testing.assert_allclose(actual, expected,
                                               rtol=1e-11, atol=1e-12)
            if size > 1:
                actual = self.checked(self.case("wiener",x,size,None))
                expected = signal.wiener(np.asarray(x),mysize=size,noise=None)
                np.testing.assert_allclose(actual, expected,
                                           rtol=1e-11, atol=1e-12)

    def test_c_abi_output_capacity_guards_and_provider(self):
        library = self.native.library
        for name in NAMES:
            encoded = name.encode("ascii")
            self.assertEqual(library.gsl_signal_stats_supported(1,encoded),1)
            self.assertEqual(library.gsl_signal_stats_supported(0,encoded),0)
            self.assertEqual(library.gsl_signal_stats_output_count(
                encoded,1,12),12)
        self.assertEqual(library.gsl_signal_stats_output_count(b"detrend",2,12),0)
        self.assertEqual(library.gsl_signal_stats_output_count(b"medfilt2d",3,3),9)
        self.assertEqual(library.gsl_signal_stats_output_count(b"medfilt",1,0),0)
        self.assertEqual(library.gsl_signal_stats_supported(1,b"lfilter"),0)
        data=(ct.c_double*3)(1.,2.,-1.)
        output=(ct.c_double*4)(11.,22.,33.,44.)
        for provider, kernel, capacity in ((0,3,4),(1,2,4),(1,3,2)):
            status = library.gsl_signal_stats_evaluate(
                provider,b"medfilt",data,1,3,1,kernel,0,0.,output,capacity)
            self.assertNotEqual(status,0)
            self.assertEqual(list(output),[11.,22.,33.,44.])
            self.assertTrue(library.gsl_signal_stats_error())
        self.assertEqual(library.gsl_signal_stats_evaluate(
            1,b"medfilt",data,1,3,1,3,0,0.,output,3),0)
        np.testing.assert_array_equal(list(output[:3]), [1.,1.,0.])
        self.assertEqual(output[3],44.)
        self.assertEqual(library.gsl_signal_stats_error(),b"")

    def test_invalid_domains_and_unsupported_shapes(self):
        for name,args in (
            ("detrend",([1.,2.], 1)),
            ("detrend",([1.,math.nan],True)),
            ("medfilt",([1.,2.],2)),
            ("medfilt",([1.,2.],-1)),
            ("medfilt",([1.,2.],128)),
            ("medfilt2d",([[1.,2.],[3.,4.]],[2,3])),
            ("wiener",([1.,2.,3.],2,.1)),
            ("wiener",([1.,2.,3.],3,-.2)),
            ("wiener",([1.,2.,3.],3,math.inf)),
        ):
            actual=self.native.evaluate(self.case(name,*args))
            with self.subTest(name=name,args=args):
                if not actual:
                    continue  # malformed Python signature is rejected early
                self.assertIn("error",actual["GSL"],actual)
        self.assertEqual(self.native.evaluate(self.case("lfilter",[1,2],3)),{})
        self.assertEqual(self.native.evaluate(
            {"api":"spatial.medfilt","args":{"arguments":[[1,2],3]}}),{})
        self.assertEqual(self.native.evaluate(
            {"api":"signal.medfilt","args":{"inputs":[[1,2],3]}}),{})
        self.assertEqual(self.native.evaluate(
            {"api":"signal.medfilt2d","args":{"arguments":[[[1,2],[3]], [3,3]]}}),{})


if __name__=="__main__":
    unittest.main()
