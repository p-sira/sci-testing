"""Independent SciPy verification of GSL complex FFT on every saved case."""

import ctypes as ct
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
from scipy import fft

from sci_testing.case_csv import read_case_file
from sci_testing.gsl_fft_simple import GSLFFTSimpleComparisons


ROOT = Path(__file__).resolve().parents[1]


def saved_cases(name):
    return read_case_file(ROOT / f"cases/generated/fft/{name}.csv.gz",
                          f"fft/{name}.csv")


def independent_scipy(case):
    z = np.asarray(case["args"]["z"], dtype=float)
    values = z[:, 0] + 1j * z[:, 1]
    output = getattr(fft, case["api"].split(".")[-1])(values)
    return np.column_stack((output.real, output.imag)).ravel()


class GSLFFTSimpleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        folder = tempfile.TemporaryDirectory()
        cls.addClassCleanup(folder.cleanup)
        path = Path(folder.name) / "libgsl_fft_simple.so"
        env = os.environ.copy()
        env.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2",
            "-Wall", "-Wextra", "-Werror", "-Wpedantic",
            "-fPIC", "-shared", str(ROOT / "native/gsl_fft_simple.cpp"),
            "-lgsl", "-lgslcblas", "-o", str(path),
        ], env=env, check=True)
        cls.adapter = GSLFFTSimpleComparisons(ct.CDLL(str(path)))

    @staticmethod
    def case(name, pairs):
        return {"api": f"fft.{name}", "args": {"z": pairs}}

    def asserted(self, case):
        answer = self.adapter.evaluate(case)
        self.assertEqual(set(answer), {"GSL"}, (case, answer))
        self.assertNotIn("error", answer["GSL"], (case, answer))
        result = np.asarray(answer["GSL"]["values"], dtype=float)
        self.assertTrue(np.all(np.isfinite(result)))
        return result

    def test_every_4080_saved_fft_and_ifft_case(self):
        verified = 0
        maxima = {}
        for name in ("fft", "ifft"):
            cases = saved_cases(name)
            self.assertEqual(len(cases), 2040, name)
            maximum = 0.
            for case in cases:
                expected = independent_scipy(case)
                with self.subTest(saved_id=case["id"]):
                    actual = self.asserted(case)
                    self.assertEqual(actual.shape, expected.shape)
                    np.testing.assert_allclose(actual, expected, rtol=3e-12,
                                               atol=5e-12, err_msg=case["id"])
                maximum = max(maximum, float(np.max(np.abs(actual-expected))))
                verified += 1
            maxima[name] = maximum
        self.assertEqual(verified, 4080)
        print(f"GSL FFT saved cases: {verified}; max abs difference: {maxima}")

    def test_arbitrary_lengths_normalization_and_roundtrip(self):
        for n in (1, 2, 3, 4, 5, 7, 8, 11, 16, 17, 31, 48, 64, 97, 127):
            rng = np.random.default_rng(n)
            x = rng.normal(size=n) + 1j * rng.normal(size=n)
            encoded = np.column_stack((x.real, x.imag)).tolist()
            with self.subTest(n=n):
                fwd = self.asserted(self.case("fft", encoded))
                backward = self.asserted(self.case("ifft", encoded))
                np.testing.assert_allclose(fwd, independent_scipy(
                    self.case("fft", encoded)), rtol=4e-13, atol=4e-13)
                np.testing.assert_allclose(backward, independent_scipy(
                    self.case("ifft", encoded)), rtol=4e-13, atol=4e-13)
                transformed = np.asarray(fwd).reshape(-1, 2).tolist()
                returned = self.asserted(self.case("ifft", transformed))
                np.testing.assert_allclose(returned.reshape(-1, 2)[:, 0]
                                           + 1j * returned.reshape(-1, 2)[:, 1],
                                           x, rtol=4e-13, atol=4e-13)

    def test_native_abi_provider_capacity_errors_and_no_partial_output(self):
        lib = self.adapter.library
        for name in (b"fft.fft", b"fft.ifft"):
            self.assertEqual(lib.gsl_fft_simple_supported(1, name), 1)
            self.assertEqual(lib.gsl_fft_simple_supported(0, name), 0)
            self.assertEqual(lib.gsl_fft_simple_output_count(name, 3), 6)
            self.assertEqual(lib.gsl_fft_simple_output_count(name, 0), 0)
            self.assertEqual(lib.gsl_fft_simple_output_count(name, 65537), 0)
        self.assertEqual(lib.gsl_fft_simple_supported(1, b"fft.fftn"), 0)
        self.assertEqual(lib.gsl_fft_simple_output_count(b"fft.fftn", 3), 0)
        data = (ct.c_double * 6)(1., 2., -3., 4., 5., -6.)
        out = (ct.c_double * 7)(*[111.] * 7)
        for provider, name, capacity, source in (
            (0, b"fft.fft", 6, data),
            (1, b"fft.fftn", 6, data),
            (1, b"fft.fft", 5, data),
            (1, b"fft.fft", 6, None),
            (1, b"fft.fft", 6, (ct.c_double*6)(math.nan, 0., 1., 2., 3., 4.)),
            (1, b"fft.ifft", 6, (ct.c_double*6)(1., 2., math.inf, 0., 1., 1.)),
        ):
            with self.subTest(provider=provider, name=name, capacity=capacity):
                self.assertEqual(lib.gsl_fft_simple_evaluate(
                    provider, name, source, 3, out, capacity), -1)
                self.assertTrue(lib.gsl_fft_simple_error())
                self.assertEqual(list(out), [111.] * 7)
        self.assertEqual(lib.gsl_fft_simple_evaluate(
            1, b"fft.fft", data, 3, out, 6), 0)
        self.assertEqual(out[6], 111.)
        self.assertEqual(lib.gsl_fft_simple_error(), b"")
        self.assertEqual(self.adapter.evaluate(self.case("fft", []))["GSL"].keys(),
                         {"error"})
        self.assertIn("error", self.adapter.evaluate(self.case(
            "ifft", [[1., math.nan]]))["GSL"])
        self.assertIn("error", self.adapter.evaluate(self.case(
            "fft", [[0., 0., 1.]]))["GSL"])
        self.assertEqual(self.adapter.evaluate(self.case(
            "fftn", [[1., 2.]])), {})
        self.assertEqual(self.adapter.evaluate({
            "api": "fft.fft", "args": {"values": [[1., 2.]]}}), {})


if __name__ == "__main__":
    unittest.main()
