"""Independent Boost.Math and GSL scalar-special comparisons on saved cases."""

import csv
import ctypes as ct
import math
import os
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy import special

from sci_testing.extra_native_comparisons import ARITIES, ExtraNativeComparisons
from tests.pruned_fixtures import removed_apis

ROOT = Path(__file__).resolve().parents[1]


class ExtraNativeComparisonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        temporary = tempfile.TemporaryDirectory()
        cls.addClassCleanup(temporary.cleanup)
        output = Path(temporary.name) / "extra-special-scalar.so"
        env = dict(os.environ)
        env.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2", "-Wall",
            "-Wextra", "-Werror", "-shared", "-fPIC",
            str(ROOT / "native/extra_special_scalar.cpp"),
            "-lgsl", "-lgslcblas", "-o", str(output),
        ], env=env, check=True)
        cls.adapter = ExtraNativeComparisons(ct.CDLL(str(output)))

    @staticmethod
    def case(api, args):
        return {"api": "special." + api,
                "args": {f"p{i}": value for i, value in enumerate(args)}}

    def check_scalar(self, name, args, rtol=3e-12, atol=3e-14):
        actual = self.adapter.evaluate(self.case(name, args))
        providers = {"Boost.Math"} if name == "cbrt" else {"Boost.Math", "GSL"}
        self.assertEqual(set(actual), providers)
        expected = float(getattr(special, name)(*args))
        self.assertTrue(math.isfinite(expected), (name, args))
        for provider in providers:
            with self.subTest(api=name, provider=provider, args=args):
                self.assertEqual(set(actual[provider]), {"values"}, actual[provider])
                self.assertEqual(len(actual[provider]["values"]), 1)
                np.testing.assert_allclose(actual[provider]["values"][0], expected,
                                           rtol=rtol, atol=atol)

    def test_all_retained_saved_fixtures(self):
        removed = removed_apis()
        measured: set[str] = set()
        for api, count in ARITIES.items():
            name = api.rsplit(".", 1)[-1]
            path = ROOT / "cases/generated/special" / (name + ".csv")
            if api in removed:
                self.assertFalse(path.exists(), f"removed fixture restored: {api}")
                continue
            with path.open(newline="") as stream:
                cases = list(csv.DictReader(stream))
            self.assertEqual(len(cases), 2001)
            measured.add(api)
            for row in cases:
                args: list[int | float] = []
                for i in range(count):
                    kind, number = row[f"p{i}"].split(":", 1)
                    args.append(int(number) if kind == "i" else float(number))
                self.check_scalar(name, args)
        self.assertEqual(measured, set(ARITIES) - removed)

    def test_distinct_domains_and_function_identities(self):
        for x in (-13., -1.25, -.001, 0., 1e-12, .75, 12.):
            self.check_scalar("cbrt", (x,))
            self.check_scalar("expm1", (x,))
            self.check_scalar("psi", (x + 13.5,))
            if x > -1:
                self.check_scalar("log1p", (x,))
        for x, degree, alpha in ((-.5, 0, .5), (2., 3, 0.),
                                 (3.75, 5, .75), (0., 12, 2.5),
                                 (0.5, 2, 17.)):
            self.check_scalar("assoc_laguerre", (x, degree, alpha))

    def test_unavailable_inputs_domains_and_output_capacity(self):
        self.assertEqual(self.adapter.evaluate(self.case("unsupported", (1.,))), {})
        self.assertEqual(self.adapter.evaluate({"api": "special.psi", "args": {}}), {})
        self.assertEqual(self.adapter.evaluate(self.case("psi", ([1., 2.],))), {})
        self.assertEqual(self.adapter.evaluate(self.case("cbrt", (float("nan"),))), {})
        for name, args in (("log1p", (-1.,)), ("log1p", (-2.,)),
                           ("psi", (0.,)), ("assoc_laguerre", (2., 2.5, .5)),
                           ("assoc_laguerre", (2., 3., -1.)),
                           ("assoc_laguerre", (2., 129., .5))):
            comparison = self.adapter.evaluate(self.case(name, args))
            for provider, result in comparison.items():
                with self.subTest(api=name, provider=provider):
                    self.assertIn("error", result)
                    self.assertTrue(result["error"])

        library = self.adapter.library
        self.assertEqual(library.extra_special_scalar_supported(0, b"cbrt"), 1)
        self.assertEqual(library.extra_special_scalar_supported(1, b"cbrt"), 0)
        arguments = (ct.c_double * 1)(2.)
        output = (ct.c_double * 1)()
        status = library.extra_special_scalar_evaluate(
            0, b"psi", arguments, 1, output, 0)
        self.assertNotEqual(status, 0)
        self.assertIn(b"buffer", library.extra_special_scalar_error())


if __name__ == "__main__":
    unittest.main()
