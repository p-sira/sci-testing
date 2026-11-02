"""Boost.Math vector logsumexp and stable softmax against public SciPy.

All checked-in generated input fixtures are run through the standalone C
adapter. No Sci internals or report/dispatch files are modified.
"""

import ctypes as ct
import csv
import gzip
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import numpy as np
from scipy import special

from sci_testing.vector_special import VectorSpecialComparisons


ROOT = Path(__file__).resolve().parent.parent
NAMES = ("logsumexp", "softmax", "log_softmax")


def _saved_vectors(name):
    directory = ROOT / "cases/generated/special"
    plain = directory / f"{name}.csv"
    compressed = directory / f"{name}.csv.gz"
    path = plain if plain.exists() else compressed
    if not path.exists():
        raise FileNotFoundError(f"No generated input fixture for {name}")
    open_file = gzip.open if path.suffix == ".gz" else open
    with open_file(path, mode="rt", newline="") as stream:
        for row in csv.DictReader(stream):
            assert row["id"].startswith(f"special.{name}/")
            assert row["p0"].startswith("v:")
            yield row["id"], [float(value) for value in row["p0"][2:].split()]


class VectorSpecialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        directory = tempfile.TemporaryDirectory()
        cls.addClassCleanup(directory.cleanup)
        shared_object = Path(directory.name) / "vector_special.so"
        environment = dict(os.environ)
        environment.pop("LD_LIBRARY_PATH", None)  # avoid the AppImage's GCC shadowing
        compiler = os.environ.get("CXX", shutil.which("clang++") or "c++")
        subprocess.run([
            compiler, "-std=c++17", "-Wall", "-Wextra", "-Werror", "-O2",
            "-shared", "-fPIC", str(ROOT / "native/vector_special.cpp"),
            "-o", str(shared_object),
        ], env=environment, check=True)
        cls.lib = ct.CDLL(str(shared_object))
        cls.native = VectorSpecialComparisons(cls.lib)

    @staticmethod
    def case(name, values, *, nested=False):
        path = "special.elementary" if nested else "special"
        return {"api": f"{path}.{name}", "args": {"p0": list(values)}}

    def assert_matches_scipy(self, name, values, *, atol=2e-12, rtol=2e-12):
        result = self.native.evaluate(self.case(name, values))
        self.assertEqual(set(result), {"Boost.Math"})
        self.assertNotIn("error", result["Boost.Math"], result)
        observed = result["Boost.Math"]["values"]
        with np.errstate(all="ignore"):
            expected = np.asarray(getattr(special, name)(
                np.array(values, dtype=np.float64))).ravel()
        self.assertEqual(len(observed), len(expected), (name, values))
        np.testing.assert_allclose(
            observed, expected, rtol=rtol, atol=atol, equal_nan=True,
            err_msg=f"{name}({values})")
        return observed

    def test_every_saved_fixture_against_scipy(self):
        for name in NAMES:
            count = 0
            for identifier, vector in _saved_vectors(name):
                with self.subTest(fixture=identifier):
                    self.assert_matches_scipy(name, vector)
                count += 1
            self.assertGreaterEqual(count, 500, name)

    def test_strongly_shifted_extreme_vectors(self):
        vectors = (
            [1.0],
            [0.0, -0.0],
            [1000.0, 1000.0, 999.0],
            [-1000.0, -1001.0, -1002.0],
            [1e308, 1e308, -1e308],
            [-1e308, -1e308, -1e308 + 1e292],
            [1e-300, -1e-300, 0.0],
            [0.0, -745.0, -1000.0, -1e100],
            [0.0, -1000.0] * 100,
            [-1200.0, -1200.0, -1200.0, -1199.0],
        )
        for values in vectors:
            for name in NAMES:
                with self.subTest(api=name, input=values[:6]):
                    observed = self.assert_matches_scipy(name, values)
                    self.assertTrue(all(math.isfinite(x) for x in observed)
                                    or name == "log_softmax" and
                                    all(math.isfinite(x) or x == -math.inf
                                        for x in observed))
                    if name == "softmax":
                        self.assertAlmostEqual(sum(observed), 1., places=13)
                    if name == "log_softmax":
                        self.assertAlmostEqual(sum(math.exp(x) for x in observed),
                                               1., places=13)

    def test_signed_infinity_and_nan_semantics(self):
        vectors = (
            [-math.inf, 0.0], [math.inf, 0.0],
            [math.inf, -math.inf], [-math.inf, -math.inf],
            [math.inf, math.inf], [math.nan, 1.0],
            [math.nan, math.inf], [0.0, -math.inf, math.inf],
            [math.inf], [-math.inf], [math.nan],
        )
        for values in vectors:
            for name in NAMES:
                with self.subTest(api=name, values=values):
                    self.assert_matches_scipy(name, values)
        self.assertEqual(self.assert_matches_scipy("logsumexp", []), [-math.inf])
        for name in ("softmax", "log_softmax"):
            record = self.native.evaluate(self.case(name, []))
            self.assertEqual(set(record), {"Boost.Math"})
            self.assertIn("error", record["Boost.Math"])

    def test_validation_and_dispatch(self):
        for name in NAMES:
            result = self.native.evaluate(self.case(name, [1., 2., 3.], nested=True))
            self.assertEqual(set(result), {"Boost.Math"})
            self.assertNotIn("error", result["Boost.Math"])
        for bad in (
            {"api": "special.foo", "args": {"p0": [1., 2.]}},
            {"api": "special.complex.softmax", "args": {"p0": [1., 2.]}},
            {"api": "special.softmax", "args": {"p0": 1.}},
            {"api": "special.softmax", "args": {"p0": [1., 2.], "p1": 3.}},
            {"api": "special.softmax", "args": {"p0": [1.+1.j, 2.]}},
            {"api": "special.logsumexp", "args": {}},
        ):
            self.assertEqual(self.native.evaluate(bad), {}, bad)

    def test_native_abi_output_sizes_and_failed_writes(self):
        for name in NAMES:
            self.assertEqual(self.lib.vector_special_supported(0, name.encode()), 1)
            self.assertEqual(self.lib.vector_special_supported(1, name.encode()), 0)
            self.assertEqual(self.lib.vector_special_output_count(
                name.encode(), 5), 1 if name == "logsumexp" else 5)
        self.assertEqual(self.lib.vector_special_output_count(b"logsumexp", 0), 1)
        self.assertEqual(self.lib.vector_special_output_count(b"softmax", 0), 0)
        self.assertEqual(self.lib.vector_special_output_count(b"other", 2), 0)

        input_values = (ct.c_double * 3)(1., 2., 3.)
        output = (ct.c_double * 3)(42., 43., 44.)
        for library, name, capacity, count, inputs in (
            (0, b"softmax", 2, 3, input_values),
            (1, b"softmax", 3, 3, input_values),
            (0, b"unknown", 3, 3, input_values),
            (0, b"softmax", 3, 3, None),
            (0, b"softmax", 3, -1, input_values),
            (0, b"logsumexp", 0, 3, input_values),
        ):
            status = self.lib.vector_special_evaluate(
                library, name, inputs, count, output, capacity)
            self.assertNotEqual(status, 0)
            self.assertEqual(list(output), [42., 43., 44.])
            self.assertTrue(self.lib.vector_special_error())


if __name__ == "__main__":
    unittest.main()
