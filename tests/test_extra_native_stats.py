"""Saved-case fidelity and C ABI tests for independent GSL-derived oracles."""

import csv
import ctypes as ct
import gzip
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
from mpmath import mp
from scipy import special, stats

from sci_testing.case_csv import decode
from sci_testing.extra_native_stats import (
    ExtraNativeStatsComparisons, SPECIAL_ARITIES, STATS_NAMES, SUPPORTED,
)


ROOT = Path(__file__).resolve().parents[1]


def saved(api):
    folder, name = api.rsplit(".", 1)
    filename = ROOT / "cases/generated" / folder / (name + ".csv")
    if not filename.exists():
        filename = Path(str(filename) + ".gz")
    with (gzip.open(filename, "rt", newline="") if filename.suffix == ".gz"
          else filename.open(newline="")) as stream:
        for row in csv.DictReader(stream):
            yield row["id"], {"api": api, "args": {
                key: decode(value) for key, value in row.items()
                if key not in ("id", "group")
            }}


def expected_stats(name, args):
    x, y, p = args["x"], args["y"], args["p"]
    if name == "entropy":
        return [stats.entropy(x, y if y else None, base=p[0] or None)]
    if name == "rankdata":
        return stats.rankdata(x, method=("average", "min", "max", "dense", "ordinal")[
            int(p[0])]).tolist()
    if name == "tiecorrect":
        return [stats.tiecorrect(x)]
    if name == "mode":
        result = stats.mode(x)
        return [result.mode, result.count]
    if name == "wasserstein_distance":
        return [stats.wasserstein_distance(x, y)]
    if name == "energy_distance":
        return [stats.energy_distance(x, y)]
    if name == "boxcox_llf":
        return [stats.boxcox_llf(p[0], x)]
    if name == "zmap":
        return np.ravel(stats.zmap(x, y, ddof=int(p[0]))).tolist()
    if name == "expectile":
        return [stats.expectile(x, p[0])]
    raise AssertionError(name)


def expected(case):
    name = case["api"].rsplit(".", 1)[-1]
    args = case["args"]
    if case["api"].startswith("stats."):
        return expected_stats(name, args)
    return [getattr(special, name)(*[args[f"p{i}"] for i in range(
        SPECIAL_ARITIES[name])])]


class ExtraNativeStatsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(temp.cleanup)
        output = Path(temp.name) / "extra-native-stats.so"
        env = dict(os.environ)
        env.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2",
            "-Wall", "-Wextra", "-Werror", "-shared", "-fPIC",
            str(ROOT / "native/extra_stats_scalar.cpp"),
            "-lgsl", "-lgslcblas", "-o", str(output),
        ], check=True, env=env)
        cls.adapter = ExtraNativeStatsComparisons(ct.CDLL(str(output)))

    def check(self, case, *, atol=3e-11, rtol=3e-10):
        values = self.adapter.evaluate(case)
        self.assertEqual(set(values), {"GSL"}, (case, values))
        self.assertEqual(set(values["GSL"]), {"values"}, (case, values))
        reference = np.ravel(expected(case)).astype(float)
        native = np.asarray(values["GSL"]["values"])
        self.assertEqual(native.shape, reference.shape, case)
        np.testing.assert_allclose(native, reference, atol=atol, rtol=rtol,
                                   err_msg=f"{case['api']} {case['args']}")

    def test_all_saved_cases(self):
        counts = {}
        for api in sorted(SUPPORTED):
            count = 0
            for case_id, case in saved(api):
                with self.subTest(id=case_id):
                    # SciPy's Box-Cox likelihood can lose ~1e-8 at lambda=1e-8;
                    # high-precision direct verification is tested below.
                    atol = 3e-8 if api == "stats.boxcox_llf" and abs(
                        case["args"]["p"][0]) < 1e-7 else 3e-11
                    self.check(case, atol=atol)
                count += 1
            counts[api] = count
            self.assertGreaterEqual(count, 2000, api)
        self.assertEqual(set(counts), SUPPORTED)

    def test_extended_statistical_domains(self):
        examples = [
            ("stats.rankdata", {"x": [1, 1, 1, 1], "y": [], "p": [method]})
            for method in range(5)
        ] + [
            ("stats.mode", {"x": [2, 1, 2, 1, 4], "y": [], "p": []}),
            ("stats.tiecorrect", {"x": [1, 1, 2, 2, 2], "y": [], "p": []}),
            ("stats.entropy", {"x": [0, 1, 3], "y": [2, 1, 4], "p": [2]}),
            ("stats.entropy", {"x": [0, 1, 3], "y": [], "p": [0]}),
            ("stats.wasserstein_distance", {"x": [-3, 0, 3], "y": [1, 2], "p": []}),
            ("stats.energy_distance", {"x": [-3, 0, 3], "y": [1, 2], "p": []}),
            ("stats.boxcox_llf", {"x": [.1, .5, 2., 7.], "y": [], "p": [0]}),
            ("stats.zmap", {"x": [2, 3], "y": [1, 2, 4, 5], "p": [0]}),
            ("stats.expectile", {"x": [-3, 0, 2, 6], "y": [], "p": [.1]}),
            ("stats.expectile", {"x": [-3, 0, 2, 6], "y": [], "p": [.9]}),
        ]
        for api, args in examples:
            with self.subTest(api=api, args=args):
                self.check({"api": api, "args": {**args, "d": "", "g": []}})
        for api, args in (
            ("special.entr", (0.,)),
            ("special.xlogy", (0., 0.)),
            ("special.xlog1py", (0., -1.)),
            ("special.rel_entr", (0., 0.)),
            ("special.kl_div", (0., 3.)),
            ("special.cosm1", (1e-8,)),
            ("special.tandg", (-25.,)),
        ):
            with self.subTest(api=api, args=args):
                self.check({"api": api, "args": {f"p{i}": v for i, v in enumerate(args)}})

    def test_boxcox_near_zero_against_high_precision_variance(self):
        with mp.workdps(80):
            for data in ([1., 2., 4., 8.], [.1, .5, 2., 7.]):
                for transform in (mp.mpf("1e-8"), mp.mpf("1e-6"), mp.mpf("0.5")):
                    xs = [mp.mpf(str(value)) for value in data]
                    values = [mp.expm1(transform * mp.log(value))/transform
                              for value in xs]
                    mean = mp.fsum(values)/len(values)
                    variance = mp.fsum((value-mean)**2 for value in values)/len(values)
                    truth = ((transform-1)*mp.fsum(mp.log(value) for value in xs)
                             - len(values)/2 * mp.log(variance))
                    result = self.adapter.evaluate({"api": "stats.boxcox_llf", "args": {
                        "d": "", "g": [], "x": list(data), "y": [], "p": [float(transform)]}})
                    with self.subTest(data=data, transform=str(transform)):
                        self.assertIn("values", result["GSL"])
                        self.assertAlmostEqual(result["GSL"]["values"][0],
                                               float(truth), delta=8e-13)

    def test_provider_ownership_invalid_domain_and_capacity(self):
        library = self.adapter.library
        for api in sorted(SUPPORTED):
            self.assertEqual(library.extra_native_stats_supported(0, api.encode()), 0)
            self.assertEqual(library.extra_native_stats_supported(1, api.encode()), 1)
        self.assertEqual(self.adapter.evaluate({"api": "stats.unknown", "args": {}}), {})
        self.assertEqual(self.adapter.evaluate({"api": "stats.rankdata", "args": {
            "x": [1, 2], "y": [], "p": []}}), {})
        for api, args in (
            ("stats.entropy", {"x": [1, -2], "y": [], "p": [0]}),
            ("stats.entropy", {"x": [1, 2], "y": [], "p": [1]}),
            ("stats.boxcox_llf", {"x": [1, -2], "y": [], "p": [0]}),
            ("stats.zmap", {"x": [1], "y": [2, 2], "p": [1]}),
            ("stats.rankdata", {"x": [1, 2], "y": [], "p": [5]}),
            ("special.entr", {"p0": -1.}),
            ("special.rel_entr", {"p0": 1., "p1": 0.}),
            ("special.xlogy", {"p0": 1., "p1": -1.}),
        ):
            result = self.adapter.evaluate({"api": api, "args": args})
            with self.subTest(api=api, args=args):
                self.assertIn("error", result["GSL"])
                self.assertTrue(result["GSL"]["error"])
        data = (ct.c_double * 2)(1., 2.)
        output = (ct.c_double * 1)()
        status = library.extra_native_stats_evaluate(
            1, b"stats.rankdata", data, 2, None, 0, (ct.c_double * 1)(0.),
            1, output, 1)
        self.assertNotEqual(status, 0)
        self.assertIn(b"buffer", library.extra_native_stats_error())


if __name__ == "__main__":
    unittest.main()
