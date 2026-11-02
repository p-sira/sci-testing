"""Saved-case fidelity and provider-ownership tests for Boost statistics."""
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
from scipy import stats

from sci_testing.boost_stats_more import BoostStatsMoreComparisons, NAMES

ROOT = Path(__file__).resolve().parent.parent


def saved(name):
    path = ROOT / f"cases/generated/stats/{name}.csv"
    if not path.is_file():
        path = Path(str(path) + ".gz")
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, mode="rt", newline="") as f:
        for row in csv.DictReader(f):
            def numbers(field):
                tag, _, value = row[field].partition(":")
                assert tag == "v", row
                return [float(v) for v in value.split()]
            yield row["id"], {
                "api": f"stats.{name}",
                "args": {field: numbers(field) for field in ("x", "y", "p")},
            }


def independent(name, a):
    x, y, p = a["x"], a["y"], a["p"]
    if name == "gmean":
        return [stats.gmean(x)]
    if name == "hmean":
        return [stats.hmean(x)]
    if name == "gstd":
        return [stats.gstd(x, ddof=int(p[0]))]
    if name == "jarque_bera":
        r = stats.jarque_bera(x)
        return [r.statistic, r.pvalue]
    if name == "mode":
        r = stats.mode(x)
        return [r.mode, r.count]
    if name == "binomtest":
        r = stats.binomtest(int(p[0]), int(p[1]), p[2],
                            alternative=("two-sided", "less", "greater")[int(p[3])])
        return [r.statistic, r.pvalue]
    if name == "fisher_exact":
        r = stats.fisher_exact(np.asarray(x).reshape(2, 2),
                               alternative=("two-sided", "less", "greater")[int(p[0])])
        return [r.statistic, r.pvalue]
    if name in ("chisquare", "power_divergence"):
        if name == "power_divergence" and p[1] == -0.5 and any(v == 0 for v in x):
            statistic = 4 * sum((math.sqrt(v)-math.sqrt(z))**2 for v, z in zip(x,y))
            return [statistic, stats.chi2.sf(statistic,len(x)-1-int(p[0]))]
        kwargs = {"ddof": int(p[0])}
        if name == "power_divergence":
            kwargs["lambda_"] = p[1]
        r = getattr(stats, name)(x, y, **kwargs)
        return [r.statistic, r.pvalue]
    raise ValueError(name)


class BoostStatsMoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(tmp.cleanup)
        target = Path(tmp.name) / "boost-stats-more.so"
        env = os.environ.copy()
        env.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2", "-shared",
            "-fPIC", str(ROOT / "native/boost_stats_more.cpp"), "-o", str(target)
        ], check=True, env=env)
        cls.adapter = BoostStatsMoreComparisons(ct.CDLL(str(target)))

    def test_all_saved_cases_against_independent_scipy(self):
        failures = []
        for name in sorted(NAMES):
            count = 0
            for case_id, case in saved(name):
                count += 1
                result = self.adapter.evaluate(case)
                if set(result) != {"Boost.Math"} or "error" in result["Boost.Math"]:
                    failures.append((case_id, result))
                    continue
                expected = np.asarray(independent(name, case["args"]), dtype=float)
                actual = np.asarray(result["Boost.Math"]["values"], dtype=float)
                if (actual.shape != expected.shape or
                    not np.allclose(actual, expected, rtol=3e-11, atol=2e-12)):
                    failures.append((case_id, actual.tolist(), expected.tolist()))
            self.assertGreaterEqual(count, 2000, name)
        self.assertFalse(failures, f"{len(failures)} numerical/ABI failures; examples {failures[:12]}")

    def test_unsupported_and_invalid(self):
        self.assertEqual(self.adapter.evaluate({"api": "stats.unknown",
                                                "args": {"x": [1.]}}), {})
        self.assertEqual(self.adapter.evaluate({"api": "special.gmean",
                                                "args": {"x": [1.]}}), {})
        for name, x, y, p in [
            ("gmean", [-1., 2.], [], []),
            ("gstd", [1., 2.], [], [2.]),
            ("chisquare", [1., 2.], [2., 3.], [0.]),
            ("binomtest", [], [], [6., 5., .5, 0.]),
            ("fisher_exact", [1., 2., 3.], [], [0.]),
        ]:
            result = self.adapter.evaluate({"api": f"stats.{name}",
                                            "args": {"x": x, "y": y, "p": p}})
            self.assertIn("error", result["Boost.Math"], result)


if __name__ == "__main__":
    unittest.main()
