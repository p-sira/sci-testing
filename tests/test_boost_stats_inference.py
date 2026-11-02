"""Full saved-case SciPy verification of Boost statistics inference adapters."""
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

from sci_testing.boost_stats_inference import BoostStatsInferenceComparisons, NAMES

ROOT = Path(__file__).resolve().parent.parent
ALT = ("two-sided", "less", "greater")


def cases(name):
    path = ROOT / "cases" / "generated" / "stats" / f"{name}.csv.gz"
    with gzip.open(path, "rt", newline="") as f:
        for row in csv.DictReader(f):
            def decode_vector(key):
                tag, separator, tail = row[key].partition(":")
                if tag != "v" or not separator:
                    raise ValueError(f"{row['id']}: malformed {key}")
                return [float(value) for value in tail.split()]

            group_tag, _, raw_groups = row["g"].partition(":")
            if group_tag == "v":
                if raw_groups.strip():
                    raise ValueError(f"{row['id']}: invalid group vector")
                groups = []
            elif group_tag == "m":
                groups = [[float(value) for value in group.split()]
                          for group in raw_groups.split(";")]
            else:
                raise ValueError(f"{row['id']}: invalid g")
            yield row["id"], {
                "api": f"stats.{name}",
                "args": {"x": decode_vector("x"), "y": decode_vector("y"),
                         "p": decode_vector("p"), "g": groups},
            }


def independent(name, args):
    x, p, g = args["x"], args["p"], args["g"]
    if name == "ttest_ind_from_stats":
        r = stats.ttest_ind_from_stats(
            *p[:6], equal_var=bool(p[6]), alternative=ALT[int(p[7])])
        _, s1, n1, _, s2, n2 = p[:6]
        a = s1*s1/n1
        b = s2*s2/n2
        df = (n1+n2-2 if p[6] else
              (a+b)**2/(a*a/(n1-1)+b*b/(n2-1)))
        return [r.statistic, r.pvalue, df]
    if name in {"f_oneway", "bartlett"}:
        result = getattr(stats, name)(*g)
        return [result.statistic, result.pvalue]
    opts = {"low": p[0], "high": p[1]}
    if name == "circstd":
        opts["normalize"] = bool(p[2])
    return [getattr(stats, name)(x, **opts)]


class BoostStatsInferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(tmp.cleanup)
        output = Path(tmp.name) / "boost-stats-inference.so"
        env = os.environ.copy()
        env.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2", "-shared", "-fPIC",
            str(ROOT / "native/boost_stats_inference.cpp"), "-o", str(output)
        ], check=True, env=env)
        cls.native = BoostStatsInferenceComparisons(ct.CDLL(str(output)))

    def test_every_saved_case_compared_to_scipy(self):
        errors = []
        counts = {}
        for name in sorted(NAMES):
            count = 0
            for case_id, case in cases(name):
                count += 1
                result = self.native.evaluate(case)
                if set(result) != {"Boost.Math"} or "error" in result["Boost.Math"]:
                    errors.append((case_id, result))
                    continue
                actual = np.asarray(result["Boost.Math"]["values"], dtype=float)
                expected = np.asarray(independent(name, case["args"]), dtype=float)
                if actual.shape != expected.shape or not np.allclose(
                        actual, expected, rtol=2e-10, atol=2e-12):
                    errors.append((case_id, actual.tolist(), expected.tolist()))
            counts[name] = count
            self.assertGreaterEqual(count, 2000, name)
        self.assertFalse(errors, f"{len(errors)} mismatches; examples {errors[:15]}")
        self.assertEqual(sum(counts.values()), 12020)

    def test_native_abi_and_invalid_cases(self):
        native = self.native.library
        self.assertEqual(native.boost_stats_inference_output_count(b"f_oneway"), 2)
        self.assertEqual(native.boost_stats_inference_output_count(b"ttest_ind_from_stats"), 3)
        self.assertEqual(native.boost_stats_inference_output_count(b"circmean"), 1)
        self.assertEqual(native.boost_stats_inference_supported(1, b"f_oneway"), 0)
        self.assertEqual(self.native.evaluate({"api": "stats.unknown", "args": {}}), {})
        for name, args in [
            ("ttest_ind_from_stats", {"p": [1,1,0,2,1,6,0,0]}),
            ("f_oneway", {"g": [[1.], [2.,3.]]}),
            ("bartlett", {"g": [[1.,1.], [2.,2.]]}),
            ("circmean", {"x": [1.,2.], "p": [1.,1.]}),
            ("circstd", {"x": [1.,2.], "p": [0.,2.,2.]}),
        ]:
            result = self.native.evaluate({"api": f"stats.{name}", "args": args})
            self.assertIn("error", result["Boost.Math"], (name,result))

    def test_nonstandard_circular_boundaries(self):
        for name in ("circmean", "circvar", "circstd"):
            for values, low, high in (([3.,4.,6.], -2., 8.),
                                      ([0.,0.,0.],0.,2*math.pi),
                                      ([5.,5.,6.], -math.pi, math.pi)):
                for normalize in ((False,True) if name=="circstd" else (False,)):
                    p=[low,high]+([float(normalize)] if name=="circstd" else [])
                    case={"api":f"stats.{name}", "args":{"x":values,"p":p}}
                    native=self.native.evaluate(case)["Boost.Math"]
                    self.assertNotIn("error",native,(name,native))
                    expected=independent(name,{"x":values,"p":p,"g":[]})
                    np.testing.assert_allclose(native["values"],expected,
                                               rtol=1e-10,atol=2e-12)


if __name__ == "__main__":
    unittest.main()
