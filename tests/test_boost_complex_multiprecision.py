"""Full saved-case audit of genuine public Boost.Multiprecision complex math.

Tests 34 complex APIs and every generated case, including signed-zero branch
banks. A native 'unavailable' status is never counted as successful coverage.
For supported samples, independently check BOTH pinned Sci f64 outputs and the
saved high-precision mathematical reference decimal, without calling NumPy.
"""

import ctypes as ct
from decimal import Decimal
import math
import os
from pathlib import Path
import subprocess
import tempfile
import time
import unittest

import mpmath as mp

from sci_testing.boost_complex_multiprecision import (
    BoostComplexMultiprecisionComparisons, NAMES, PAIR_NAMES, SCALAR_NAMES,
    PROVIDER,
)
from sci_testing.case_csv import (
    read_case_file, read_measurement_content, csv_text,
)


ROOT = Path(__file__).resolve().parents[1]
NUMERIC_REL = Decimal("2e-10")
NUMERIC_ABS = Decimal("2e-13")


class BoostComplexMultiprecisionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        temporary = tempfile.TemporaryDirectory()
        cls.addClassCleanup(temporary.cleanup)
        shared = Path(temporary.name) / "boost-complex-mp.so"
        environment = os.environ.copy()
        environment.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2",
            "-Wall", "-Wextra", "-Werror", "-shared", "-fPIC",
            str(ROOT / "native/boost_complex_multiprecision.cpp"),
            "-o", str(shared),
        ], env=environment, check=True)
        cls.adapter = BoostComplexMultiprecisionComparisons(ct.CDLL(str(shared)))

    def test_exact_registration_and_signature_shapes(self):
        self.assertEqual(len(NAMES), 34)
        self.assertEqual(len(SCALAR_NAMES), 5)
        self.assertEqual(len(PAIR_NAMES), 29)
        self.assertEqual(NAMES, {
            p.name.split(".csv", 1)[0] for p in
            (ROOT / "cases/generated/complex").glob("*.csv*")
        })
        for name in NAMES:
            count = self.adapter.library.boost_complex_mp_supported(name.encode())
            self.assertEqual(count, 1 if name in SCALAR_NAMES else 2, name)
        self.assertEqual(self.adapter.library.boost_complex_mp_supported(
            b"not-an-api"), 0)
        self.assertEqual(self.adapter.evaluate({"api": "special.erf", "args": {}}), {})
        self.assertEqual(self.adapter.evaluate({
            "api": "complex.add", "args": {"p0": [1, 2]}
        }), {})

    def test_full_70919_saved_sci_and_independent_references(self):
        totals = {
            "total": 0, "supported": 0,
            "unavailable": 0, "unexpected_native_error": 0,
            "failed_numeric": 0, "failed_signed_zero": 0,
        }
        per_api = {}
        started = time.perf_counter()
        for path in sorted((ROOT / "cases/generated/complex").glob("*.csv*")):
            name = path.name.split(".csv", 1)[0]
            cases = read_case_file(path, f"complex/{name}.csv")
            measurement = ROOT / "results/measurements/complex" / f"{name}.csv"
            if not measurement.exists():
                measurement = measurement.with_suffix(".csv.gz")
            records = read_measurement_content(csv_text(measurement))
            self.assertEqual({case["id"] for case in cases}, set(records), name)
            counts = {"total": 0, "supported": 0, "unavailable": 0}
            for case in cases:
                counts["total"] += 1
                record = records[case["id"]]
                actual = self.adapter.evaluate(case)
                if not actual:
                    counts["unavailable"] += 1
                    continue
                self.assertEqual(set(actual), {PROVIDER}, case["id"])
                response = actual[PROVIDER]
                self.assertNotIn("error", response, (case["id"], response))
                self.assertIn("values", response, case["id"])
                observed = response["values"]
                sci = record["libraries"]["Sci"]
                self.assertIn("values", sci, (case["id"], "saved Sci failed"))
                baseline = sci["values"]
                reference = record["reference"]
                self.assertEqual(len(observed), len(baseline), case["id"])
                self.assertEqual(len(observed), len(reference), case["id"])
                for i, (value, target, oracle) in enumerate(
                    zip(observed, baseline, reference)
                ):
                    self.assertTrue(math.isfinite(value), (case["id"], i))
                    self.assertTrue(math.isclose(value, target, rel_tol=2e-10,
                                                 abs_tol=2e-13),
                                    (case["id"], i, value, target))
                    # The reference is a 100-165-digit decimal from independent
                    # mpmath mathematics, NOT a rounded Sci double.
                    reference_decimal = Decimal(oracle)
                    self.assertTrue(reference_decimal.is_finite(), case["id"])
                    observed_decimal = Decimal.from_float(value)
                    delta = abs(observed_decimal-reference_decimal)
                    threshold = max(
                        NUMERIC_ABS,
                        NUMERIC_REL*abs(reference_decimal),
                    )
                    self.assertLessEqual(
                        delta, threshold,
                        (case["id"], i, str(observed_decimal),
                         str(reference_decimal), str(delta)),
                    )
                    if value == target == 0:
                        self.assertEqual(
                            math.copysign(1., value), math.copysign(1., target),
                            (case["id"], i, "signed-zero mismatch"),
                        )
                counts["supported"] += 1
            self.assertEqual(counts["total"], len(cases))
            per_api[name] = counts
            for label in counts:
                totals[label] += counts[label]
        elapsed = time.perf_counter()-started
        self.assertEqual(len(per_api), 34)
        self.assertEqual(totals["total"], 70_919)
        self.assertEqual(
            totals["supported"] + totals["unavailable"], totals["total"])
        self.assertGreater(totals["supported"], 60_000)
        self.assertLess(elapsed, 120.)
        print(
            f"\nBoost cpp_complex_50 full saved audit: "
            f"{totals['supported']}/{totals['total']} supported, "
            f"{totals['unavailable']} unavailable; "
            f"runtime {elapsed:.3f}s "
            f"({totals['total']/elapsed:.0f} case checks/s)"
        )
        for name, counts in per_api.items():
            print(
                f"  complex.{name}: {counts['supported']}/"
                f"{counts['total']} supported"
            )

    def test_independent_high_precision_transcendental_spot_checks(self):
        samples = [
            (.499, .2495), (-.499, .2495),
            (.12, -.45), (2.5, .3), (-2.3, -.8),
        ]
        methods = {
            "acos": mp.acos, "acosh": mp.acosh,
            "asin": mp.asin, "asinh": mp.asinh,
            "atan": mp.atan, "atanh": mp.atanh,
            "cos": mp.cos, "cosh": mp.cosh,
            "exp": mp.exp, "ln": mp.log,
            "sin": mp.sin, "sinh": mp.sinh,
            "sqrt": mp.sqrt, "tan": mp.tan, "tanh": mp.tanh,
            "cbrt": lambda z: mp.root(z, 3),
            "reciprocal": lambda z: 1/z,
        }
        with mp.workdps(100):
            for name, method in methods.items():
                for real, imag in samples:
                    case = {
                        "api": "complex." + name,
                        "args": {"x": real, "y": imag} if name in ("exp","sqrt")
                                else {"p0": [real, imag]},
                    }
                    response = self.adapter.evaluate(case)
                    self.assertIn(PROVIDER, response, (name, real, imag))
                    self.assertIn("values", response[PROVIDER])
                    result = response[PROVIDER]["values"]
                    target = method(mp.mpc(real, imag))
                    for got, oracle in zip(result, (target.real, target.imag)):
                        self.assertTrue(math.isclose(
                            got, float(oracle), rel_tol=2e-11,
                            abs_tol=2e-13), (name, real, imag, got, oracle))

    def test_unsupported_domains_and_error_no_output_mutation(self):
        examples = [
            ("ln", (-1., -0.)),
            ("sqrt", (-1., 0.)),
            ("atan", (1., -5e-324)),
            ("asinh", (0., -2.)),
            ("acos", (1., 0.)),
            ("reciprocal", (0., 0.)),
            ("divide", (1., 1., 0., 0.)),
            ("powu", (1., 1., -2.)),
            ("powi", (1., 1., 1.2)),
            ("exp", (1e13, 0.)),
            ("sin", (float("inf"), 1.)),
        ]
        for name, args in examples:
            with self.subTest(name=name, args=args):
                inp = (ct.c_double * len(args))(*args)
                output = (ct.c_double * 2)(123.5, -99.25)
                status = self.adapter.library.boost_complex_mp_evaluate(
                    name.encode(), inp, len(args), output, 2)
                self.assertEqual(status, 1)
                self.assertEqual(list(output), [123.5, -99.25])
        inp = (ct.c_double * 2)(1., 2.)
        output = (ct.c_double * 2)(123.5, -99.25)
        self.assertEqual(self.adapter.library.boost_complex_mp_evaluate(
            b"acos", inp, 1, output, 2), -1)
        self.assertEqual(list(output), [123.5, -99.25])


if __name__ == "__main__":
    unittest.main()
