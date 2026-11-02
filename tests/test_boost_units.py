"""Full saved-case and independent SciPy tests for Boost.Units conversions.

Only public Boost.Units (Boost >= 1.83) temperature/unit/physical-constant
functionality is used by the C++ implementation. Four saved CODATA lookup
APIs mostly remain unsupported because Boost.Units' older reference tables
are NOT the same numerical source as SciPy CODATA 2022. The exact SI speed
of light is an independently verified shared value and uncertainty.
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

from scipy import constants as scipy_constants

from sci_testing.boost_units import (
    BoostUnitsComparisons, PROVIDER, SUPPORTED_APIS, NUMERIC_APIS, KEYED_APIS,
)
from sci_testing.case_csv import (
    read_case_file, read_measurement_content, csv_text,
)


ROOT = Path(__file__).resolve().parents[1]
ALL_SAVED = {
    "constants.convert_temperature",
    "constants.lambda2nu", "constants.nu2lambda",
    "constants.value", "constants.unit", "constants.precision",
    "constants.find",
}
EXPLICITLY_UNSUPPORTED = ALL_SAVED - SUPPORTED_APIS


class BoostUnitsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        temporary = tempfile.TemporaryDirectory()
        cls.addClassCleanup(temporary.cleanup)
        shared = Path(temporary.name) / "boost-units.so"
        env = dict(os.environ)
        env.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2",
            "-Wall", "-Wextra", "-Werror", "-shared", "-fPIC",
            str(ROOT / "native/boost_units.cpp"), "-o", str(shared),
        ], check=True, env=env)
        cls.adapter = BoostUnitsComparisons(ct.CDLL(str(shared)))

    @staticmethod
    def case(api, x):
        return {"api": api, "args": {"x": x}}

    def value(self, api, x):
        result = self.adapter.evaluate(self.case(api, x))
        self.assertEqual(set(result), {PROVIDER}, (api, x, result))
        self.assertNotIn("error", result[PROVIDER], (api, x, result))
        values = result[PROVIDER]["values"]
        self.assertEqual(len(values), 1)
        self.assertTrue(math.isfinite(values[0]), (api, x, values))
        return values[0]

    def test_provenance_and_complete_seven_api_inventory(self):
        saved = {
            f"constants.{p.name.split('.csv', 1)[0]}"
            for p in (ROOT / "cases/generated/constants").glob("*.csv*")
        }
        self.assertEqual(saved, ALL_SAVED)
        self.assertEqual(NUMERIC_APIS, {
            "constants.convert_temperature", "constants.lambda2nu",
            "constants.nu2lambda",
        })
        self.assertEqual(KEYED_APIS, {
            "constants.value", "constants.precision",
        })
        self.assertEqual(EXPLICITLY_UNSUPPORTED, {
            "constants.find", "constants.unit",
        })
        for api in NUMERIC_APIS:
            self.assertEqual(self.adapter.library.boost_units_supported(
                api.rsplit(".", 1)[1].encode()), 1)
        for api in KEYED_APIS:
            name = api.rsplit(".", 1)[1].encode()
            self.assertEqual(
                self.adapter.library.boost_units_key_supported(
                    name, b"speed of light in vacuum"), 1)
            self.assertEqual(self.adapter.library.boost_units_key_supported(
                name, b"Planck constant"), 0)
            self.assertEqual(self.adapter.evaluate({
                "api": api,
                "args": {"key": "Planck constant"},
            }), {})
        for api in EXPLICITLY_UNSUPPORTED:
            self.assertEqual(self.adapter.evaluate(
                {"api": api, "args": {"key": "Boltzmann constant"}}), {})
            self.assertEqual(self.adapter.library.boost_units_supported(
                api.rsplit(".", 1)[1].encode()), 0)

    def test_all_6252_saved_cases_against_independent_scipy_and_reference(self):
        start = time.perf_counter()
        count = 0
        per_api = {}
        for api in sorted(NUMERIC_APIS):
            name = api.rsplit(".", 1)[1]
            cases = read_case_file(
                ROOT / f"cases/generated/constants/{name}.csv",
                f"constants/{name}.csv",
            )
            path = ROOT / f"results/measurements/constants/{name}.csv"
            if not path.exists():
                path = path.with_suffix(".csv.gz")
            records = read_measurement_content(csv_text(path))
            self.assertEqual({case["id"] for case in cases}, set(records), api)
            max_absolute = 0.
            max_scaled = 0.
            for case in cases:
                with self.subTest(api=api, id=case["id"]):
                    x = float(case["args"]["x"])
                    actual = self.value(api, x)
                    if name == "convert_temperature":
                        independent_scipy = float(scipy_constants.convert_temperature(
                            x, "Celsius", "Fahrenheit"))
                    elif name == "lambda2nu":
                        independent_scipy = float(scipy_constants.lambda2nu(x))
                    else:
                        independent_scipy = float(scipy_constants.nu2lambda(x))
                    saved = records[case["id"]]
                    saved_sci = saved["libraries"]["Sci"]["values"]
                    self.assertEqual(len(saved_sci), 1)
                    self.assertEqual(len(saved["reference"]), 1)
                    baseline = saved_sci[0]
                    # Sci is the project's Rust implementation, and SciPy
                    # independently chooses its own floating expression order
                    # for affine temperatures. They may differ by one ULP.
                    self.assertTrue(math.isclose(
                        independent_scipy, baseline,
                        rel_tol=3e-15, abs_tol=1e-13,
                    ), (case["id"], independent_scipy, baseline))
                    error = abs(actual - independent_scipy)
                    scale = max(1., abs(independent_scipy))
                    max_absolute = max(max_absolute, error)
                    max_scaled = max(max_scaled, error/scale)
                    self.assertTrue(math.isclose(
                        actual, independent_scipy, rel_tol=4e-14,
                        abs_tol=2e-12,
                    ), (case["id"], x, actual, independent_scipy))
                    self.assertTrue(math.isclose(
                        actual, baseline, rel_tol=4e-14, abs_tol=2e-12,
                    ), (case["id"], x, actual, baseline))
                    # Separately verify against the independently generated
                    # arbitrary-precision mathematical reference as decimal.
                    oracle = Decimal(saved["reference"][0])
                    self.assertTrue(oracle.is_finite())
                    delta = abs(Decimal.from_float(actual) - oracle)
                    self.assertLessEqual(
                        delta,
                        max(Decimal("2e-12"), abs(oracle)*Decimal("4e-14")),
                        (case["id"], actual, str(oracle), str(delta)),
                    )
            per_api[name] = (len(cases), max_absolute, max_scaled)
            count += len(cases)
        elapsed = time.perf_counter() - start
        self.assertEqual(per_api["convert_temperature"][0], 2000)
        self.assertEqual(per_api["lambda2nu"][0], 2126)
        self.assertEqual(per_api["nu2lambda"][0], 2126)
        self.assertEqual(count, 6252)
        self.assertEqual(per_api["lambda2nu"][1], 0.)
        self.assertEqual(per_api["nu2lambda"][1], 0.)
        self.assertLess(elapsed, 30.)
        print(f"\nBoost.Units: {count} saved+independent SciPy checks in "
              f"{elapsed:.3f}s (~{count/elapsed:.0f}/s)")
        for name, (n, maxabs, maxscaled) in per_api.items():
            print(f"  constants.{name}: {n} passed; "
                  f"maxabs={maxabs:.3g}, maxscaled={maxscaled:.3g}")

    def test_exact_codataset_speed_of_light_keyed_subset(self):
        supported_count = 0
        for api in sorted(KEYED_APIS):
            name = api.rsplit(".", 1)[1]
            records_path = ROOT / f"results/measurements/constants/{name}.csv"
            if not records_path.exists():
                records_path = records_path.with_suffix(".csv.gz")
            records = read_measurement_content(csv_text(records_path))
            cases = read_case_file(
                ROOT / f"cases/generated/constants/{name}.csv",
                f"constants/{name}.csv",
            )
            self.assertEqual(len(cases), 2000)
            for case in cases:
                result = self.adapter.evaluate(case)
                key = case["args"]["key"]
                if key != "speed of light in vacuum":
                    self.assertEqual(result, {}, case["id"])
                    continue
                with self.subTest(api=api, id=case["id"]):
                    supported_count += 1
                    self.assertEqual(set(result), {PROVIDER})
                    value = result[PROVIDER]["values"][0]
                    expected = getattr(scipy_constants, name)(key)
                    saved = records[case["id"]]
                    self.assertEqual(value, expected)
                    self.assertEqual(value, saved["libraries"]["Sci"]["values"][0])
                    self.assertEqual(
                        Decimal.from_float(value),
                        Decimal(saved["reference"][0]),
                    )
        # The key appears four times in value and twice in precision fixtures.
        self.assertEqual(supported_count, 6)

    def test_absolute_temperature_reference_points(self):
        for celsius, fahrenheit in (
            (-273.15, -459.67), (-40., -40.), (0., 32.),
            (37., 98.6), (100., 212.), (1e4, 18032.),
        ):
            with self.subTest(celsius=celsius):
                value = self.value("constants.convert_temperature", celsius)
                self.assertTrue(math.isclose(
                    value, fahrenheit, rel_tol=4e-14, abs_tol=3e-12))

    def test_frequency_wavelength_dimensional_roundtrips(self):
        for x in (-1e150, -1., -1e-150, 1e-150, .1, 1., 1e150):
            with self.subTest(x=x):
                f = self.value("constants.lambda2nu", x)
                result = self.value("constants.nu2lambda", f)
                self.assertTrue(math.isclose(
                    result, x, rel_tol=2e-15, abs_tol=0.))
        # Speed of light is an exact physical constant in the Boost.Units
        # public SI CODATA catalogue; both conversion directions reuse it.
        self.assertEqual(self.value("constants.lambda2nu", 1.), 299792458.)
        self.assertEqual(self.value("constants.nu2lambda", 1.), 299792458.)

    def test_unsupported_domains_and_c_abi_does_not_mutate_output(self):
        for name, number in (
            ("lambda2nu", 0.), ("nu2lambda", -0.),
            ("convert_temperature", float("nan")),
            ("nu2lambda", float("inf")),
        ):
            with self.subTest(name=name, number=number):
                input_array = (ct.c_double * 1)(number)
                output_array = (ct.c_double * 1)(12345.5)
                status = self.adapter.library.boost_units_evaluate(
                    name.encode(), input_array, 1, output_array, 1)
                self.assertEqual(status, 1)
                self.assertEqual(output_array[0], 12345.5)
                self.assertEqual(self.adapter.evaluate(
                    self.case("constants."+name, number)), {})
        input_array = (ct.c_double * 1)(1.)
        output_array = (ct.c_double * 1)(12345.5)
        for name, size, capacity in (
            (b"find", 1, 1), (b"lambda2nu", 2, 1),
            (b"nu2lambda", 1, 0),
        ):
            with self.subTest(name=name, size=size, capacity=capacity):
                status = self.adapter.library.boost_units_evaluate(
                    name, input_array, size, output_array, capacity)
                self.assertEqual(status, -1)
                self.assertEqual(output_array[0], 12345.5)
        self.assertEqual(self.adapter.library.boost_units_evaluate_key(
            b"value", b"Planck constant", output_array, 1), -1)
        self.assertEqual(output_array[0], 12345.5)


if __name__ == "__main__":
    unittest.main()
