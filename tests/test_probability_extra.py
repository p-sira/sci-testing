"""SciPy cross-validation for additional native Boost.Math/GSL probability APIs."""

import ctypes
import csv
import gzip
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
import warnings

from scipy import special

from sci_testing.probability_extra import ProbabilityExtraComparisons


ROOT = Path(__file__).resolve().parents[1]


class ProbabilityExtraTests(unittest.TestCase):
    ALL = {
        "nbdtr": 3, "nbdtrc": 3, "nbdtri": 3,
        "chndtr": 3, "chndtrix": 3, "chndtridf": 3, "chndtrinc": 3,
        "ncfdtr": 4, "ncfdtri": 4, "ncfdtrinc": 4,
        "ncfdtridfd": 4, "ncfdtridfn": 4,
        "nctdtr": 3, "nctdtrit": 3, "nctdtrinc": 3,
        "nctdtridf": 3,
        "btdtria": 3, "btdtrib": 3, "bdtrin": 3,
        "nbdtrik": 3, "nbdtrin": 3, "pdtrik": 2,
        "gdtrib": 3, "chdtriv": 2,
        "kolmogorov": 1, "kolmogi": 1,
    }

    @classmethod
    def setUpClass(cls):
        cls.tempdir = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tempdir.cleanup)
        library = Path(cls.tempdir.name) / "probability_extra.so"
        env = dict(os.environ)
        env.pop("LD_LIBRARY_PATH", None)  # AppImage loader can break system C++ compiler.
        subprocess.run(["c++", "-std=c++17", "-O1", "-shared", "-fPIC",
                        str(ROOT / "native/probability_extra.cpp"),
                        "-lgsl", "-lgslcblas", "-o", str(library)],
                       check=True, env=env)
        cls.comparison = ProbabilityExtraComparisons(ctypes.CDLL(str(library)))

    @staticmethod
    def case(name, args):
        return {"api": f"special.{name}",
                "args": {f"p{i}": value for i, value in enumerate(args)}}

    def test_exact_coverage_and_exclusions(self):
        bridge = self.comparison.library
        for name, arity in self.ALL.items():
            self.assertEqual(bridge.probability_extra_supported(0, name.encode()), arity)
            self.assertEqual(bridge.probability_extra_supported(1, name.encode()),
                             arity if name in ("nbdtr", "nbdtrc", "nbdtri") else 0)
        for name in ("bad_name",):
            for library in (0, 1):
                self.assertEqual(bridge.probability_extra_supported(library, name.encode()), 0)
        self.assertEqual(self.comparison.evaluate(self.case("bad_name", [1])), {})

    def test_degree_of_freedom_inversion_uses_portable_boost_cdf_only(self):
        # Ubuntu's older Boost distributions do not have the newer static
        # df inverses. Keep this a hard regression guard even when the test
        # runner itself happens to be built with recent Boost headers.
        source = (ROOT / "native/probability_extra.cpp").read_text()
        for method in ("::find_v1(", "::find_v2(", "::find_degrees_of_freedom("):
            self.assertNotIn(method, source)
        self.assertIn("invert_degree_of_freedom", source)

    def test_direct_api_values_against_scipy(self):
        values = {
            "nbdtr": [(2, 5, .35), (7, 2.5, .72)],
            "nbdtrc": [(2, 5, .35), (7, 2.5, .72)],
            "nbdtri": [(2, 5, .35), (7, 2.5, .72)],
            "chndtr": [(3.5, 5, .7), (12, 3.5, 4)],
            "chndtrix": [(.35, 5, .7), (.7, 3.5, 4)],
            "ncfdtr": [(3, 8, .7, 1.5), (10, 5, 3, 2.1)],
            "ncfdtri": [(3, 8, .7, .35), (10, 5, 3, .75)],
            "ncfdtridfd": [(3, .25637506749643246, .5, .5)],
            "ncfdtridfn": [(.25637506749643246, 5, .5, .5)],
            "nctdtr": [(3, .5, .3), (8, -2.1, -1.5)],
            "nctdtrit": [(3, .5, .35), (8, -2.1, .75)],
            "nctdtridf": [(.4120716320853016, .5, .3)],
            "btdtria": [(.35, 2.5, .5), (.75, 5, .2)],
            "btdtrib": [(2.5, .35, .5), (5, .75, .2)],
            "bdtrin": [(2, .35, .5), (5, .6, .25)],
            "nbdtrik": [(.35, 5, .5), (.65, 3, .8)],
            "nbdtrin": [(2, .35, .5), (5, .6, .25)],
            "pdtrik": [(.35, 5), (.6, 10)],
            "gdtrib": [(1.5, .35, 2.5), (.5, .6, 10)],
            "chdtriv": [(.35, 2.5), (.6, 10)],
            "kolmogorov": [(.25,), (1.4,), (2.9,)],
            "kolmogi": [(.3,), (.65,), (.95,)],
        }
        for name, samples in values.items():
            for args in samples:
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", message="floating point number truncated")
                    expected = float(getattr(special, name)(*args))
                self.assertTrue(math.isfinite(expected), (name, args, expected))
                actual = self.comparison.evaluate(self.case(name, args))
                self.assertIn("Boost.Math", actual)
                for library, value in actual.items():
                    with self.subTest(name=name, args=args, library=library):
                        self.assertNotIn("error", value, value)
                        self.assertTrue(math.isclose(value["values"][0], expected,
                                                      rel_tol=(5e-5 if name == "ncfdtridfd" else
                                                               5e-6 if name == "ncfdtridfn" else 4e-7),
                                                      abs_tol=4e-9),
                                        (name, args, library, value, expected))

    def test_generated_fixture_parameters(self):
        # Check case wiring and representative real generated values in addition
        # to the independent hand-picked inputs above. Deliberately inspect
        # several positions in each sweep, rather than a single first row.
        for name, arity in self.ALL.items():
            folder = ROOT / "cases/generated/special"
            path = folder / f"{name}.csv"
            if not path.exists():
                path = folder / f"{name}.csv.gz"
            self.assertTrue(path.exists(), name)
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "rt", newline="") as stream:
                rows = list(csv.DictReader(stream))
            self.assertGreater(len(rows), 10, name)
            for row in (rows[0], rows[len(rows)//4], rows[len(rows)//2], rows[-1]):
                args = tuple(float(row[f"p{i}"].split(":", 1)[1]) for i in range(arity))
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", message="floating point number truncated")
                    expected = float(getattr(special, name)(*args))
                if not math.isfinite(expected):
                    continue
                observed = self.comparison.evaluate(self.case(name, args))
                with self.subTest(name=name, args=args):
                    self.assertIn("Boost.Math", observed)
                    for library, result in observed.items():
                        self.assertNotIn("error", result, (name, args, library, result))
                        self.assertTrue(math.isclose(result["values"][0], expected,
                                                      rel_tol=5e-5, abs_tol=1e-6),
                                        (name, args, library, result, expected))

    def test_all_compressed_degree_of_freedom_cases(self):
        # SciPy inverse accuracy is lower than forward CDF accuracy: e.g. at
        # (3, .25637506749643246, .5, .5), SciPy reports 5.00021997982672
        # while Boost finds exactly 5.0. Assert the same branch at SciPy's
        # precision, and independently enforce a tight forward CDF residual.
        expectations = {
            "ncfdtridfd": (4.5e-5, 2e-13),
            "ncfdtridfn": (4e-6, 2e-13),
            "nctdtridf": (1e-7, 2e-13),
        }
        for name, (relative_tolerance, residual_tolerance) in expectations.items():
            with gzip.open(ROOT / f"cases/generated/special/{name}.csv.gz", "rt") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 2001, name)
            for row in rows:
                argc = self.ALL[name]
                args = tuple(float(row[f"p{i}"].split(":", 1)[1]) for i in range(argc))
                actual = self.comparison.evaluate(self.case(name, args))
                with self.subTest(api=name, id=row["id"]):
                    self.assertIn("Boost.Math", actual)
                    result = actual["Boost.Math"]
                    self.assertNotIn("error", result, result)
                    value = result["values"][0]
                    reference = float(getattr(special, name)(*args))
                    self.assertTrue(math.isfinite(reference))
                    self.assertTrue(math.isclose(value, reference,
                                                  rel_tol=relative_tolerance, abs_tol=1e-8),
                                    (name, row["id"], value, reference))
                    if name == "ncfdtridfd":
                        target, forward = args[1], special.ncfdtr(args[0], value, args[2], args[3])
                    elif name == "ncfdtridfn":
                        target, forward = args[0], special.ncfdtr(value, args[1], args[2], args[3])
                    else:
                        target, forward = args[0], special.nctdtr(value, args[1], args[2])
                    self.assertLessEqual(abs(float(forward)-target), residual_tolerance)

    def test_inverted_noncentral_parameters_against_scipy(self):
        # Construct valid inversion targets from the forward functions to avoid
        # asking any solver to infer a non-existent or ambiguous parameter.
        cases = [
            ("chndtridf", (4.2, special.chndtr(4.2, 3.3, .8), .8), 3.3),
            ("chndtridf", (10, special.chndtr(10, 5, 2), 2), 5),
            ("chndtrinc", (4.2, 3.3, special.chndtr(4.2, 3.3, .8)), .8),
            ("chndtrinc", (10, 5, special.chndtr(10, 5, 2)), 2),
            ("ncfdtrinc", (3, 8, special.ncfdtr(3, 8, .7, 1.5), 1.5), .7),
            ("ncfdtrinc", (10, 5, special.ncfdtr(10, 5, 3, 2.1), 2.1), 3),
            ("nctdtrinc", (3, special.nctdtr(3, .5, .3), .3), .5),
            ("nctdtrinc", (8, special.nctdtr(8, -2.1, -1.5), -1.5), -2.1),
        ]
        for name, args, expected_parameter in cases:
            scipy_inverse = float(getattr(special, name)(*args))
            actual = self.comparison.evaluate(self.case(name, args))
            with self.subTest(name=name, args=args):
                self.assertIn("Boost.Math", actual)
                self.assertNotIn("error", actual["Boost.Math"], actual)
                value = actual["Boost.Math"]["values"][0]
                self.assertAlmostEqual(value, expected_parameter, delta=1e-5)
                self.assertAlmostEqual(value, scipy_inverse, delta=3e-4)

    def test_errors_are_explicit(self):
        for name, args in (("nbdtr", (2, -3, .4)),
                           ("ncfdtr", (4, -1, 1, 2)),
                           ("btdtria", (.3, -2, .4)),
                           ("chndtrinc", (1, 3, .999999999)),
                           ("kolmogi", (1.2,))):
            results = self.comparison.evaluate(self.case(name, args))
            for library, result in results.items():
                with self.subTest(name=name, library=library):
                    self.assertIn("error", result)
        self.assertIn("error", self.comparison.evaluate(
            {"api": "special.nbdtri", "args": {"p0": 2}})["Boost.Math"])

    def test_discrete_and_kolmogorov_boundaries(self):
        for name, args, expected in (
            ("nbdtr", (2, 5, 0), 0),
            ("nbdtr", (2, 5, 1), 1),
            ("nbdtrc", (2, 5, 0), 1),
            ("nbdtrc", (2, 5, 1), 0),
            ("nbdtri", (2, 5, 0), 0),
            ("nbdtri", (2, 5, 1), 1),
            ("nbdtrik", (.01, 5, .5), 0),
            ("nbdtrik", (.001, 5, .5), 0),
            ("kolmogorov", (-1,), 1),
            ("kolmogorov", (0,), 1),
            ("kolmogi", (1,), 0),
        ):
            values = self.comparison.evaluate(self.case(name, args))
            for label, result in values.items():
                with self.subTest(name=name, args=args, library=label):
                    self.assertEqual(result, {"values": [expected]})

    def test_ambiguous_noncentral_f_df_inverse_is_explicit_error(self):
        # SciPy documents non-monotonicity in *either* noncentral-F degree
        # of freedom; both constructed targets cross the forward CDF twice.
        cases = (
            ("ncfdtridfd", (5, .01, 0, .05),
             lambda v: special.ncfdtr(5, v, 0, .05), .01),
            ("ncfdtridfn", (.14, 1, 3, .5),
             lambda v: special.ncfdtr(v, 1, 3, .5), .14),
        )
        for name, args, forward, target in cases:
            with self.subTest(name=name):
                scanned = [float(forward(v)) for v in (1e-6, .01, .1, 1, 10, 1e5)]
                signs = [v > target for v in scanned]
                self.assertGreaterEqual(sum(a != b for a, b in zip(signs, signs[1:])), 2)
                result = self.comparison.evaluate(self.case(name, args))
                self.assertIn("ambiguous", result["Boost.Math"]["error"])

    def test_nonidentifiable_t_df_returns_error_not_arbitrary_value(self):
        # For nc=t=0, CDF_t(0)=.5 regardless of degrees of freedom. SciPy
        # may return an arbitrary positive number; no unique inverse exists.
        args = (.5, 0, 0)
        for df in (.001, 1, 10, 1e5):
            self.assertAlmostEqual(float(special.nctdtr(df, 0, 0)), .5)
        result = self.comparison.evaluate(self.case("nctdtridf", args))
        self.assertIn("error", result["Boost.Math"])


if __name__ == "__main__":
    unittest.main()
