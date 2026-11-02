"""Direct Boost.Math and GSL spherical-special ABI tests against SciPy."""

import csv
import ctypes as ct
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

from scipy import special

from sci_testing.spherical_extra import SphericalExtraComparisons


ROOT = Path(__file__).resolve().parents[1]
LIBRARIES = {"Boost.Math", "GSL"}


class SphericalExtraTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        directory = tempfile.TemporaryDirectory()
        cls.addClassCleanup(directory.cleanup)
        shared = Path(directory.name) / "spherical-extra.so"
        environment = dict(os.environ)
        environment.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2", "-Wall",
            "-Wextra", "-Werror", "-shared", "-fPIC",
            str(ROOT / "native/spherical_extra.cpp"),
            "-lgsl", "-lgslcblas", "-o", str(shared),
        ], check=True, env=environment)
        cls.library = ct.CDLL(str(shared))
        cls.native = SphericalExtraComparisons(cls.library)

    @staticmethod
    def case(name, args):
        return {"api": f"special.{name}",
                "args": {f"p{i}": value for i, value in enumerate(args)}}

    def result(self, name, args):
        result = self.native.evaluate(self.case(name, args))
        self.assertEqual(set(result), LIBRARIES, (name, args, result))
        for label, entry in result.items():
            self.assertNotIn("error", entry, (name, args, label, entry))
            self.assertEqual(len(entry["values"]), 2 if name == "sph_harm_y" else 1)
        return result

    def assert_spherical_close(self, name, args, result, *, tolerance=1e-12):
        expected = (complex(special.sph_harm_y(*args)) if name == "sph_harm_y"
                    else float(special.sph_legendre_p(*args)[0]))
        for label, entry in result.items():
            actual = (complex(*entry["values"]) if name == "sph_harm_y"
                      else entry["values"][0])
            self.assertLessEqual(abs(actual - expected), tolerance,
                                 (name, args, label, actual, expected))

    def test_public_symbols_arity_and_output_shape(self):
        bridge = self.library
        for lib in (0, 1):
            for name, arity, shape in ((b"sph_harm_y", 4, 2),
                                       (b"sph_legendre_p", 3, 1)):
                self.assertEqual(bridge.spherical_extra_supported(lib, name), arity)
                self.assertEqual(bridge.spherical_extra_output_count(name), shape)
            self.assertEqual(bridge.spherical_extra_supported(lib, b"unknown"), 0)
        self.assertEqual(bridge.spherical_extra_supported(2, b"sph_harm_y"), 0)
        self.assertEqual(bridge.spherical_extra_output_count(b"unknown"), 0)
        self.assertEqual(self.native.evaluate(self.case("unknown", [1])), {})

    def test_all_2001_saved_spherical_cases_per_api(self):
        for name, arity in (("sph_harm_y", 4), ("sph_legendre_p", 3)):
            with (ROOT / f"cases/generated/special/{name}.csv").open(
                    newline="") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 2001, name)
            for row in rows:
                args = tuple(
                    int(row[f"p{i}"][2:]) if row[f"p{i}"].startswith("i:")
                    else float(row[f"p{i}"][2:]) for i in range(arity))
                with self.subTest(case=row["id"]):
                    self.assert_spherical_close(name, args, self.result(name, args),
                                                tolerance=2e-12)

    def test_independent_degree_order_and_angle_grid(self):
        for n in range(7):
            for m in range(-n-1, n+2):
                for theta, phi in ((0., .3), (.2, 1.4), (.5, .3),
                                   (math.pi/2, 3.5), (math.pi, 2*math.pi),
                                   (2.7, -1.7)):
                    with self.subTest(n=n, m=m, theta=theta, phi=phi):
                        harmonic = self.result("sph_harm_y", (n, m, theta, phi))
                        legendre = self.result("sph_legendre_p", (n, m, theta))
                        self.assert_spherical_close("sph_harm_y", (n,m,theta,phi), harmonic)
                        self.assert_spherical_close("sph_legendre_p", (n,m,theta), legendre)
                        # The real normalized Legendre amplitude times the
                        # complex azimuthal phase reconstructs Y_n^m.
                        phase = complex(math.cos(m*phi), math.sin(m*phi))
                        for label in LIBRARIES:
                            y = complex(*harmonic[label]["values"])
                            amp = legendre[label]["values"][0]
                            self.assertLessEqual(abs(y - amp*phase), 3e-12)

    def test_negative_order_condon_shortley_and_conjugation(self):
        for n in range(1, 8):
            for m in range(1, n+1):
                for theta, phi in ((.3, 1.2), (1.5, 5.8)):
                    y_pos = self.result("sph_harm_y", (n,m,theta,phi))
                    y_neg = self.result("sph_harm_y", (n,-m,theta,phi))
                    p_pos = self.result("sph_legendre_p", (n,m,theta))
                    p_neg = self.result("sph_legendre_p", (n,-m,theta))
                    for label in LIBRARIES:
                        with self.subTest(n=n, m=m, theta=theta, library=label):
                            self.assertAlmostEqual(p_neg[label]["values"][0],
                                                   (-1)**m*p_pos[label]["values"][0],
                                                   delta=1e-12)
                            self.assertAlmostEqual(
                                abs(complex(*y_neg[label]["values"]) -
                                    (-1)**m*complex(*y_pos[label]["values"]).conjugate()),
                                0, delta=1e-12)

    def test_high_degrees_and_polar_precision_boundary(self):
        # Large |m| exercises the factorial ratio (n-|m|)!/(n+|m|)!,
        # while the two libraries remain independent for moderate angles.
        for n, m, theta in ((12,-5,.1), (30,10,2.), (50,-25,1.),
                            (100,10,.3), (128,128,1.5), (128,-128,1.5),
                            (200,30,1.3)):
            for name, args in (("sph_harm_y", (n,m,theta,.6)),
                               ("sph_legendre_p", (n,m,theta))):
                with self.subTest(name=name, n=n, m=m):
                    self.assert_spherical_close(name, args, self.result(name, args),
                                                tolerance=2e-11)

        # Using x=cos(theta) can round small *nonzero* theta to the exact
        # pole. GSL's P_l^m(x) would then silently output 0 for m>0,
        # whereas Boost can use theta directly and agree with SciPy.
        for name, args in (("sph_harm_y", (3,1,1e-8,.6)),
                           ("sph_legendre_p", (3,-1,1e-8))):
            with self.subTest(name=name):
                result = self.native.evaluate(self.case(name,args))
                self.assertIn("polar-angle precision", result["GSL"]["error"])
                self.assertNotIn("error", result["Boost.Math"])
                expected = (complex(special.sph_harm_y(*args)) if name == "sph_harm_y"
                            else float(special.sph_legendre_p(*args)[0]))
                actual = (complex(*result["Boost.Math"]["values"])
                          if name == "sph_harm_y" else result["Boost.Math"]["values"][0])
                self.assertAlmostEqual(actual.real if isinstance(actual,complex) else actual,
                                       expected.real if isinstance(expected,complex) else expected,
                                       delta=1e-13)

    def test_rejects_bad_arguments_and_output_capacity(self):
        for name, args in (
            ("sph_harm_y", (-1, 0, .5, .2)),
            ("sph_harm_y", (3.5, 1, .5, .2)),
            ("sph_harm_y", (3, 1.5, .5, .2)),
            ("sph_harm_y", (3, 1, float("nan"), .2)),
            ("sph_harm_y", (3, 1, .5, float("inf"))),
            ("sph_harm_y", (3, 1, math.pi+.1, .2)),
            ("sph_legendre_p", (2, 1, -.5)),
            ("sph_legendre_p", (2000, 1, .5)),
        ):
            actual = self.native.evaluate(self.case(name, args))
            for label in LIBRARIES:
                with self.subTest(name=name, args=args, library=label):
                    self.assertTrue(actual[label]["error"])

        # A missing required argument is reported as an explicit error.
        actual = self.native.evaluate(self.case("sph_harm_y", (3,1,.5)))
        for label in LIBRARIES:
            self.assertIn("expected 4", actual[label]["error"])

        parameters = (ct.c_double * 4)(3, 1, .5, .3)
        too_small = (ct.c_double * 1)()
        error = self.library.spherical_extra_evaluate(
            0, b"sph_harm_y", parameters, 4, too_small, 1)
        self.assertNotEqual(error, 0)
        self.assertIn(b"capacity", self.library.spherical_extra_error())


if __name__ == "__main__":
    unittest.main()
