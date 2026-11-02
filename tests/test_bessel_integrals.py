"""Independent numerical tests for the Bessel-integrals extension."""
import csv
import ctypes as ct
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
from scipy import special

from sci_testing.bessel_integrals import BesselIntegralComparisons

ROOT = Path(__file__).resolve().parent.parent
NAMES = ("besselpoly", "it2i0k0", "it2j0y0", "iti0k0", "itj0y0")


class BesselIntegralTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(tmp.cleanup)
        output = Path(tmp.name) / "bessel-integrals.so"
        env = os.environ.copy()
        env.pop("LD_LIBRARY_PATH", None)
        subprocess.run(
            [os.environ.get("CXX", "c++"), "-std=c++17", "-O2",
             "-fPIC", "-shared", str(ROOT / "native/bessel_integrals.cpp"),
             "-lgsl", "-lgslcblas", "-o", str(output)],
            env=env, check=True)
        cls.native = BesselIntegralComparisons(ct.CDLL(str(output)))

    @staticmethod
    def case(name, *args):
        return {"api": f"special.{name}",
                "args": {f"p{i}": x for i, x in enumerate(args)}}

    def check(self, name, *args, rtol=3e-10, atol=2e-12):
        expected = np.asarray(getattr(special, name)(*args), dtype=float).ravel()
        results = self.native.evaluate(self.case(name, *args))
        self.assertEqual(set(results), {"Boost.Math", "GSL"})
        for lib, result in results.items():
            with self.subTest(name=name, args=args, library=lib):
                self.assertNotIn("error", result, result)
                np.testing.assert_allclose(result["values"], expected,
                                           rtol=rtol, atol=atol)

    def test_all_saved_inputs(self):
        for name in NAMES:
            with (ROOT / f"cases/generated/special/{name}.csv").open(
                    newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2001, name)
            for row in rows:
                args = tuple(float(row[f"p{i}"].split(":", 1)[1])
                             for i in range(3 if name == "besselpoly" else 1))
                self.check(name, *args)

    def test_additional_domains(self):
        self.check("itj0y0", 0.)
        self.check("iti0k0", 0.)
        for x in (0.02, 0.1, 0.5, 1., 3., 7.):
            for name in NAMES[1:]:
                self.check(name, x, rtol=2e-9, atol=4e-12)
        for args in ((0., .5, 2.), (0., 0., 0.), (1., .5, 2.),
                     (3., 1., 1.), (2., 2., 3.)):
            self.check("besselpoly", *args, rtol=2e-9)
        for name, args in (
                ("it2i0k0", (0.,)), ("it2j0y0", (11.,)),
                ("besselpoly", (-1., .5, 2.)),
                ("besselpoly", (1., .5, -2.))):
            for result in self.native.evaluate(self.case(name, *args)).values():
                self.assertIn("error", result)
        self.assertEqual(self.native.evaluate(self.case("unknown", 1.)), {})
        self.assertEqual(self.native.evaluate(self.case("itj0y0", 1+2j)), {})


if __name__ == "__main__":
    unittest.main()
