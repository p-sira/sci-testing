"""Verify GSL 2.7 public fixed-Gauss-Legendre oracle across saved Sci cases."""

import ctypes as ct
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
from scipy import integrate

from sci_testing.case_csv import case_paths, read_case_file
from sci_testing.gsl_numerical_core import GSLNumericalCoreComparisons


class GSLNumericalCoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root = Path(__file__).resolve().parent.parent
        cls.tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tmp.cleanup)
        target = Path(cls.tmp.name) / 'gsl-numerical-core.so'
        env = os.environ.copy()
        env.pop('LD_LIBRARY_PATH', None)
        subprocess.run([os.environ.get('CXX', 'c++'), '-std=c++17', '-O2',
                        '-shared', '-fPIC', str(root / 'native/gsl_numerical_core.cpp'),
                        '-lgsl', '-lgslcblas', '-o', str(target)],
                       env=env, check=True)
        cls.adapter = GSLNumericalCoreComparisons(ct.CDLL(str(target)))
        paths = case_paths(root / 'cases/generated')
        cls.cases = read_case_file(paths['integrate/fixed_quad.csv'],
                                   'integrate/fixed_quad.csv')

    def test_all_saved_fixed_quadrature_cases(self):
        self.assertEqual(len(self.cases), 2068)
        for case in self.cases:
            self.assertEqual(set(case['args']), {'a', 'b'})
            a, b = case['args']['a'], case['args']['b']
            result = self.adapter.evaluate(case)
            self.assertEqual(set(result), {'GSL'})
            self.assertNotIn('error', result['GSL'], case)
            actual = result['GSL']['values']
            self.assertEqual(len(actual), 1)
            reference = integrate.fixed_quad(np.exp, a, b, n=16)[0]
            self.assertTrue(math.isclose(actual[0], reference, rel_tol=2e-11,
                                         abs_tol=1e-11), (case['id'],actual,reference))

    def test_signature_and_capacity_fail_closed(self):
        library = self.adapter.library
        self.assertEqual(library.gsl_numerical_core_supported(1,b'fixed_quad'),2)
        self.assertEqual(library.gsl_numerical_core_supported(0,b'fixed_quad'),0)
        self.assertEqual(self.adapter.evaluate({'api':'integrate.fixed_quad','args':{'a':1}}),{})
        values=(ct.c_double*2)(0,1)
        guard=(ct.c_double*1)(987654.)
        self.assertNotEqual(library.gsl_numerical_core_evaluate(1,b'fixed_quad',values,2,guard,0),0)
        self.assertEqual(guard[0],987654.)
        self.assertNotEqual(library.gsl_numerical_core_evaluate(0,b'fixed_quad',values,2,guard,1),0)
        self.assertEqual(guard[0],987654.)


if __name__ == '__main__':
    unittest.main()
