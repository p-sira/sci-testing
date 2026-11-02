"""Independent high-precision references and scientific complex comparisons."""

import math
import unittest

import mpmath as mp

from sci_testing import report
from sci_testing.complex_accuracy import numpy_comparison, reference
from sci_testing.case_csv import case_paths, read_case_file
from sci_testing import ROOT


class ComplexAccuracyTests(unittest.TestCase):
    def test_cancellation_and_extreme_magnitudes_have_high_precision_oracles(self):
        examples = [
            ('multiply', {'p0': [1.2e-162, 1.2e-162],
                          'p1': [1.0000000000000002, 0.9999999999999998]}),
            ('powi', {'p0': [0., 1.], 'p1': 2147483647}),
            ('powf', {'p0': [1.7e308, 1.7e308], 'p1': -0.5}),
        ]
        with mp.workdps(160):
            for name, args in examples:
                case = {'api': f'complex.{name}', 'args': args}
                result = report.high_precision(case)
                self.assertEqual(result, reference(case))
                if name == 'multiply':
                    x, y = [mp.mpf(v) for v in args['p0']]
                    a, b = [mp.mpf(v) for v in args['p1']]
                    self.assertEqual(result, [x*a-y*b, x*b+y*a])
                    self.assertGreater(abs(result[0]), 0)
                elif name == 'powi':
                    self.assertEqual(result, [0, -1])
                    self.assertNotEqual(numpy_comparison(case)['NumPy']['values'][0], 0)
                else:
                    z = mp.mpc(*[mp.mpf(v) for v in args['p0']])
                    self.assertLess(abs(mp.mpc(*result) - 1/mp.sqrt(z)), mp.mpf('1e-168'))

    def test_principal_branch_distinguishes_negative_imaginary_zero(self):
        with mp.workdps(160):
            for name, args in [('ln', {'p0': [-1.25, -0.]}),
                               ('sqrt', {'x': -1., 'y': -0.}),
                               ('powf', {'p0': [-2., -0.], 'p1': .5})]:
                real, imag = reference({'api': f'complex.{name}', 'args': args})
                self.assertLess(imag, 0, name)
                if name == 'ln':self.assertEqual(imag, -mp.pi)
                if name == 'sqrt':self.assertLess(abs(real), mp.mpf('1e-350'))

    def test_inverse_function_branch_axes_match_the_correct_one_sided_limit(self):
        examples = [
            ('acos', [2., 0.], -1), ('acos', [2., -0.], 1),
            ('asin', [2., 0.], 1), ('asin', [2., -0.], -1),
            ('acosh', [0.5, -0.], -1),
            ('atanh', [2., 0.], 1), ('atanh', [2., -0.], -1),
        ]
        with mp.workdps(200):
            for name, value, side in examples:
                with self.subTest(name=name, value=value):
                    result=reference({'api':f'complex.{name}','args':{'p0':value}})
                    self.assertEqual(mp.sign(result[1]), side)
            for name in ('asinh', 'atan'):
                for real in (0., -0.):
                    for imaginary in (2., -2.):
                        observed=reference({'api':f'complex.{name}',
                                            'args':{'p0':[real,imaginary]}})
                        self.assertEqual(mp.sign(observed[0]),
                                         -1 if math.copysign(1.,real)<0 else 1)

    def test_subnormal_competitors_use_sufficient_oracle_precision(self):
        case={'id':'complex.atan/large','api':'complex.atan',
              'args':{'p0':[0.,1e308]}}
        self.assertEqual(report.reference_precisions(case)[0],(400,600))
        actual=report.stable_reference(case)
        self.assertGreater(actual[1],0)
        with mp.workdps(600):
            self.assertLess(abs(actual[1]-mp.mpf(1)/mp.mpf(1e308)),mp.mpf('1e-600'))

    def test_all_public_complex_numpy_competitors_return_correct_shape(self):
        for name, path in sorted(case_paths(ROOT/'cases/generated').items()):
            if not name.startswith('complex/'):
                continue
            case = read_case_file(path, name)[0]
            compared = report.competitors(case)
            with self.subTest(case=case['id']):
                self.assertEqual(set(compared), {'NumPy'})
                self.assertNotIn('error', compared['NumPy'])
                length = 1 if case['api'].split('.')[-1] in (
                    'abs', 'arg', 'imag', 'real', 'norm_sqr') else 2
                self.assertEqual(len(compared['NumPy']['values']), length)
                with mp.workdps(100):
                    self.assertEqual(len(report.high_precision(case)), length)

    def test_python_stdlib_never_supplies_comparison_results(self):
        for api, args in [
            ('special.erf', {'x': 0.5}),
            ('stats.normal_cdf', {'x': -1.0}),
            ('stats.variance', {'values': [1., 2., 3.]}),
        ]:
            with self.subTest(api=api):
                result = report.competitors({'api': api, 'args': args})
                self.assertNotIn('stdlib', result)
                self.assertIn('SciPy', result)


if __name__ == '__main__':
    unittest.main()
