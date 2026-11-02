"""Public GSL complex special comparisons through the checked native ABI."""
import ctypes as ct
import math
import unittest

import mpmath as mp

from sci_testing.native_compare import NativeComparisons


class NativeComplexTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.native = NativeComparisons()

    def case(self, name, real, imag):
        return {'api': f'special.complex.{name}',
                'args': {'special_complex03_case': True, 'z': [real, imag]}}

    def test_matches_independent_complex_references(self):
        with mp.workdps(80):
            for real, imag in ((.25, .5), (-1.25, .25), (-1.25, -.25),
                               (10., 10.)):
                z = mp.mpc(real, imag)
                for name, reference in (('digamma', mp.digamma(z)),
                                        ('spence', mp.polylog(2, 1 - z))):
                    with self.subTest(name=name, z=(real, imag)):
                        observed = self.native.evaluate_complex(self.case(name, real, imag))
                        self.assertEqual(set(observed), {'GSL'})
                        self.assertNotIn('error', observed['GSL'])
                        actual = observed['GSL']['values']
                        self.assertEqual(len(actual), 2)
                        for x, y in zip(actual, (float(mp.re(reference)),
                                                 float(mp.im(reference)))):
                            self.assertTrue(math.isclose(x, y, rel_tol=2e-11,
                                                         abs_tol=2e-12),
                                            f'{name}({z}): {actual} vs {reference}')

    def test_spence_negative_real_cut_exposes_gsl_signed_zero_convention(self):
        above = self.native.evaluate_complex(self.case('spence', -1.25, 0.))['GSL']['values']
        below = self.native.evaluate_complex(self.case('spence', -1.25, -0.))['GSL']['values']
        self.assertEqual(above, below)
        self.assertAlmostEqual(above[1], -math.pi * math.log(2.25), places=12)
        # The mathematical lower-bank limit has the opposite imaginary part.
        self.assertAlmostEqual(-below[1], math.pi * math.log(2.25), places=12)
        lower_limit = self.native.evaluate_complex(self.case('spence', -1.25, -1e-12))
        self.assertGreater(lower_limit['GSL']['values'][1], 0.)

    def test_nonfinite_input_and_gsl_domain_error_are_reported(self):
        nonfinite = self.native.evaluate_complex(self.case('digamma', math.inf, 1.))
        self.assertIn('nonfinite complex input', nonfinite['GSL']['error'])
        pole = self.native.evaluate_complex(self.case('digamma', 0., 0.))
        self.assertIn('error', pole['GSL'])
        self.assertNotIn('values', pole['GSL'])

    def test_unsupported_or_malformed_input(self):
        # Gamma is now handled by the dedicated checked GSL complex bridge.
        gamma = self.native.evaluate_complex(self.case('gamma', 1., 1.))
        self.assertEqual(set(gamma), {'GSL'})
        self.assertNotIn('error', gamma['GSL'])
        self.assertEqual(len(gamma['GSL']['values']), 2)
        self.assertEqual(self.native.evaluate_complex(self.case('wofz', 1., 1.)), {})
        self.assertEqual(self.native.evaluate_complex({'api': 'special.digamma',
                                                       'args': {'z': [.5, .5]}}), {})
        with self.assertRaisesRegex(ValueError, 'expected a complex'):
            self.native.evaluate_complex({'api': 'special.complex.digamma',
                                          'args': {'z': [.5]}})

    def test_abi_rejects_short_output_without_overwriting_buffer(self):
        self.assertEqual(self.native.library.comparison_complex_supported(b'digamma'), 2)
        self.assertEqual(self.native.library.comparison_complex_supported(b'spence'), 2)
        self.assertEqual(self.native.library.comparison_complex_supported(b'gamma'), 0)
        inputs = (ct.c_double * 2)(.25, .5)
        output = (ct.c_double * 2)(42., 43.)
        status = self.native.library.comparison_complex_evaluate(
            1, b'digamma', inputs, 2, output, 1)
        self.assertEqual(status, -1)
        self.assertEqual(list(output), [42., 43.])
        self.assertIn(b'output capacity', self.native.library.comparison_error())
        status = self.native.library.comparison_complex_evaluate(
            0, b'digamma', inputs, 2, output, 2)
        self.assertEqual(status, -1)
        self.assertEqual(list(output), [42., 43.])


if __name__ == '__main__':
    unittest.main()
