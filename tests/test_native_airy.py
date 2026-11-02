"""Independent Boost.Math/GSL Airy comparisons and checked vector ABI."""
import ctypes as ct
import math
import unittest

import mpmath as mp

from sci_testing.native_compare import NativeComparisons


def airy_reference(x, *, scaled=False):
    with mp.workdps(100):
        value = mp.mpf(x)
        result = [mp.airyai(value), mp.airyai(value, 1),
                  mp.airybi(value), mp.airybi(value, 1)]
        if scaled:
            zeta = mp.mpf(2) / 3 * value * mp.sqrt(value)
            result = [result[0] * mp.exp(zeta), result[1] * mp.exp(zeta),
                      result[2] * mp.exp(-zeta), result[3] * mp.exp(-zeta)]
        return [float(item) for item in result]


class NativeAiryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.native = NativeComparisons()

    def check_values(self, name, x):
        actual = self.native.evaluate({'api': f'special.{name}', 'args': {'p0': x}})
        expected = airy_reference(x, scaled=name == 'airye')
        self.assertEqual(set(actual), {'Boost.Math', 'GSL'})
        for library, result in actual.items():
            with self.subTest(api=name, x=x, library=library):
                self.assertNotIn('error', result)
                self.assertEqual(len(result['values']), 4)
                for index, (observed, reference) in enumerate(zip(result['values'], expected)):
                    self.assertTrue(math.isclose(observed, reference, rel_tol=2e-12,
                                                  abs_tol=2e-14),
                                    f'{name}({x}) {library} output {index}: '
                                    f'{observed} != {reference}')

    def test_unscaled_airy_preserves_scipy_component_order(self):
        for x in (-8., -1., 0., .5, 2., 20.):
            self.check_values('airy', x)

    def test_scaled_airy_uses_unscaled_derivatives_times_exponential(self):
        for x in (0., .5, 2., 20., 100., 150.):
            self.check_values('airye', x)

    def test_negative_real_scaled_airy_is_not_a_gsl_equivalent(self):
        # GSL returns unscaled real Ai for x < 0; Sci's real airye does not.
        for result in self.native.evaluate({
            'api': 'special.airye', 'args': {'p0': -1.}
        }).values():
            self.assertIn('real airye requires x >= 0', result['error'])
            self.assertNotIn('values', result)

    def test_complex_airy_is_outside_the_real_native_abi(self):
        self.assertEqual(self.native.evaluate({
            'api': 'special.complex.airy', 'args': {'p0': [.5, .25],
                                                    'special_complex03_case': True}
        }), {})

    def test_scalar_native_api_retains_single_output(self):
        result = self.native.evaluate({'api': 'special.erf', 'args': {'p0': .5}})
        for entry in result.values():
            self.assertNotIn('error', entry)
            self.assertEqual(len(entry['values']), 1)
            self.assertAlmostEqual(entry['values'][0], math.erf(.5))

    def test_checked_native_capacity_prevents_short_array_writes(self):
        values = (ct.c_double * 1)(.5)
        output = (ct.c_double * 3)(10., 20., 30.)
        status = self.native.library.comparison_evaluate_many(1, b'airy', values,
                                                               1, output, 3)
        self.assertEqual(status, -1)
        self.assertEqual(list(output), [10., 20., 30.])
        self.assertIn(b'output capacity', self.native.library.comparison_error())


if __name__ == '__main__':
    unittest.main()
