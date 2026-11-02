"""Independent numerical regressions for adapted Boost.Math/GSL conventions."""
import math
import unittest

import mpmath as mp
from scipy import special

from sci_testing.native_compare import NativeComparisons


class NativeSpecialAccuracyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.native = NativeComparisons()

    def compare(self, name, arguments, expected, *, libraries=('Boost.Math', 'GSL'),
                rtol=2e-11, atol=2e-13):
        case = {'api': f'special.{name}',
                'args': {f'p{i}': arg for i, arg in enumerate(arguments)}}
        results = self.native.evaluate(case)
        self.assertEqual(set(results), set(libraries), case)
        if not isinstance(expected, (list, tuple)):
            expected = [expected]
        for provider, record in results.items():
            with self.subTest(function=name, args=arguments, provider=provider):
                self.assertNotIn('error', record, record)
                self.assertEqual(len(record['values']), len(expected))
                for observed, reference in zip(record['values'], expected):
                    self.assertTrue(math.isclose(observed, float(reference),
                                                  rel_tol=rtol, abs_tol=atol),
                                    f'{name} {provider}: {observed} != {reference}')

    def test_carlson_and_elliptic_characteristic_conventions(self):
        with mp.workdps(65):
            for m in (0.1, 0.6):
                for n in (-0.2, 0.35):
                    self.compare('ellippi', (n, m), mp.ellippi(n, m))
                    for phi in (-0.7, 0.2, 1.1):
                        expected = mp.ellippi(n, phi, m)
                        self.compare('ellippiinc', (phi, n, m), expected)
                        self.compare('ellippiinc_bulirsch', (phi, n, m), expected)
            self.compare('elliprc', (0.2, 1.3), mp.elliprc(.2, 1.3))
            self.compare('elliprf', (0., 0.2, 1.), mp.elliprf(0., .2, 1.))
            self.compare('elliprd', (0., 0.2, 1.), mp.elliprd(0., .2, 1.))
            self.compare('elliprj', (0., .2, 1., 1.7), mp.elliprj(0., .2, 1., 1.7))

    def test_ellipkm1_preserves_tiny_complement(self):
        for complement in (1e-200, 1e-50, 1e-10, .25):
            self.compare('ellipkm1', (complement,), special.ellipkm1(complement))

    def test_jacobi_amplitude_keeps_periods(self):
        m = 0.45
        k = special.ellipk(m)
        for argument in (-5. * k, -2.3 * k, .2 * k, 2.3 * k, 5. * k):
            expected = special.ellipj(argument, m)
            self.compare('ellipj', (argument, m), expected, atol=3e-12)

    def test_real_gegenbauer_laguerre_and_spherical_derivatives(self):
        # Boost's ordinary associated Laguerre takes integer alpha. Sci's
        # generalised polynomial allows noninteger alpha.
        self.compare('eval_genlaguerre', (3, .5, .2), special.eval_genlaguerre(3, .5, .2))
        for name, scipy_function in (('spherical_jn', special.spherical_jn),
                                     ('spherical_yn', special.spherical_yn)):
            for n in (0, 1, 3):
                for derivative in (False, True):
                    self.compare(name, (n, 1.75, int(derivative)),
                                 scipy_function(n, 1.75, derivative=derivative))

    def test_inverse_incomplete_functions_and_new_sci_03_edges(self):
        self.compare('betaincinv', (2.5, 3.5, .4), special.betaincinv(2.5, 3.5, .4))
        self.compare('betainccinv', (2.5, 3.5, .4), special.betainccinv(2.5, 3.5, .4))
        self.compare('gammaincinv', (2.5, .4), special.gammaincinv(2.5, .4))
        self.compare('gammainccinv', (2.5, .4), special.gammainccinv(2.5, .4))
        with mp.workdps(70):
            a, b = mp.mpf('5000.1'), mp.mpf('5000.2')
            reference = mp.exp(mp.loggamma(a) - mp.loggamma(b))
            self.compare('gamma_ratio', (float(a), float(b)),
                         mp.exp(mp.loggamma(mp.mpf(float(a))) -
                                mp.loggamma(mp.mpf(float(b)))), rtol=5e-11)
        native = self.native.evaluate({
            'api': 'special.integrals.scaled_exp1',
            'args': {'special03_case': True, 'x': 0.0}})
        self.assertEqual(native, {'GSL': {'values': [0.0]}})


if __name__ == '__main__':
    unittest.main()
