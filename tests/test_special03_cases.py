"""Independent references and deterministic coverage for Sci 0.3 special APIs."""

from collections import Counter
import math
import unittest

import mpmath as mp
from scipy import special as scipy_special

from sci_testing import special03_cases as special03


def case(api, **kwargs):
    return {'api': api, 'args': dict(special03_case=True, **kwargs)}


class Special03Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.batch = list(special03.cases(1000, 20261102))
        cls.by_api = {api: [c for c in cls.batch if c['api'] == api]
                      for api in special03.APIS}

    def test_coverage_determinism_and_selection(self):
        counts = Counter(c['api'] for c in self.batch)
        self.assertEqual(len(special03.APIS), 9)
        self.assertEqual(set(counts), set(special03.APIS))
        self.assertEqual(len({c['id'] for c in self.batch}), len(self.batch))
        for api, batch in self.by_api.items():
            self.assertGreaterEqual(counts[api], 2005)
            self.assertEqual(sum('/linspace-' in c['id'] for c in batch), 1000)
            self.assertEqual(sum('/random-' in c['id'] for c in batch), 1000)
            self.assertTrue(all(c['args']['special03_case'] for c in batch))
            self.assertEqual(batch, list(special03.cases(1000, 20261102, [api])))
        api = special03.APIS[0]
        self.assertNotEqual(list(special03.cases(5, 1, [api])),
                            list(special03.cases(5, 2, [api])))
        with self.assertRaises(ValueError):
            list(special03.cases(0, 1))
        with self.assertRaises(ValueError):
            list(special03.cases(1, 1, ['special.not_public']))

    def test_typed_domain_profiles(self):
        for c in self.batch:
            api, a = c['api'], c['args']
            if api.startswith('special.bessel.'):
                self.assertGreaterEqual(a['v'], .5)
                if api.endswith('iv_ratioinv'):
                    self.assertTrue(0 <= a['r'] < 1)
                else:
                    self.assertGreaterEqual(a['x'], 0)
            if api.endswith('normalized_gen_harmonic'):
                self.assertTrue(1 <= a['j'] <= a['k'] <= a['n'] <= 1_000_000)
                self.assertTrue(math.isfinite(a['a']))
            if api.endswith('hurwitz_zeta'):
                self.assertTrue(len(a['s']) == len(a['q']) == 2)
                self.assertTrue(0 <= a['s'][0] <= 30)
                self.assertTrue(0 < a['q'][0] <= 64)
                self.assertNotEqual(a['s'], [1, 0])
            if api.endswith('log1pmx'):
                self.assertGreater(a['x'], -1)
            if api.endswith('gamma_ratio'):
                self.assertTrue(a['a'] > 0 and a['b'] > 0)
            if api.endswith('scaled_exp1'):
                self.assertGreaterEqual(a['x'], 0)

    def test_mathematical_identities(self):
        with mp.workdps(100):
            ref = lambda api, **kw: special03.reference(case(api, **kw))
            for x in (1e-30, 0.1, 1., 30.):
                r = ref('special.bessel.iv_ratio', v=.5, x=x)[0]
                c = ref('special.bessel.iv_ratio_c', v=.5, x=x)[0]
                self.assertTrue(mp.almosteq(r, mp.tanh(mp.mpf(x))))
                self.assertTrue(mp.almosteq(r + c, 1))
            for v, r in ((1., .01), (1., .5), (2., .9), (8., .3), (128., .7)):
                x = ref('special.bessel.iv_ratioinv', v=v, r=r)[0]
                self.assertTrue(mp.almosteq(mp.besseli(v, x) / mp.besseli(v-1, x),
                                            mp.mpf(r), abs_eps=mp.mpf('1e-88')))
            ng = 'special.combinatorics.normalized_gen_harmonic'
            self.assertEqual(ref(ng, j=2, k=3, n=4, a=0.), [mp.mpf('.5')])
            self.assertEqual(ref(ng, j=1, k=4, n=4, a=-200.), [mp.mpf(1)])
            hz = 'special.complex.hurwitz_zeta'
            self.assertTrue(mp.almosteq(ref(hz, s=[2., 0.], q=[1., 0.])[0], mp.pi**2/6))
            self.assertEqual(ref(hz, s=[0., 0.], q=[1.5, .5]), [mp.mpf(-1), mp.mpf('-.5')])
            lm = 'special.elementary.log1pmx'
            small = ref(lm, x=1e-150)[0]
            self.assertTrue(small < 0)
            self.assertTrue(mp.almosteq(small, -mp.mpf(1e-150)**2/2))
            self.assertEqual(ref(lm, x=0.), [mp.mpf(0)])
            # The alternating series has a positive cubic term. Both signs
            # expose a missing power-series sign that tiny x cannot detect.
            for x in (-0.07242022334607315, 0.08008987436254666):
                self.assertTrue(mp.almosteq(ref(lm, x=x)[0],
                                            mp.log1p(mp.mpf(x)) - mp.mpf(x),
                                            rel_eps=mp.mpf('1e-90')))
            for y in (-40., -3., -0.5772156649015329, 0., 10., 40.):
                inv = ref('special.gamma.digammainv', y=y)[0]
                self.assertTrue(inv > 0)
                self.assertTrue(mp.almosteq(mp.digamma(inv), mp.mpf(y),
                                            abs_eps=mp.mpf('1e-85')))
            gr = 'special.gamma.gamma_ratio'
            self.assertTrue(mp.almosteq(ref(gr, a=201., b=200.)[0], mp.mpf(200),
                                        rel_eps=mp.mpf('1e-90')))
            self.assertTrue(mp.almosteq(ref(gr, a=.5, b=1.)[0], mp.sqrt(mp.pi)))
            se = 'special.integrals.scaled_exp1'
            self.assertEqual(ref(se, x=0.), [mp.mpf(0)])
            for x in (1e-10, .5, 2., 1000.):
                expected = mp.quad(lambda t: mp.exp(-t) * x/(x+t), [0, mp.inf])
                self.assertTrue(mp.almosteq(ref(se, x=x)[0], expected,
                                            rel_eps=mp.mpf('1e-85')))

    def test_selected_public_scipy_crosschecks(self):
        # SciPy provides the components, but not most of these composite
        # Sci exports. Calling independent public functions tests conventions.
        with mp.workdps(100):
            ref = lambda api, **kw: float(special03.reference(case(api, **kw))[0])
            for v, x in ((.5, .1), (1., 1.), (2., 8.), (8., 20.)):
                ratio = scipy_special.ive(v, x) / scipy_special.ive(v - 1, x)
                self.assertAlmostEqual(ref('special.bessel.iv_ratio', v=v, x=x),
                                       ratio, delta=2e-13)
                self.assertAlmostEqual(ref('special.bessel.iv_ratio_c', v=v, x=x),
                                       1 - ratio, delta=2e-13)
            h = 'special.combinatorics.normalized_gen_harmonic'
            expected = ((scipy_special.zeta(2., 2.) - scipy_special.zeta(2., 6.)) /
                        (scipy_special.zeta(2., 1.) - scipy_special.zeta(2., 11.)))
            self.assertAlmostEqual(ref(h, j=2, k=5, n=10, a=2.), expected, delta=2e-14)
            z = special03.reference(case('special.complex.hurwitz_zeta',
                                         s=[2., 0.], q=[2., 0.]))
            self.assertAlmostEqual(float(z[0]), scipy_special.zeta(2., 2.), delta=1e-14)
            self.assertEqual(float(z[1]), 0.)
            self.assertAlmostEqual(ref('special.elementary.log1pmx', x=-.3),
                                   math.log1p(-.3) + .3, delta=2e-16)
            x = ref('special.gamma.digammainv', y=-2.22)
            self.assertAlmostEqual(scipy_special.digamma(x), -2.22, delta=3e-14)
            a, b = 20., 20.25
            self.assertAlmostEqual(ref('special.gamma.gamma_ratio', a=a, b=b),
                                   math.exp(scipy_special.gammaln(a) -
                                            scipy_special.gammaln(b)), delta=2e-13)
            for x in (1e-6, .5, 1., 10., 50.):
                self.assertAlmostEqual(ref('special.integrals.scaled_exp1', x=x),
                                       x * math.exp(x) * scipy_special.exp1(x), delta=2e-14)

    def test_100_and_160_digit_references_agree(self):
        for api, batch in self.by_api.items():
            for c in (batch[0], batch[2], batch[5], batch[-1]):
                with self.subTest(case=c['id']):
                    with mp.workdps(100):
                        low = special03.reference(c)
                    with mp.workdps(160):
                        high = special03.reference(c)
                        self.assertEqual(len(low), len(high))
                        self.assertTrue(all(mp.isfinite(x) for x in high))
                        scale = max([mp.mpf(2)**-1022] + [abs(x) for x in high])
                        self.assertTrue(all(abs(x-y) <= mp.mpf('1e-90') * scale
                                            for x,y in zip(low,high)), c['id'])


if __name__ == '__main__':
    unittest.main()
