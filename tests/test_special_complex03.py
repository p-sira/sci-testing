"""Independent tests for every Sci 0.3 complex-special public overload."""

from collections import Counter
import json
import math
import os
import subprocess
import unittest

import mpmath as mp
from scipy import special

from sci_testing import special_complex03_cases as complex_cases

# Recorded on Sci 0.3, SciPy 1.18.1, z=-1.5-0j, v=0. These eight Sci
# values are the mathematically correct lower-bank limits; public SciPy
# returns the upper-bank values in the second column for -0j. The comparison
# deliberately uses both independent mpmath limits and Rust results.
_LOWER_CUT_RECORDS = {
    'hankel1': (1.5354830152077545 + 0.38244892379775886j,
                -0.5118276717359181 + 0.38244892379775886j),
    'hankel1e': (-0.272875112416297 + 1.5586899763356903j,
                 -0.417696141355075 - 0.483492169911526j),
    'hankel2': (-0.5118276717359181 - 0.38244892379775886j,
                 1.5354830152077545 - 0.38244892379775886j),
    'hankel2e': (-0.417696141355075 + 0.483492169911526j,
                  -0.272875112416297 - 1.5586899763356903j),
    'kv': (0.21380556264752565 + 5.173333475486464j,
           0.21380556264752562 - 5.173333475486464j),
    'kve': (0.04770646943416754 + 1.1543267268865278j,
            0.04770646943416754 - 1.1543267268865278j),
    'yv': (0.38244892379775886 - 1.0236553434718363j,
           0.38244892379775886 + 1.0236553434718363j),
    'yve': (0.38244892379775886 - 1.0236553434718363j,
            0.38244892379775886 + 1.0236553434718363j),
}


def _rust_responses(batch):
    command = os.environ['SCI_TESTING_COMPLEX_RUNNER']
    run = subprocess.run([command],
                         input=''.join(json.dumps(c, allow_nan=False) + '\n' for c in batch),
                         text=True, capture_output=True, check=True)
    results = [json.loads(line) for line in run.stdout.splitlines()]
    if [r['id'] for r in results] != [c['id'] for c in batch]:
        raise AssertionError('Rust runner dropped, duplicated or reordered cases')
    return results


class ComplexSpecial03Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.samples = list(complex_cases.cases(4, 20261102))

    def test_ci_does_not_skip_compiled_complex_regressions(self):
        # The signed-zero Airy regression was previously skipped by CI even
        # though the workflow was green. Keep local oracle-only runs optional,
        # but require the actual compiled Sci runner in GitHub Actions.
        if os.environ.get('CI') == 'true':
            self.assertTrue(os.environ.get('SCI_TESTING_COMPLEX_RUNNER'),
                            'CI must set SCI_TESTING_COMPLEX_RUNNER to run '
                            'the Rust-backed Airy/branch-cut regressions')

    def find(self, name, edge_index):
        api = 'special.complex.' + name
        return next(case for case in self.samples
                    if case['id'] == f'{api}/edge-{edge_index:03d}')

    def test_complete_unique_namespace_and_reproducible_sweeps(self):
        counts = Counter(case['api'] for case in self.samples)
        self.assertEqual(len(complex_cases.APIS), 32)
        self.assertEqual(set(counts), set(complex_cases.APIS))
        self.assertEqual(len(self.samples), len({case['id'] for case in self.samples}))
        self.assertEqual(self.samples, list(complex_cases.cases(4, 20261102)))
        for case in self.samples:
            self.assertTrue(case['api'].startswith('special.complex.'))
            self.assertIs(case['args']['special_complex03_case'], True)
            self.assertGreaterEqual(counts[case['api']], 13)
            self.assertNotIn('p0', case['args'])
        self.assertNotEqual(
            list(complex_cases.cases(3, 1, ['special.complex.gamma'])),
            list(complex_cases.cases(3, 2, ['special.complex.gamma'])))

    def test_validation_and_selected_api_order(self):
        with self.assertRaises(ValueError):
            list(complex_cases.cases(0, 1))
        with self.assertRaises(ValueError):
            list(complex_cases.cases(1, 1, ['special.gamma']))
        with self.assertRaises(ValueError):
            complex_cases.reference({'api': 'special.gamma',
                                     'args': {'special_complex03_case': True}})
        sub = list(complex_cases.cases(2, 8, ['special.complex.zeta',
                                             'special.complex.airy']))
        self.assertEqual(list(dict.fromkeys(case['api'] for case in sub)),
                         ['special.complex.airy', 'special.complex.zeta'])

    def test_mpmath_reference_is_finite_and_has_exact_output_shape(self):
        for case in self.samples:
            with self.subTest(case=case['id']), mp.workdps(65):
                values = complex_cases.reference(case)
                name = case['api'].split('.')[-1]
                expected = 8 if name in ('airy', 'airye') else 4 if name == 'fresnel' else 2
                self.assertEqual(len(values), expected)
                self.assertTrue(all(isinstance(x, mp.mpf) and mp.isfinite(x) for x in values))

    def test_known_values_and_result_orders(self):
        with mp.workdps(80):
            airy = complex_cases.reference(self.find('airy', 1))
            scipy_airy = special.airy(1. + 0.j)
            for expected, observed in zip(scipy_airy, (complex(*airy[i:i + 2])
                                                       for i in range(0, 8, 2))):
                self.assertAlmostEqual(observed.real, expected.real, places=12)
                self.assertAlmostEqual(observed.imag, expected.imag, places=12)
            s, c = special.fresnel(1. + 0.j)
            fresnel = complex_cases.reference(self.find('fresnel', 1))
            self.assertAlmostEqual(float(fresnel[0]), s.real, places=12)
            self.assertAlmostEqual(float(fresnel[2]), c.real, places=12)
            self.assertEqual(complex_cases.reference(self.find('zeta', 0)),
                             [mp.mpf(-.5), mp.mpf(0)])
            self.assertTrue(mp.almosteq(complex_cases.reference(self.find('zeta', 1))[0],
                                        mp.pi**2 / 6))
            self.assertTrue(mp.almosteq(complex_cases.reference(self.find('lambertw', 1))[0],
                                        mp.lambertw(1).real))

    def test_scaled_airy_and_bessel_are_mathematical_scalings(self):
        with mp.workdps(80):
            z = [1.25, -.75]
            values = [complex_cases.reference({
                'api': 'special.complex.' + name,
                'args': {'special_complex03_case': True, 'z': z, **({'v': 0.5}
                                                                  if name in ('jve', 'ive', 'kve', 'hankel1e') else {})}
            }) for name in ('airy', 'airye', 'jve', 'ive', 'kve', 'hankel1e')]
            z_mp = mp.mpc(*z)
            phase = mp.mpf(2) / 3 * z_mp * mp.sqrt(z_mp)
            self.assertTrue(mp.almosteq(mp.mpc(*values[1][0:2]),
                                        mp.mpc(*values[0][0:2]) * mp.exp(phase)))
            self.assertTrue(mp.almosteq(mp.mpc(*values[1][4:6]),
                                        mp.mpc(*values[0][4:6]) * mp.exp(-abs(mp.re(phase)))))
            for name, scaled, scale in (
                ('jv', values[2], mp.exp(-abs(z_mp.imag))),
                ('iv', values[3], mp.exp(-abs(z_mp.real))),
                ('kv', values[4], mp.exp(z_mp)),
                ('hankel1', values[5], mp.exp(-mp.j*z_mp)),
            ):
                unscaled = complex_cases.reference({
                    'api': 'special.complex.' + name,
                    'args': {'special_complex03_case': True, 'z': z, 'v': 0.5}})
                self.assertTrue(mp.almosteq(mp.mpc(*scaled), mp.mpc(*unscaled)*scale), name)

    def test_independent_scipy_agrees_with_mpmath_at_interior_points(self):
        # Check formula selection, argument order, derivatives and exponential
        # scaling using two independent public numerical implementations.
        with mp.workdps(100):
            for api in complex_cases.APIS:
                name = api.rsplit('.', 1)[-1]
                args = {'special_complex03_case': True, 'z': [1.25, .75]}
                parameters = [complex(1.25, .75)]
                if name in complex_cases._BESSEL:
                    args['v'] = 0.5
                    parameters.insert(0, 0.5)
                elif name == 'lambertw':
                    args['k'] = 0
                    parameters.append(0)
                elif name in ('xlog1py', 'xlogy'):
                    args = {'special_complex03_case': True,
                            'x': [.8, .3], 'y': [1.25, .75]}
                    parameters = [complex(.8, .3), complex(1.25, .75)]
                observed = complex_cases.reference({'api': api, 'args': args})
                expected = getattr(special, name)(*parameters)
                expected = list(expected) if name in ('airy', 'airye', 'fresnel') else [expected]
                with self.subTest(api=api):
                    for i, value in enumerate(expected):
                        actual = complex(float(observed[2*i]), float(observed[2*i+1]))
                        self.assertLessEqual(abs(actual - value), 2e-12 * max(1., abs(value)))

    def test_precision_stability_at_100_and_160_digits(self):
        # References must be reproducible under the report's two mpmath dps.
        for api in complex_cases.APIS:
            case = next(case for case in self.samples if case['api'] == api
                        and '/random-' in case['id'])
            with mp.workdps(100):
                r100 = list(map(float, complex_cases.reference(case)))
            with mp.workdps(160):
                r160 = list(map(float, complex_cases.reference(case)))
            with self.subTest(api=api):
                self.assertEqual(r100, r160)

    def test_every_edge_is_stable_at_100_and_160_digits(self):
        # A previous cut-oracle used a dps-dependent imaginary displacement.
        # Integer-order J/I then had artificial 1e-50 imaginary components
        # at 100 dps and 1e-80 at 160 dps. Check all 313 edge oracles.
        for api in complex_cases.APIS:
            for case in complex_cases.cases(1, 20261102, [api]):
                if '/edge-' not in case['id']:
                    continue
                with mp.workdps(100):
                    r100 = list(map(float, complex_cases.reference(case)))
                with mp.workdps(160):
                    r160 = list(map(float, complex_cases.reference(case)))
                with self.subTest(case=case['id']):
                    self.assertEqual(r100, r160)

    def test_recorded_lower_cut_scipy_discrepancies(self):
        # These are observed differences in *public* SciPy, not grounds to
        # replace the independent reference with the upper-cut values.
        self.assertEqual(set(_LOWER_CUT_RECORDS),
                         complex_cases._SCIPY_IGNORES_LOWER_ZERO)
        for name, (lower_expected, scipy_expected) in _LOWER_CUT_RECORDS.items():
            args = {'special_complex03_case': True, 'v': 0., 'z': [-1.5, -0.]}
            case = {'api': 'special.complex.' + name, 'args': args}
            with self.subTest(api=name), mp.workdps(100):
                oracle = complex(*complex_cases.reference(case))
                scipy_value = complex(getattr(special, name)(0., complex(-1.5, -0.)))
                self.assertLess(abs(oracle-lower_expected), 1e-12)
                self.assertLess(abs(scipy_value-scipy_expected), 1e-12)
                self.assertGreater(abs(oracle-scipy_value), .25)

    def test_recorded_exact_axis_conventions_vs_analytic_limits(self):
        # Sci's Ei(x±i0) is the real principal value rather than the analytic
        # continuation Ei(x)±iπ. At z=-1.25 the actual Rust result was
        # (-0.14641337252591016, -0.0) on both banks.
        with mp.workdps(100):
            for side in (0., -0.):
                case = {'api': 'special.complex.expi',
                        'args': {'special_complex03_case': True,
                                 'z': [-1.25, side]}}
                exact = complex(*complex_cases.reference(case))
                limit = mp.ei(mp.mpc(-1.25, mp.mpf('1e-45') *
                                      (-1 if math.copysign(1., side) < 0 else 1)))
                with self.subTest(side=math.copysign(1., side)):
                    self.assertLess(abs(exact - complex(-.14641337252591016, 0)), 1e-14)
                    self.assertAlmostEqual(abs(complex(limit)-exact), math.pi, places=12)

            # Sci's W_-1(-0.3-0j) is the real W_-1(-0.3), while analytic
            # approach from below has a distinct value. Observed Rust:
            # -1.7813370234216275-0j; independent lower-bank limit:
            # -3.3002378364383755-7.436294411632747j.
            case = {'api': 'special.complex.lambertw',
                    'args': {'special_complex03_case': True,
                             'z': [-.3, -0.], 'k': -1}}
            exact = complex(*complex_cases.reference(case))
            limit = complex(mp.lambertw(mp.mpc(-.3, -mp.mpf('1e-45')), -1))
            self.assertLess(abs(exact - complex(-1.7813370234216275, 0)), 1e-14)
            self.assertLess(abs(limit - complex(-3.3002378364383755,
                                                -7.436294411632747)), 1e-12)
            self.assertGreater(abs(exact-limit), 7.)

    def test_airye_lower_cut_oracle_is_independent_of_scipy_axis_convention(self):
        # SciPy 1.18.1 still differs from the independent lower-bank limit for
        # the *exact* -0j input. This diagnostic must never alter the oracle
        # or excuse a mismatch in the current Rust implementation.
        api = 'special.complex.airye'
        lower_case = {'api': api, 'args': {
            'special_complex03_case': True, 'z': [-1.25, -0.]}}
        upper_case = {'api': api, 'args': {
            'special_complex03_case': True, 'z': [-1.25, 0.]}}
        with mp.workdps(100):
            lower_math = complex(*complex_cases.reference(lower_case)[:2])
            upper_math = complex(*complex_cases.reference(upper_case)[:2])
            scipy_exact_lower = complex(special.airye(complex(-1.25, -0.))[0])
            scipy_approach_lower = complex(special.airye(complex(-1.25, -1e-15))[0])
            self.assertLess(abs(lower_math - complex(0.3101938067761736,
                                                     0.41740519982220414)), 1e-12)
            self.assertLess(abs(upper_math - complex(0.3101938067761736,
                                                     -0.41740519982220414)), 1e-12)
            self.assertLess(abs(scipy_exact_lower - complex(-0.17350425503349404,
                                                           -0.19502321479086732)), 1e-12)
            self.assertLess(abs(scipy_approach_lower - lower_math), 1e-12)
            self.assertGreater(abs(scipy_exact_lower-lower_math), .5)

    def test_lower_cut_bessel_oracle_uses_mathematical_continuation(self):
        # SciPy 1.18.1 treats -0j like +0j for these AMOS overloads.
        # On the lower side mpmath's limiting continuation is different.
        for name in complex_cases._SCIPY_IGNORES_LOWER_ZERO:
            case = {'api': 'special.complex.' + name, 'args': {
                'special_complex03_case': True, 'v': 0., 'z': [-1.5, -0.]}}
            with self.subTest(api=name), mp.workdps(100):
                ref = complex(*complex_cases.reference(case))
                if name.startswith('hankel1'):
                    raw = mp.hankel1(0, mp.mpc(-1.5, -mp.mpf('1e-35')))
                elif name.startswith('hankel2'):
                    raw = mp.hankel2(0, mp.mpc(-1.5, -mp.mpf('1e-35')))
                elif name.startswith('kv'):
                    raw = mp.besselk(0, mp.mpc(-1.5, -mp.mpf('1e-35')))
                else:
                    raw = mp.bessely(0, mp.mpc(-1.5, -mp.mpf('1e-35')))
                z = mp.mpc(-1.5)
                factor = (mp.exp(-mp.j*z) if name == 'hankel1e' else
                          mp.exp(mp.j*z) if name == 'hankel2e' else
                          mp.exp(z) if name == 'kve' else 1)
                lower = complex(raw * factor)
                upper = complex(getattr(special, name)(0., complex(-1.5, +0.)))
                self.assertLess(abs(ref-lower), 1e-12)
                self.assertGreater(abs(ref-upper), .25)

    def test_cut_bank_sign_and_json_transport(self):
        for name in ('exp1', 'spence', 'lambertw', 'hankel1', 'yv', 'xlogy'):
            api = 'special.complex.' + name
            bank_key = 'y' if name == 'xlogy' else 'z'
            cases = [case for case in self.samples if case['api'] == api and
                     case['args'][bank_key][0] < 0 and
                     case['args'][bank_key][1] == 0.]
            self.assertGreaterEqual(len(cases), 2)
            self.assertEqual({math.copysign(1., case['args'][bank_key][1]) for case in cases},
                             {-1., 1.})
            for case in cases:
                with self.subTest(case=case['id']):
                    sent = json.loads(json.dumps(case, allow_nan=False))
                    self.assertEqual(math.copysign(1., sent['args'][bank_key][1]),
                                     math.copysign(1., case['args'][bank_key][1]))
                    self.assertTrue(all(mp.isfinite(x) for x in complex_cases.reference(sent)))

    @unittest.skipUnless(os.environ.get('SCI_TESTING_COMPLEX_RUNNER'),
                         'set SCI_TESTING_COMPLEX_RUNNER to the compiled Rust runner')
    def test_public_rust_overloads_produce_correct_shape_and_finite_results(self):
        # One moderate interior point per API; include arrays, tuples, order and
        # branch parameters in a real JSON->Rust->JSON transport round-trip.
        batch = [next(case for case in self.samples if case['api'] == api and
                      '/random-' in case['id']) for api in complex_cases.APIS]
        responses = _rust_responses(batch)
        self.assertEqual([row['id'] for row in responses], [case['id'] for case in batch])
        for case, result in zip(batch, responses):
            with self.subTest(case=case['id']):
                self.assertNotIn('error', result)
                reference = complex_cases.reference(case)
                self.assertEqual(len(result['values']), len(reference))
                self.assertTrue(all(math.isfinite(x) for x in result['values']))
                for i in range(0, len(reference), 2):
                    actual = complex(*result['values'][i:i+2])
                    expected = complex(float(reference[i]), float(reference[i+1]))
                    self.assertLessEqual(abs(actual-expected), 2e-10*max(1.,abs(expected)))

    @unittest.skipUnless(os.environ.get('SCI_TESTING_COMPLEX_RUNNER'),
                         'set SCI_TESTING_COMPLEX_RUNNER to the compiled Rust runner')
    def test_rust_cut_values_match_independent_oracle_and_record_conventions(self):
        batch = [{'id': f'special.complex.{name}/canonical-lower',
                  'api': f'special.complex.{name}',
                  'args': {'special_complex03_case': True, 'v': 0.,
                           'z': [-1.5, -0.]}}
                 for name in sorted(_LOWER_CUT_RECORDS)]
        batch.extend({'id': f'special.complex.expi/axis-{i}',
                      'api': 'special.complex.expi',
                      'args': {'special_complex03_case': True,
                               'z': [-1.25, sign]}}
                     for i, sign in enumerate((0., -0.)))
        batch.append({'id': 'special.complex.lambertw/axis--1',
                      'api': 'special.complex.lambertw',
                      'args': {'special_complex03_case': True,
                               'z': [-.3, -0.], 'k': -1}})
        for case, result in zip(batch, _rust_responses(batch)):
            with self.subTest(case=case['id']), mp.workdps(100):
                self.assertNotIn('error', result)
                oracle = complex(*complex_cases.reference(case))
                rust = complex(*result['values'])
                self.assertLessEqual(abs(rust-oracle), 2e-12*max(1.,abs(oracle)))

    @unittest.skipUnless(os.environ.get('SCI_TESTING_COMPLEX_RUNNER'),
                         'set SCI_TESTING_COMPLEX_RUNNER to the compiled Rust runner')
    def test_actual_rust_airye_matches_independent_lower_cut_reference(self):
        # The former Sci revision returned an incorrect lower-bank result for
        # -1.25-0j. This is now a positive regression: every component of
        # Ai, Ai', Bi and Bi' must agree with the independent 100-digit limit
        # on BOTH signed-zero banks. Preserve the sign through JSON transport.
        batch = [{'id': f'special.complex.airye/exact-{i}',
                  'api': 'special.complex.airye',
                  'args': {'special_complex03_case': True,
                           'z': [-1.25, side]}}
                 for i, side in enumerate((0., -0.))]
        for case, row in zip(batch, _rust_responses(batch)):
            with self.subTest(case=case['id']), mp.workdps(100):
                self.assertNotIn('error', row)
                mathematical = complex_cases.reference(case)
                self.assertEqual(len(row['values']), 8)
                self.assertEqual(len(mathematical), 8)
                for i, label in enumerate(('Ai', "Ai'", 'Bi', "Bi'")):
                    actual = complex(*row['values'][2*i:2*i+2])
                    expected = complex(*mathematical[2*i:2*i+2])
                    with self.subTest(component=label):
                        self.assertTrue(math.isfinite(actual.real) and
                                        math.isfinite(actual.imag))
                        self.assertLessEqual(abs(actual-expected),
                                             2e-12*max(1., abs(expected)))


if __name__ == '__main__':
    unittest.main()
