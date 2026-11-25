"""WolframEngine reference generation: signatures, precision and failure modes."""

# These regression tests deliberately inspect private parser and codegen helpers.
# pyright: reportPrivateUsage=false

import json
import math
import shutil
import subprocess
import unittest
from decimal import Decimal, localcontext
from typing import cast
from unittest.mock import patch

from sci_testing import wolfram_reference as wolfram


def case(name: str, *values: object) -> dict[str, object]:
    return {'api': name, 'args': {f'p{i}': value for i, value in enumerate(values)}}


class ReferenceValidationTests(unittest.TestCase):
    def test_curated_api_signatures(self):
        for name, args in [
            ('special.erf', {'x': .125}),
            ('special.gamma', {'p0': .125}),
            ('special.ellippiinc', {'p0': .7, 'p1': -.25, 'p2': -.5}),
            ('special.elliprj', {'p0': 0, 'p1': 1, 'p2': 2, 'p3': 3}),
            ('special.elliprg', {'p0': 0, 'p1': 1, 'p2': 2}),
            ('stats.normal_cdf', {'x': 0}),
            ('stats.variance', {'values': [1, 2, 3]}),
            ('integrate.fixed_quad', {'a': 0, 'b': 1}),
        ]:
            with self.subTest(name=name):
                self.assertTrue(wolfram.supports({'api': name, 'args': args}))
        for bad in [
            case('special.unreviewed', .5),
            case('special.erf', .5, .1),
            {'api': 'special.erf', 'args': {'p1': .1}},
            case('special.gamma', -0.0),
            case('special.gamma', '-0'),
            case('special.erf', float('nan')),
            case('special.gamma', 'Erf[0];Quit[]'),
            {'api': 'special.ellippiinc', 'args': {'phi': .5, 'n': .2, 'm': .5}},
            {'api': 'stats.variance', 'args': {'values': [1]}},
            {'api': 'special.complex.exp1', 'args': {'z': [1, 0]}},
        ]:
            with self.subTest(bad=bad):
                self.assertFalse(wolfram.supports(bad))
                with self.assertRaises(wolfram.UnsupportedReferenceError):
                    _ = wolfram.reference(bad)

    def test_definition_provenance_is_a_static_formula(self):
        self.assertEqual(wolfram.definition(case('special.ellippiinc', .5, -.2, -.3)),
                         'EllipticPi[{1},{0},{2}]')
        self.assertEqual(wolfram.definition(case('special.elliprf', .5, .25, .75)),
                         'CarlsonRF[{0},{1},{2}]')
        self.assertEqual(wolfram.definition({'api': 'stats.variance',
                                             'args': {'values': [1, 2, 3]}}),
                         'Variance[{values...}]')
        self.assertEqual(wolfram.definition({'api': 'integrate.fixed_quad',
                                             'args': {'a': 0, 'b': 1}}),
                         'Exp[{1}]-Exp[{0}]  (arguments: a, b)')

    def test_binary64_is_transmitted_as_exact_rational(self):
        value = .1
        numerator, denominator = value.as_integer_ratio()
        self.assertEqual(wolfram._expression(case('special.erf', value)),
                         f'Erf[Rational[{numerator},{denominator}]]')
        self.assertEqual(wolfram._expression(case('special.erf', '0.1')),
                         'Erf[Rational[1,10]]')
        self.assertEqual(wolfram._expression(case('special.gamma', 2**75)),
                         f'Gamma[Rational[{2**75},1]]')

    def test_wolfram_number_format(self):
        self.assertEqual(wolfram._number('1.125`160.'), '1.125')
        self.assertEqual(wolfram._number('-3.125`100.*^-450'), '-3.125e-450')
        self.assertEqual(wolfram._number('0'), '0')
        for bad in ['Indeterminate', 'ComplexInfinity', '1+I', '2/3',
                    '1.`MachinePrecision', 'NaN', '5`10*^foo']:
            with self.assertRaises(wolfram.WolframEvaluationError):
                _ = wolfram._number(bad)
        result = wolfram._normalize(
            ['OK', '1.23456789012345678987654321`20.', '0', '20'],
            digits=20, label='precision-test')
        self.assertEqual(result, ['1.2345678901234567899'])

    def test_imaginary_absolute_accuracy_zero_does_not_weaken_precision(self):
        # Im[N[Hypergeometric1F1[...], 180]] can produce 0``160.01:
        # zero absolute error bounded at 160 digits, NOT 160-digit precision.
        row = ['OK', '1.25`160.', '0``160.0167016254448', '160.']
        self.assertEqual(wolfram._normalize(row, digits=160, label='test'), ['1.25'])
        # Requiring 157 relative digits for a small real value also requires
        # more absolute accuracy in the spurious imaginary component.
        for bad in [
            ['OK', '1.25`160.*^-100', '0``160.0167016254448', '160.'],
            ['OK', '1.25`160.', '0``20.', '160.'],
            ['OK', '1.25`160.', '1``160.0167016254448', '160.'],
            ['OK', '0``160.0167016254448', '0', '160.'],
            ['OK', '0', '0``160.0167016254448', '160.'],
            ['OK', '1.25`30.', '0``160.0167016254448', '30.'],
        ]:
            with self.subTest(row=bad), self.assertRaises(wolfram.WolframEvaluationError):
                _ = wolfram._normalize(bad, digits=160, label='test')
        with self.assertRaises(wolfram.WolframEvaluationError):
            _ = wolfram._number('0``160.0167016254448')

    def test_batch_parser_and_no_shell(self):
        output = ('banner\n' + wolfram._MARKER +
                  json.dumps([['OK', '0.5`60.', '0', '60'],
                              ['OK', '2`60.', '0', '60']]) + '\nNull\n')
        fake = subprocess.CompletedProcess(args=[], returncode=0, stdout=output, stderr='')
        with patch.object(wolfram, '_executable', return_value='/usr/bin/wolframscript'), \
             patch.object(subprocess, 'run', return_value=fake) as runner:
            result = wolfram.references([case('special.erf', .5),
                                         case('special.gamma', 2)], digits=60)
        self.assertEqual(result, [['0.5'], ['2']])
        positional, keyword = cast(tuple[tuple[list[str]], dict[str, object]],
                                   runner.call_args)
        self.assertEqual(positional[0][0], '/usr/bin/wolframscript')
        self.assertEqual(positional[0][1], '-code')
        self.assertIn('Rational[1,2]', positional[0][2])
        self.assertNotIn('shell', keyword)
        self.assertEqual(keyword['timeout'], 90.0)

    def test_failure_is_not_a_fabricated_number(self):
        for value in [
            ['ERROR', '$Failed', '0', '0'],
            ['PRECISION', '16', '0', '0'],
            ['OK', '0.5`60.', '1.*^-100', '60'],
            ['OK', '2.0`16.', '0', '16'],
        ]:
            with self.subTest(value=value), self.assertRaises(wolfram.WolframEvaluationError):
                _ = wolfram._normalize(value, digits=60, label='test-case')
        self.assertEqual(wolfram.references([]), [])
        with self.assertRaises(ValueError):
            _ = wolfram.references([], digits=4)

    def test_batch_auto_splits_to_respect_argument_limit(self):
        calls: list[int] = []

        def fake_run(argv: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
            del kwargs
            count = argv[2].count('Module[{r=Quiet[')
            calls.append(count)
            data: list[list[str]] = [['OK', '1.25`60.', '0', '60']] * count
            return subprocess.CompletedProcess(
                args=argv, returncode=0,
                stdout=wolfram._MARKER + json.dumps(data) + '\nNull\n', stderr='')

        with patch.object(wolfram, '_MAX_ARGV', 800), \
             patch.object(wolfram, '_executable', return_value='/usr/bin/wolframscript'), \
             patch.object(subprocess, 'run', side_effect=fake_run):
            result = wolfram.references([case('special.erf', k) for k in range(9)],
                                        digits=60, batch_size=9)
        self.assertEqual(result, [['1.25']] * 9)
        self.assertGreater(len(calls), 1)
        self.assertEqual(sum(calls), 9)


@unittest.skipUnless(shutil.which('wolframscript'), 'WolframEngine is optional')
class WolframEngineIntegrationTests(unittest.TestCase):
    def test_hypergeometric_zero_imaginary_accuracy(self):
        value = wolfram.reference(case('special.hyp1f1', 1.5, 2.5, .4995), digits=160)
        self.assertEqual(len(value), 1)
        self.assertGreater(Decimal(value[0]), 1)

    def test_all_curated_special_names_have_numeric_wolfram_implementations(self):
        samples: list[dict[str, object]] = []
        for name, (_, count) in wolfram._SPECS.items():
            values = [.25] * count
            if name == 'lambertwm1':
                values = [-.2]
            elif name == 'polygamma':
                values = [2, 1.5]
            elif name == 'expn':
                values = [2, .5]
            elif name.startswith('eval_'):
                values[0] = 2
            samples.append(case('special.' + name, *values))
        results = wolfram.references(samples, digits=40)
        self.assertEqual(len(results), len(wolfram._SPECS))
        self.assertTrue(all(Decimal(row[0]).is_finite() for row in results))

    def test_engine_version_and_provenance(self):
        identity = wolfram.version()
        self.assertIn('Linux', identity)
        record = wolfram.provenance(digits=80)
        self.assertEqual(record['version'], identity)
        self.assertEqual(record['digits'], 80)
        self.assertIn('binary64', str(record['float_input']))

    def test_precision_binary64_and_real_recipes(self):
        samples = [
            case('special.erf', .1),
            case('special.gamma', .5),
            case('special.j0', 0),
            case('special.ellipk', 0),
            case('special.ellipe', 0),
            case('special.elliprf', 0, 1, 1),
            case('special.elliprd', 0, 1, 1),
            case('special.ellippiinc', .5, 0, 0),
            {'api': 'stats.normal_cdf', 'args': {'x': 0}},
            {'api': 'stats.variance', 'args': {'values': [1, 2, 3]}},
        ]
        results = wolfram.references(samples, digits=80)
        self.assertEqual(len(results), len(samples))
        self.assertTrue(all(len(row) == 1 for row in results))
        numbers = [Decimal(row[0]) for row in results]
        self.assertTrue(all(n.is_finite() for n in numbers))
        with localcontext() as ctx:
            ctx.prec = 60
            self.assertAlmostEqual(float(numbers[0]), math.erf(.1), places=15)
            self.assertAlmostEqual(float(numbers[1]), math.sqrt(math.pi), places=15)
            self.assertEqual(numbers[2], 1)
            for idx in (3, 4, 5):
                self.assertAlmostEqual(float(numbers[idx]), math.pi/2, places=14)
            self.assertAlmostEqual(float(numbers[6]), 3*math.pi/4, places=14)
            self.assertEqual(numbers[7], Decimal('0.5'))
            self.assertEqual(numbers[8], Decimal('0.5'))
            self.assertEqual(numbers[9], 1)
            # Exact binary64 .1 differs from exact decimal 0.1.  This result
            # also demonstrates that the final output carries many more than
            # binary64's ~16 significant digits.
            self.assertGreater(len(results[0][0]), 65)
            decimal_input = wolfram.reference(case('special.erf', '0.1'), digits=80)
            self.assertNotEqual(decimal_input, results[0])

    def test_elliptic_negative_parameter_and_branch_errors(self):
        case_ = case('special.ellippiinc', .5, -.4, -1.)
        output = wolfram.reference(case_, digits=65)
        self.assertGreater(Decimal(output[0]), Decimal('0.45'))
        self.assertLess(Decimal(output[0]), Decimal('0.5'))
        with self.assertRaises(wolfram.WolframEvaluationError):
            _ = wolfram.reference(case('special.ellipk', 2.), digits=60)
        with self.assertRaises(wolfram.WolframEvaluationError):
            _ = wolfram.reference(case('special.gamma', -1), digits=60)


if __name__ == '__main__':
    _ = unittest.main()
