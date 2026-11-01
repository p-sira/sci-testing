"""Independent SciPy/NumPy coverage of newly exported Sci 0.3 structures."""

from collections import Counter
import json
import math
from pathlib import Path
import tempfile
import unittest
import warnings

import numpy as np
from scipy import signal

from sci_testing import structural03_cases as suite
from sci_testing.case_csv import read_cases, write_cases


class Structural03CasesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.batch = list(suite.cases(1000, 20261020))
        cls.by_id = {case['id']: case for case in cls.batch}

    def fixture(self, api, label):
        return self.by_id[f'{api}/{label}']

    def test_all_six_public_functions_have_distinct_sampled_sweeps(self):
        counts = Counter(case['api'] for case in self.batch)
        self.assertEqual(set(counts), set(suite.APIS))
        self.assertEqual(len(self.batch), len(self.by_id))
        for api in suite.APIS:
            sampled = [case for case in self.batch
                       if case['api'] == api and '/sampled-' in case['id']]
            self.assertEqual(len(sampled), 1000)
            self.assertEqual(len({json.dumps(case['args'], sort_keys=True)
                                  for case in sampled}), len(sampled))
            self.assertTrue(all(case['args']['structural03_case'] is True
                                for case in sampled))
            self.assertTrue(any(case['args']['operation'] == 'reject'
                                for case in self.batch if case['api'] == api))

    def test_seed_selection_and_order_independence(self):
        self.assertEqual(self.batch, list(suite.cases(1000, 20261020)))
        one_api = 'sparse.array.csr_indices_i64'
        subset = list(suite.cases(10, 20261020, apis=[one_api]))
        original = [c for c in self.batch if c['api'] == one_api]
        self.assertEqual(subset, original[:len(subset)])
        self.assertNotEqual(list(suite.cases(10, 1, [one_api]))[-1]['args'],
                            subset[-1]['args'])
        with self.assertRaisesRegex(ValueError, 'unknown structural03'):
            list(suite.cases(1, 1, ['sparse.array.unknown']))
        with self.assertRaisesRegex(ValueError, 'nonnegative'):
            list(suite.cases(-1, 1))

    def test_structural_storage_and_signal_cases_survive_csv_roundtrip(self):
        # These cases contain nested triples, empty matrices, explicit -0.0,
        # integer pointers and a boolean scaled flag. Verify the actual
        # generator CSV transport without touching historical fixtures.
        originals = [self.fixture(api, name) for api, name in [
            ('core.csr.from_coo', 'single-negative-zero'),
            ('core.csr.from_coo', 'bad-row'),
            ('core.csr.to_dense', 'empty-0x3'),
            ('core.csr.validate', 'unsorted-columns'),
            ('sparse.array.csr_indices_i32', 'duplicate-columns'),
            ('sparse.array.csr_indices_i64', 'subnormal'),
            ('signal.utilities.closest_stft_dual_window', 'nonconstant-desired'),
            ('signal.utilities.closest_stft_dual_window', 'zero-hop'),
        ]]
        with tempfile.TemporaryDirectory() as temporary:
            write_cases(Path(temporary), originals)
            restored = sorted(read_cases(Path(temporary)), key=lambda c: c['id'])
        self.assertEqual(sorted(originals, key=lambda c: c['id']), restored)
        for original, loaded in zip(sorted(originals, key=lambda c: c['id']), restored):
            self.assertEqual(suite.reference(original), suite.reference(loaded))

    def test_every_independent_oracle_output_is_flat_and_finite(self):
        for case in self.batch:
            with self.subTest(case=case['id']):
                reference = suite.reference(case)
                self.assertIsInstance(reference, list)
                self.assertGreater(len(reference), 0)
                self.assertTrue(all(type(v) is float and math.isfinite(v)
                                    for v in reference))
                if case['args']['operation'] == 'reject':
                    self.assertEqual(reference, [1.])

    def test_coo_duplicates_signed_zeros_and_canonical_order(self):
        prefix = 'core.csr.from_coo'
        cancel = suite.reference(self.fixture(prefix, 'duplicates-cancel'))
        # rows=2, cols=3, indptr=[0,1,2], indices=[0,2], data=[0,0]
        self.assertEqual(cancel, [2., 3., 3., 0., 1., 2.,
                                  2., 0., 2., 2., 0., 0., 0., 0.])
        negative_zero = suite.reference(self.fixture(prefix, 'single-negative-zero'))
        self.assertEqual(negative_zero[-2:], [0., 1.])
        subnormal = suite.reference(self.fixture(prefix, 'subnormal'))
        self.assertIn(2. ** -1074, subnormal)
        self.assertIn(-(2. ** -1074), subnormal)
        unsorted = suite.reference(self.fixture(prefix, 'unsorted-mixed'))
        self.assertEqual(unsorted[:8], [3., 4., 4., 0., 2., 3., 4., 4.])

    def test_dense_retains_explicit_negative_zero_and_correct_shape(self):
        prefix = 'core.csr.to_dense'
        negative_zero = suite.reference(self.fixture(prefix, 'single-negative-zero'))
        self.assertEqual(negative_zero, [1., 1., 0., 1.])
        zeros = suite.reference(self.fixture(prefix, 'empty-3x0'))
        self.assertEqual(zeros, [3., 0.])
        self.assertEqual(suite.reference(self.fixture(prefix, 'empty-0x3')), [1.])
        self.assertEqual(suite.reference(self.fixture(prefix, 'dense-work-limit')), [1.])

    def test_indices_conserve_nnz_rows_and_integer_exactness(self):
        for width in ('i32', 'i64'):
            api = f'sparse.array.csr_indices_{width}'
            for case in [c for c in self.batch if c['api'] == api
                         and c['args']['operation'] == 'success']:
                data = case['args']
                values = suite.reference(case)
                pointers = values[1:1 + data['rows'] + 1]
                indices = values[2 + data['rows'] + 1:]
                self.assertEqual(pointers, [float(x) for x in data['indptr']])
                self.assertEqual(indices, [float(x) for x in data['indices']])
                self.assertEqual(pointers[-1], len(data['data']))

    def test_stft_dual_constraint_including_scaled_target(self):
        api = 'signal.utilities.closest_stft_dual_window'
        for case in [c for c in self.batch if c['api'] == api
                     and c['args']['operation'] == 'success']:
            args = case['args']
            output = suite.reference(case)
            n, alpha, *dual = output
            self.assertEqual(int(n), len(args['win']))
            self.assertEqual(len(dual), len(args['desired']))
            if not args['scaled']:
                self.assertEqual(alpha, 1.)
            win = np.asarray(args['win'])
            dual = np.asarray(dual)
            for residue in range(args['hop']):
                self.assertAlmostEqual(float(np.dot(win[residue::args['hop']],
                                                     dual[residue::args['hop']])),
                                       1., delta=3e-12)

    def test_stft_errors_match_public_scipy(self):
        api = 'signal.utilities.closest_stft_dual_window'
        for case in [c for c in self.batch if c['api'] == api
                     and c['args']['operation'] == 'reject']:
            args = case['args']
            with self.subTest(case=case['id']):
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', RuntimeWarning)
                    with self.assertRaises(ValueError):
                        signal.closest_STFT_dual_window(
                            np.asarray(args['win']), args['hop'],
                            np.asarray(args['desired']), scaled=args['scaled'])


if __name__ == '__main__':
    unittest.main()
