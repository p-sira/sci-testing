"""Independent numerical I/O fixtures exercised through pinned public SciPy."""

from collections import Counter
import io
import json
import math
import os
import subprocess
import unittest
import warnings

from scipy import io as sio

from sci_testing import io_cases


class IoCasesTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.batch = list(io_cases.cases(1000, 20261020))
        cls.index = {case['id']: case for case in cls.batch}

    def find(self, api, name):
        return next(case for case in self.batch
                    if case['api'] == api and case['id'].split('/')[1] == name)

    def test_every_public_io_function_has_comparative_cases(self):
        counts = Counter(case['api'] for case in self.batch)
        self.assertEqual(set(counts), set(io_cases.APIS))
        self.assertEqual(len(self.index), len(self.batch))
        self.assertTrue(all(count >= 8 for count in counts.values()))
        for case in self.batch:
            self.assertTrue(case['args']['io_case'])
            self.assertTrue(bytes.fromhex(case['args']['payload_hex']))
            self.assertIn(case['args']['operation'], ('success', 'reject'))

    def test_selection_and_source_bytes_are_deterministic(self):
        self.assertEqual(self.batch, list(io_cases.cases(1000, 20261020)))
        one = list(io_cases.cases(2, 1, ['io.matlab.loadmat']))
        self.assertEqual(len([c for c in one if c['args']['operation'] == 'success']), 4)
        self.assertTrue(all(c['api'] == 'io.matlab.loadmat' for c in one))
        self.assertTrue(all('/double_' in c['id'] for c in one
                            if c['args']['operation'] == 'success'))

    def test_independent_scipy_references_produce_finite_numeric_outputs(self):
        for case in self.batch:
            with self.subTest(case=case['id']):
                output = io_cases.reference(case)
                self.assertGreater(len(output), 0)
                self.assertTrue(all(isinstance(v, (float, int)) and math.isfinite(v)
                                    for v in output))
                if case['args']['operation'] == 'reject':
                    self.assertEqual(output, [1.])

    def test_shapes_classes_and_boundary_samples(self):
        mminfo = io_cases.reference(self.find('io.matrix_market.mminfo',
                                              'symmetric-success'))
        self.assertEqual(mminfo, [3., 3., 9., 0., 0., 1.])
        mmread = io_cases.reference(self.find('io.matrix_market.mmread',
                                              'integer-success'))
        self.assertEqual(mmread[:3], [1., 2., 3.])
        self.assertIn(9007199254740992., mmread)
        pcm24 = io_cases.reference(self.find('io.wavfile.read', 'i24-success'))
        self.assertEqual(pcm24[:5], [48000., 2., 24., 5., 6.])
        self.assertIn(-2147483648., pcm24)
        self.assertEqual(io_cases.reference(self.find('io.wavfile.signal', 'u8-0-success')),
                         [3., 0., -1., 1., 0., 0., 1., 0., 0.9921875, 1.])
        self.assertEqual(io_cases.reference(self.find('io.matlab.matfile_version',
                                                     'double_plain-success')), [1., 0.])
        metadata = io_cases.reference(self.find('io.matlab.whosmat',
                                                'double_compressed-success'))
        self.assertEqual(metadata, [1., 1., 120., 2., 2., 3., 1., 0.])
        # 64-bit integer values are split into two exact 32-bit words; a
        # binary64 cast alone could not distinguish 2**53 from 2**53+1.
        large = io_cases.reference(self.find('io.matlab.loadmat',
                                             'uint64_plain-success'))
        self.assertEqual(large[:9], [1., 1., 120., 2., 1., 3., 0., 8., 3.])
        self.assertIn(4294967295., large)

    @unittest.skipUnless(os.environ.get('SCI_TESTING_IO_RUNNER'),
                         'set SCI_TESTING_IO_RUNNER to the compiled Rust test runner')
    def test_actual_rust_writer_bytes_decode_with_public_scipy(self):
        """Read *Rust-produced* bytes with SciPy; detect correlated Rust codec bugs.

        The ordinary case runner emits numerical values only. A focused byte
        probe exposes the three writer streams to independent SciPy consumers;
        its payloads are not counted as additional numerical comparison cases.
        """
        batch = [case for case in self.batch if case['api'] in (
            'io.matrix_market.mmwrite', 'io.wavfile.write', 'io.matlab.savemat')]
        requests = [{**case, 'args': {**case['args'], 'operation': 'emit_bytes'}}
                    for case in batch]
        run = subprocess.run([os.environ['SCI_TESTING_IO_RUNNER']],
                             input=''.join(json.dumps(case, separators=(',', ':'))+'\n'
                                           for case in requests),
                             text=True, capture_output=True, check=True)
        answers = [json.loads(line) for line in run.stdout.splitlines()]
        self.assertEqual(len(answers), len(batch))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            for case, answer in zip(batch, answers):
                with self.subTest(case=case['id']):
                    self.assertEqual(answer['id'], case['id'])
                    self.assertNotIn('error', answer)
                    byte_values = answer['values']
                    self.assertTrue(all(int(v) == v and 0 <= v <= 255
                                        for v in byte_values))
                    data = bytes(int(value) for value in byte_values)
                    self.assertTrue(data)
                    if case['api'] == 'io.matrix_market.mmwrite':
                        loaded = sio.mmread(io.BytesIO(data), spmatrix=False)
                        observed = io_cases._market_matrix(loaded)
                    elif case['api'] == 'io.wavfile.write':
                        observed = io_cases._wav(data)
                    else:
                        observed = io_cases._mat_values(data)
                    expected = io_cases.reference(case)
                    self.assertEqual(len(observed), len(expected))
                    for i, (actual, target) in enumerate(zip(observed, expected)):
                        self.assertTrue(math.isclose(actual, target, rel_tol=1e-12,
                                                     abs_tol=1e-12),
                                        f'{case["id"]}[{i}]: {actual} != {target}')


if __name__ == '__main__':
    unittest.main()
