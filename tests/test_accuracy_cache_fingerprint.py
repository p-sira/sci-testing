"""Numerical-cache invalidation must be independent of report prose."""

from pathlib import Path
import tempfile
import unittest

from sci_testing import ROOT
from sci_testing.accuracy_cache import metric_source_sha256


class MetricFingerprintTests(unittest.TestCase):
    def setUp(self):
        self.original = (ROOT / 'sci_testing/report.py').read_text()
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / 'report.py'

    def fingerprint(self, source):
        self.path.write_text(source)
        return metric_source_sha256(self.path)

    def test_prose_changes_preserve_numerical_cache(self):
        baseline = self.fingerprint(self.original)
        # Neither explanatory notes nor Markdown presentation affect the
        # numerical per-API rows; both still affect digest_sources().
        self.assertNotEqual(self.original.count('None of these adapters implies'), 0)
        self.assertEqual(baseline, self.fingerprint(self.original.replace(
            'None of these adapters implies', 'None of the listed adapters implies')))
        self.assertIn('Measured binary64 (`f64`) results', self.original)
        self.assertEqual(baseline, self.fingerprint(self.original.replace(
            'Measured binary64 (`f64`) results',
            'Measured binary64 (`f64`) numerical results')))

    def test_metric_logic_and_precision_changes_invalidate_cache(self):
        baseline = self.fingerprint(self.original)
        self.assertIn('reference_sum+=mp.fsum(errors)', self.original)
        self.assertNotEqual(baseline, self.fingerprint(self.original.replace(
            'reference_sum+=mp.fsum(errors)', 'reference_sum+=mp.fsum(errors)*2')))
        self.assertIn('METRIC_PRECISION_DPS = 160', self.original)
        self.assertNotEqual(baseline, self.fingerprint(self.original.replace(
            'METRIC_PRECISION_DPS = 160', 'METRIC_PRECISION_DPS = 80')))

    def test_missing_metric_function_fails_closed(self):
        with self.assertRaisesRegex(ValueError, 'Missing numerical metric definition'):
            self.fingerprint(self.original.replace('def summarize_api(', 'def renamed_api('))


if __name__ == '__main__':
    unittest.main()
