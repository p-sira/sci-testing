"""Source provenance must include the organized Python modules and tests."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from sci_testing import report


class SourceDigestTests(unittest.TestCase):
    def test_package_and_test_changes_invalidate_source_digest(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in (
                'Cargo.toml', 'Cargo.lock', 'pyproject.toml', 'uv.lock',
                'sci-revision.txt', 'api_inventory.csv', 'README.md',
                'sci_testing/report.py', 'tests/test_report.py',
            ):
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text('original')
            with patch.object(report, 'ROOT', root):
                original = report.digest_sources()
                for name in ('sci_testing/report.py', 'tests/test_report.py'):
                    with self.subTest(path=name):
                        path = root / name
                        path.write_text('changed')
                        self.assertNotEqual(report.digest_sources(), original)
                        path.write_text('original')
                cache = root / 'sci_testing/__pycache__/report.pyc'
                cache.parent.mkdir()
                cache.write_bytes(b'cache')
                self.assertEqual(report.digest_sources(), original)
