"""Source provenance must include the organized Python modules and tests."""
import gzip
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from sci_testing import report
from sci_testing.case_csv import (case_directory_digest, case_paths, directory_digest,
                                  compact_existing_cases, finalize_extended_cache,
                                  measurement_files, read_cases, stage_case_extension,
                                  write_cases, write_measurements)
from sci_testing.case_set import cached_cases, ExtendedCaseSet, historical_cases


class SourceDigestTests(unittest.TestCase):
    def test_package_and_test_changes_invalidate_source_digest(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in (
                '.gitmodules', 'Cargo.toml', 'Cargo.lock', 'pyproject.toml', 'uv.lock',
                'sci-revision.txt', 'api_inventory.csv', 'README.md',
                'sci_testing/report.py', 'tests/test_report.py',
            ):
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text('original')
            with patch.object(report, 'ROOT', root):
                original = report.digest_sources()
                for name in ('.gitmodules', 'sci_testing/report.py', 'tests/test_report.py'):
                    with self.subTest(path=name):
                        path = root / name
                        path.write_text('changed')
                        self.assertNotEqual(report.digest_sources(), original)
                        path.write_text('original')
                cache = root / 'sci_testing/__pycache__/report.pyc'
                cache.parent.mkdir()
                cache.write_bytes(b'cache')
                self.assertEqual(report.digest_sources(), original)


class CaseStorageTests(unittest.TestCase):
    def setUp(self):
        self.directory=tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root=Path(self.directory.name)
        self.cases=[
            dict(id='signal.spline/a',api='signal.spline',group='UTF-8 ± and comma, newline\nnext',
                 args={'x':[-0.0,1e-300,1e300],'p':{'provenance':'SciPy', 'precision':160}}),
            dict(id='signal.spline/b',api='signal.spline',group='reference recipe',
                 args={'x':[-0.0,1e-300,1e300],'p':{'provenance':'independent', 'precision':100}}),
            dict(id='stats.normal_cdf/a',api='stats.normal_cdf',group='edge',args={'x':-0.0}),
        ]

    def test_compressed_inputs_preserve_case_ids_values_and_original_digest(self):
        normal=self.root/'normal'
        compact=self.root/'compact'
        write_cases(normal,self.cases,compress_above=1000000)
        expected=directory_digest(normal)
        write_cases(compact,self.cases,compress_above=1,signature='recipe-1')
        self.assertTrue(all(path.suffix=='.gz' for path in case_paths(compact).values()))
        self.assertEqual(case_directory_digest(compact),expected)
        self.assertEqual(read_cases(compact),self.cases)
        self.assertEqual(math.copysign(1,read_cases(compact)[0]['args']['x'][0]),-1.)
        self.assertEqual(len(cached_cases(compact,'recipe-1')),len(self.cases))
        self.assertEqual(list(cached_cases(compact,'recipe-1')),self.cases)
        self.assertEqual(cached_cases(compact,'recipe-1').module_counts('signal'),(1,2))
        self.assertIsNone(cached_cases(compact,'changed-recipe'))
        archived=(compact/'signal/spline.csv.gz').read_bytes()
        self.assertEqual(gzip.decompress(archived),(normal/'signal/spline.csv').read_bytes())
        write_cases(compact,self.cases,compress_above=1,signature='recipe-1')
        self.assertEqual((compact/'signal/spline.csv.gz').read_bytes(),archived)

    def test_reject_modified_or_ambiguous_cached_inputs(self):
        write_cases(self.root,self.cases,compress_above=1,signature='recipe-2')
        path=self.root/'signal/spline.csv.gz'
        data=path.read_bytes()
        path.write_bytes(data[:-1]+bytes([data[-1]^1]))
        with self.assertRaisesRegex(ValueError,'Stale case CSV checksum'):
            cached_cases(self.root,'recipe-2')
        path.write_bytes(data)
        (self.root/'signal/spline.csv').write_bytes(gzip.decompress(data))
        with self.assertRaisesRegex(ValueError,'Duplicate case input files'):
            read_cases(self.root)

    def test_small_csv_is_kept_readable_and_stale_gzip_removed(self):
        write_cases(self.root,self.cases,compress_above=1)
        self.assertTrue((self.root/'signal/spline.csv.gz').exists())
        write_cases(self.root,self.cases)
        self.assertFalse((self.root/'signal/spline.csv.gz').exists())
        self.assertTrue((self.root/'signal/spline.csv').exists())
        self.assertEqual(read_cases(self.root),self.cases)
        self.assertFalse((self.root/'.cache.json').exists())

    def test_migrate_existing_input_snapshot_only_after_provenance_match(self):
        write_cases(self.root,self.cases,compress_above=1000000)
        original=directory_digest(self.root)
        with self.assertRaisesRegex(ValueError,'saved snapshot digest'):
            compact_existing_cases(self.root,signature='recipe',expected_digest='wrong',compress_above=1)
        self.assertFalse((self.root/'.cache.json').exists())
        self.assertTrue((self.root/'signal/spline.csv').exists())
        compact_existing_cases(self.root,signature='recipe',expected_digest=original,compress_above=1)
        self.assertEqual(case_directory_digest(self.root),original)
        self.assertEqual(list(cached_cases(self.root,'recipe')),self.cases)
        self.assertEqual(cached_cases(self.root,'recipe').module_counts('signal'),(1,2))


class HistoricalExtensionTests(unittest.TestCase):
    def setUp(self):
        directory=tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root=Path(directory.name)
        self.input_dir=self.root/'historical'
        self.measurements=self.root/'measurements'
        self.staged=self.root/'staged-inputs'
        self.staged_measurements=self.root/'staged-measurements'
        self.original=[
            dict(id='special.erf/a',api='special.erf',group='original',args={'x':-0.0}),
            dict(id='special.erf/b',api='special.erf',group='original',args={'x':1e-300}),
        ]
        self.new=[
            dict(id='io.wavfile.read/fixture-a',api='io.wavfile.read',group='pinned SciPy',
                 args={'io_case':True,'operation':'success','payload_hex':'0001'}),
            dict(id='special.new_function/linspace-0000',api='special.new_function',group='new',
                 args={'special03_case':True,'value':-0.0}),
        ]
        self.records={case['id']:dict(id=case['id'],reference=['0.0000000000000000000001'],
                       libraries={'Sci':{'values':[-0.0]},'SciPy':{'values':[0.0]}})
                      for case in [*self.original,*self.new]}
        write_cases(self.input_dir,self.original,signature='old-generator',compress_above=1)
        write_measurements(self.measurements,measurement_files(self.original,self.records),compress_above=0)
        self.old_digest=case_directory_digest(self.input_dir)
        self.measurement_digest=directory_digest(self.measurements)

    def historical(self):
        return historical_cases(self.input_dir,expected_cases_sha256=self.old_digest,
                                measurements_directory=self.measurements,
                                expected_measurements_sha256=self.measurement_digest)

    def test_extension_preserves_historical_input_bytes_and_waits_for_measurements(self):
        original={name:path.read_bytes() for name,path in case_paths(self.input_dir).items()}
        historic=self.historical()
        self.assertIsNone(cached_cases(self.input_dir,'new-generator'))
        view=ExtendedCaseSet(historic,self.new)
        self.assertEqual(len(view),4)
        self.assertEqual(view.module_counts('special'),(2,3))
        self.assertEqual(view.module_counts('io'),(1,1))
        self.assertEqual(view.function('special.erf'),self.original)
        self.assertEqual(list(view)[-1]['id'],'special.new_function/linspace-0000')
        counts=stage_case_extension(historic,self.new,self.staged,compress_above=1)
        self.assertEqual(counts,{'io.wavfile/read.csv':1,'special/new_function.csv':1})
        self.assertFalse((self.staged/'.cache.json').exists())
        for name,data in original.items():
            self.assertEqual(case_paths(self.staged)[name].read_bytes(),data)
        self.staged_measurements.mkdir()
        for path in self.measurements.rglob('*.csv.gz'):
            destination=self.staged_measurements/path.relative_to(self.measurements)
            destination.parent.mkdir(parents=True,exist_ok=True)
            destination.write_bytes(path.read_bytes())
        arguments=dict(signature='new-generator',historical_cases_sha256=self.old_digest,
                       historical_measurements_directory=self.measurements,
                       historical_measurements_sha256=self.measurement_digest,
                       measurements_directory=self.staged_measurements)
        with self.assertRaisesRegex(ValueError,'Extended measurements'):
            finalize_extended_cache(self.staged,historic,self.new,**arguments)
        self.assertFalse((self.staged/'.cache.json').exists())
        write_measurements(self.staged_measurements,measurement_files(self.new,self.records),compress_above=0)
        extended_digest=finalize_extended_cache(self.staged,historic,self.new,**arguments)
        self.assertNotEqual(extended_digest,self.old_digest)
        self.assertEqual(case_directory_digest(self.input_dir),self.old_digest)
        self.assertEqual(cached_cases(self.staged,'new-generator').provenance['generator_equivalence_verified'],False)
        self.assertEqual(list(cached_cases(self.staged,'new-generator')),list(view))
        self.assertEqual({name:path.read_bytes() for name,path in case_paths(self.input_dir).items()},original)

    def test_historical_integrity_and_edited_old_apis_are_rejected(self):
        with self.assertRaisesRegex(ValueError,'Historical case input digest'):
            historical_cases(self.input_dir,expected_cases_sha256='wrong',
                             measurements_directory=self.measurements,
                             expected_measurements_sha256=self.measurement_digest)
        with self.assertRaisesRegex(ValueError,'Historical measurement digest'):
            historical_cases(self.input_dir,expected_cases_sha256=self.old_digest,
                             measurements_directory=self.measurements,
                             expected_measurements_sha256='wrong')
        with self.assertRaisesRegex(ValueError,'historical function'):
            ExtendedCaseSet(self.historical(),[dict(id='special.erf/new',api='special.erf',args={})])
        path=next(iter(case_paths(self.input_dir).values()))
        data=path.read_bytes()
        path.write_bytes(data[:-1]+bytes([data[-1]^1]))
        with self.assertRaisesRegex(ValueError,'Stale case CSV checksum'):
            self.historical()

    def test_extension_rejects_changed_historical_stage_and_incomplete_new_measurements(self):
        historic=self.historical()
        stage_case_extension(historic,self.new,self.staged,compress_above=1)
        write_measurements(self.staged_measurements,measurement_files(self.original,self.records),compress_above=0)
        write_measurements(self.staged_measurements,measurement_files(self.new,self.records),compress_above=0)
        arguments=dict(signature='new',historical_cases_sha256=self.old_digest,
                       historical_measurements_directory=self.measurements,
                       historical_measurements_sha256=self.measurement_digest,
                       measurements_directory=self.staged_measurements)
        staged_original=case_paths(self.staged)['special/erf.csv']
        data=staged_original.read_bytes()
        staged_original.write_bytes(data+b'changed')
        with self.assertRaisesRegex(ValueError,'Changed historical input bytes'):
            finalize_extended_cache(self.staged,historic,self.new,**arguments)
        self.assertFalse((self.staged/'.cache.json').exists())
        staged_original.write_bytes(data)
        missing=self.staged_measurements/'io.wavfile/read.csv.gz'
        missing.unlink()
        with self.assertRaisesRegex(ValueError,'Extended measurements'):
            finalize_extended_cache(self.staged,historic,self.new,**arguments)
        self.assertFalse((self.staged/'.cache.json').exists())

    def test_overlay_avoids_generating_old_cases_after_signature_changes(self):
        from sci_testing import generate_cases
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            provider=root/'sci_testing/io_cases.py'
            provider.parent.mkdir()
            provider.write_text('# indicator for installed provider\n')
            with patch.object(generate_cases,'ROOT',root), \
                 patch('sci_testing.io_cases.cases',return_value=self.new), \
                 patch.object(generate_cases,'build_cases',side_effect=AssertionError('regenerated old inputs')):
                overlay=generate_cases.extend_historical_cases(self.historical())
            self.assertEqual(len(overlay),len(self.original)+len(self.new))
            self.assertEqual([case['id'] for case in overlay.function('special.erf')],
                             [case['id'] for case in self.original])
            self.assertEqual(set(overlay.added),{'io.wavfile.read','special.new_function'})
