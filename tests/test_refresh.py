"""Sci refreshes retain the saved comparison snapshot and support namespaces."""
from contextlib import ExitStack, redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from sci_testing import report
from sci_testing.case_csv import (directory_digest, measurement_files,
                                 read_measurements, write_cases, write_measurements)


class Cases(list):
    @property
    def apis(self):
        return sorted({case['api'] for case in self})

    def function(self, api):
        return [case for case in self if case['api'] == api]

    def module_counts(self, module):
        apis=[api for api in self.apis if api.startswith(module+'.')]
        return len(apis),sum(len(self.function(api)) for api in apis)


class RefreshTests(unittest.TestCase):
    def setUp(self):
        self.directory=tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root=Path(self.directory.name)
        self.output=self.root/'pages'
        self.cases=Cases(dict(id=f'{api}/{group}-{i:04d}',api=api,group=group,args={'x':.5})
                         for api in ['sparse.linalg.norm','special.erf']
                         for group in ['linspace','random'] for i in range(1000))
        self.records={case['id']:dict(id=case['id'],reference=['0.50000000000000000000000000000000000000001'],
                      libraries={'Sci':{'error':'old failure'} if case['id'].endswith('-0000') else {'values':[.5]},
                                 'SciPy':{'values':[-0.], 'warnings':['first\nsecond']},
                                 'GSL':{'values':[.25],'error':'saved error'}})
                      for case in self.cases}
        write_cases(self.root/'cases/generated',self.cases)
        write_measurements(self.root/'results/measurements',measurement_files(self.cases,self.records),compress_above=0)
        self.meta=dict(sci_commit='old',sci_source_dirty=False,testing_sources_sha256='old harness',
                       cases_sha256=directory_digest(self.root/'cases/generated'),
                       measurements_sha256=directory_digest(self.root/'results/measurements'),
                       native_libraries={'GSL':'saved'},cxx='saved compiler',native_build_flags='saved flags',
                       module_libraries={module:['Sci','SciPy','GSL'] for module in ['sparse','special']},
                       numpy='2.5.4',scipy='1.18.1',filterpy='1.4.5',mpmath='1.3.0',
                       rust='rustc',python='3.13',platform='test',processor='test')
        (self.root/'results/environment.json').write_text(json.dumps(self.meta))
        (self.root/'sci-revision.txt').write_text('new\n')
        self.stack=ExitStack()
        self.addCleanup(self.stack.close)
        for name,value in [('ROOT',self.root),('MODULES',['sparse','special'])]:
            self.stack.enter_context(patch.object(report,name,value))
        self.stack.enter_context(patch.object(report,'build_cases',return_value=self.cases))
        self.stack.enter_context(patch.object(report,'command',return_value='new'))
        self.stack.enter_context(patch.object(report,'digest_sources',return_value='new harness'))
        self.stack.enter_context(patch.object(report,'environment',return_value={**self.meta,'sci_commit':'new','testing_sources_sha256':'new harness'}))
        self.stack.enter_context(patch.object(report.subprocess,'run'))
        self.measure=self.stack.enter_context(patch.object(report,'measure_sci',side_effect=lambda batch,executable:
                          [dict(id=case['id'],values=[.75]) for case in batch]))
        self.native=self.stack.enter_context(patch.object(report,'NativeComparisons',side_effect=AssertionError('native evaluation')))
        self.competitors=self.stack.enter_context(patch.object(report,'competitors',side_effect=AssertionError('competitor evaluation')))
        self.reference=self.stack.enter_context(patch.object(report,'high_precision',side_effect=AssertionError('reference evaluation')))
        self.stack.enter_context(patch('sci_testing.coverage_audit.inventory',return_value=[
            dict(api=api,canonical=api) for api in self.cases.apis]))

    def run_report(self, *arguments):
        with patch.object(report.sys,'argv',['report','--output',str(self.output),*arguments]), redirect_stdout(io.StringIO()):
            report.main()

    def assert_preserved(self, selected):
        actual=read_measurements(self.root/'results/measurements')
        for case in self.cases:
            before=self.records[case['id']]
            after=actual[case['id']]
            self.assertEqual(after['reference'],before['reference'])
            self.assertEqual(after['libraries']['SciPy'],before['libraries']['SciPy'])
            self.assertEqual(after['libraries']['GSL'],before['libraries']['GSL'])
            self.assertEqual(after['libraries']['Sci'],{'values':[.75]} if case['api'] in selected else before['libraries']['Sci'])
        self.native.assert_not_called()
        self.competitors.assert_not_called()
        self.reference.assert_not_called()
        return json.loads((self.root/'results/environment.json').read_text())

    def test_default_refreshes_only_sci_and_passes_check(self):
        self.run_report()
        meta=self.assert_preserved(set(self.cases.apis))
        self.assertEqual(meta['reference_environment'],self.meta)
        self.assertEqual(self.measure.call_count,2)
        self.run_report('--check')

    def test_full_refresh_recomputes_references_and_competitors(self):
        self.native.side_effect=None
        self.native.return_value.versions={'GSL':'fresh'}
        self.native.return_value.evaluate.return_value={}
        self.reference.side_effect=None
        self.reference.return_value=[.5]
        self.competitors.side_effect=None
        self.competitors.return_value={'SciPy':{'values':[.25]}}
        self.run_report('--full')
        actual=read_measurements(self.root/'results/measurements')
        for record in actual.values():
            self.assertEqual(record['reference'],['0.5'])
            self.assertEqual(record['libraries'],{'Sci':{'values':[.75]},'SciPy':{'values':[.25]}})
        self.assertEqual(self.reference.call_count,2*len(self.cases))
        self.assertEqual(self.competitors.call_count,len(self.cases))
        meta=json.loads((self.root/'results/environment.json').read_text())
        self.assertNotIn('reference_environment',meta)
        self.run_report('--check')

    def test_nested_module_preserves_unselected_file_and_revision(self):
        path=self.root/'results/measurements/special/erf.csv.gz'
        before=path.read_bytes()
        self.run_report('--module','sparse.linalg')
        meta=self.assert_preserved({'sparse.linalg.norm'})
        self.assertEqual(path.read_bytes(),before)
        self.assertEqual(meta['sci_measurement_commits'],{'sparse.linalg.norm':'new','special.erf':'old'})
        self.assertEqual(self.measure.call_count,1)
        self.run_report('--check')
        self.run_report('--module','special','--module','sparse')
        meta=self.assert_preserved(set(self.cases.apis))
        self.assertEqual(meta['reference_environment'],self.meta)
        self.assertEqual(set(meta['sci_measurement_commits'].values()),{'new'})

    def test_invalid_selection_is_rejected_before_measurement(self):
        for arguments in [('--module','unknown'),('--full','--module','special'),('--check','--module','special')]:
            with self.subTest(arguments=arguments), redirect_stdout(io.StringIO()), patch('sys.stderr',new_callable=io.StringIO):
                with self.assertRaises(SystemExit):self.run_report(*arguments)
        self.measure.assert_not_called()

    def test_inputs_and_measurements_must_match_before_refresh(self):
        for path in [self.root/'cases/generated/special/erf.csv',self.root/'results/measurements/special/erf.csv.gz']:
            before=path.read_bytes()
            path.write_bytes(b'changed')
            with self.assertRaises(ValueError):self.run_report()
            path.write_bytes(before)
        self.measure.assert_not_called()

    def test_failed_coverage_keeps_original_snapshot(self):
        before={p:p.read_bytes() for p in (self.root/'results').rglob('*') if p.is_file()}
        self.measure.side_effect=lambda batch,executable:[dict(id=c['id'],error='new failure') for c in batch]
        with self.assertRaisesRegex(ValueError,'without a successful comparison'):self.run_report()
        self.assertEqual(before,{p:p.read_bytes() for p in before})

    def test_refresh_applies_existing_shape_and_assignment_checks(self):
        cases=[dict(id='cluster.vq.vq/a',api='cluster.vq.vq',args={'obs':[[1.]]}),
               dict(id='special.erf/a',api='special.erf',args={'x':.5})]
        records={c['id']:dict(id=c['id'],reference=['0'],libraries={'Sci':{'values':[0.]}}) for c in cases}
        self.measure.side_effect=lambda batch,executable:[dict(id=c['id'],values=[1.] if c['api']=='cluster.vq.vq' else [1.,2.]) for c in batch]
        report.refresh_sci(cases,records,Path('runner'),self.root/'staged',{c['api'] for c in cases})
        actual=read_measurements(self.root/'staged')
        self.assertEqual(actual['cluster.vq.vq/a']['libraries']['Sci']['error'],'incorrect code assignment')
        self.assertIn('shape mismatch',actual['special.erf/a']['libraries']['Sci']['error'])
