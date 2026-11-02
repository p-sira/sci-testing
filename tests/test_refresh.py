"""Sci refreshes retain the saved comparison snapshot and support namespaces."""
from contextlib import ExitStack, redirect_stdout
import io
import json
import os
from pathlib import Path
import shutil
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
        from sci_testing.generate_cases import case_signature
        write_cases(self.root/'cases/generated',self.cases,signature=case_signature())
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
        for name,value in [('ROOT',self.root),('SCI',self.root/'sci'),('MODULES',['sparse','special'])]:
            self.stack.enter_context(patch.object(report,name,value))
        self.stack.enter_context(patch.object(report,'build_cases',return_value=self.cases))
        self.stack.enter_context(patch.object(report,'command',return_value='new'))
        self.stack.enter_context(patch.object(report,'digest_sources',return_value='new harness'))
        self.stack.enter_context(patch.object(report,'environment',return_value={**self.meta,'sci_commit':'new','testing_sources_sha256':'new harness'}))
        self.stack.enter_context(patch.object(report.subprocess,'run'))
        self.stack.enter_context(patch.object(report.subprocess,'check_output',
            return_value=json.dumps({'target_directory':str(self.root/'target')})))
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

    def test_render_rebuilds_tables_without_running_sci_or_references(self):
        before=read_measurements(self.root/'results/measurements')
        stored=json.loads((self.root/'results/environment.json').read_text())
        self.run_report('--render')
        self.assertEqual(before,read_measurements(self.root/'results/measurements'))
        self.measure.assert_not_called()
        self.native.assert_not_called()
        self.competitors.assert_not_called()
        self.reference.assert_not_called()
        meta=json.loads((self.root/'results/environment.json').read_text())
        self.assertEqual(meta['sci_commit'],stored['sci_commit'])
        self.assertEqual(meta['reference_environment'],stored)
        self.assertEqual(meta['testing_sources_sha256'],'new harness')
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

    def test_exact_api_refresh_preserves_other_api_measurements(self):
        untouched=self.root/'results/measurements/sparse.linalg/norm.csv.gz'
        before=untouched.read_bytes()
        extension={'apis':['special.erf'],'sci_commit':'new-api-original-reference'}
        self.meta['reference_extensions']=[extension]
        (self.root/'results/environment.json').write_text(json.dumps(self.meta))
        self.run_report('--module','special.erf')
        meta=self.assert_preserved({'special.erf'})
        self.assertEqual(untouched.read_bytes(),before)
        self.assertEqual(meta['reference_extensions'],[extension])
        self.assertEqual(meta['sci_measurement_commits'],
                         {'special.erf':'new','sparse.linalg.norm':'old'})
        self.assertEqual(self.measure.call_count,1)
        self.run_report('--check')

    def test_sci_refresh_retains_comparator_backfill_provenance(self):
        history=[{'apis':['special.erf'],'scipy':'1.18.1',
                  'testing_sources_sha256':'original-comparator-code'}]
        self.meta['comparison_refreshes']=history
        (self.root/'results/environment.json').write_text(json.dumps(self.meta))
        self.run_report('--module','special.erf')
        meta=self.assert_preserved({'special.erf'})
        self.assertEqual(meta['comparison_refreshes'],history)
        self.run_report('--check')

    def test_comparator_backfill_keeps_sci_references_and_unselected_measurements(self):
        untouched=self.root/'results/measurements/sparse.linalg/norm.csv.gz'
        original=untouched.read_bytes()
        self.native.side_effect=None
        self.native.return_value.versions={'Boost.Math':'test','GSL':'test'}
        self.native.return_value.evaluate.return_value={'Boost.Math':{'values':[.5]}}
        self.competitors.side_effect=None
        self.competitors.return_value={'SciPy':{'values':[.5]}}
        self.run_report('--comparators-only','--module','special.erf')
        self.assertEqual(untouched.read_bytes(),original)
        updated=read_measurements(self.root/'results/measurements')
        for case in self.cases:
            saved=self.records[case['id']]
            actual=updated[case['id']]
            self.assertEqual(actual['reference'],saved['reference'])
            self.assertEqual(actual['libraries']['Sci'],saved['libraries']['Sci'])
            if case['api']=='special.erf':
                self.assertEqual(actual['libraries']['Boost.Math'],{'values':[.5]})
                self.assertEqual(actual['libraries']['SciPy'],{'values':[.5]})
            else:self.assertEqual(actual,saved)
        self.measure.assert_not_called()
        self.reference.assert_not_called()
        meta=json.loads((self.root/'results/environment.json').read_text())
        self.assertEqual(meta['reference_environment'],self.meta)
        self.assertEqual(meta['comparison_refreshes'][-1]['apis'],['special.erf'])
        self.assertIn('Boost.Math',meta['module_libraries']['special'])
        self.run_report('--check')

    def test_incremental_backfill_reuses_coverage_and_untouched_report_page(self):
        self.run_report('--render')
        untouched=(self.output/'sparse.md').read_bytes()
        baseline_coverage=(self.root/'results/coverage.csv').read_bytes()
        self.native.side_effect=None
        self.native.return_value.versions={'Boost.Math':'fresh','GSL':'fresh'}
        self.native.return_value.evaluate.return_value={'GSL':{'values':[.5]}}
        self.competitors.side_effect=None
        self.competitors.return_value={'SciPy':{'values':[.5]}}
        with patch.object(report,'render',wraps=report.render) as rendered, \
             patch('sci_testing.coverage_audit.audit',side_effect=AssertionError('unnecessary audit')):
            self.run_report('--comparators-only','--module','special.erf')
            self.assertEqual(rendered.call_args.kwargs['modules'],['special'])
        self.assertEqual((self.output/'sparse.md').read_bytes(),untouched)
        self.assertEqual((self.root/'results/coverage.csv').read_bytes(),baseline_coverage)
        self.run_report('--check')

    def test_targeted_render_keeps_unrelated_module_page(self):
        self.run_report('--render')
        untouched=(self.output/'special.md').read_bytes()
        with patch.object(report,'render',wraps=report.render) as rendered:
            self.run_report('--render','--module','sparse.linalg.norm')
            self.assertEqual(rendered.call_args.kwargs['modules'],['sparse'])
        self.assertEqual((self.output/'special.md').read_bytes(),untouched)
        self.run_report('--check')

    def test_fast_check_reuses_metrics_but_exhaustive_check_recomputes_them(self):
        self.run_report('--render')
        cache=self.root/'results/accuracy-summaries.json'
        self.assertEqual(set(json.loads(cache.read_text())['entries']),set(self.cases.apis))
        with patch.object(report,'component_error',side_effect=AssertionError('recomputed metric')), \
             patch('sci_testing.coverage_audit.audit',side_effect=AssertionError('recomputed coverage')):
            self.run_report('--check')
        with patch.object(report,'component_error',wraps=report.component_error) as metric, \
             patch('sci_testing.coverage_audit.audit',wraps=__import__(
                 'sci_testing.coverage_audit',fromlist=['audit']).audit) as coverage:
            self.run_report('--check','--exhaustive')
            self.assertGreater(metric.call_count,0)
            coverage.assert_called_once()

    def test_accuracy_summary_digest_and_per_file_keys_fail_closed(self):
        self.run_report('--render')
        cache=self.root/'results/accuracy-summaries.json'
        original=cache.read_bytes()
        cache.write_bytes(original+b'\n')
        with self.assertRaisesRegex(ValueError,'Accuracy summary digest mismatch'):
            self.run_report('--check')
        cache.write_bytes(original)
        data=json.loads(original)
        api=self.cases.apis[0]
        data['entries'][api]['key']['measurement']='0'*64
        altered=(json.dumps(data,sort_keys=True,separators=(',',':'),ensure_ascii=False)+'\n').encode()
        cache.write_bytes(altered)
        meta=json.loads((self.root/'results/environment.json').read_text())
        import hashlib
        meta['accuracy_summaries_sha256']=hashlib.sha256(altered).hexdigest()
        (self.root/'results/environment.json').write_text(json.dumps(meta))
        with self.assertRaisesRegex(ValueError,'Stale or missing accuracy summary'):
            self.run_report('--check')
        cache.write_bytes(original)

    def test_exhaustive_check_detects_authentic_but_falsified_summary(self):
        self.run_report('--render')
        cache=self.root/'results/accuracy-summaries.json'
        saved=json.loads(cache.read_text())
        saved['entries'][self.cases.apis[0]]['summary']['row']='| falsified |'
        serialized=(json.dumps(saved,sort_keys=True,separators=(',',':'),ensure_ascii=False)+'\n').encode()
        cache.write_bytes(serialized)
        meta=json.loads((self.root/'results/environment.json').read_text())
        import hashlib
        meta['accuracy_summaries_sha256']=hashlib.sha256(serialized).hexdigest()
        (self.root/'results/environment.json').write_text(json.dumps(meta))
        with self.assertRaisesRegex(ValueError,'Stale accuracy summaries'):
            self.run_report('--check','--exhaustive')

    def test_comparator_backfill_is_atomic_if_reporting_fails(self):
        before={p:p.read_bytes() for p in (self.root/'results').rglob('*') if p.is_file()}
        self.native.side_effect=None
        self.native.return_value.versions={'Boost.Math':'test','GSL':'test'}
        self.native.return_value.evaluate.return_value={'Boost.Math':{'values':[.5]}}
        self.competitors.side_effect=None
        self.competitors.return_value={}
        with patch.object(report,'render',side_effect=ValueError('report failed')):
            with self.assertRaisesRegex(ValueError,'report failed'):
                self.run_report('--comparators-only','--module','special.erf')
        self.assertEqual(before,{p:p.read_bytes() for p in before})


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

    def test_failed_render_keeps_original_measurements_and_metadata(self):
        before={p:p.read_bytes() for p in (self.root/'results').rglob('*') if p.is_file()}
        with patch.object(report,'render',side_effect=ValueError('page generation failed')):
            with self.assertRaisesRegex(ValueError,'page generation failed'):self.run_report()
        self.assertEqual(before,{p:p.read_bytes() for p in before})

    def test_staged_hardlinks_preserve_saved_references_and_unselected_files(self):
        source=self.root/'results/measurements'
        untouched=source/'special/erf.csv.gz'
        original=untouched.read_bytes()
        selected=source/'sparse.linalg/norm.csv.gz'
        before=selected.read_bytes()
        with tempfile.TemporaryDirectory(dir=self.root) as directory:
            staged=Path(directory)/'measurements'
            shutil.copytree(source,staged,copy_function=os.link)
            self.assertEqual(os.stat(selected).st_ino,os.stat(staged/'sparse.linalg/norm.csv.gz').st_ino)
            changed=report.refresh_sci(self.cases,self.records,Path('runner'),staged,{'sparse.linalg.norm'})
            self.assertEqual(changed,{'sparse.linalg/norm.csv.gz'})
            self.assertNotEqual(os.stat(selected).st_ino,os.stat(staged/'sparse.linalg/norm.csv.gz').st_ino)
            self.assertEqual((staged/'special/erf.csv.gz').read_bytes(),original)
            self.assertEqual(selected.read_bytes(),before)
            updated=read_measurements(staged)
            self.assertEqual(updated[self.cases[0]['id']]['reference'],self.records[self.cases[0]['id']]['reference'])

    def test_verified_cached_inputs_skip_regeneration_even_when_compressed(self):
        from sci_testing.generate_cases import case_signature
        from sci_testing.case_csv import case_directory_digest
        inputs=self.root/'cases/generated'
        write_cases(inputs,self.cases,compress_above=0,signature=case_signature())
        self.assertTrue((inputs/'special/erf.csv.gz').exists())
        self.assertFalse((inputs/'special/erf.csv').exists())
        self.assertEqual(case_directory_digest(inputs),self.meta['cases_sha256'])
        with patch.object(report,'build_cases',side_effect=AssertionError('regenerated inputs')):
            self.run_report('--module','special')
            self.run_report('--check')
        self.assert_preserved({'special.erf'})

    def test_extend_measures_only_new_io_cases_and_preserves_saved_oracles(self):
        from sci_testing.io_cases import cases as io_cases
        extra=Cases(io_cases(1000,20261020,apis=['io.matlab.matfile_version']))
        combined=Cases(sorted([*self.cases,*extra],key=lambda case:(case['api'],case['id'])))
        old_results=read_measurements(self.root/'results/measurements')
        self.reference.side_effect=None
        self.reference.return_value=[.5]
        self.competitors.side_effect=None
        self.competitors.return_value={'SciPy':{'values':[.5]}}
        from sci_testing.case_set import ExtendedCaseSet
        with patch('sci_testing.generate_cases.extend_historical_cases',
                   side_effect=lambda historical: ExtendedCaseSet(historical,extra)), \
             patch.object(report,'build_cases',side_effect=AssertionError('historical case regeneration')), \
             patch.object(report,'MODULES',['io','sparse','special']), \
             patch('sci_testing.coverage_audit.inventory',return_value=[
                 dict(api=api,canonical=api,status='covered') for api in combined.apis]):
            self.run_report('--extend')
            updated=read_measurements(self.root/'results/measurements')
            for identifier,record in old_results.items():
                self.assertEqual(updated[identifier],record)
            self.assertEqual(len(updated)-len(old_results),len(extra))
            self.assertEqual(self.measure.call_count,1)
            self.assertEqual(self.reference.call_count,2*len(extra))
            self.assertEqual(self.competitors.call_count,len(extra))
            self.native.assert_not_called()
            meta=json.loads((self.root/'results/environment.json').read_text())
            self.assertEqual(meta['reference_environment'],self.meta)
            self.assertEqual(meta['reference_extensions'][0]['apis'],extra.apis)
            self.assertEqual(meta['sci_measurement_commits']['special.erf'],'old')
            self.assertEqual(meta['sci_measurement_commits']['io.matlab.matfile_version'],'new')
            self.run_report('--check')

    def test_extend_rejects_overlapping_historical_api_without_publishing(self):
        from sci_testing.case_set import ExtendedCaseSet
        prior={p:p.read_bytes() for p in (self.root/'results').rglob('*') if p.is_file()}
        with patch('sci_testing.generate_cases.extend_historical_cases',
                   side_effect=lambda historical:ExtendedCaseSet(historical,[
                       {**self.cases[0],'args':{'x':.75}}])):
            with self.assertRaisesRegex(ValueError,'cannot change a historical function'):
                self.run_report('--extend')
        self.assertEqual(prior,{p:p.read_bytes() for p in prior})
        self.measure.assert_not_called()

    def test_extend_rejects_new_api_with_no_success_without_publishing(self):
        from sci_testing.case_set import ExtendedCaseSet
        from sci_testing.io_cases import cases as io_cases
        extra=list(io_cases(1000,20261020,apis=['io.matlab.matfile_version']))
        prior={p:p.read_bytes() for p in (self.root/'results').rglob('*') if p.is_file()}
        inputs_before={p:p.read_bytes() for p in (self.root/'cases/generated').rglob('*') if p.is_file()}
        self.reference.side_effect=None
        self.reference.return_value=[.5]
        self.competitors.side_effect=None
        self.competitors.return_value={}
        self.measure.side_effect=lambda batch,executable:[dict(id=c['id'],error='unimplemented') for c in batch]
        with patch('sci_testing.generate_cases.extend_historical_cases',
                   side_effect=lambda historical:ExtendedCaseSet(historical,extra)), \
             patch.object(report,'MODULES',['io','sparse','special']), \
             patch('sci_testing.coverage_audit.inventory',return_value=[
                 dict(api=api,canonical=api,status='new_api_untested' if api.startswith('io.') else 'covered')
                 for api in [*self.cases.apis,'io.matlab.matfile_version']]):
            with self.assertRaisesRegex(ValueError,'New public APIs without successful comparisons'):
                self.run_report('--extend')
        self.assertEqual(prior,{p:p.read_bytes() for p in prior})
        self.assertEqual(inputs_before,{p:p.read_bytes() for p in inputs_before})

    def test_reordered_runner_output_is_rejected_without_publishing(self):
        before={p:p.read_bytes() for p in (self.root/'results').rglob('*') if p.is_file()}
        self.measure.side_effect=lambda batch,executable:[dict(id=c['id']+'-wrong',values=[.75]) for c in batch]
        with self.assertRaisesRegex(ValueError,'lost, duplicated or reordered'):
            self.run_report('--module','special')
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
