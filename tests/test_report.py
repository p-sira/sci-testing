"""Regression checks for metrics, reference definitions and failure reporting."""
from pathlib import Path
import tempfile
import json
from collections import Counter
from unittest.mock import patch
import unittest
import warnings
import mpmath as mp
from sci_testing import report
from sci_testing.case_csv import (case_files, decode, disk_files, encode, files_digest,
                      measurement_files, read_cases, read_measurements, write_files)
from sci_testing.generate_cases import CASES_PER_METHOD


class AccuracyTests(unittest.TestCase):
    def test_sparse_reference_recipes_and_nonfinite_tags(self):
        import math
        from sci_testing.sparse_cases import recipes, reference, tagged
        cases=recipes()
        self.assertEqual(len({c['api'] for c in cases}),98)
        for case in cases:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                values=reference(case)
            self.assertTrue(all(map(math.isfinite,values)),case['id'])
            self.assertEqual(len(values)%2,0,case['id'])
        self.assertEqual(tagged([float('inf'),-float('inf'),float('nan'),-0.]),
                         [1.,0.,-1.,0.,2.,0.,0.,-0.])
        self.assertLess(math.copysign(1.,tagged([-0.])[-1]),0.)
        disconnected=next(c for c in cases if c['args']['raw']['api']=='shortest_path')
        self.assertIn(1.,reference(disconnected)[::2])

    def test_sparse_dense_sweeps_preserve_valid_domains(self):
        from sci_testing.sparse_cases import recipes
        from sci_testing.registered_cases import expanded_registered,reference
        from sci_testing.generate_cases import SEED
        import math,warnings
        apis={c['api'] for c in recipes()}
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            for case in expanded_registered(4,SEED,apis):
                self.assertTrue(all(map(math.isfinite,reference(case))),case['id'])
                raw=case['args']['raw']
                if raw['api']=='maximum_flow':
                    self.assertTrue(all(v>=0 and int(v)==v for row in raw['a'] for v in row))

    def test_binary64_equality_shortcut_keeps_high_precision_differences(self):
        with mp.workdps(160):
            exact=mp.nstr(mp.mpf(.1),165)
            for actual in [.1,.10000000000000002,-.1,0.]:
                self.assertEqual(report.component_error(actual,exact,binary64=True),
                                 report.symmetric_error(actual,mp.mpf(exact)))
            self.assertGreater(report.component_error(.1,'0.1',binary64=False),0)
            with self.assertRaises(ValueError):report.component_error(float('inf'),'inf',binary64=True)

    def test_large_measurements_use_deterministic_lossless_gzip(self):
        import gzip
        from sci_testing.case_csv import write_measurements, MeasurementStore, directory_digest
        case={'id':'stats.variance/a','api':'stats.variance','args':{}}
        record={'id':case['id'],'reference':['0.5'],'libraries':{'Sci':{'values':[.5]}}}
        files=measurement_files([case],{case['id']:record})
        with tempfile.TemporaryDirectory() as directory:
            names=write_measurements(directory,files,compress_above=0)
            self.assertEqual(names,{'stats/variance.csv.gz'})
            path=Path(directory)/next(iter(names));original=path.read_bytes()
            self.assertEqual(gzip.decompress(original).decode(),next(iter(files.values())))
            digest=directory_digest(directory)
            write_measurements(directory,files,compress_above=0)
            self.assertEqual(path.read_bytes(),original)
            self.assertEqual(directory_digest(directory),digest)
            self.assertEqual(MeasurementStore(directory)[case['id']],record)
            self.assertEqual(read_measurements(directory)[case['id']],record)
            write_measurements(directory,files,compress_above=10**9)
            self.assertFalse(path.exists())
            self.assertEqual(MeasurementStore(directory)[case['id']],record)

    def test_density_sweeps_preserve_support_and_discrete_arguments(self):
        import math
        from sci_testing.registered_cases import expanded_registered, reference
        from sci_testing.generate_cases import SEED
        for c in expanded_registered(1000,SEED,{'stats.pdf','stats.logpdf'}):
            self.assertTrue(all(math.isfinite(x) for x in reference(c)),c['id'])

    def test_duplicate_inputs_preserve_case_ids_signed_zero_and_policy_calls(self):
        from types import SimpleNamespace
        sent=[]
        def run(command,**kwargs):
            batch=[json.loads(line) for line in kwargs['input'].splitlines()]
            sent.extend(batch)
            return SimpleNamespace(stdout='\n'.join(json.dumps({'id':c['id'],'values':[c['args']['x']]}) for c in batch))
        cases=[{'id':f'special.erf/{i}','api':'special.erf','args':{'x':x}} for i,x in enumerate([1.,1.,-0.,0.])]
        with patch.object(report.subprocess,'run',side_effect=run):
            rows=report.measure_sci(cases,Path('runner'))
            self.assertEqual(len(sent),3)
            self.assertEqual([r['id'] for r in rows],[c['id'] for c in cases])
            self.assertEqual([json.dumps(r['values']) for r in rows],['[1.0]','[1.0]','[-0.0]','[0.0]'])
            sent.clear()
            for c in cases[:2]:c['api']='special.seterr'
            report.measure_sci(cases[:2],Path('runner'))
            self.assertEqual(len(sent),2)

    def test_large_flattened_reference_csv_cells(self):
        case={'id':'signal.qspline2d/large','api':'signal.qspline2d','args':{}}
        reference=['1.'+'0'*160]*4096
        record={'id':case['id'],'reference':reference,'libraries':{'Sci':{'values':[1.]*4096}}}
        with tempfile.TemporaryDirectory() as directory:
            write_files(directory,measurement_files([case],{case['id']:record}))
            self.assertEqual(read_measurements(directory)[case['id']],record)

    def test_dense_sweeps_preserve_constrained_domains(self):
        import math
        import numpy as np
        from sci_testing.registered_cases import expanded_registered, reference
        from sci_testing.generate_cases import SEED
        apis={'linalg.cossin','stats.fisher_exact','stats.barnard_exact',
              'stats.boschloo_exact','stats.multinomial_pmf',
              'stats.multinomial_logpmf','stats.pointbiserialr','stats.truncate'}
        for case in expanded_registered(4,SEED,apis):
            if '/recipe-' in case['id']:continue
            api,args=case['api'],case['args']
            if api=='linalg.cossin':
                matrix=np.asarray(args['raw']['args'][0])
                np.testing.assert_allclose(matrix.T@matrix,np.eye(len(matrix)),atol=1e-14)
            elif api in {'stats.fisher_exact','stats.barnard_exact','stats.boschloo_exact'}:
                self.assertTrue(all(x>=0 and int(x)==x for x in args['x']))
                self.assertLessEqual(sum(args['x']),64)
            elif api.startswith('stats.multinomial_'):
                self.assertTrue(all(x>=0 and int(x)==x for x in args['x']))
                self.assertAlmostEqual(sum(args['y']),1.)
            elif api=='stats.pointbiserialr':self.assertLessEqual(set(args['x']),{0,1})
            elif api=='stats.truncate':
                self.assertTrue(all(args['p'][0]<x<args['p'][1] for x in args['x']))
            self.assertTrue(all(math.isfinite(x) for x in reference(case)),case['id'])

    def test_case_recipes_are_reproducible_and_preserve_base_cases(self):
        from sci_testing.case_set import CaseSet, profiles
        from sci_testing.generate_cases import dense_cases, SEED
        full=report.build_cases()
        base=read_cases(report.ROOT/'cases/base')
        original=base+list(dense_cases(sorted({c['api'] for c in base}),count=4))
        cases=CaseSet(original,4,SEED)
        self.assertEqual(case_files(cases),case_files(CaseSet(original,4,SEED)))
        by_id={c['id']:c for c in cases}
        self.assertEqual([by_id[c['id']] for c in base],base)
        ids=[c['id'] for c in cases]
        self.assertEqual(len(ids),len(set(ids)))
        self.assertEqual({c['api'].split('.')[0] for c in cases},set(report.MODULES))
        self.assertEqual(full.count,1000)
        for api in cases.apis:
            batch=cases.function(api)
            for method in ['linspace','random']:
                selected=[c for c in batch if c['id'].startswith(f'{api}/{method}-')]
                self.assertEqual(len(selected),4,api)
                if api in profiles() and method=='linspace':
                    self.assertGreater(len({json.dumps(c['args'],sort_keys=True) for c in selected}),1,api)
            json.dumps(batch,allow_nan=False)
        # Check the complete per-function counts on representative scalar,
        # matrix and structured adapters without retaining the full suite.
        for api in ['special.erf','linalg.solve','signal.savgol_filter']:
            batch=full.function(api)
            for method in ['linspace','random']:
                self.assertEqual(sum(c['id'].startswith(f'{api}/{method}-') for c in batch),1000)
        with tempfile.TemporaryDirectory() as directory:
            write_files(directory,case_files(cases))
            self.assertEqual(read_cases(directory),list(cases))
        grid=[by_id[f'special.erf/linspace-{i:04d}']['args']['x'] for i in range(4)]
        self.assertEqual((grid[0],grid[-1]),(-8.,8.))
        self.assertTrue(all(abs((b-a)-16/3)<4e-15 for a,b in zip(grid,grid[1:])))

    def test_generated_check_rejects_stale_or_tampered_files(self):
        from sci_testing.case_csv import write_measurements, directory_digest
        cases=[{'id':'stats.variance/a','api':'stats.variance','group':'ordinary','args':{'values':[1.,2.]}}]
        measurements=[{'id':'stats.variance/a','reference':['.5'],'libraries':{'Sci':{'values':[.5]}}}]
        raw=measurement_files(cases,{'stats.variance/a':measurements[0]})
        inputs=case_files(cases)
        meta={'sci_commit':'revision','sci_source_dirty':False,
              'testing_sources_sha256':report.digest_sources(),
              'cases_sha256':files_digest(inputs),
              'measurements_sha256':files_digest(raw),
              'numpy':'2.5.4','scipy':'1.18.1','filterpy':'1.4.5','mpmath':'1.3.0',
              'rust':'rustc','python':'3.13','platform':'test','processor':'test'}
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); output=root/'pages'
            (root/'results').mkdir()
            write_files(root/'cases/generated',inputs)
            write_measurements(root/'results/measurements',raw,compress_above=0)
            meta['measurements_sha256']=directory_digest(root/'results/measurements')
            (root/'results/environment.json').write_text(json.dumps(meta))
            report.render(cases,{'stats.variance/a':measurements[0]},meta,output)
            before={p:p.read_bytes() for p in root.rglob('*') if p.is_file()}
            report.check_generated(cases,output,root=root)
            self.assertEqual(before,{p:p.read_bytes() for p in before})
            with patch('sci_testing.report.digest_sources',return_value='changed'):
                with self.assertRaisesRegex(ValueError,'Testing sources changed'):
                    report.check_generated(cases,output,root=root)
            for path,bad,expected in [
                (root/'cases/generated/stats/variance.csv',b'bad','Stale case CSV'),
                (root/'results/measurements/stats/variance.csv.gz',b'bad','measurement digest'),
                (output/'stats.md',b'edited table','Stale accuracy page'),
            ]:
                original=path.read_bytes(); path.write_bytes(bad)
                with self.assertRaisesRegex(ValueError,expected):
                    report.check_generated(cases,output,root=root)
                path.write_bytes(original)
            (output/'extra.md').write_text('extra')
            with self.assertRaisesRegex(ValueError,'page set differs'):
                report.check_generated(cases,output,root=root)

    def test_csv_preserves_numbers_shapes_errors_and_multiline_warnings(self):
        for value in [-0., 5e-324, 1e308, [1., -0.], [[1., 2.], [3., 4.]]]:
            self.assertEqual(encode(decode(encode(value))),encode(value))
        for text in ['s:nan', 'v:1 inf', 'm:1 2;3 nan', 'bad:1']:
            with self.assertRaises(ValueError): decode(text)
        cases=[dict(id='special.erf/test',api='special.erf',group='a, "quoted" group',args={'x':-0.})]
        records={'special.erf/test':dict(id='special.erf/test',reference=['0.12345678901234567890123456789'],
            libraries={'Sci':{'values':[-0.,5e-324],'warnings':['line one\nline two','comma, quote "']},
                       'failed':{'error':'ValueError: failure, with details','values':[1.]}})}
        with tempfile.TemporaryDirectory() as directory:
            write_files(directory,measurement_files(cases,records))
            self.assertEqual(read_measurements(directory),records)
            write_files(directory,{'extra/old.csv':'stale'})
            self.assertEqual(disk_files(directory),{'extra/old.csv':'stale'})

    def test_random_streams_do_not_depend_on_other_apis(self):
        from sci_testing.generate_cases import dense_cases
        self.assertEqual(list(dense_cases(['special.erf'])),
                         [c for c in dense_cases(['special.gamma','special.erf']) if c['api']=='special.erf'])

    def test_export_inventory_has_no_missing_function_recipes(self):
        import csv
        from io import StringIO
        rows=list(csv.DictReader(StringIO((report.ROOT/'api_inventory.csv').read_text())))
        aliases={r['api']:r['canonical'] for r in rows}
        apis=set(report.build_cases().apis)
        self.assertFalse(apis-set(aliases))
        self.assertEqual({aliases[a] for a in apis},set(aliases.values()))

    def test_generated_rust_matches_the_pinned_dispatch_sources(self):
        from sci_testing.api_registry import rust_registrations
        from sci_testing.generate_adapters import generated as extracted
        from sci_testing.generate_linalg_adapter import generated as linalg
        from sci_testing.generate_interpolate_adapter import generated as interpolation
        from sci_testing.generate_sparse_adapter import generated as sparse
        from sci_testing.rust_source import formatted
        for filename,generator in [('generated.rs',rust_registrations),('extracted.rs',extracted),
                                   ('linalg_adapter.rs',linalg),('interpolate_adapter.rs',interpolation),('sparse_adapter.rs',sparse)]:
            self.assertEqual((report.ROOT/'src'/filename).read_text(),formatted(generator()),filename)

    def test_native_comparisons_and_individual_library_failures(self):
        from sci_testing.native_compare import NativeComparisons
        native=NativeComparisons()
        self.assertEqual(set(native.versions),{'Boost.Math','GSL'})
        for api,args,want in [('special.erf',{'x':.5},float(mp.erf(.5))),
                             ('special.betainc',{'p0':2.,'p1':3.,'p2':.25},.26171875)]:
            for result in native.evaluate({'api':api,'args':args}).values():
                self.assertNotIn('error',result)
                self.assertAlmostEqual(result['values'][0],want,delta=2e-15)
        results=native.evaluate({'api':'special.iv','args':{'p0':-.5,'p1':1.}})
        self.assertIn('values',results['Boost.Math'])
        self.assertNotIn('error',results['Boost.Math'])
        self.assertIn('error',results['GSL'])

    def test_mixed_successes_keep_metrics_and_count_only_failed_cases(self):
        cases=[dict(id=f'special.erf/{i}',api='special.erf',group='sample',args={'x':.5}) for i in range(2)]
        records={c['id']:dict(reference=['.5'],libraries={'Sci':{'values':[.5]},'other':{'values':[.5]}}) for c in cases}
        records[cases[0]['id']]['libraries']['other']={'error':'domain error'}
        meta=dict(sci_commit='test',sci_source_dirty=False,testing_sources_sha256='test',numpy='test',scipy='test',filterpy='test',mpmath='test',rust='test',python='test',platform='test',processor='test',cases_sha256='test')
        with tempfile.TemporaryDirectory() as directory:
            report.render(cases,records,meta,Path(directory))
            text=(Path(directory)/'special.md').read_text()
            row=next(line for line in text.splitlines() if line.startswith('| `erf`'))
            self.assertIn('0.0 / 0.0; failed 1/2',row)
            self.assertEqual(row.count('failed'),1)

    def test_tables_compare_libraries_with_sci_and_report_missing_baselines(self):
        cases=[dict(id=f'special.erf/{i}',api='special.erf',group='sample',args={'x':.5}) for i in range(4)]
        records={c['id']:dict(reference=['2.0'],libraries={
            'Sci':{'values':[1.]},'SciPy':{'values':[1.+2**-52]},
            'Boost.Math':{'values':[1.]}}) for c in cases}
        records[cases[1]['id']]['libraries']['Sci']={'error':'baseline error'}
        records[cases[2]['id']]['libraries']['SciPy']={'error':'competitor error'}
        records[cases[3]['id']]['libraries']['Sci']={'error':'baseline error'}
        records[cases[3]['id']]['libraries']['SciPy']={'error':'competitor error'}
        meta=dict(sci_commit='test',sci_source_dirty=False,testing_sources_sha256='test',numpy='test',scipy='test',filterpy='test',mpmath='test',rust='test',python='test',platform='test',processor='test',cases_sha256='test')
        for batch in [cases, type('Partition', (), {
                'apis':['special.erf'], 'function':lambda _, api:cases,
                'module_counts':lambda _, module:(1,4) if module=='special' else (0,0),
                '__len__':lambda _:4})()]:
            meta['module_libraries']={m:['Sci','SciPy','Boost.Math'] if m=='special' else [] for m in report.MODULES}
            with self.subTest(partitioned=hasattr(batch,'function')), tempfile.TemporaryDirectory() as directory:
                report.render(batch,records,meta,Path(directory))
                text=(Path(directory)/'special.md').read_text()
                self.assertIn('| Function | Cases | SciPy | Boost.Math |',text)
                self.assertIn('| `erf` | 4; Sci failed 2/4 | 1.0 / 1.0; failed 2/4; Sci unavailable 1/4 | 0.0 / 0.0; Sci unavailable 2/4 |',text)
                self.assertIn('**Sci-special as the reference**',text)
                self.assertIn('(Sci): 2 failed; baseline error',text)
                self.assertIn('(SciPy): 2 failed; competitor error',text)

    def test_symmetric_denominator(self):
        with mp.workdps(160):
            self.assertEqual(report.symmetric_error(1,2),1/report.EPS)
            self.assertEqual(report.symmetric_error(2,1),1/report.EPS)
            self.assertEqual(report.symmetric_error(-1,-2),1/report.EPS)
            self.assertEqual(report.symmetric_error(-1,1),2/report.EPS)
            self.assertEqual(report.symmetric_error(1.,1.+2**-52),1)

    def test_zeros_subnormals_and_nonfinite(self):
        self.assertEqual(report.symmetric_error(0.,-0.),0)
        self.assertEqual(report.symmetric_error(0.,5e-324),0)
        self.assertEqual(report.symmetric_error(0.,1.),1/report.EPS)
        self.assertEqual(report.symmetric_error(float.fromhex('0x1p-1022'),0.),1/report.EPS)
        for x in [float('nan'),float('inf'),-float('inf')]:
            with self.assertRaises(ValueError): report.symmetric_error(x,1)

    def test_reference_keeps_fractional_epsilon(self):
        with mp.workdps(160):
            error=report.symmetric_error(1,1+report.EPS/4)
            self.assertEqual(error,mp.mpf('.25'))

    def test_exact_input_conversion(self):
        with mp.workdps(160):
            r=report.high_precision({'api':'constants.lambda2nu','args':{'x':.1}})[0]
            self.assertEqual(r,mp.mpf(299792458)/mp.mpf(.1))
            self.assertNotEqual(r,mp.mpf(299792458)/mp.mpf('0.1'))

    def test_independent_dft_normalization_and_order(self):
        case={'api':'fft.fft','args':{'z':[[1.,0.],[0.,0.],[0.,0.],[0.,0.]]}}
        self.assertEqual(report.high_precision(case),[1,0,1,0,1,0,1,0])
        case['api']='fft.ifft'
        self.assertEqual(report.high_precision(case),[.25,0,.25,0,.25,0,.25,0])

    def test_discrete_integral_and_correlation(self):
        r=report.high_precision({'api':'integrate.trapezoid','args':{'coordinates':[0.,1.,2.],'y':[0.,1.,4.]}})
        self.assertEqual(r,[3])
        r=report.high_precision({'api':'signal.correlate','args':{'x':[1.,2.],'h':[3.,4.]}})
        self.assertEqual(r,[4,11,6])

    def test_nearly_parallel_cosine_reference_is_stable(self):
        case={'api':'spatial.cosine','args':{'u':[1.,2.,3.],'v':[1.000001,2.,3.]}}
        with mp.workdps(100): low=report.high_precision(case)[0]
        with mp.workdps(160):
            high=report.high_precision(case)[0]
            self.assertLess(abs(low-high)/high,mp.mpf('1e-90'))

    def test_all_comparator_adapters_evaluate_cases(self):
        from sci_testing.case_set import profiles
        cases=read_cases(report.ROOT/'cases/base')+[batch[0] for batch in profiles().values()]
        for case in cases:
            with self.subTest(case=case['id']):
                results=report.competitors(case)
                if case['api'] in ['special.seterr','special.geterr']:continue
                self.assertTrue(results)
                for result in results.values():
                    self.assertTrue('error' in result or 'values' in result)

    def test_failures_are_visible(self):
        cases=[{'id':'a','api':'stats.variance','group':'ordinary','args':{'values':[1.,2.]}}]
        records={'a':{'reference':['.5'],'libraries':{'Sci':{'error':'failed'},'SciPy':{'values':[.5]},'Additional library':{'values':[.5]}}}}
        meta={'sci_commit':'revision','sci_source_dirty':False,'testing_sources_sha256':'digest','numpy':'2.5.4','scipy':'1.18.1','filterpy':'1.4.5','mpmath':'1.3.0','rust':'rustc','python':'3.13','platform':'test','processor':'test','cases_sha256':'inputs'}
        with tempfile.TemporaryDirectory() as d:
            report.render(cases,records,meta,Path(d))
            page=(Path(d)/'stats.md').read_text()
            self.assertIn('1; Sci failed 1/1',page)
            self.assertIn('—; Sci unavailable 1/1',page)
            self.assertIn('## Evaluation failures',page)
            self.assertIn('`stats.variance` (Sci): 1 failed;',page)
            self.assertIn('Example: `a`.',page)
            self.assertNotIn('0.0 / 0.0',page)
            self.assertIn('| Additional library |',page)
            rows=[row for row in page.splitlines() if row.startswith('|')]
            self.assertEqual(len({len(row.split('|')) for row in rows}),1)
            records['a']['libraries']['SciPy']['values']=[]
            with self.assertRaisesRegex(ValueError,'shape mismatch'):
                report.render(cases,records,meta,Path(d))


if __name__ == '__main__': unittest.main()
