"""Regression checks for metrics, reference definitions and failure reporting."""
from pathlib import Path
import tempfile
import json
import unittest
import mpmath as mp
import report


class AccuracyTests(unittest.TestCase):
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
        cases=json.loads((report.ROOT/'cases.json').read_text())
        for case in cases:
            with self.subTest(case=case['id']):
                results=report.competitors(case)
                self.assertTrue(results)
                for library,result in results.items():
                    self.assertNotIn('error',result,(library,result))
                    self.assertTrue(result['values'])

    def test_failures_are_visible(self):
        cases=[{'id':'a','api':'stats.variance','group':'ordinary','args':{'values':[1.,2.]}}]
        records={'a':{'reference':['.5'],'libraries':{'Sci':{'error':'failed'},'SciPy':{'values':[.5]}}}}
        meta={'sci_commit':'revision','sci_source_dirty':False,'testing_sources_sha256':'digest','numpy':'2.5.4','scipy':'1.18.1','filterpy':'1.4.5','mpmath':'1.3.0','rust':'rustc','python':'3.13','platform':'test','processor':'test','cases_sha256':'inputs'}
        with tempfile.TemporaryDirectory() as d:
            report.render(cases,records,meta,Path(d))
            page=(Path(d)/'stats.md').read_text()
            self.assertIn('No finite results; failures = 1/1',page)
            self.assertIn('failures = 0/1',page)
            rows=[row for row in page.splitlines() if row.startswith('|')]
            self.assertEqual(len({len(row.split('|')) for row in rows}),1)
            records['a']['libraries']['SciPy']['values']=[]
            with self.assertRaisesRegex(ValueError,'shape mismatch'):
                report.render(cases,records,meta,Path(d))


if __name__ == '__main__': unittest.main()
