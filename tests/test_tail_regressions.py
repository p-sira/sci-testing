"""Independent tail recipes and append-only measurement preservation."""
import csv
import io
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

import mpmath as mp

from sci_testing import report, tail_regressions
from sci_testing.append_tail_regressions import append_csv, append_files
from sci_testing.case_csv import csv_text, read_measurement_content, write_measurements
from sci_testing.reference_generation import recipe_source


class TailRegressionTests(unittest.TestCase):
    def test_reference_dispatch_and_relative_stability(self):
        cases = list(tail_regressions.cases())
        self.assertEqual(len({c['id'] for c in cases}), len(cases))
        case = next(c for c in cases if c['api']=='special.ndtr' and c['args']['p0']==-38.48)
        self.assertEqual(float(report.stable_reference(case)[0]).hex(), '0x0.0000000000001p-1022')
        self.assertIn('tail_regressions', recipe_source(case)[1])
        historical = {**case, 'id': 'special.ndtr/random-0000'}
        self.assertIn('registered_cases', recipe_source(historical)[1])
        # A false zero is not accepted even far below binary64's subnormal range.
        with patch.object(report, 'high_precision', side_effect=[[mp.mpf(0)], [mp.mpf('1e-600')]]):
            with self.assertRaisesRegex(ValueError, 'Unstable independent reference'):
                report.stable_reference(case)

    def test_positive_mixture_against_independent_reference(self):
        case = next(c for c in tail_regressions.cases() if c['api']=='special.chndtr')
        with mp.workdps(100):
            actual = tail_regressions.reference(case)[0]
            self.assertLess(abs(actual/mp.mpf('9.35193547998677e-80')-1), mp.mpf('2e-15'))

    def test_csv_append_preserves_old_reference_comparators_and_warnings(self):
        old = ('id,library,reference,values,error,warning_0\n'
               'a,Sci,0.10000000000000000000000000001,0.1,,\n'
               'a,SciPy,,-0.0,,"first\nsecond"\n')
        new = ('id,library,reference,values,error,warning_0,warning_1\n'
               'b,Sci,1.2e-320,1.2e-320,,,\n'
               'b,GSL,,,underflow,one,two\n')
        merged = append_csv(old, new, measurements=True)
        result = read_measurement_content(merged)
        self.assertEqual(result['a'], read_measurement_content(old)['a'])
        self.assertEqual(result['b'], read_measurement_content(new)['b'])
        self.assertEqual(result['a']['libraries']['SciPy']['values'][0].hex(), '-0x0.0p+0')
        with self.assertRaisesRegex(ValueError,'Duplicate historical'):
            append_csv(old, old, measurements=True)
        with self.assertRaisesRegex(ValueError,'schema differs'):
            append_csv('id,x\na,1\n','id,y\nb,2\n')

    def test_append_breaks_hardlinks_and_preserves_unselected_files(self):
        with tempfile.TemporaryDirectory() as directory:
            src, dest = Path(directory)/'source', Path(directory)/'staged'
            old='id,group,p0\nspecial.ndtr/random-0000,old,s:-38.0\n'
            new='id,group,p0\nspecial.ndtr/tail-regression-000,new,s:-38.48\n'
            write_measurements(src, {'special/ndtr.csv':old, 'special/erfc.csv':old}, compress_above=0)
            shutil.copytree(src,dest,copy_function=__import__('os').link)
            untouched=(src/'special/erfc.csv.gz').read_bytes()
            append_files(src,dest,{'special/ndtr.csv':new})
            self.assertEqual(csv_text(src/'special/ndtr.csv.gz'),old)
            self.assertEqual(csv_text(dest/'special/ndtr.csv'),old+new.partition('\n')[2])
            self.assertEqual((dest/'special/erfc.csv.gz').read_bytes(),untouched)


if __name__ == '__main__':
    unittest.main()
