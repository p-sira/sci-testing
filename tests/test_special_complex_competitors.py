"""Public SciPy competitor values for all Sci 0.3 complex special overloads."""
import unittest

from scipy import special

from sci_testing.special_complex03_cases import APIS, cases, scipy_comparison


class ComplexCompetitorTests(unittest.TestCase):
    def test_every_overload_has_finite_scipy_values_for_edge_inputs(self):
        for api in APIS:
            batch=[c for c in cases(1,20261020,apis=[api]) if '/edge-' in c['id']]
            for case in batch:
                with self.subTest(case=case['id']):
                    row=scipy_comparison(case)['SciPy']
                    self.assertNotIn('error',row)
                    self.assertEqual(len(row['values']) % 2,0)
                    self.assertEqual(len(row['values']),8 if api.endswith(('.airy','.airye'))
                                     else 4 if api.endswith('.fresnel') else 2)

    def test_hankel_lower_bank_discrepancy_is_visible(self):
        case=dict(api='special.complex.hankel1',args={
            'special_complex03_case':True,'v':0.,'z':[-1.5,-0.]})
        row=scipy_comparison(case)['SciPy']
        expected=special.hankel1(0.,complex(-1.5,-0.))
        self.assertEqual(row['values'],[expected.real,expected.imag])
        self.assertLess(row['values'][0],0)

    def test_reference_oracle_does_not_depend_on_scipy_comparison(self):
        from sci_testing.special_complex03_cases import reference
        case=dict(api='special.complex.hankel1',args={
            'special_complex03_case':True,'v':0.,'z':[-1.5,-0.]})
        actual=reference(case)
        compared=scipy_comparison(case)['SciPy']['values']
        self.assertGreater(float(actual[0]),0)
        self.assertLess(compared[0],0)


if __name__=='__main__':
    unittest.main()
