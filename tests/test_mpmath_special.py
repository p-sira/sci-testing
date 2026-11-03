"""Public mpmath comparisons for special cases, shapes and conventions."""

import math
import unittest

import mpmath as mp

from sci_testing import mpmath_special as adapter
from sci_testing import special03_cases, special_complex03_cases


def historical(name, *args):
    return {'api': f'special.{name}', 'args': {f'p{i}': x for i,x in enumerate(args)}}


def typed(name, **kwargs):
    return {'api':name, 'args':{'special03_case':True,**kwargs}}


def complex_case(name, **kwargs):
    return {'api':f'special.complex.{name}',
            'args':{'special_complex03_case':True,**kwargs}}


def values(case):
    result=adapter.evaluate(case)
    assert set(result)=={'mpmath'}, (case,result)
    nums=result['mpmath']['values']
    assert nums and all(isinstance(x,float) and math.isfinite(x) for x in nums)
    return nums


class MpmathSpecialTests(unittest.TestCase):
    def assertVectorClose(self, actual, expected, *, tol=2e-13):
        self.assertEqual(len(actual),len(expected))
        for x,y in zip(actual,expected):
            self.assertTrue(math.isclose(x,float(y),rel_tol=tol,abs_tol=tol), (x,y))

    def test_unavailable_and_provider_contract(self):
        self.assertEqual(adapter.evaluate({'api':'stats.ndtr','args':{'p0':1}}),{})
        self.assertEqual(adapter.evaluate(historical('unknown',2)),{})
        self.assertEqual(adapter.evaluate(historical('erfc',2,3)),{})
        self.assertEqual(adapter.evaluate({'api':'special.gamma', 'args':{'p1':2}}),{})
        self.assertEqual(adapter.evaluate(historical('roots_hermite_unavailable',4)),{})
        self.assertEqual(adapter.evaluate({'api':'special.complex.erf','args':{'z':[.5,.4]}}),{})
        self.assertEqual(adapter.evaluate(historical('ellipk',2)),{})

    def test_all_typed_03_apis_and_shapes(self):
        for name in special03_cases.APIS:
            edges=list(special03_cases.cases(1,543,[name]))
            for case in (edges[0],edges[-1]):
                with self.subTest(id=case['id']):
                    result=values(case)
                    self.assertEqual(len(result),2 if name.endswith('hurwitz_zeta') else 1)

        self.assertVectorClose(values(typed('special.bessel.iv_ratio',v=.5,x=1)),[math.tanh(1)])
        self.assertVectorClose(values(typed('special.bessel.iv_ratio_c',v=.5,x=1)),[1-math.tanh(1)])
        self.assertVectorClose(values(typed('special.bessel.iv_ratioinv',v=.5,r=math.tanh(1))),[1])
        self.assertVectorClose(values(typed('special.combinatorics.normalized_gen_harmonic',
                                             j=2,k=3,n=4,a=0)),[.5])
        self.assertVectorClose(values(typed('special.complex.hurwitz_zeta',
                                             s=[2,0],q=[1,0])),[math.pi**2/6,0])
        self.assertVectorClose(values(typed('special.elementary.log1pmx',x=.25)),[math.log1p(.25)-.25])
        self.assertVectorClose(values(typed('special.gamma.digammainv',
                                            y=-.5772156649015329)),[1],tol=1e-10)
        self.assertVectorClose(values(typed('special.gamma.gamma_ratio',a=6,b=5)),[5])
        with mp.workdps(50):
            self.assertVectorClose(values(typed('special.integrals.scaled_exp1',x=1)),
                                   [mp.e*mp.e1(1)])

    def test_all_complex_03_apis_and_branch_conventions(self):
        for name in special_complex03_cases.APIS:
            cases=list(special_complex03_cases.cases(1,2811,[name]))
            for case in (cases[0],cases[-1]):
                with self.subTest(id=case['id']):
                    nums=values(case)
                    expected_len=8 if name.endswith(('airy','airye')) else 4 if name.endswith('fresnel') else 2
                    self.assertEqual(len(nums),expected_len)
        self.assertVectorClose(values(complex_case('exp1',z=[-1,-0.])),
                               [-mp.ei(1),math.pi],tol=1e-12)
        self.assertVectorClose(values(complex_case('expi',z=[-1,-0.])),[mp.ei(-1),0])
        self.assertVectorClose(values(complex_case('lambertw',z=[-1/4,-0.],k=-1)),
                               [mp.lambertw(-mp.mpf(1)/4,-1),0])
        self.assertEqual(values(complex_case('xlogy',x=[0,0],y=[0,0])),[0.,0.])

    def test_elliptic_and_polynomial_identities(self):
        self.assertVectorClose(values(historical('ellipke',.5)),[mp.ellipk(.5),mp.ellipe(.5)])
        self.assertVectorClose(values(historical('ellipd',0)),[math.pi/4])
        self.assertVectorClose(values(historical('ellipdinc',0.7,0)),
                               [(0.7-math.sin(.7)*math.cos(.7))/2])
        self.assertVectorClose(values(historical('ellippi',.2,.5)),[mp.ellippi(.2,.5)])
        self.assertVectorClose(values(historical('ellipj',3.,0)),
                               [math.sin(3),math.cos(3),1,3])
        self.assertVectorClose(values(historical('eval_chebyc',3,.2)),
                               [2*mp.chebyt(3,.1)])
        self.assertVectorClose(values(historical('eval_hermitenorm',3,.2)),
                               [.2**3-3*.2])
        self.assertVectorClose(values(historical('eval_jacobi',2,.5,.3,.2)),
                               [mp.jacobi(2,.5,.3,.2)])
        # Saved alias signature is (x, degree, alpha); Sci's modern
        # eval_genlaguerre order is (degree, alpha, x).
        self.assertVectorClose(values(historical('assoc_laguerre',.3,3,.5)),[1.028])

    def test_sequence_and_complex_shapes(self):
        self.assertVectorClose(values(historical('legendre_p_all',3,.3)),
                               [1,.3,-.365,-.3825])
        with mp.workdps(60):
            expected=[[mp.legenp(n,m,mp.mpf('.3')) if abs(m)<=n else 0
                       for m in (0,1,2,-2,-1)] for n in range(4)]
            self.assertVectorClose(values(historical('assoc_legendre_p_all',3,2,.3)),
                                   [float(y) for row in expected for y in row])
            q0=mp.log(mp.mpf('1.3')/mp.mpf('.7'))/2
            q1=mp.mpf('.3')*q0-1
            self.assertVectorClose(values(historical('lqn',1,.3)),
                                   [q0,q1,1/(1-mp.mpf('.3')**2),q0+mp.mpf('.3')/(1-mp.mpf('.3')**2)])
        self.assertEqual(len(values(historical('lqmn',2,3,.3))),24)
        for x in (2,-2):
            q0=math.log((x+1)/(x-1))/2
            q1=x*q0-1
            self.assertVectorClose(values(historical('lqn',1,x))[:2],[q0,q1])
            self.assertVectorClose(values(historical('lqmn',1,0,x))[:2],
                                   [q0,-1/math.sqrt(3)])
        self.assertEqual(len(values(historical('sph_harm_y_all',2,1,.3,.5))),18)
        self.assertEqual(len(values(historical('sph_harm_y',3,1,.5,.3))),2)
        self.assertEqual(len(values(historical('airy',.5))),4)
        self.assertEqual(len(values(historical('fresnel',.5))),2)
        self.assertEqual(len(values(historical('wofz',[.5,-.25]))),2)

    def test_associated_legendre_recurrences_at_80_digits(self):
        """Every matrix cell agrees with the public Ferrers P / real Q."""
        from sci_testing.mpmath_special import (
            _associated_legendre_p_all,
            _associated_legendre_q_all,
        )
        with mp.workdps(80):
            for n, m, x in ((3, 2, .3), (4, 3, -.3),
                            (5, 4, .999), (3, 3, .02)):
                t=mp.mpf(x)
                actual=_associated_legendre_p_all(n,m,t)
                orders=[*range(m+1), *range(-m,0)]
                self.assertEqual(len(actual), n+1)
                for degree, row in enumerate(actual):
                    self.assertEqual(len(row),2*m+1)
                    for order, value in zip(orders,row):
                        expected=(mp.legenp(degree,order,t)
                                  if abs(order)<=degree else mp.mpf(0))
                        self.assertLess(abs(value-expected),
                                        mp.mpf('1e-68')*max(1,abs(expected)))
                self.assertEqual(len(values(historical('assoc_legendre_p_all',
                                                       n,m,x))), (n+1)*(2*m+1))
            self.assertEqual(len(values(historical('assoc_legendre_p_all',
                                                   3,3,0.))), 28)
            for m,n,x in ((2,3,.3),(2,3,-.3),(2,3,2.),(2,3,-2.)):
                t=mp.mpf(x)
                q,dq=_associated_legendre_q_all(m,n,t)
                self.assertEqual([len(row) for row in q],[n+1]*(m+1))
                self.assertEqual([len(row) for row in dq],[n+1]*(m+1))
                for order in range(m+1):
                    for degree in range(n+1):
                        expected=(mp.legenq(degree,order,t,type=2) if abs(t)<1
                                  else (-1 if t<0 and order%2 else 1)
                                  *mp.re(mp.legenq(degree,order,t,type=3)))
                        self.assertLess(abs(q[order][degree]-expected),
                                        mp.mpf('1e-68')*max(1,abs(expected)))
                # Derivatives are checked against independent public
                # mpmath differentiation for both inner and outer domains.
                for order,degree in ((0,0),(1,2),(2,3)):
                    f=lambda u,degree=degree,order=order: (mp.legenq(degree,order,u,type=2)
                                 if abs(u)<1 else
                                 (-1 if u<0 and order%2 else 1)
                                 *mp.re(mp.legenq(degree,order,u,type=3)))
                    expected=mp.diff(f,t)
                    self.assertLess(abs(dq[order][degree]-expected),
                                    mp.mpf('1e-68')*max(1,abs(expected)))
                self.assertEqual(len(values(historical('lqmn',m,n,x))),
                                 2*(n+1)*(m+1))
            # Forward degree recurrences outside [-1,1] require guard
            # digits: Q_l is the decaying solution, but rounding also
            # excites a growing one.  Degree 50 at x=10 exposes this.
            q,_=_associated_legendre_q_all(1,50,mp.mpf(10))
            expected=mp.re(mp.legenq(50,0,10,type=3))
            self.assertLess(abs(q[0][50]-expected),
                            mp.mpf('1e-68')*abs(expected))

    def test_bessel_derivative_order_shift_at_80_digits(self):
        from sci_testing.mpmath_special import _bessel_derivative
        funcs={'jvp':mp.besselj,'yvp':mp.bessely,'ivp':mp.besseli,
               'kvp':mp.besselk,'h1vp':mp.hankel1,'h2vp':mp.hankel2}
        with mp.workdps(80):
            order=mp.mpf('.5')
            x=mp.mpf('2')
            for name,fn in funcs.items():
                for derivative_order in (0,1,2,3):
                    expected=mp.diff(lambda t,fn=fn:fn(order,t),x,derivative_order)
                    actual=_bessel_derivative(name,order,x,derivative_order)
                    self.assertLess(abs(actual-expected),
                                    mp.mpf('1e-68')*max(1,abs(expected)))
                actual=values(historical(name,.5,2.,2))
                self.assertEqual(len(actual),2 if name in ('h1vp','h2vp') else 1)
                self.assertVectorClose(actual,
                    ([mp.re(_bessel_derivative(name,order,x,2)),
                      mp.im(_bessel_derivative(name,order,x,2))]
                     if name in ('h1vp','h2vp') else
                     [_bessel_derivative(name,order,x,2)]))

    def test_legacy_named_special_inputs(self):
        self.assertVectorClose(values({'api':'special.gamma','args':{'x':2.5}}),
                               [math.gamma(2.5)])
        self.assertVectorClose(values({'api':'special.erf','args':{'x':-.1}}),
                               [math.erf(-.1)])
        self.assertEqual(adapter.evaluate({'api':'special.geterr',
                                           'args':{'explicit':True}}),{})
        self.assertEqual(adapter.evaluate({'api':'special.gamma',
                                           'args':{'x':2.5,'p0':2.5}}),{})

    def test_public_mpmath_numeric_identities(self):
        self.assertVectorClose(values(historical('erfi',.5)),[mp.erfi(.5)])
        self.assertVectorClose(values(historical('jvp',1,2,1)),
                               [(mp.besselj(0,2)-mp.besselj(2,2))/2])
        self.assertVectorClose(values(historical('hankel1e',.5,2)),
                               [mp.re(mp.hankel1(.5,2)*mp.exp(-2*mp.j)),
                                mp.im(mp.hankel1(.5,2)*mp.exp(-2*mp.j))])
        self.assertVectorClose(values(historical('gammainc',2,3)),
                               [1-4*math.exp(-3)])
        self.assertVectorClose(values(historical('betaincc',2,3,.5)),[.3125])
        self.assertVectorClose(values(historical('softmax',[0,0,0])),[1/3]*3)
        self.assertVectorClose(values(historical('logsumexp',[1000,1000,999])),
                               [1000+math.log(2+1/math.e)])
        self.assertVectorClose(values(historical('ln_multinomial',5,[2,3])),[math.log(10)])
        self.assertVectorClose(values(historical('voigt_profile',0,1,0)),
                               [1/math.sqrt(2*math.pi)])

    def test_zero_sequences_and_gaussian_quadrature(self):
        roots=values(historical('roots_legendre',3))
        self.assertVectorClose(roots,[-math.sqrt(3/5),0,math.sqrt(3/5),
                                      5/9,8/9,5/9])
        self.assertVectorClose(values(historical('roots_chebyc',3))[:3],
                               [-math.sqrt(3),0,math.sqrt(3)])
        self.assertVectorClose(values(historical('roots_sh_legendre',3))[:3],
                               [(1-math.sqrt(3/5))/2,.5,(1+math.sqrt(3/5))/2])
        self.assertVectorClose(values(historical('roots_hermitenorm',3))[:3],
                               [-math.sqrt(3),0,math.sqrt(3)])
        self.assertVectorClose(values(historical('roots_genlaguerre',3,.5))[:3],
                               [.6663259077023709,2.8007750541502565,7.032899038147373])
        self.assertVectorClose(values(historical('roots_laguerre',3))[:3],
                               [.4157745567834791,2.2942803602790414,6.28994508293748])
        self.assertVectorClose(values(historical('roots_jacobi',3,.5,.3))[:3],
                               [-.7441752557062968,-.03758275496411322,.6935227165527629])
        self.assertVectorClose(values(historical('jn_zeros',0,2)),
                               [2.4048255576957728,5.5200781102863106])
        self.assertVectorClose(values(historical('jnp_zeros',0,1)),
                               [3.8317059702075125])
        self.assertEqual(len(values(historical('jnyn_zeros',2,3))),12)
        self.assertVectorClose(values(historical('ai_zeros',1))[:2],
                               [-2.338107410459767,-1.018792971647471])
        self.assertVectorClose(values(historical('euler',4)),[1,0,-1,0,5])
        self.assertVectorClose(values(historical('stirling2',8,3)),[966])

    def test_repeated_root_sequences_reuse_fixed_precision_results(self):
        from sci_testing.mpmath_special import _discrete_sequence
        previous = _discrete_sequence.cache_info().hits
        first = values(historical('ai_zeros', 2))
        second = values(historical('ai_zeros', 2))
        self.assertEqual(first, second)
        self.assertGreater(_discrete_sequence.cache_info().hits, previous)

    def test_kelvin_values_are_bessel_components(self):
        with mp.workdps(50):
            x=mp.mpf('.5')
            expected=[mp.ber(0,x),mp.bei(0,x),mp.ker(0,x),mp.kei(0,x)]
            actual=values(historical('kelvin',.5))
            self.assertEqual(len(actual),8)
            self.assertVectorClose(actual[:4],expected)
            self.assertVectorClose(values(historical('berp',.5)),
                                   [mp.diff(lambda t:mp.ber(0,t),x)])


if __name__=='__main__':
    unittest.main()
