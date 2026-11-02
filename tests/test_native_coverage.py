"""Guard the public Boost.Math/GSL scalar comparison coverage.

This inventory is intentionally independent of the C++ registration table:
otherwise removing a native adapter would silently remove its own test.
Entries were checked against the installed Boost.Math special_functions and
distributions headers and GSL gsl_sf_*, gsl_cdf and gsl_randist headers.

Each entry names a Sci public API, one regular finite sample in Sci argument
order, the libraries with a matching public mathematical operation, and the
number of resulting doubles. A provider with no matching API is allowed to
report an explicit unsupported-library error. Shared names alone do not imply
that both providers can evaluate them.

An assertion failure means either an adapter is missing or a claimed mapping
fails even at an ordinary point; it does not imply the functions agree over
all domains, parameter conventions, or branch cuts.
"""

import csv
import hashlib
import math
import re
import unittest
from pathlib import Path

from sci_testing import ROOT
from sci_testing.native_compare import NativeComparisons


BOOST = 'Boost.Math'
GSL = 'GSL'
BOTH = frozenset((BOOST, GSL))


def _inventory():
    """Explicit independent expectation, including each provider's coverage."""
    entries = {}

    def add(name, args, providers=BOTH, outputs=1, module='special'):
        api = f'{module}.{name}'
        if api in entries:
            raise AssertionError(f'duplicate native expectation: {api}')
        assert providers and providers <= BOTH
        assert outputs >= 1
        entries[api] = (tuple(float(value) for value in args), providers, outputs)

    # Error, gamma, beta and their normalized/inverse forms. gsl_cdf_*inv
    # supplies the inverse normalized gamma/beta operations where needed.
    for name, args in (
        ('erf', (.5,)), ('erfc', (.5,)), ('gamma', (2.5,)),
        ('gammaln', (2.5,)), ('loggamma', (2.5,)),
        ('gammasgn', (-.5,)), ('rgamma', (2.5,)),
        ('digamma', (2.5,)), ('trigamma', (2.5,)),
        ('polygamma', (2., 2.5)),
        ('beta', (2.5, 3.5)), ('betaln', (2.5, 3.5)),
        ('gammainc', (2.5, 1.2)), ('gammaincc', (2.5, 1.2)),
        ('gammaincinv', (2.5, .5)), ('gammainccinv', (2.5, .5)),
        ('betainc', (2.5, 3.5, .4)), ('betaincc', (2.5, 3.5, .4)),
        ('betaincinv', (2.5, 3.5, .4)),
        ('betainccinv', (2.5, 3.5, .4)),
    ):
        add(name, args)

    # GSL exposes no dedicated inverse erf; Boost.Math does. The normal CDF
    # inverse is separately audited, so a future GSL conversion is welcome.
    add('erfinv', (.3,), frozenset((BOOST,)))
    add('erfcinv', (.3,), frozenset((BOOST,)))

    # Bessel variants: GSL offers fixed/integer/general orders and scaled
    # variants. Boost.Math offers general orders and exponential rescaling.
    for name, args in (
        ('j0', (1.5,)), ('j1', (1.5,)), ('y0', (1.5,)), ('y1', (1.5,)),
        ('i0', (1.5,)), ('i1', (1.5,)), ('k0', (1.5,)), ('k1', (1.5,)),
        ('i0e', (1.5,)), ('i1e', (1.5,)),
        ('k0e', (1.5,)), ('k1e', (1.5,)),
        ('jn', (2., 1.5)), ('yn', (2., 1.5)), ('kn', (2., 1.5)),
        ('jv', (1.5, 2.)), ('yv', (1.5, 2.)),
        ('iv', (1.5, 2.)), ('kv', (1.5, 2.)),
        ('jve', (1.5, 2.)), ('yve', (1.5, 2.)),
        ('ive', (1.5, 2.)), ('kve', (1.5, 2.)),
        ('spherical_jn', (2., 1.5, 0.)), ('spherical_yn', (2., 1.5, 0.)),
    ):
        add(name, args)

    # Sci's elliptic 'm' is a parameter, while both providers use modulus k.
    # GSL's third-kind characteristic additionally has the opposite sign.
    for name, args in (
        ('ellipk', (.4,)), ('ellipe', (.4,)), ('ellipkm1', (.4,)),
        ('ellipd', (.4,)), ('ellipf', (.5, .4)),
        ('ellipkinc', (.5, .4)), ('ellipeinc', (.5, .4)),
        ('ellippi', (.2, .4)), ('ellippiinc', (.5, .2, .4)),
        ('ellippiinc_bulirsch', (.5, .2, .4)),
        ('ellipdinc', (.5, .4)), ('elliprc', (1., 2.)),
        ('elliprd', (1., 2., 3.)), ('elliprf', (1., 2., 3.)),
        ('elliprj', (1., 2., 3., 4.)),
    ):
        add(name, args)
    add('ellipke', (.4,), outputs=2)
    add('ellipj', (.5, .4), outputs=4)
    # GSL exports RC, RD, RF and RJ, but does not export symmetric RG.
    add('elliprg', (1., 2., 3.), frozenset((BOOST,)))

    # Hypergeometric functions, orthogonal polynomial evaluations, and other
    # separately exported real scalar functions.
    for name, args in (
        ('hyp0f1', (1.5, .2)), ('hyp1f1', (1.5, 2.5, .2)),
        ('hyp2f1', (.5, .75, 1.25, .2)),
        ('eval_legendre', (3., .25)),
        ('eval_gegenbauer', (3., .8, .25)),
        ('eval_genlaguerre', (3., .8, .25)),
        ('eval_hermite', (3., .25)),
        ('eval_hermitenorm', (3., .25)),
        ('exp1', (1.5,)), ('expi', (1.5,)), ('expn', (2., 1.5)),
        ('zetac', (2.5,)),
        ('ndtr', (.5,)), ('ndtri', (.7,)),
        ('lambertw', (.5,)), ('lambertwm1', (-.1,)),
    ):
        add(name, args)

    add('airy', (.5,), outputs=4)
    add('airye', (.5,), outputs=4)

    # Integer combinatorics and generalized Pochhammer are exposed directly
    # by gsl_sf_gamma.h and Boost.Math's factorials/gamma/binomial headers.
    add('factorial', (5.,))
    add('factorial2', (7.,))
    add('combinatorics.ln_factorial', (7.,))
    add('ln_binomial', (7., 3.))
    add('binom', (7., 3.))
    add('comb', (7., 3.))
    add('poch', (2.5, 1.5))

    # Orthogonal polynomial and elliptic relatives with direct public
    # providers. These include shifted argument and parameter conversions.
    add('eval_sh_legendre', (3., .25))
    add('eval_laguerre', (3., .25))
    for name, args in (
        ('eval_chebyt', (3., .25)), ('eval_chebyu', (3., .25)),
        ('eval_sh_chebyt', (3., .25)), ('eval_sh_chebyu', (3., .25)),
        ('eval_jacobi', (3., .5, .75, .25)),
        ('jacobi_zeta', (.5, .4)), ('heuman_lambda', (.5, .4)),
    ):
        add(name, args, frozenset((BOOST,)))

    # GSL's rising factorial evaluates Gamma(a)/Gamma(b) at (b,a-b), too.
    add('gamma.gamma_ratio', (2.5, 3.5))

    # Hurwitz zeta has two arguments, unlike the one-argument Riemann zeta.
    # GSL's gsl_sf_hzeta provides the public two-argument operation.
    add('zeta', (2.5, 1.25), frozenset((GSL,)))

    # Public GSL functions without a direct Boost.Math counterpart in the
    # installed headers. Two pairs return (Si,Ci) and (Shi,Chi).
    for name, args in (
        ('dawsn', (.5,)), ('spence', (.5,)),
        ('hyperu', (.5, 1.25, 2.)),
        ('integrals.scaled_exp1', (1.5,)),
        ('mathieu_a', (2., .25)), ('mathieu_b', (2., .25)),
    ):
        add(name, args, frozenset((GSL,)))
    add('sici', (1.5,), frozenset((GSL,)), outputs=2)
    add('shichi', (1.5,), frozenset((GSL,)), outputs=2)

    # Public Boost.Math functions without a matching GSL special function.
    for name, args in (
        ('owens_t', (.7, .9)), ('powm1', (1.25, .2)),
    ):
        add(name, args, frozenset((BOOST,)))

    # The logistic distribution's CDF/inverse CDF exactly match expit/logit.
    add('expit', (.5,))
    add('logit', (.4,))
    # Boost.Math log1p.hpp also exports log1pmx directly, independently of GSL.
    add('elementary.log1pmx', (.1,))

    # These three Sci statistics APIs are genuinely scalar standard-normal
    # wrappers. The separate stats bridges handle vector/distribution inputs.
    add('normal_cdf', (.5,), module='stats')
    add('normal_pdf', (.5,), module='stats')
    add('normal_ppf', (.7,), module='stats')

    # These Sci special-function CDFs are scalar, with public Boost.Math
    # distribution operations and corresponding GSL CDFs.
    for name, args in (
        ('chdtr', (3., 1.5)), ('chdtrc', (3., 1.5)),
        ('chdtri', (3., .5)), ('stdtr', (3., 1.5)),
        ('stdtrit', (3., .75)), ('fdtr', (3., 4., 1.5)),
        ('fdtrc', (3., 4., 1.5)), ('fdtri', (3., 4., .75)),
        ('gdtr', (2., 3., 1.5)), ('gdtrc', (2., 3., 1.5)),
        ('gdtrix', (2., 3., .75)),
        ('pdtr', (2., 1.5)), ('pdtrc', (2., 1.5)),
        ('pdtri', (2., .75)), ('bdtr', (2., 8., .4)),
        ('bdtrc', (2., 8., .4)), ('bdtri', (2., 8., .75)),
    ):
        add(name, args)
    return entries


EQUIVALENTS = _inventory()

# Every single-provider expectation must explain the missing library API.
# These omissions are intentional. An implementation that computes the same
# function from other public calls is allowed and will be tested if present.
PROVIDER_GAPS = {
    'special.erfinv': (GSL, 'no public GSL inverse-error special function'),
    'special.erfcinv': (GSL, 'no public GSL inverse-complementary-error special function'),
    'special.elliprg': (GSL, 'GSL offers RC/RD/RF/RJ but has no symmetric RG'),
    'special.eval_chebyt': (GSL, 'no GSL direct Chebyshev T polynomial evaluator'),
    'special.eval_chebyu': (GSL, 'no GSL direct Chebyshev U polynomial evaluator'),
    'special.eval_sh_chebyt': (GSL, 'no GSL direct shifted Chebyshev T evaluator'),
    'special.eval_sh_chebyu': (GSL, 'no GSL direct shifted Chebyshev U evaluator'),
    'special.eval_jacobi': (GSL, 'no GSL Jacobi-polynomial evaluator'),
    'special.jacobi_zeta': (GSL, 'GSL elliptic APIs do not include Jacobi zeta'),
    'special.heuman_lambda': (GSL, 'GSL elliptic APIs do not include Heuman lambda'),
    'special.zeta': (BOOST, 'Boost zeta is Riemann zeta; Sci zeta is Hurwitz zeta'),
    'special.dawsn': (BOOST, 'no direct Dawson integral in Boost.Math headers'),
    'special.spence': (BOOST, 'no direct Spence dilogarithm in Boost.Math headers'),
    'special.hyperu': (BOOST, 'no direct confluent hypergeometric U in Boost.Math headers'),
    'special.integrals.scaled_exp1': (BOOST, 'Boost expint is unscaled'),
    'special.mathieu_a': (BOOST, 'Boost.Math has no Mathieu characteristic value function'),
    'special.mathieu_b': (BOOST, 'Boost.Math has no Mathieu characteristic value function'),
    'special.sici': (BOOST, 'Boost.Math lacks Si/Ci integral functions'),
    'special.shichi': (BOOST, 'Boost.Math lacks Shi/Chi integral functions'),
    'special.owens_t': (GSL, 'GSL has no Owen T special function'),
    'special.powm1': (GSL, 'GSL integer power does not cover arbitrary real powm1'),
}

# Deliberate complex-ABI boundaries, not unimplemented statistics coverage.
# Statistics with vector or distribution inputs are audited separately below.
NONSCALAR_GAPS = {
    'special.complex.gamma': 'complex input/output needs the separate complex ABI',
    'special.complex.airy': 'complex input and four complex outputs',
}


# Match Sci's named distribution parameters, not the function registration
# table in native/stats_distributions.cpp. The audit protects each of the
# five operations for every family, including library-specific omissions.
# (parameters, interior sample observation)
STATS_DISTRIBUTIONS = {
    'Normal': ((1.5, 2.25), 1.),
    'Beta': ((2.5, 3.5), .42),
    'Gamma': ((2.5, 1.5), 3.1),
    'Erlang': ((3., 2.), 4.1),
    'ChiSquared': ((3.5,), 1.2),
    'Chi': ((3.,), 1.2),
    'StudentT': ((4.5,), .6),
    'F': ((3., 7.), 1.2),
    'Exponential': ((2.,), 1.2),
    'Cauchy': ((1., 2.), .6),
    'Laplace': ((1., 2.), .6),
    'LogNormal': ((.3, .7), 1.2),
    'Uniform': ((-2., 3.), 1.2),
    'Weibull': ((1.5, 2.), 1.2),
    'Pareto': ((1., 3.), 1.2),
    'InverseGamma': ((3., 2.), 1.2),
    'Gumbel': ((1., 2.), .6),
    'Triangular': ((-2., 3., 1.), 1.2),
    'Bernoulli': ((.3,), 0.),
    'Binomial': ((10., .3), 3.),
    'Poisson': ((2.,), 3.),
    'Geometric': ((.3,), 3.),
    'NegativeBinomial': ((3., .3), 3.),
    'Hypergeometric': ((30., 10., 8.), 3.),
    'DiscreteUniform': ((1., 6.), 3.),
    'Levy': ((0., 1.), 1.),
}
STATS_OPERATIONS = ('pdf', 'logpdf', 'cdf', 'sf', 'ppf')
BOOST_DISTRIBUTION_GAPS = frozenset(('Chi', 'DiscreteUniform', 'Levy'))
GSL_DISTRIBUTION_GAPS = frozenset(('InverseGamma', 'Triangular', 'DiscreteUniform', 'Levy'))
GSL_DISCRETE_PPF_GAPS = frozenset((
    'Bernoulli', 'Binomial', 'Poisson', 'Geometric',
    'NegativeBinomial', 'Hypergeometric',
))


# Sci vector statistics use x, an optional paired sample y, and parameters p.
# Unlike the real scalar ABI, the output may contain several components.
# Samples have variation and contain no ties to exercise ordinary domains.
_SAMPLE = (1., 2., 3., 4., 7., 9.)
_PAIR = (4., 1., 6., 3., 9., 5.)
STATS_SUMMARY_CASES = {
    'mean': (_SAMPLE, (), (), 1),
    'variance': (_SAMPLE, (), (1.,), 1),
    'std': (_SAMPLE, (), (1.,), 1),
    'sem': (_SAMPLE, (), (1.,), 1),
    'variation': (_SAMPLE, (), (1.,), 1),
    'skew': (_SAMPLE, (), (1.,), 1),
    'kurtosis': (_SAMPLE, (), (1., 1.), 1),
    'median': (_SAMPLE, (), (), 1),
    'median_abs_deviation': (_SAMPLE, (), (1.,), 1),
    'quantile': (_SAMPLE, (), (.4,), 1),
    'scoreatpercentile': (_SAMPLE, (), (40.,), 1),
    'iqr': (_SAMPLE, (), (), 1),
    'describe': (_SAMPLE, (), (1., 1.), 7),
    'zscore': (_SAMPLE, (), (1.,), len(_SAMPLE)),
    'pearsonr': (_SAMPLE, _PAIR, (0.,), 2),
    'spearmanr': (_SAMPLE, _PAIR, (0.,), 2),
    'spearmanrho': (_SAMPLE, _PAIR, (0.,), 2),
    'pointbiserialr': ((0., 1., 0., 1., 1., 0.), _PAIR, (), 2),
}
BOOST_SUMMARY_GAPS = frozenset((
    'quantile', 'scoreatpercentile', 'iqr', 'spearmanr', 'spearmanrho',
))


# Registrations in the independent native translation units are deliberately
# written here, rather than imported from the modules being audited. A removed
# C++ registration therefore fails the coverage test.
SPECIAL_EXTRA_CASES = {
    'jvp': ((.5, 2., 1.), BOTH, 1),
    'yvp': ((.5, 2., 1.), BOTH, 1),
    'ivp': ((.5, 2., 1.), BOTH, 1),
    'kvp': ((.5, 2., 1.), BOTH, 1),
    'assoc_legendre_p': ((3., 1., .25), BOTH, 1),
    'spherical_in': ((2., 1.5, 0.), BOTH, 1),
    'spherical_kn': ((2., 1.5, 0.), BOTH, 1),
    'eval_sh_jacobi': ((3., 2.5, 1.5, .3), frozenset((BOOST,)), 1),
    'eval_chebyc': ((3., .3), frozenset((BOOST,)), 1),
    'eval_chebys': ((3., .3), frozenset((BOOST,)), 1),
    'perm': ((7., 3.), frozenset((BOOST,)), 1),
    'jn_zeros': ((2., 3.), BOTH, 3),
    'yn_zeros': ((2., 3.), frozenset((BOOST,)), 3),
    'legendre_p': ((3., .3), BOTH, 1),
    'sinc': ((.3,), BOTH, 1),
    'exprel': ((.3,), frozenset((GSL,)), 1),
}

PROBABILITY_EXTRA_CASES = {
    'nbdtr': ((2., 5., .35), BOTH),
    'nbdtrc': ((2., 5., .35), BOTH),
    'nbdtri': ((2., 5., .35), BOTH),
    'chndtr': ((3.5, 5., .7), frozenset((BOOST,))),
    'chndtrix': ((.35, 5., .7), frozenset((BOOST,))),
    'chndtridf': ((4.2, .3, .8), frozenset((BOOST,))),
    'chndtrinc': ((4.2, 3.3, .3), frozenset((BOOST,))),
    'ncfdtr': ((3., 8., .7, 1.5), frozenset((BOOST,))),
    'ncfdtri': ((3., 8., .7, .35), frozenset((BOOST,))),
    'ncfdtrinc': ((3., 8., .35, 1.5), frozenset((BOOST,))),
    'ncfdtridfd': ((3., .25637506749643246, .5, .5), frozenset((BOOST,))),
    'ncfdtridfn': ((.25637506749643246, 5., .5, .5), frozenset((BOOST,))),
    'nctdtr': ((3., .5, .3), frozenset((BOOST,))),
    'nctdtrit': ((3., .5, .35), frozenset((BOOST,))),
    'nctdtrinc': ((3., .4, .3), frozenset((BOOST,))),
    'nctdtridf': ((.4120716320853016, .5, .3), frozenset((BOOST,))),
    'btdtria': ((.35, 2.5, .5), frozenset((BOOST,))),
    'btdtrib': ((2.5, .35, .5), frozenset((BOOST,))),
    'bdtrin': ((2., .35, .5), frozenset((BOOST,))),
    'nbdtrik': ((.35, 5., .5), frozenset((BOOST,))),
    'nbdtrin': ((2., .35, .5), frozenset((BOOST,))),
    'pdtrik': ((.35, 5.), frozenset((BOOST,))),
    'gdtrib': ((1.5, .35, 2.5), frozenset((BOOST,))),
    'chdtriv': ((.35, 2.5), frozenset((BOOST,))),
    'kolmogorov': ((1.4,), frozenset((BOOST,))),
    'kolmogi': ((.65,), frozenset((BOOST,))),
}

STATS_EXTRA_CASES = {
    'trim_mean': (_SAMPLE, (), (.2,), 1, frozenset((GSL,))),
    'tmean': (_SAMPLE, (), (0., 10., 1., 1., 1.), 1, BOTH),
    'tvar': (_SAMPLE, (), (0., 10., 1., 1., 1.), 1, BOTH),
    'tstd': (_SAMPLE, (), (0., 10., 1., 1., 1.), 1, BOTH),
    'tsem': (_SAMPLE, (), (0., 10., 1., 1., 1.), 1, BOTH),
    'tmin': (_SAMPLE, (), (1., 1.), 1, frozenset((GSL,))),
    'tmax': (_SAMPLE, (), (10., 1.), 1, frozenset((GSL,))),
    'ttest_1samp': (_SAMPLE, (), (2., 0.), 3, BOTH),
    'ttest_rel': (_SAMPLE, _PAIR, (0.,), 3, BOTH),
    'ttest_ind': (_SAMPLE, _PAIR, (1., 0.), 3, BOTH),
    'linregress': (_SAMPLE, _PAIR, (0.,), 6, BOTH),
    'dirichlet_pdf': ((.3, .7), (2., 3.), (), 1, frozenset((GSL,))),
    'dirichlet_logpdf': ((.3, .7), (2., 3.), (), 1, frozenset((GSL,))),
    'multinomial_pmf': ((1., 2.), (.4, .6), (), 1, frozenset((GSL,))),
    'multinomial_logpmf': ((1., 2.), (.4, .6), (), 1, frozenset((GSL,))),
    'multivariate_normal_pdf': ((0., 1.), (1., 2.), (), 1, frozenset((GSL,))),
    'multivariate_normal_logpdf': ((0., 1.), (1., 2.), (), 1, frozenset((GSL,))),
}
STATS_EXTRA_GAUSSIAN_COVARIANCE = ((2., .3), (.3, 1.))

GSL_COMPLEX_CASES = frozenset(('digamma', 'spence', 'gamma', 'loggamma', 'rgamma'))
VECTOR_SPECIAL_CASES = {
    'logsumexp': ((-1., 0., 1., 2.), 1),
    'softmax': ((-1., 0., 1., 2.), 4),
    'log_softmax': ((-1., 0., 1., 2.), 4),
}
COMBINATORICS_EXTRA_CASES = {
    'bernoulli': ((8.,), frozenset((BOOST,)), 9),
    'harmonic': ((10.,), BOTH, 1),
    'gen_harmonic': ((10., 2.), frozenset((GSL,)), 1),
}
SPHERICAL_EXTRA_CASES = {
    'sph_harm_y': ((3., 1., .4, .6), BOTH, 2),
    'sph_legendre_p': ((3., 1., .3), BOTH, 1),
}


def _saved_report_names(module):
    """The saved 405/169 API universe, including historical recipe aliases.

    The published accuracy tables contain the exact sampled names, including
    the nested `special.complex.*` and `stats.contingency.*` fixtures. The
    public-inventory CSV also includes *unsampled* public aliases, so its
    cardinality cannot replace these report rows.
    """
    path = ROOT / 'docs' / 'accuracy' / f'{module}.md'
    return set(re.findall(r'^\| `([^`]+)` \| [0-9]+ \|',
                          path.read_text(), re.MULTILINE))


# Classification is intentionally conservative: a name not found in a bridge
# does not establish that Boost.Math or GSL lacks the mathematical operation.
# A reason of "unconfirmed" calls out the still-open investigation explicitly.
# These strings are selected *per API*, not used as global unsupported flags.
_LEGACY_SPECIAL = frozenset('''
    assoc_laguerre cbrt exp10 exp2 expm1 geterr log1p psi round seterr
'''.split())
_ARRAY_SPECIAL = frozenset('''
    ai_zeros bi_zeros erf_zeros fresnel_zeros fresnelc_zeros fresnels_zeros
    bei_zeros beip_zeros ber_zeros berp_zeros kei_zeros keip_zeros
    ker_zeros kerp_zeros kelvin_zeros y0_zeros y1_zeros y1p_zeros
    ynp_zeros jnjnp_zeros jnp_zeros jnyn_zeros
    assoc_legendre_p_all legendre_p_all sph_legendre_p_all sph_harm_y_all
    pbdn_seq pbdv_seq pbvv_seq pro_cv_seq obl_cv_seq
    chebyc chebys chebyt chebyu gegenbauer genlaguerre hermite hermitenorm
    jacobi laguerre legendre sh_chebyt sh_chebyu sh_jacobi sh_legendre
    roots_chebyc roots_chebys roots_chebyt roots_chebyu roots_gegenbauer
    roots_genlaguerre roots_hermite roots_hermitenorm roots_jacobi
    roots_laguerre roots_legendre roots_sh_chebyt roots_sh_chebyu
    roots_sh_jacobi roots_sh_legendre
    bernoulli euler diric lmbda lqmn lqn riccati_jn riccati_yn
    multinomial ln_multinomial softmax log_softmax
    shichi sici airy airye fresnel itairy kelvin
'''.split())
_COMPLEX_SPECIAL = frozenset('''
    h1vp h2vp hankel1 hankel1e hankel2 hankel2e wofz
'''.split())
_MISMATCH_SPECIAL = frozenset('''
    bdtrik fdtridfd gdtria
    nrdtrimn nrdtrisd stdtridf smirnov smirnovi ndtri_exp log_ndtr
    sph_legendre_p
    pbdv pbvv pbwa mathieu_cem mathieu_sem mathieu_modcem1
    mathieu_modsem1 mathieu_modcem2 mathieu_modsem2
    pro_ang1 pro_ang1_cv pro_rad1 pro_rad1_cv pro_rad2 pro_rad2_cv
    obl_ang1 obl_ang1_cv obl_rad1 obl_rad1_cv obl_rad2 obl_rad2_cv
'''.split())
# Header review scope for the following classifications:
# /usr/include/boost/math/{special_functions,distributions,statistics}/
# /usr/include/gsl/{gsl_sf_*.h,gsl_cdf.h,gsl_randist.h,gsl_statistics_double.h}
# These categories mean the *whole Sci callable* was not matched. They do
# not deny that the libraries expose useful component functions.
_DERIVED_SPECIAL = frozenset('''
    agm boxcox boxcox1p cosdg cosm1 cotdg entr huber inv_boxcox
    inv_boxcox1p kl_div log_expit pseudo_huber radian rel_entr sindg
    softplus tandg xlog1py xlogy
'''.split())
_BESSEL_UNCONFIRMED = frozenset('''
    bei beip ber berp bessel.iv_ratio bessel.iv_ratio_c bessel.iv_ratioinv
    besselpoly it2i0k0 it2j0y0 it2struve0 iti0k0 itj0y0 itmodstruve0
    itstruve0 kei keip ker kerp modstruve struve voigt_profile
'''.split())
_ELLIPTIC_UNCONFIRMED = frozenset('''
    cel cel1 cel2 cel3 el1 el2 el3 ellip_harm ellip_harm_2 ellip_normal
'''.split())
_OTHER_SPECIAL_UNCONFIRMED = frozenset('''
    combinatorics.normalized_gen_harmonic erfcx erfi factorialk
    gamma.digammainv gen_harmonic log_wright_bessel lpmv mathieu_even_coef
    mathieu_odd_coef modfresnelm modfresnelp multigammaln obl_cv pro_cv
    stirling2 tklmbda wright_bessel wrightomega
'''.split())

_INFERENCE_STATS = frozenset('''
    alexandergovern anderson anderson_ksamp ansari barnard_exact bartlett
    binomtest boschloo_exact brunnermunzel chatterjeexi chi2_contingency
    chisquare combine_pvalues cramervonmises cramervonmises_2samp dunnett
    epps_singleton_2samp f_oneway fisher_exact fligner friedmanchisquare
    goodness_of_fit jarque_bera kendalltau kruskal ks_1samp ks_2samp
    kurtosistest levene mannwhitneyu median_test mood normaltest
    obrientransform page_trend_test poisson_means_test power_divergence
    quantile_test ranksums shapiro skewtest somersd tukey_hsd weightedtau
    wilcoxon
'''.split())
_DESCRIPTIVE_STATS = frozenset('''
    circmean circstd circvar differential_entropy directional_stats
    entropy expectile gmean gstd gzscore hmean kstat kstatvar lmoment
    mode moment percentileofscore pmean rankdata sigmaclip tiecorrect
    trim1 trimboth zmap
'''.split())
_RESAMPLING_STATS = frozenset('''
    bootstrap bws_test goodness_of_fit monte_carlo_test permutation_test
    sobol_indices false_discovery_control power
'''.split())
_REGRESSION_STATS = frozenset('''
    linregress siegelslopes theilslopes logrank multiscale_graphcorr
    wasserstein_distance wasserstein_distance_nd energy_distance
'''.split())
_FITTING_STATS = frozenset('''
    boxcox boxcox_llf boxcox_normmax boxcox_normplot
    yeojohnson yeojohnson_llf yeojohnson_normmax yeojohnson_normplot
    ppcc_max ppcc_plot probplot
'''.split())
_STRUCTURED_STATS = frozenset('''
    binned_statistic binned_statistic_2d binned_statistic_dd cumfreq relfreq
    categorical_cdf categorical_pmf bayes_mvs mvsdist ecdf estimated_cdf
    make_distribution fit order_statistic truncate
    ttest_ind_from_stats
'''.split())
_CONFIRMED_UNADAPTED = {
}


def _known_coverage_names(module):
    """Explicit tested adapter names; no introspection of native registries."""
    if module == 'special':
        return ({name.removeprefix('special.') for name in EQUIVALENTS
                 if name.startswith('special.')}
                | set(SPECIAL_EXTRA_CASES) | set(PROBABILITY_EXTRA_CASES)
                | set(VECTOR_SPECIAL_CASES) | set(COMBINATORICS_EXTRA_CASES)
                | set(SPHERICAL_EXTRA_CASES)
                | {'complex.' + name for name in GSL_COMPLEX_CASES}
                | {'ln_factorial'})  # historical flat recipe for canonical combinatorics.ln_factorial
    return ({name.removeprefix('stats.') for name in EQUIVALENTS
             if name.startswith('stats.')}
            | set(STATS_SUMMARY_CASES) | set(STATS_EXTRA_CASES)
            | set(STATS_OPERATIONS))


def remaining_native_coverage(module):
    """Every unmapped saved API paired with a scoped, reviewable reason.

    `unconfirmed` is a research status, not a claim that public equivalents
    do not exist. This result can be printed or inspected to get the complete
    concrete list, with no hidden implicit omissions or wildcard exclusions.
    """
    remaining = _saved_report_names(module) - _known_coverage_names(module)
    outcome = {}
    for name in sorted(remaining):
        qualified = f'{module}.{name}'
        if qualified in _CONFIRMED_UNADAPTED:
            reason = _CONFIRMED_UNADAPTED[qualified]
        elif module == 'special' and name.startswith('complex.'):
            reason = ('complex input/output or branch conventions; no verified '
                      'matching full-function public Boost/GSL adapter')
        elif module == 'special' and name in _LEGACY_SPECIAL:
            reason = ('historical recipe alias; no current independent Sci 0.3 '
                      'public callable confirmed with this saved signature')
        elif module == 'special' and name in _ARRAY_SPECIAL:
            reason = ('vector, coefficient array or multi-output signature; '
                      'matching complete result needs a separate shaped ABI')
        elif module == 'special' and name in _COMPLEX_SPECIAL:
            reason = ('complex mathematical result or nontrivial branch convention; '
                      'current real-scalar bridge does not cover full Sci contract')
        elif module == 'special' and name in _MISMATCH_SPECIAL:
            reason = ('nontrivial inverse, special normalization or multi-output '
                      'contract; direct full-function equivalence unconfirmed')
        elif module == 'special' and name in _DERIVED_SPECIAL:
            reason = ('Sci-specific elementary/limit-preserving transform; '
                      'no confirmed direct named public Boost/GSL full-function API')
        elif module == 'special' and name in _BESSEL_UNCONFIRMED:
            reason = ('Kelvin, Struve, integrated Bessel, Bessel-ratio or Voigt '
                      'contract; no confirmed corresponding full public function')
        elif module == 'special' and name in _ELLIPTIC_UNCONFIRMED:
            reason = ('Bulirsch or ellipsoidal integral contract; Carlson/Legendre '
                      'public primitives exist but complete Sci function unmatched')
        elif module == 'special' and name in _OTHER_SPECIAL_UNCONFIRMED:
            reason = ('specialized gamma, hypergeometric, polynomial or inverse '
                      'contract; matching complete public function unconfirmed')
        elif module == 'stats' and name in {'abs', 'exp', 'log', 'kstest'}:
            reason = ('historical stats compatibility recipe; no matching '
                      'full-function stats distribution or summary comparator registered')
        elif module == 'stats' and name in {
                'make_distribution', 'fit', 'truncate', 'order_statistic',
                'ecdf', 'estimated_cdf', 'bayes_mvs', 'mvsdist', 'probplot',
                'boxcox_normplot', 'yeojohnson_normplot', 'ppcc_plot'}:
            reason = ('structured distribution, fitting or diagnostic result; '
                      'no verified matching public full-operation ABI')
        elif module == 'stats' and name.startswith(('contingency.', 'qmc.')):
            reason = ('structured table, contingency or quasi-random vectors; '
                      'no verified matching full-output public library API')
        elif module == 'stats' and name in {
                'multivariate_normal_pdf', 'multivariate_normal_logpdf',
                'multivariate_t_pdf', 'multivariate_t_logpdf',
                'dirichlet_pdf', 'dirichlet_logpdf',
                'multinomial_pmf', 'multinomial_logpmf'}:
            reason = ('matrix/vector parameters and multivariate density; '
                      'separate multivariate ABI and parameter contract needed')
        elif module == 'stats' and name in _INFERENCE_STATS:
            reason = ('full hypothesis-test statistic and tail/p-value contract; '
                      'public Boost/GSL components exist but an exact equivalent '
                      'with all fixture alternatives remains unverified')
        elif module == 'stats' and name in _RESAMPLING_STATS:
            reason = ('resampling scheme, explicit draws or distribution of test '
                      'statistics; no verified matching full-operation public ABI')
        elif module == 'stats' and name in _DESCRIPTIVE_STATS:
            reason = ('vector statistic with Sci bias, rank, weight or boundary '
                      'conventions; no confirmed named public full-function equivalent')
        elif module == 'stats' and name in _REGRESSION_STATS:
            reason = ('paired vector regression, rank inference or distance '
                      'algorithm; complete public operation unverified')
        elif module == 'stats' and name in _FITTING_STATS:
            reason = ('transformation optimization or diagnostic curve output; '
                      'full Sci fitting objective/return contract not verified')
        elif module == 'stats' and name in _STRUCTURED_STATS:
            reason = ('vector, structured binning, fit or transformed-result '
                      'contract; no verified matching full-function public API')
        else:
            reason = ('unconfirmed: no verified matching full Sci operation '
                      'among the public Boost.Math/GSL interfaces inspected; '
                      'further signature/semantic review required')
        outcome[qualified] = reason
    return outcome


class NativeCoverageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.native = NativeComparisons()

    def test_expected_functions_are_public_exports_and_gaps_are_explained(self):
        with (ROOT / 'api_inventory.csv').open(newline='') as stream:
            public = {row['api']: row for row in csv.DictReader(stream)}
        self.assertFalse(EQUIVALENTS.keys() & NONSCALAR_GAPS.keys())
        expected_gaps = {
            (api, provider)
            for api, (_, supported, _) in EQUIVALENTS.items()
            for provider in BOTH - supported
        }
        self.assertEqual(expected_gaps, {
            (api, provider) for api, (provider, _) in PROVIDER_GAPS.items()
        })
        for api, (_, explanation) in PROVIDER_GAPS.items():
            self.assertGreater(len(explanation), 20, api)
        for api in (*EQUIVALENTS, *NONSCALAR_GAPS):
            with self.subTest(api=api):
                self.assertTrue(api in public, f'{api}: missing from public API inventory')
                self.assertTrue(public[api]['canonical'], api)
        for api, explanation in NONSCALAR_GAPS.items():
            self.assertGreater(len(explanation), 20, api)
        self.assertTrue(any(api.startswith('special.') for api in EQUIVALENTS))
        self.assertTrue(any(api.startswith('stats.') for api in EQUIVALENTS))

    def test_known_native_equivalents_have_working_library_adapters(self):
        """Report *all* missing adapters together so a mapping PR is reviewable."""
        failures = []
        for api, (arguments, required, output_count) in sorted(EQUIVALENTS.items()):
            name = api.rsplit('.', 1)[-1]
            encoded = name.encode('ascii')
            registered_arity = self.native.library.comparison_supported(encoded)
            if registered_arity != len(arguments):
                failures.append(f'{api}: shared arity {registered_arity}, expected {len(arguments)}')
                continue
            registered_outputs = self.native.library.comparison_output_count(encoded)
            if registered_outputs != output_count:
                failures.append(f'{api}: outputs {registered_outputs}, expected {output_count}')
                continue
            case = {'api': api, 'args': {f'p{i}': value for i, value in enumerate(arguments)}}
            try:
                observed = self.native.evaluate(case)
            except (TypeError, ValueError, OSError) as exc:
                failures.append(f'{api}: native wrapper raised {exc!r}')
                continue
            for provider in (BOOST, GSL):
                row = observed.get(provider)
                if provider in required:
                    if row is None or 'error' in row:
                        failures.append(f'{api} / {provider}: {row!r}')
                    elif len(row.get('values', ())) != output_count or not all(
                            isinstance(value, (float, int)) and math.isfinite(value)
                            for value in row['values']):
                        failures.append(f'{api} / {provider}: invalid outputs {row!r}')
                elif row is not None:
                    # Provider-specific API absence must be visible and cannot
                    # masquerade as a runtime domain/overflow failure.
                    if 'values' in row:
                        self.assertTrue(all(math.isfinite(v) for v in row['values']),
                                        f'{api} / {provider}: nonfinite optional adapter')
                    elif not any(word in row.get('error', '').lower() for word in
                                 ('unavailable', 'unsupported', 'not supported', 'not implemented')):
                        failures.append(f'{api} / {provider}: unclear provider gap {row!r}')
        self.assertFalse(failures, '\n' + '\n'.join(failures))

    def test_stats_distribution_families_have_per_library_operation_coverage(self):
        """Every expected family/operation must be connected to the live dispatcher."""
        with (ROOT / 'api_inventory.csv').open(newline='') as stream:
            public = {row['api'] for row in csv.DictReader(stream)}
        self.assertLessEqual({'stats.' + op for op in STATS_OPERATIONS}, public)
        self.assertEqual(len(STATS_DISTRIBUTIONS), 26)
        self.assertEqual(len(STATS_DISTRIBUTIONS) - len(BOOST_DISTRIBUTION_GAPS), 23)
        self.assertEqual(len(STATS_DISTRIBUTIONS) - len(GSL_DISTRIBUTION_GAPS), 22)

        failures = []
        for family, (parameters, observation) in sorted(STATS_DISTRIBUTIONS.items()):
            for operation in STATS_OPERATIONS:
                expected = {
                    lib for lib, excluded in (
                        (BOOST, BOOST_DISTRIBUTION_GAPS),
                        (GSL, GSL_DISTRIBUTION_GAPS),
                    ) if family not in excluded and not (
                        lib == GSL and operation == 'ppf'
                        and family in GSL_DISCRETE_PPF_GAPS
                    )
                }
                for index, library in enumerate((BOOST, GSL)):
                    registered = self.native.library.stats_distribution_supported(
                        index, family.encode('ascii'), operation.encode('ascii'))
                    expected_arity = len(parameters) if library in expected else 0
                    if registered != expected_arity:
                        failures.append(f'{family} {operation} / {library}: '
                                        f'arity {registered}, expected {expected_arity}')

                case = {
                    'api': 'stats.' + operation,
                    'args': {'d': family, 'p': list(parameters),
                             'x': [.42 if operation == 'ppf' else observation]},
                }
                observed = self.native.evaluate(case)
                if set(observed) != expected:
                    failures.append(f'{family} {operation}: returned {set(observed)}, '
                                    f'expected {expected}')
                for library in expected:
                    result = observed.get(library, {})
                    if ('error' in result or len(result.get('values', [])) != 1
                            or not all(math.isfinite(v) for v in result.get('values', []))):
                        failures.append(f'{family} {operation} / {library}: {result!r}')
        self.assertFalse(failures, '\n' + '\n'.join(failures))

    def test_stats_summary_vector_and_pair_adapters_are_reachable(self):
        """The summary ABI must be linked, registered and used by native.evaluate."""
        with (ROOT / 'api_inventory.csv').open(newline='') as stream:
            public = {row['api'] for row in csv.DictReader(stream)}
        self.assertLessEqual({'stats.' + name for name in STATS_SUMMARY_CASES}, public)

        # Importing only here keeps the independent scalar/distribution tests
        # available if a summary build is temporarily incomplete.
        from sci_testing.stats_summary import StatsSummaryComparisons

        summary = StatsSummaryComparisons(self.native)
        failures = []
        for name, (x, y, p, count) in sorted(STATS_SUMMARY_CASES.items()):
            expected = {GSL} if name in BOOST_SUMMARY_GAPS else {BOOST, GSL}
            for index, library in enumerate((BOOST, GSL)):
                registered = summary.library.stats_summary_supported(index, name.encode('ascii'))
                if bool(registered) != (library in expected):
                    failures.append(f'stats.{name} / {library}: registration {registered}, '
                                    f'expected {library in expected}')
            output_count = summary.library.stats_summary_output_count(name.encode('ascii'), len(x))
            if output_count != count:
                failures.append(f'stats.{name}: output count {output_count}, expected {count}')
            case = {'api': f'stats.{name}',
                    'args': {'x': list(x), 'y': list(y), 'p': list(p)}}
            observed = self.native.evaluate(case)
            if set(observed) != expected:
                failures.append(f'stats.{name}: routed {set(observed)}, expected {expected}')
            for library in expected:
                result = observed.get(library, {})
                if ('error' in result or len(result.get('values', ())) != count
                        or not all(math.isfinite(v) for v in result.get('values', ()))):
                    failures.append(f'stats.{name} / {library}: {result!r}')
        self.assertFalse(failures, '\n' + '\n'.join(failures))

    def test_saved_api_universe_is_exhaustively_audited(self):
        """Enumerate each saved name and retain a reason for each open mapping."""
        expected = {
            'special': (405, '4d88d0726eb7dc0d3de954eeae10ffdc4a895906e18be35b2307ea349a16d5dd'),
            'stats': (169, 'f334ac796c2fb74098b79fccf86d67d07cef5c893ca660e23d515af052c805e3'),
        }
        for module, (count, signature) in expected.items():
            saved = _saved_report_names(module)
            self.assertEqual(len(saved), count)
            self.assertEqual(hashlib.sha256('\n'.join(sorted(saved)).encode()).hexdigest(),
                             signature, f'{module}: saved API names changed; re-audit gaps')
            covered = _known_coverage_names(module)
            remaining = remaining_native_coverage(module)
            self.assertFalse(set(remaining) & {module + '.' + x for x in covered})
            self.assertEqual({module + '.' + name for name in saved},
                             {module + '.' + name for name in saved & covered} | set(remaining))
            for api, reason in remaining.items():
                with self.subTest(api=api):
                    self.assertGreater(len(reason), 35)
                    self.assertNotIn('unsupported by both libraries', reason.lower())
            unresolved = {api for api, reason in remaining.items()
                          if reason.startswith('unconfirmed:')}
            self.assertFalse(unresolved, f'{module}: new API needs a specific audit reason')

    def test_native_registration_catalog_has_no_unaudited_scalar_or_probability_names(self):
        """A newly registered native API must also receive an independent test."""
        scalar_source = (ROOT / 'native/comparisons.cpp').read_text()
        scalar_names = set(re.findall(r'^\s*X\((\w+),\s*\d+\s*,',
                                      scalar_source, re.MULTILINE))
        self.assertGreaterEqual(len(scalar_names), 120)
        expected_scalar = {name.rsplit('.', 1)[-1] for name in EQUIVALENTS}
        self.assertEqual(scalar_names - expected_scalar, set(),
                         'new scalar bridge registration requires a declared test case')

        # These worker registries have a separate stable ABI. Do not treat an
        # extra export as covered merely because the C++ source includes it.
        probability_source = (ROOT / 'native/probability_extra.cpp').read_text()
        table = re.search(r'constexpr Function functions\[\] = \{(.*?)\};',
                          probability_source, re.DOTALL)
        self.assertIsNotNone(table, 'probability registration table changed; review audit')
        native_names = set(re.findall(r'\{"([a-zA-Z0-9_]+)"\s*,', table.group(1)))
        self.assertEqual(native_names, set(PROBABILITY_EXTRA_CASES),
                         'probability registration/test expectation mismatch')

        from sci_testing.stats_extra import NAMES as stats_extra_names
        from sci_testing.stats_summary import SUPPORTED_NAMES as stats_summary_names
        from sci_testing.complex_extra import SUPPORTED_APIS as complex_extra_names
        from sci_testing.vector_special import SUPPORTED_NAMES as vector_special_names
        from sci_testing.spherical_extra import APIS as spherical_extra_names

        self.assertEqual(set(stats_extra_names), set(STATS_EXTRA_CASES))
        self.assertEqual(set(stats_summary_names), set(STATS_SUMMARY_CASES))
        self.assertEqual(set(complex_extra_names),
                         {'special.complex.' + name
                          for name in GSL_COMPLEX_CASES - {'digamma', 'spence'}})
        self.assertEqual(set(vector_special_names), set(VECTOR_SPECIAL_CASES))
        self.assertEqual(set(spherical_extra_names),
                         {'special.' + name for name in SPHERICAL_EXTRA_CASES})

    def test_special_extra_adapters_are_not_lost(self):
        failures = []
        bridge = getattr(self.native, 'special_extra', None)
        if bridge is None:
            self.fail('NativeComparisons has not integrated special_extra adapters')
        for name, (args, providers, count) in sorted(SPECIAL_EXTRA_CASES.items()):
            for index, provider in enumerate((BOOST, GSL)):
                arity = bridge.library.special_extra_supported(index, name.encode())
                if arity != (len(args) if provider in providers else 0):
                    failures.append(f'special.{name} / {provider}: unexpected arity {arity}')
            case = {'api': 'special.' + name,
                    'args': {f'p{i}': v for i, v in enumerate(args)}}
            observed = self.native.evaluate(case)
            if set(observed) != providers:
                failures.append(f'special.{name}: returned providers {set(observed)} vs {providers}')
            for provider in providers:
                value = observed.get(provider, {})
                if ('error' in value or len(value.get('values', ())) != count or
                        not all(math.isfinite(x) for x in value.get('values', ()))):
                    failures.append(f'special.{name} / {provider}: {value!r}')
        self.assertFalse(failures, '\n' + '\n'.join(failures))

    def test_probability_extra_api_registrations_and_dispatch(self):
        failures = []
        bridge = getattr(self.native, 'probability_extra', None)
        if bridge is None:
            self.fail('NativeComparisons has not integrated probability_extra adapters')
        # Some parameter-inverse functions have restricted attainable targets.
        # Registration is checked for every name; smoke calls cover ordinary
        # interior mathematical arguments, permitting explicit domain errors.
        for name, (args, providers) in sorted(PROBABILITY_EXTRA_CASES.items()):
            for index, provider in enumerate((BOOST, GSL)):
                arity = bridge.library.probability_extra_supported(index, name.encode())
                if arity != (len(args) if provider in providers else 0):
                    failures.append(f'special.{name} / {provider}: unexpected arity {arity}')
            case = {'api': 'special.' + name,
                    'args': {f'p{i}': v for i, v in enumerate(args)}}
            observed = self.native.evaluate(case)
            if set(observed) != providers:
                failures.append(f'special.{name}: providers {set(observed)} vs {providers}')
            for provider in providers:
                value = observed.get(provider, {})
                if 'values' in value and (len(value['values']) != 1 or
                                          not all(math.isfinite(x) for x in value['values'])):
                    failures.append(f'special.{name} / {provider}: {value!r}')
                if not value.get('values') and not value.get('error'):
                    failures.append(f'special.{name} / {provider}: no result or explicit error')
        self.assertFalse(failures, '\n' + '\n'.join(failures))

    def test_stats_extra_registration_and_dispatch(self):
        failures = []
        bridge = getattr(self.native, 'stats_extra', None)
        if bridge is None:
            self.fail('NativeComparisons has not integrated stats_extra adapters')
        for name, (x, y, p, output_count, providers) in sorted(STATS_EXTRA_CASES.items()):
            for index, provider in enumerate((BOOST, GSL)):
                support = bridge.library.stats_extra_supported(index, name.encode())
                if bool(support) != (provider in providers):
                    failures.append(f'stats.{name} / {provider}: registration {support}')
            if bridge.library.stats_extra_output_count(name.encode()) != output_count:
                failures.append(f'stats.{name}: output shape changed')
            case = {'api': 'stats.' + name,
                    'args': {'x': list(x), 'y': list(y), 'p': list(p)}}
            if name in ('multivariate_normal_pdf', 'multivariate_normal_logpdf'):
                case['args']['g'] = [list(row) for row in STATS_EXTRA_GAUSSIAN_COVARIANCE]
            observed = self.native.evaluate(case)
            if set(observed) != providers:
                failures.append(f'stats.{name}: providers {set(observed)} vs {providers}')
            for provider in providers:
                value = observed.get(provider, {})
                if ('error' in value or len(value.get('values', ())) != output_count or
                        not all(math.isfinite(v) for v in value.get('values', ()))):
                    failures.append(f'stats.{name} / {provider}: {value!r}')
        self.assertFalse(failures, '\n' + '\n'.join(failures))

    def test_gsl_complex_gamma_family_including_extra_bridge(self):
        failures = []
        for name in sorted(GSL_COMPLEX_CASES):
            case = {'api': 'special.complex.' + name, 'args': {'z': [.5, .25]}}
            observed = self.native.evaluate_complex(case)
            value = observed.get(GSL, {})
            if ('error' in value or len(value.get('values', ())) != 2 or
                    not all(math.isfinite(v) for v in value.get('values', ()))):
                failures.append(f'special.complex.{name}: {observed!r}')
        self.assertFalse(failures, '\n' + '\n'.join(failures))

    def test_public_vector_special_adapter_coverage(self):
        bridge = getattr(self.native, 'vector_special', None)
        self.assertIsNotNone(bridge, 'missing Boost vector_special bridge')
        failures = []
        for name, (values, output_count) in VECTOR_SPECIAL_CASES.items():
            encoded = name.encode()
            if bridge.library.vector_special_supported(0, encoded) != 1:
                failures.append(f'{name}: no Boost vector registration')
            if bridge.library.vector_special_supported(1, encoded):
                failures.append(f'{name}: unverified GSL vector registration')
            if bridge.library.vector_special_output_count(encoded, len(values)) != output_count:
                failures.append(f'{name}: output count mismatch')
            observed = self.native.evaluate({'api': 'special.' + name,
                                             'args': {'p0': list(values)}})
            if set(observed) != {BOOST}:
                failures.append(f'{name}: providers {set(observed)}')
            result = observed.get(BOOST, {})
            if ('error' in result or len(result.get('values', ())) != output_count or
                    not all(math.isfinite(v) for v in result.get('values', ()))):
                failures.append(f'{name}: {result!r}')
        self.assertFalse(failures, '\n' + '\n'.join(failures))

    def test_combinatorics_extra_registration_and_vector_shape(self):
        bridge = getattr(self.native, 'combinatorics_extra', None)
        self.assertIsNotNone(bridge, 'missing Bernoulli/harmonic bridge')
        failures = []
        for name, (args, providers, output_count) in COMBINATORICS_EXTRA_CASES.items():
            for index, library in enumerate((BOOST, GSL)):
                arity = bridge.library.combinatorics_extra_supported(index, name.encode())
                if arity != (len(args) if library in providers else 0):
                    failures.append(f'{name}/{library}: arity {arity}')
            case = {'api': 'special.' + name,
                    'args': {f'p{i}': value for i, value in enumerate(args)}}
            actual = self.native.evaluate(case)
            if set(actual) != providers:
                failures.append(f'{name}: providers {set(actual)} vs {providers}')
            for library in providers:
                result = actual.get(library, {})
                if ('error' in result or len(result.get('values', ())) != output_count or
                        not all(math.isfinite(v) for v in result.get('values', ()))):
                    failures.append(f'{name}/{library}: {result!r}')
        self.assertFalse(failures, '\n' + '\n'.join(failures))

    def test_spherical_complex_and_real_abi_coverage(self):
        bridge = getattr(self.native, 'spherical_extra', None)
        self.assertIsNotNone(bridge, 'missing spherical special bridge')
        failures = []
        for name, (args, providers, output_count) in SPHERICAL_EXTRA_CASES.items():
            for index, library in enumerate((BOOST, GSL)):
                arity = bridge.library.spherical_extra_supported(index, name.encode())
                if arity != (len(args) if library in providers else 0):
                    failures.append(f'{name}/{library}: arity {arity}')
            if bridge.library.spherical_extra_output_count(name.encode()) != output_count:
                failures.append(f'{name}: wrong output count')
            case = {'api': 'special.' + name,
                    'args': {f'p{i}': value for i, value in enumerate(args)}}
            actual = self.native.evaluate(case)
            if set(actual) != providers:
                failures.append(f'{name}: providers {set(actual)} vs {providers}')
            for library in providers:
                result = actual.get(library, {})
                if ('error' in result or len(result.get('values', ())) != output_count or
                        not all(math.isfinite(v) for v in result.get('values', ()))):
                    failures.append(f'{name}/{library}: {result!r}')
        self.assertFalse(failures, '\n' + '\n'.join(failures))


if __name__ == '__main__':
    unittest.main()
