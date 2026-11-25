"""Independent high-precision probability regression fixtures.

Arguments use existing runner signatures; all floats are the exact parsed
binary64 inputs. References return unrounded mpmath numbers at caller precision.
"""
import math
import sys

import mpmath as mp


_GROUP = 'high-precision tail regression'


def is_case(case):
    return case.get('id', '').startswith(case.get('api', '') + '/tail-regression-')


def cases(count=None, seed=None):
    """Small deterministic fixture set (independent of sweep count and seed)."""
    def positional(api, rows):
        for i, row in enumerate(rows):
            yield dict(id=f'{api}/tail-regression-{i:03d}', api=api, group=_GROUP,
                       args={f'p{j}': value for j, value in enumerate(row)})

    yield from positional('special.chndtr', [(2., 5., 400.), (10., 5., 400.)])
    yield from positional('special.ncfdtr', [(5., 30., 400., .5), (5., 30., 400., 2.)])
    yield from positional('special.nctdtr', [
        (30., 20., .5), (30., 20., 2.), (10000., 30., 2.),
        (30., 5., -2.), (30., 5., -4.), (30., 10., -2.), (30., 20., -2.),
        (10000., 30., -2.), (.01, 5., -1e100), (3., -5., -1000.),
        (.01, 40., -2.), (30., .5, -5.141358682658483),
        (30., .5, -5.148029306436635),
    ])
    normal = [-37., math.nextafter(-37., -math.inf), -37.6, -37.6996996996997,
              -37.7, -38., -38.4, -38.48, -38.49, -38.5, -38.51]
    yield from positional('special.ndtr', [(x,) for x in normal])
    yield from positional('special.erfc', [(x,) for x in
                          (26.5, 26.7, 27., 27.2, 27.22, 27.23)])
    for i, x in enumerate(normal):
        yield dict(id=f'stats.normal_cdf/tail-regression-{i:03d}',
                   api='stats.normal_cdf', group=_GROUP, args={'x': x})
    # Equal, separated groups give a closed-form rank-sum statistic and
    # permit one-rounding tests at the minimum-subnormal two-sided p-value.
    for i, n in enumerate((987, 988, 989)):
        yield dict(id=f'stats.ranksums/tail-regression-{i:03d}', api='stats.ranksums',
                   group=_GROUP, args=dict(p=[0.], x=list(range(n)),
                                          y=list(range(n, 2*n)), g=[], d=''))
    for api in ('stats.cdf', 'stats.sf', 'stats.pdf', 'stats.logpdf'):
        for i, (location, scale, x) in enumerate([
            (10., .01, 10.01483966966967),
            (-sys.float_info.max, sys.float_info.max, sys.float_info.max),
            (0., 1e-320, 10.),
            (-sys.float_info.max, 1e-320, sys.float_info.max),
        ]):
            yield dict(id=f'{api}/tail-regression-{i:03d}', api=api,
                       group=_GROUP, args=dict(p=[location, scale], x=[x],
                                              y=[], g=[], d='Levy'))


def _poisson_mixture(mean, conditional):
    # Sum from j=0, with a rigorous remainder bound from decreasing
    # conditional CDFs and a geometric upper bound on future Poisson weights.
    weight = mp.exp(-mean)
    total = mp.mpf(0)
    for j in range(4096):
        cdf = conditional(j)
        total += weight * cdf
        ratio = mean / (j + 1)
        if ratio < 1 and weight * cdf * ratio / (1 - ratio) <= mp.eps * total:
            return total
        weight *= ratio
    raise ArithmeticError('Poisson mixture did not converge')


def _normal_cdf(x):
    return mp.erfc(-x/mp.sqrt(2))/2


def _noncentral_t(df, nc, t):
    # T=(Z+nc)/sqrt(V/df), with independent standard normal Z and
    # chi-square V. Integrate the conditional CDF in log(v), v=sqrt(V/df).
    # Mode scaling preserves tiny probabilities and avoids absolute-error
    # quadrature on their original scale. This does not use a beta series.
    def log_density(x):
        v = mp.exp(x)
        return df*x - df*v*v/2 + mp.log(_normal_cdf(t*v - nc))

    def derivative(x):
        v = mp.exp(x)
        z = t*v - nc
        mills = mp.exp(-z*z/2) / (mp.sqrt(2*mp.pi)*_normal_cdf(z))
        return df - df*v*v + t*v*mills

    if df <= 0 or not all(map(mp.isfinite, (df, nc, t))):
        raise ValueError('expected finite noncentral t parameters and positive df')

    origin = -mp.log(max(mp.mpf(1), abs(t)))

    def bracket(sign):
        point = origin
        step = mp.mpf(1)
        for _ in range(32):
            if sign*derivative(point) <= 0:
                return point
            point += sign*step
            step *= 2
        raise ArithmeticError('Conditional normal integration mode not bracketed')

    lo, hi = bracket(-1), bracket(1)
    for _ in range(100):
        mid = (lo + hi)/2
        if derivative(mid) > 0:
            lo = mid
        else:
            hi = mid
    mode = (lo + hi)/2
    anchor = log_density(mode)
    width = 1/mp.sqrt(df)
    # The extra 25 decimal orders accommodate integration over wide log-v
    # ranges at small df. The report independently repeats the calculation at
    # 100 and 160 digits to detect any residual truncation or quadrature error.
    threshold = -(mp.mp.dps + 25)*mp.log(10)
    bounds = [mode]
    for sign in (-1, 1):
        step = width
        for _ in range(32):
            point = mode + sign*step
            bounds.append(point)
            if log_density(point) - anchor < threshold:
                break
            step *= 2
        else:
            raise ArithmeticError('Conditional normal integration range exceeded')
    bounds.sort()
    area = mp.quad(lambda x: mp.exp(log_density(x)-anchor), bounds)
    return 2*mp.exp((df/2)*mp.log(df/2)-mp.loggamma(df/2)+anchor)*area


def reference(case):
    """Evaluate a new regression independently at active mpmath precision."""
    if not is_case(case):
        raise ValueError('not a probability tail regression case')
    api, a = case['api'], case['args']
    if api == 'stats.normal_cdf':
        return [_normal_cdf(mp.mpf(a['x']))]
    if api == 'stats.ranksums':
        n, m = map(mp.mpf, (len(a['x']), len(a['y'])))
        if a['p'] != [0.] or max(a['x']) >= min(a['y']):
            raise ValueError('expected disjoint rank-sum samples and two-sided test')
        z = -n*m/2 / mp.sqrt(n*m*(n+m+1)/12)
        return [z, mp.erfc(abs(z)/mp.sqrt(2))]
    if api in ('stats.cdf', 'stats.sf', 'stats.pdf', 'stats.logpdf'):
        location, scale = map(mp.mpf, a['p'])
        delta = mp.mpf(a['x'][0]) - location
        if a['d'] != 'Levy' or delta <= 0 or scale <= 0:
            raise ValueError('expected interior Levy distribution case')
        r = scale/(2*delta)
        if api == 'stats.cdf':
            return [mp.erfc(mp.sqrt(r))]
        if api == 'stats.sf':
            return [mp.erf(mp.sqrt(r))]
        log_density = (mp.log(scale)-mp.log(2*mp.pi)-3*mp.log(delta))/2-r
        return [log_density if api == 'stats.logpdf' else mp.exp(log_density)]
    v = [mp.mpf(a[f'p{i}']) for i in range(len(a))]
    if api == 'special.ndtr':
        return [_normal_cdf(v[0])]
    if api == 'special.erfc':
        return [mp.erfc(v[0])]
    if api == 'special.chndtr':
        x, df, nc = v
        return [_poisson_mixture(nc/2, lambda j: mp.gammainc(df/2+j, 0, x/2,
                                                             regularized=True))]
    if api == 'special.ncfdtr':
        dfn, dfd, nc, f = v
        x = dfn*f/(dfd+dfn*f)
        return [_poisson_mixture(nc/2, lambda j: mp.betainc(dfn/2+j, dfd/2, 0, x,
                                                           regularized=True))]
    if api == 'special.nctdtr':
        return [_noncentral_t(*v)]
    raise ValueError(api)
