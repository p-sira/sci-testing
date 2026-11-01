"""Independent, typed comparisons for nine new Sci 0.3 special exports.

References use mpmath mathematical definitions, not Sci or its private
numerical backends. ``reference`` respects the caller's mpmath precision so
the report can check stability at 100 and 160 decimal digits. Case arguments
are binary64 JSON numbers; mpmath converts those inputs exactly.

The regular profile supplies ``count`` endpoint-inclusive linspace cases and
``count`` independently seeded samples per API. Small explicit edge profiles
cover limits and identities that a random sweep would almost never exercise.
Rejected domains and nonfinite limits belong to the ordinary Sci API tests:
the comparison report expects a finite mathematical reference for each case.
"""

import hashlib
import math
import random

import mpmath as mp


APIS = (
    'special.bessel.iv_ratio',
    'special.bessel.iv_ratio_c',
    'special.bessel.iv_ratioinv',
    'special.combinatorics.normalized_gen_harmonic',
    'special.complex.hurwitz_zeta',
    'special.elementary.log1pmx',
    'special.gamma.digammainv',
    'special.gamma.gamma_ratio',
    'special.integrals.scaled_exp1',
)


def _arguments(api, i, u, v):
    """Choose valid, varied domains without relying on computed Sci outputs."""
    if api.startswith('special.bessel.'):
        orders = (0.5, 1., 1.5, 2., 3., 8., 32., 128.)
        order = orders[i % len(orders)]
        if api.endswith('iv_ratioinv'):
            # The inverse becomes badly conditioned as r approaches one.
            return dict(v=order, r=0.0001 + 0.9699 * u)
        return dict(v=order, x=10. ** (-8. + 11. * u))
    if api.endswith('normalized_gen_harmonic'):
        n = 1 + 3 * (i % 32)
        j = 1 + min(n - 1, int(n * u))
        k = j + min(n - j, int((n - j + 1) * v))
        exponent = (0., 1., -1., -25. + 50. * u, -8. + 16. * v)[i % 5]
        return dict(j=j, k=k, n=n, a=exponent)
    if api.endswith('hurwitz_zeta'):
        # Stay in the reliable, non-pole Euler–Maclaurin range of Sci 0.3.
        return dict(s=[2. + 8. * u, 8. * (v - 0.5)],
                    q=[1. + 4. * v, 2. * (u - 0.5)])
    if api.endswith('log1pmx'):
        mode = i % 6
        if mode == 0:
            x = -1. + 10. ** (-1. - 14. * u)
        elif mode == 1:
            x = -10. ** (-1. - 140. * u)
        elif mode == 2:
            x = 10. ** (-1. - 140. * u)
        elif mode == 3:
            x = -0.9 + 1.8 * u
        elif mode == 4:
            x = 10. ** (-2. + 8. * u)
        else:
            x = 10. ** (2. + 10. * u)
        return dict(x=x)
    if api.endswith('digammainv'):
        return dict(y=(-50. + 100. * u) if i % 5 == 0 else -10. + 20. * u)
    if api.endswith('gamma_ratio'):
        mode = i % 4
        if mode == 0:
            return dict(a=0.05 + 60. * u, b=0.05 + 60. * v)
        if mode == 1:
            a = 1. + 149. * u
            return dict(a=a, b=a + 0.5 * (2. * v - 1.))
        if mode == 2:
            a = 2. + 100. * u
            return dict(a=a, b=a + (20. * v - 10.)) if a + 20. * v - 10. > 0 else dict(a=a, b=a + 0.5)
        a = 500. + 5000. * u
        return dict(a=a, b=a + (v - 0.5))
    if api.endswith('scaled_exp1'):
        return dict(x=10. ** (-12. + 18. * u))
    raise ValueError(f'unsupported Sci 0.3 special API: {api}')


_EDGES = {
    'special.bessel.iv_ratio': [
        dict(v=v, x=x) for v, x in ((0.5, 0.), (0.5, 1e-12), (0.5, 20.),
                                  (1., 1.), (1., 1e-200), (1., 1e8),
                                  (2., 1e4), (128., 0.125))],
    'special.bessel.iv_ratio_c': [
        dict(v=v, x=x) for v, x in ((0.5, 0.), (0.5, 1e-12), (0.5, 20.),
                                  (0.5, 200.), (1., 1.), (1., 1e8),
                                  (2., 1e5), (128., 0.125))],
    'special.bessel.iv_ratioinv': [
        dict(v=v, r=r) for v, r in ((0.5, 0.), (0.5, 0.5), (0.5, 0.99999),
                                  (1., 0.), (1., 0.5), (1., 1e-12),
                                  (2., 0.999), (128., 0.3))],
    'special.combinatorics.normalized_gen_harmonic': [
        dict(j=j, k=k, n=n, a=a) for j, k, n, a in
        ((1, 1, 1, 0.), (1, 4, 4, -200.), (2, 3, 4, 0.),
         (1, 10, 10, 1.), (2, 10, 10, -8.), (100, 100, 100, 300.),
         (999, 1000, 1000, -200.), (1, 10, 1000, 2.))],
    'special.complex.hurwitz_zeta': [
        dict(s=list(s), q=list(q)) for s, q in
        (((0., 0.), (1., 0.)), ((0., 0.), (1.5, 0.5)),
         ((2., 0.), (1., 0.)), ((2., 0.), (0.5, 0.)),
         ((2., 1.), (1., 0.)), ((3., 2.), (2., -1.)),
         ((8., 0.), (8., 0.)), ((12., 4.), (3., 1.)))],
    'special.elementary.log1pmx': [
        dict(x=x) for x in (0., -0.0, 1e-150, -1e-150, 1e-8, -1e-8,
                            -0.5, math.nextafter(-1., 0.), 0.5, 1e8)],
    'special.gamma.digammainv': [
        dict(y=y) for y in (-100., -40., -10., -2.22, -0.5772156649015329,
                            0., 1., 8., 36., 40.)],
    'special.gamma.gamma_ratio': [
        dict(a=a, b=b) for a, b in ((1., 1.), (201., 200.), (200., 201.),
                                     (0.5, 0.5), (0.5, 1.), (1., 0.5),
                                     (1e-10, 2e-10), (5000., 5000.25),
                                     (1000., 999.5), (2., 160.))],
    'special.integrals.scaled_exp1': [
        dict(x=x) for x in (0., 1e-300, 1e-50, 1e-12, 0.5, 1.,
                            2., 30., 1000., 1e10)],
}


def _case(api, name, args, group):
    return dict(id=f'{api}/{name}', api=api, group=group,
                args=dict(special03_case=True, **args))


def cases(count, seed, apis=None):
    """Yield reproducible 2*count sweeps plus a few distinctive limits per API."""
    if count < 1:
        raise ValueError('count must be positive')
    selected = set(APIS if apis is None else apis)
    if selected.difference(APIS):
        raise ValueError(f'unknown Sci 0.3 special APIs: {sorted(selected.difference(APIS))}')
    for api in APIS:
        if api not in selected:
            continue
        for i, args in enumerate(_EDGES[api]):
            yield _case(api, f'edge-{i:03d}', args, 'independent mathematical edge / identity')
        digest = hashlib.sha256(f'{seed}:{api}'.encode()).digest()
        rng = random.Random(int.from_bytes(digest, 'big'))
        for method in ('linspace', 'random'):
            for i in range(count):
                if method == 'linspace':
                    u = i / (count - 1) if count > 1 else 0.5
                    v = ((i + 37 * (count // 2)) % count) / (count - 1) if count > 1 else 0.5
                else:
                    u, v = rng.random(), rng.random()
                yield _case(api, f'{method}-{i:04d}', _arguments(api, i, u, v),
                            f'{method}, {count} cases; independent mathematical reference')


def _bessel_ratio(v, x):
    if not x:
        return mp.mpf(0)
    if v == mp.mpf('0.5'):
        return mp.tanh(x)
    return mp.besseli(v, x) / mp.besseli(v - 1, x)


def _bessel_inverse(v, r):
    if not r:
        return mp.mpf(0)
    if v == mp.mpf('0.5'):
        return mp.atanh(r)
    # A positive seed, followed by Newton steps on the *independent* Bessel
    # definition. The Riccati identity gives an analytic ratio derivative.
    with mp.extradps(12):
        x = 2 * v * r / (1 - r * r)
        for _ in range(90):
            value = _bessel_ratio(v, x)
            if abs(value - r) <= 100 * mp.eps:
                return +x
            derivative = 1 - value * value - (2 * v - 1) * value / x
            delta = (value - r) / derivative
            next_x = x - delta
            if next_x <= 0 or not mp.isfinite(next_x):
                next_x = x / 2
            if abs(next_x - x) <= 4 * mp.eps * max(1, abs(x)):
                return +next_x
            x = next_x
        raise ArithmeticError(f'independent Bessel inverse did not converge: v={v}, r={r}')


def reference(case):
    """Return an mpmath vector at the active caller-selected precision.

    The complex result is [real, imaginary], matching Rust's flattened format.
    """
    api, a = case['api'], case['args']
    if api not in APIS or not a.get('special03_case'):
        raise ValueError('not a Sci 0.3 special comparison case')
    number = mp.mpf  # Exact binary64-to-mpmath conversion.
    if api.startswith('special.bessel.'):
        v = number(a['v'])
        if api.endswith('iv_ratioinv'):
            return [_bessel_inverse(v, number(a['r']))]
        x = number(a['x'])
        if api.endswith('iv_ratio_c'):
            if v == mp.mpf('0.5'):
                return [2 / (1 + mp.exp(2 * x))]
            return [1 - _bessel_ratio(v, x)]
        return [_bessel_ratio(v, x)]
    if api.endswith('normalized_gen_harmonic'):
        j, k, n, exponent = a['j'], a['k'], a['n'], number(a['a'])
        # Scale by the largest weight; neither sum can overflow even at a=-200.
        anchor = n if exponent < 0 else 1
        terms = [(number(i) / anchor) ** (-exponent) for i in range(1, n + 1)]
        return [mp.fsum(terms[j - 1:k]) / mp.fsum(terms)]
    if api.endswith('hurwitz_zeta'):
        s = mp.mpc(*map(number, a['s']))
        q = mp.mpc(*map(number, a['q']))
        z = mp.zeta(s, q)
        return [z.real, z.imag]
    if api.endswith('log1pmx'):
        x = number(a['x'])
        # Subtracting log1p(x)-x at finite precision loses small remainders.
        # The mathematical power series preserves their relative accuracy.
        if abs(x) < mp.mpf('0.125'):
            power = -x * x
            total = power / 2
            for k in range(3, 512):
                power *= -x
                term = power / k
                next_total = total + term
                if next_total == total:
                    break
                total = next_total
            return [total]
        return [mp.log1p(x) - x]
    if api.endswith('digammainv'):
        y = number(a['y'])
        with mp.extradps(12):
            x = mp.exp(y) + mp.mpf('0.5') if y >= -mp.mpf('2.22') else -1 / (y + mp.euler)
            for _ in range(90):
                delta = (mp.digamma(x) - y) / mp.polygamma(1, x)
                next_x = x - delta
                if next_x <= 0:
                    next_x = x / 2
                if abs(next_x - x) <= 4 * mp.eps * max(1, abs(x)):
                    return [+next_x]
                x = next_x
            raise ArithmeticError(f'independent inverse digamma did not converge: {y}')
    if api.endswith('gamma_ratio'):
        return [mp.exp(mp.loggamma(number(a['a'])) - mp.loggamma(number(a['b'])))]
    if api.endswith('scaled_exp1'):
        x = number(a['x'])
        return [x * mp.exp(x) * mp.e1(x)] if x else [mp.mpf(0)]
    raise ValueError(f'unknown Sci 0.3 special API: {api}')
