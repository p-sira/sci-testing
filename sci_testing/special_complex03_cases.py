"""Sci 0.3 public complex-special overloads, with independent complex oracles.

All 32 names refer to ``sci::special::complex``. They are distinct from the
real facade ``sci::special::{gamma,bessel,elementary,...}``. Every complex
argument is a JSON `[re, im]` binary64 pair; the returned vector interleaves
real and imaginary components for each complex result. Airy outputs are
`Ai, Ai', Bi, Bi'`; Fresnel outputs are `S, C`.

mpmath provides high-precision independent definitions at the caller's dps.
An infinitesimal imaginary perturbation selects the correct one-sided limit
at a negative-real branch cut without requiring SciPy reference values.
Exceptions are documented exact-axis conventions: ``expi(x ± 0j)`` returns
real principal-value Ei(x) for x<0; Lambert W branch -1 returns the real
branch on [-1/e,0) for either sign of imaginary zero. These exact-axis
conventions differ from selected complex analytic one-sided limits.

SciPy 1.18.1 evaluates negative-imaginary-zero Hankel, K and Y inputs on the
upper bank; unlike the mathematical lower-bank limit, this is a known public
SciPy comparison discrepancy. Never use SciPy values to construct the oracle.
The sampling avoids singularities and known bounded-work limits in Sci 0.3.
"""

import hashlib
import math
import random
import warnings

import mpmath as mp


_NAMES = (
    'airy', 'airye', 'dawsn', 'digamma', 'erf', 'erfc', 'erfcx', 'erfi',
    'exp1', 'expi', 'fresnel', 'gamma', 'hankel1', 'hankel1e', 'hankel2',
    'hankel2e', 'iv', 'ive', 'jv', 'jve', 'kv', 'kve', 'lambertw',
    'loggamma', 'rgamma', 'spence', 'wofz', 'xlog1py', 'xlogy', 'yv',
    'yve', 'zeta',
)
APIS = tuple('special.complex.' + name for name in _NAMES)
_BESSEL = frozenset(('hankel1', 'hankel1e', 'hankel2', 'hankel2e',
                     'iv', 'ive', 'jv', 'jve', 'kv', 'kve', 'yv', 'yve'))
_BESSEL_ORDERS = (0., 0.5, -0.5, 1., -1., 2.5, 4.)
_CUT = frozenset(_BESSEL | {'airye', 'exp1', 'expi', 'loggamma', 'spence',
                            'lambertw', 'xlogy', 'xlog1py'})
# Diagnostic only: these public SciPy 1.18.1 functions collapse -0j to +0j
# on the negative real axis. The independent reference retains the lower bank.
_SCIPY_IGNORES_LOWER_ZERO = frozenset(('hankel1', 'hankel1e',
                                      'hankel2', 'hankel2e',
                                      'kv', 'kve', 'yv', 'yve'))
_ENTIRE = frozenset(('airy', 'airye', 'dawsn', 'erf', 'erfc', 'erfcx',
                      'erfi', 'fresnel', 'wofz', 'rgamma'))


def scipy_comparison(case):
    """Compare against public SciPy complex overloads independently of mpmath.

    Use the input's signed imaginary zero unchanged. Some SciPy functions select
    the upper bank at -0j whereas Sci and the stored mpmath oracle select the
    lower bank; those differences must be visible in the comparison table.
    """
    from scipy import special

    name=case['api'].split('.')[-1]
    if name not in _NAMES:
        raise ValueError(f'No registered complex SciPy comparison: {name}')
    a=case['args']
    try:
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter('always')
            if name in ('xlogy','xlog1py'):
                result=getattr(special,name)(complex(*a['x']),complex(*a['y']))
            elif name=='lambertw':
                result=special.lambertw(complex(*a['z']),k=int(a['k']))
            elif name in _BESSEL:
                result=getattr(special,name)(a['v'],complex(*a['z']))
            else:
                result=getattr(special,name)(complex(*a['z']))
            parts=result if isinstance(result,(list,tuple)) else (result,)
            values=[v for part in parts for v in (float(complex(part).real),
                                                  float(complex(part).imag))]
            if not all(math.isfinite(v) for v in values):
                raise ValueError('nonfinite output')
            row={'values':values}
            if captured:row['warnings']=[str(w.message) for w in captured]
            return {'SciPy':row}
    except (ValueError,TypeError,ArithmeticError) as error:
        return {'SciPy':{'error':f'{type(error).__name__}: {error}'}}


def _case(api, name, args, description):
    return {'id': f'{api}/{name}', 'api': api, 'group': description,
            'args': {'special_complex03_case': True, **args}}


def _edges(name):
    """Distinct identities, conjugate quadrants, cut banks, and small limits."""
    if name == 'zeta':
        return [dict(z=list(z)) for z in (
            (0., 0.), (2., 0.), (4., 0.), (0.5, 14.),
            (2., 0.5), (3., -3.), (5., 7.), (12., -1.))]
    if name == 'lambertw':
        return [dict(z=list(z), k=k) for z, k in (
            ((0., 0.), 0), ((1., 0.), 0), ((-0.3, 0.), 0),
            ((-0.3, -0.), 0), ((-0.3, 0.), -1),
            ((-0.3, -0.), -1), ((1., 0.5), 1),
            ((0.25, -1.), -1), ((2., -0.25), 0))]
    if name in ('xlogy', 'xlog1py'):
        ys = ((0., 0.), (1., 0.), (-0.5, 0.),
              (-2., 1e-9), (-2., -1e-9), (1e-10, 1e-10))
        examples = [dict(x=[0., 0.], y=[0., 0.]),
                    dict(x=[0., -0.], y=[-1., 0.] if name == 'xlog1py' else [0., 0.])]
        examples += [dict(x=[x[0], x[1]], y=list(y)) for x, y in zip(
            ((1., 0.), (-2., 0.5), (0.25, -0.7),
             (-1.25, 2.), (1e-10, 0.), (0.5, 0.5)), ys[1:])]
        examples.append(dict(x=[1., 0.], y=[-2., -0.]))
        return examples
    if name in _BESSEL:
        return [dict(v=_BESSEL_ORDERS[i % len(_BESSEL_ORDERS)], z=list(z))
                for i, z in enumerate(((1., 0.), (0.5, 0.5),
                                       (2., -0.75), (0.2, 3.),
                                       (-1.5, 0.25), (-1.5, -0.25),
                                       (-1.5, 0.), (-1.5, -0.),
                                       (0.0001, -0.0002)))]
    values = [(.25, .5), (1., 0.), (1., -0.),
              (-.75, .25), (-.75, -.25), (2., 1.),
              (.125, -2.), (-1.25, 1e-9), (-1.25, -1e-9)]
    if name in _ENTIRE:
        values += [(0., 0.), (0., -0.)]
    if name in _CUT:
        values += [(-1.25, 0.), (-1.25, -0.)]
    return [dict(z=list(z)) for z in values]


def _arguments(name, i, u, t, w):
    # U/V parameter sweeps stay within the binary64 AMOS and analytic
    # continuation working ranges, with nonzero imaginary cuts except edges.
    re = -1.8 + 5.0 * u
    im = -2.5 + 5.0 * t
    # Independent small and moderate scales, including all four quadrants.
    if i % 7 == 0:
        re, im = re * .1, im * .1
    elif i % 11 == 0:
        re, im = re * 1.5, im * 1.25
    if name == 'zeta':
        return dict(z=[0.2 + 13. * u, -10. + 20. * t])
    if name == 'lambertw':
        return dict(z=[re, im], k=(-1, 0, 1)[i % 3])
    if name in ('xlogy', 'xlog1py'):
        x = [.3 + 2.5 * u, -.8 + 1.6 * w]
        # Include small complex y to exercise xlog1py cancellation.
        y = [re, im] if i % 3 else [re * 1e-9, im * 1e-9]
        return dict(x=x, y=y)
    result = dict(z=[re, im])
    if name in _BESSEL:
        result['v'] = _BESSEL_ORDERS[i % len(_BESSEL_ORDERS)]
    return result


def cases(count, seed, apis=None):
    """Yield edge cases plus ``count`` seeded and ``count`` linspace samples/API."""
    if count < 1:
        raise ValueError('count must be positive')
    selected = set(APIS if apis is None else apis)
    unknown = selected.difference(APIS)
    if unknown:
        raise ValueError(f'unrecognized public Sci complex APIs: {sorted(unknown)}')
    for api in APIS:
        if api not in selected:
            continue
        name = api.rsplit('.', 1)[1]
        for i, args in enumerate(_edges(name)):
            yield _case(api, f'edge-{i:03d}', args,
                        'public complex special: edge / branch / identity')
        digest = hashlib.sha256(f'{seed}:{api}'.encode()).digest()
        rng = random.Random(int.from_bytes(digest, 'big'))
        for mode in ('linspace', 'random'):
            for i in range(count):
                if mode == 'linspace':
                    u = i / (count - 1) if count > 1 else 0.5
                    t = ((i * 37 + count // 2) % count) / (count - 1) if count > 1 else .5
                    w = ((i * 43 + count // 3) % count) / (count - 1) if count > 1 else .5
                else:
                    u, t, w = rng.random(), rng.random(), rng.random()
                yield _case(api, f'{mode}-{i:04d}',
                            _arguments(name, i, u, t, w),
                            f'{mode}: independent complex special reference')


def _z(value):
    return mp.mpc(mp.mpf(value[0]), mp.mpf(value[1]))


def _complexes(values):
    return [part for z in values for part in (mp.re(z), mp.im(z))]


def reference(case):
    """Independent, precision-aware mathematical oracle for each public overload."""
    api, a = case['api'], case['args']
    if api not in APIS or a.get('special_complex03_case') is not True:
        raise ValueError('not a public Sci 0.3 complex special case')
    name = api.rsplit('.', 1)[1]
    branch_key = 'y' if name in ('xlogy', 'xlog1py') else 'z'
    branch_value = a.get(branch_key)
    on_negative_cut = (name in _CUT and branch_value is not None
                       and branch_value[0] < 0. and branch_value[1] == 0.)
    if on_negative_cut:
        real = mp.mpf(branch_value[0])
        if name == 'expi':
            # The explicitly real-axis PV convention is distinct from the
            # analytic complex limit Ei(x ± i0) = Ei(x) ± i*pi (x<0).
            return [mp.re(mp.ei(real)), mp.mpf(0)]
        if name == 'lambertw' and a['k'] == -1 and -1 / mp.e <= real < 0:
            # Sci's exact negative-real branch -1 returns real W_-1(x) even
            # for -0j; the lower complex limit is on another W sheet.
            return _complexes([mp.lambertw(real, -1)])
        # mpmath mpc() erases signed zero. An arbitrarily small nonzero
        # imaginary component identifies the branch bank; mpmath's exponents
        # are unbounded, so this can be smaller than binary64 subnormals even
        # with a 100/160-digit mantissa. The limiting binary64 value, including
        # components mathematically equal to zero, then stays dps-independent.
        side = -1 if math.copysign(1., branch_value[1]) < 0 else 1
        tiny = side * mp.power(10, -max(400, mp.mp.dps * 4))
        a = {**a, branch_key: [branch_value[0], tiny]}
    z = _z(a['z']) if 'z' in a else None
    if name in _BESSEL:
        v = mp.mpf(a['v'])
        base = name.rstrip('e') if name.endswith('e') else name
        if base == 'hankel1':
            result = mp.hankel1(v, z)
        elif base == 'hankel2':
            result = mp.hankel2(v, z)
        else:
            result = {'iv': mp.besseli, 'jv': mp.besselj,
                      'kv': mp.besselk, 'yv': mp.bessely}[base](v, z)
        scale = {
            'hankel1e': lambda: mp.exp(-mp.j * z),
            'hankel2e': lambda: mp.exp(mp.j * z),
            'ive': lambda: mp.exp(-abs(mp.re(z))),
            'jve': lambda: mp.exp(-abs(mp.im(z))),
            'kve': lambda: mp.exp(z),
            'yve': lambda: mp.exp(-abs(mp.im(z))),
        }
        if name in scale:
            result *= scale[name]()
        if on_negative_cut and name in ('jv', 'jve', 'iv', 'ive') and v == int(v):
            # Integer-order J and I (including these real scale factors) are
            # real on the real axis. Numerical continuation through ±i eps
            # can leave a tiny, precision-dependent imaginary rounding tail.
            result = mp.mpc(mp.re(result), mp.mpf(0))
        return _complexes([result])
    if name in ('airy', 'airye'):
        result = [mp.airyai(z), mp.airyai(z, 1), mp.airybi(z), mp.airybi(z, 1)]
        if name == 'airye':
            phase = mp.mpf(2) / 3 * z * mp.sqrt(z)
            result = [result[0] * mp.exp(phase), result[1] * mp.exp(phase),
                      result[2] * mp.exp(-abs(mp.re(phase))),
                      result[3] * mp.exp(-abs(mp.re(phase)))]
            if on_negative_cut:
                # Bi and Bi' are real on the negative real axis. Their scale
                # exp(-abs(Re phase)) is also real, so perturbation tails must
                # vanish in the exact boundary limit.
                result[2:] = [mp.mpc(mp.re(value), mp.mpf(0))
                              for value in result[2:]]
        return _complexes(result)
    if name == 'fresnel':
        return _complexes([mp.fresnels(z), mp.fresnelc(z)])
    if name == 'lambertw':
        return _complexes([mp.lambertw(z, a['k'])])
    if name == 'zeta':
        return _complexes([mp.zeta(z)])
    if name == 'spence':
        return _complexes([mp.polylog(2, 1 - z)])
    if name in ('xlogy', 'xlog1py'):
        x, y = _z(a['x']), _z(a['y'])
        result = mp.mpf(0) if x == 0 else x * (mp.log1p(y) if name == 'xlog1py' else mp.log(y))
        return _complexes([result])
    special = {
        'dawsn': lambda: mp.sqrt(mp.pi) / 2 * mp.exp(-z*z) * mp.erfi(z),
        'digamma': lambda: mp.digamma(z),
        'erf': lambda: mp.erf(z),
        'erfc': lambda: mp.erfc(z),
        'erfcx': lambda: mp.exp(z*z) * mp.erfc(z),
        'erfi': lambda: mp.erfi(z),
        'exp1': lambda: mp.e1(z),
        'expi': lambda: mp.ei(z),
        'gamma': lambda: mp.gamma(z),
        'loggamma': lambda: mp.loggamma(z),
        'rgamma': lambda: mp.rgamma(z),
        'wofz': lambda: mp.exp(-z*z) * mp.erfc(-mp.j * z),
    }
    return _complexes([special[name]()])
