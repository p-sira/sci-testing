"""Independent complex oracles and NumPy comparisons for Sci's complex facade.

SciPy does not export basic complex arithmetic/transcendental functions; its
NumPy dependency provides the matching public binary64 operations. The oracle
uses mpmath with the caller's precision and preserves signed-zero cut banks.
"""

import math
import warnings

import mpmath as mp
import numpy as np


def _mpz(pair):
    return mp.mpc(mp.mpf(pair[0]), mp.mpf(pair[1]))


def reference(case):
    """Evaluate a complex API without using binary64 complex arithmetic."""
    name = case['api'].split('.')[-1]
    a = case['args']
    if 'x' in a and 'y' in a:
        z = _mpz((a['x'], a['y']))
    elif name == 'from_polar':
        radius, theta = mp.mpf(a['p0']), mp.mpf(a['p1'])
        value = radius * mp.exp(mp.j * theta)
        return [mp.re(value), mp.im(value)]
    else:
        z = _mpz(a['p0'])
    w = a.get('p1')
    other = _mpz(w) if isinstance(w, list) else mp.mpf(w) if isinstance(w, float) else w

    if name == 'real': return [mp.re(z)]
    if name == 'imag': return [mp.im(z)]
    if name == 'arg':
        if mp.im(z) == 0 and mp.re(z) < 0:
            return [-mp.pi if math.copysign(1., a['p0'][1]) < 0 else mp.pi]
        return [mp.arg(z)]
    if name == 'abs': return [abs(z)]
    if name == 'norm_sqr': return [mp.re(z)**2 + mp.im(z)**2]
    if name == 'to_polar':
        angle = (-mp.pi if mp.im(z) == 0 and mp.re(z) < 0 and
                 math.copysign(1., a['p0'][1]) < 0 else mp.arg(z))
        return [abs(z), angle]

    # mpmath drops signed zeros. Correct exact-axis branch selection analytically
    # below; an infinitesimal perturbation is unreliable at branch points.
    pair = (a['x'], a['y']) if 'x' in a else a['p0']

    if name == 'add': result = z + other
    elif name == 'subtract': result = z - other
    elif name == 'multiply': result = z * other
    elif name == 'divide': result = z / other
    elif name == 'negate': result = -z
    elif name == 'conjugate': result = mp.conj(z)
    elif name == 'reciprocal': result = 1 / z
    elif name == 'powi' or name == 'powu': result = z ** other
    elif name == 'powf' or name == 'powc':
        # Integer exponents give exact algebraic powers on the axes.
        if mp.im(other) == 0 and mp.re(other) == int(mp.re(other)):
            result = z ** int(mp.re(other))
        else:
            result = mp.power(z, other)
    elif name == 'cbrt': result = mp.root(z, 3)
    elif name == 'ln': result = mp.log(z)
    else:
        functions = {
            'acos': mp.acos, 'acosh': mp.acosh, 'asin': mp.asin,
            'asinh': mp.asinh, 'atan': mp.atan, 'atanh': mp.atanh,
            'cos': mp.cos, 'cosh': mp.cosh, 'exp': mp.exp,
            'sin': mp.sin, 'sinh': mp.sinh, 'sqrt': mp.sqrt,
            'tan': mp.tan, 'tanh': mp.tanh,
        }
        result = functions[name](z)
    if pair[1] == 0.:
        lower = math.copysign(1., pair[1]) < 0
        if name in ('asin', 'atanh') and abs(pair[0]) > 1:
            result = mp.mpc(mp.re(result), -abs(mp.im(result)) if lower else abs(mp.im(result)))
        elif name == 'acos' and abs(pair[0]) > 1:
            result = mp.mpc(mp.re(result), abs(mp.im(result)) if lower else -abs(mp.im(result)))
        elif name == 'acosh' and pair[0] < 1:
            result = mp.mpc(mp.re(result), -abs(mp.im(result)) if lower else abs(mp.im(result)))
        elif lower and pair[0] < 0:
            if name in ('ln', 'sqrt', 'cbrt'):
                result = mp.conj(result)
            elif name in ('powf', 'powc') and mp.im(other) != 0:
                result = mp.exp(other * (mp.log(-mp.re(z)) - mp.pi * mp.j))
            elif name in ('powf', 'powc'):
                result = mp.conj(result)
        if name in ('powf', 'powc') and pair[0] < 0 and mp.im(other) == 0:
            doubled = 2 * mp.re(other)
            if doubled == int(doubled) and int(doubled) % 2:
                result = mp.mpc(0, mp.im(result))
    elif pair[0] == 0. and name in ('asinh', 'atan') and abs(pair[1]) > 1:
        real = abs(mp.re(result))
        result = mp.mpc(-real if math.copysign(1., pair[0]) < 0 else real, mp.im(result))
    return [mp.re(result), mp.im(result)]


def numpy_comparison(case):
    """Public NumPy complex operations; never present a formula as SciPy."""
    name, a = case['api'].split('.')[-1], case['args']
    try:
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter('always')
            if name == 'from_polar':
                result = a['p0'] * (np.cos(a['p1']) + 1j * np.sin(a['p1']))
            else:
                pair = (a['x'], a['y']) if 'x' in a else a['p0']
                z = np.complex128(complex(*pair))
                w = a.get('p1')
                other = np.complex128(complex(*w)) if isinstance(w, list) else w
                binary = {'add': np.add, 'subtract': np.subtract,
                          'multiply': np.multiply, 'divide': np.divide,
                          'powf': np.power, 'powc': np.power,
                          'powi': np.power, 'powu': np.power}
                unary = {'acos': np.arccos, 'acosh': np.arccosh,
                         'asin': np.arcsin, 'asinh': np.arcsinh,
                         'atan': np.arctan, 'atanh': np.arctanh,
                         'ln': np.log, 'conjugate': np.conjugate,
                         'negate': np.negative, 'reciprocal': np.reciprocal,
                         'abs': np.abs, 'arg': np.angle, 'real': np.real,
                         'imag': np.imag, 'sin': np.sin, 'sinh': np.sinh,
                         'cos': np.cos, 'cosh': np.cosh, 'tan': np.tan,
                         'tanh': np.tanh, 'sqrt': np.sqrt, 'exp': np.exp}
                if name in binary: result = binary[name](z, other)
                elif name in unary: result = unary[name](z)
                elif name == 'cbrt': result = np.exp(np.log(z) / 3)
                elif name == 'norm_sqr': result = np.real(np.multiply(z, np.conjugate(z)))
                elif name == 'to_polar': result = (np.abs(z), np.angle(z))
                else: raise ValueError(f'No NumPy complex comparison: {name}')
            parts = result if isinstance(result, tuple) else (result,)
            values = [component for part in parts for component in
                      ([float(np.real(part)), float(np.imag(part))] if np.iscomplexobj(part)
                       else [float(part)])]
            if not all(math.isfinite(value) for value in values):
                raise ValueError('nonfinite output')
            output = {'values': values}
            if captured: output['warnings'] = [str(item.message) for item in captured]
            return {'NumPy': output}
    except (ArithmeticError, ValueError, TypeError) as error:
        return {'NumPy': {'error': f'{type(error).__name__}: {error}'}}
