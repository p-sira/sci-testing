"""Independent GSL fixed-Gauss-Legendre integration for Sci's saved callback.

Only the explicitly saved `integrate.fixed_quad` 16-point integral of exp(x)
is supported. The GSL rule is direct, not interchangeable with adaptive QAG.
"""

import ctypes as ct
import math


class GSLNumericalCoreComparisons:
    def __init__(self, library):
        self.library = library
        library.gsl_numerical_core_supported.argtypes = (ct.c_int, ct.c_char_p)
        library.gsl_numerical_core_supported.restype = ct.c_int
        library.gsl_numerical_core_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
        )
        library.gsl_numerical_core_evaluate.restype = ct.c_int
        library.gsl_numerical_core_error.argtypes = ()
        library.gsl_numerical_core_error.restype = ct.c_char_p

    def evaluate(self, case):
        if case.get('api') != 'integrate.fixed_quad':
            return {}
        args = case.get('args', {})
        if set(args) != {'a', 'b'}:
            return {}
        try:
            a, b = float(args['a']), float(args['b'])
        except (ValueError, TypeError, OverflowError):
            return {}
        if not math.isfinite(a) or not math.isfinite(b):
            return {}
        values = (ct.c_double * 2)(a, b)
        output = (ct.c_double * 1)()
        status = self.library.gsl_numerical_core_evaluate(
            1, b'fixed_quad', values, 2, output, 1)
        if status:
            return {'GSL': {'error': self.library.gsl_numerical_core_error().decode()}}
        return {'GSL': {'values': [output[0]]}}
