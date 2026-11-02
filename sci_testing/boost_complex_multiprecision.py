"""Bounded independent complex comparisons using Boost.Multiprecision.

The C++ backend is Boost's public cpp_complex_50 (Boost >= 1.83),
not std::complex and not a standalone Boost.Math special function.
native/boost_complex_multiprecision.cpp must be linked in the native
shared object. A return of {} means unavailable for a saved input domain;
this is intentionally NOT an unsupported-case numerical result.
"""

import ctypes as ct
import math


SCALAR_NAMES = frozenset({"abs", "arg", "norm_sqr", "real", "imag"})
PAIR_NAMES = frozenset({
    "acos", "acosh", "add", "asin", "asinh", "atan", "atanh",
    "cbrt", "conjugate", "cos", "cosh", "divide", "exp",
    "from_polar", "ln", "multiply", "negate", "powc", "powf",
    "powi", "powu", "reciprocal", "sin", "sinh", "sqrt",
    "subtract", "tan", "tanh", "to_polar",
})
NAMES = SCALAR_NAMES | PAIR_NAMES
PROVIDER = "Boost.Multiprecision"
_COMPLEX_SECOND = frozenset({"add", "divide", "multiply", "powc", "subtract"})
_REAL_SECOND = frozenset({"powf", "powi", "powu"})


def _components(case, name):
    values = case.get("args", {})
    if name in ("exp", "sqrt"):
        return (values["x"], values["y"])
    if name == "from_polar":
        return (values["p0"], values["p1"])
    z = values["p0"]
    if not isinstance(z, (list, tuple)) or len(z) != 2:
        raise ValueError("expected p0 complex [real,imag]")
    args = list(z)
    if name in _COMPLEX_SECOND:
        w = values["p1"]
        if not isinstance(w, (list, tuple)) or len(w) != 2:
            raise ValueError("expected p1 complex [real,imag]")
        args.extend(w)
    elif name in _REAL_SECOND:
        args.append(values["p1"])
    return args


class BoostComplexMultiprecisionComparisons:
    """Public Boost complex C ABI, with per-input unavailable result semantics."""

    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.boost_complex_mp_supported.argtypes = (ct.c_char_p,)
        self.library.boost_complex_mp_supported.restype = ct.c_int
        self.library.boost_complex_mp_evaluate.argtypes = (
            ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
        )
        self.library.boost_complex_mp_evaluate.restype = ct.c_int
        self.library.boost_complex_mp_error.argtypes = ()
        self.library.boost_complex_mp_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("complex."):
            return {}
        name = api.split(".", 1)[1]
        if name not in NAMES:
            return {}
        encoded = name.encode("ascii")
        count = self.library.boost_complex_mp_supported(encoded)
        if count not in (1, 2):
            return {}
        try:
            args = tuple(float(x) for x in _components(case, name))
        except (KeyError, ValueError, TypeError, OverflowError):
            return {}
        if not all(math.isfinite(x) for x in args):
            return {}
        values = (ct.c_double * len(args))(*args)
        result = (ct.c_double * count)()
        status = self.library.boost_complex_mp_evaluate(
            encoded, values, len(args), result, count,
        )
        if status == 1:
            return {}
        if status != 0:
            error = self.library.boost_complex_mp_error()
            return {PROVIDER: {"error": (error or b"unknown native error").decode(
                "utf-8", errors="replace")}}
        return {PROVIDER: {"values": list(result)}}
