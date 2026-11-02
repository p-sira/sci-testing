"""Full real Legendre arrays from Boost.Math's public polynomial functions.

The independent C ABI exports `boost_special_arrays_*`. Its output concatenates
the two SciPy vectors for `special.roots_legendre` and `special.lqn`.
Boost.Math directly provides the roots, Legendre derivatives, and interior
Legendre Q; standard Legendre recurrences complete the two array contracts.
"""

import ctypes as ct
import math


NAMES = frozenset(("roots_legendre", "lqn"))


class BoostSpecialArraysComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.boost_special_arrays_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.boost_special_arrays_supported.restype = ct.c_int
        self.library.boost_special_arrays_output_count.argtypes = (
            ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int)
        self.library.boost_special_arrays_output_count.restype = ct.c_int
        self.library.boost_special_arrays_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int)
        self.library.boost_special_arrays_evaluate.restype = ct.c_int
        self.library.boost_special_arrays_error.argtypes = ()
        self.library.boost_special_arrays_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("special."):
            return {}
        name = api.rsplit(".", 1)[-1]
        if name not in NAMES:
            return {}
        args = case.get("args", {})
        if set(args) != ({"p0"} if name == "roots_legendre" else {"p0", "p1"}):
            return {}
        try:
            n = float(args["p0"])
            x = float(args["p1"]) if name == "lqn" else None
        except (TypeError, ValueError, OverflowError):
            return {}
        if not math.isfinite(n) or n != int(n) or n < (1 if x is None else 0) or n > 128:
            return {}
        if x is not None and (not math.isfinite(x) or abs(x) == 1 or
                              (abs(x) > 1 and
                               (abs(x) > 1e12 or (abs(x) < 1.0001 and n > 8)))):
            return {}
        encoded = name.encode("ascii")
        values = [n] if x is None else [n, x]
        params = (ct.c_double * len(values))(*values)
        size = self.library.boost_special_arrays_output_count(encoded, params, len(values))
        if size <= 0 or self.library.boost_special_arrays_supported(0, encoded) != 2:
            return {}
        output = (ct.c_double * size)()
        status = self.library.boost_special_arrays_evaluate(
            0, encoded, params, len(values), output, size)
        if status:
            error = self.library.boost_special_arrays_error()
            return {"Boost.Math": {"error": (error or b"Boost array error").decode(
                "utf-8", errors="replace")}}
        if not all(math.isfinite(value) for value in output):
            return {"Boost.Math": {"error": "nonfinite Legendre array output"}}
        return {"Boost.Math": {"values": list(output)}}
