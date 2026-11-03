"""Additional public Boost.Math/GSL oracles for saved real special APIs.

``assoc_laguerre`` follows SciPy's ``(x, n, k)`` convention, and Boost's
noninteger-parameter case uses its public confluent hypergeometric function.
Only actual native library calls are attributed to Boost.Math and GSL.
"""

import ctypes as ct
import math

ARITIES = {
    "special.psi": 1,
    "special.assoc_laguerre": 3,
    "special.expm1": 1,
    "special.log1p": 1,
    "special.cbrt": 1,
}


class ExtraNativeComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.extra_special_scalar_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.extra_special_scalar_supported.restype = ct.c_int
        self.library.extra_special_scalar_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
        )
        self.library.extra_special_scalar_evaluate.restype = ct.c_int
        self.library.extra_special_scalar_error.argtypes = ()
        self.library.extra_special_scalar_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        required = ARITIES.get(api)
        if required is None:
            return {}
        arguments = case.get("args")
        if not isinstance(arguments, dict) or set(arguments) != {
                f"p{i}" for i in range(required)}:
            return {}
        try:
            values = [float(arguments[f"p{i}"]) for i in range(required)]
        except (ValueError, TypeError, OverflowError):
            return {}
        if not all(math.isfinite(value) for value in values):
            return {}
        name = api.rsplit(".", 1)[-1].encode("ascii")
        inputs = (ct.c_double * required)(*values)
        result = {}
        for index, provider in enumerate(("Boost.Math", "GSL")):
            if self.library.extra_special_scalar_supported(index, name) != required:
                continue
            output = (ct.c_double * 1)()
            status = self.library.extra_special_scalar_evaluate(
                index, name, inputs, required, output, 1)
            if status:
                error = self.library.extra_special_scalar_error()
                result[provider] = {"error": (error or b"native scalar special error").decode(
                    "utf-8", errors="replace")}
            elif math.isfinite(output[0]):
                result[provider] = {"values": [output[0]]}
            else:
                result[provider] = {"error": "nonfinite native scalar special result"}
        return result
