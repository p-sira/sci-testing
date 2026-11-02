"""Public Boost.Math pFq-derived Bessel, Kelvin and Wright comparisons.

This is a standalone native C ABI with explicitly restricted *real* domains.
The Boost implementation calls only public gamma/hypergeometric_pFq APIs;
it does not claim full SciPy complex, arbitrary-order or asymptotic support.

The native translation unit is ``native/bessel_hypergeom.cpp``. The prime
agent may link it into NativeComparisons' shared object and evaluate routing.
"""

import ctypes as ct
import math


BOOST = "Boost.Math"
NAMES = frozenset((
    "struve", "modstruve", "ber", "bei", "berp", "beip",
    "wright_bessel", "log_wright_bessel",
))


class BesselHypergeomComparisons:
    """Use C ABI functions bessel_hypergeom_supported/evaluate/error."""

    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.bessel_hypergeom_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.bessel_hypergeom_supported.restype = ct.c_int
        self.library.bessel_hypergeom_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
        )
        self.library.bessel_hypergeom_evaluate.restype = ct.c_int
        self.library.bessel_hypergeom_error.argtypes = ()
        self.library.bessel_hypergeom_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("special.") or api.startswith("special.complex."):
            return {}
        name = api.rsplit(".", 1)[-1]
        if name not in NAMES:
            return {}
        encoded = name.encode("ascii")
        count = self.library.bessel_hypergeom_supported(0, encoded)
        if count == 0:
            return {}
        params = case.get("args", {})
        if any(f"p{i}" not in params for i in range(count)):
            return {}
        try:
            values = [float(params[f"p{i}"]) for i in range(count)]
        except (ValueError, TypeError, OverflowError):
            return {}
        if any(not math.isfinite(v) for v in values):
            return {BOOST: {"error": "nonfinite real input is unsupported"}}
        data = (ct.c_double * count)(*values)
        result = (ct.c_double * 1)()
        status = self.library.bessel_hypergeom_evaluate(
            0, encoded, data, count, result, 1,
        )
        if status:
            error = self.library.bessel_hypergeom_error()
            return {BOOST: {"error": (error or b"Boost hypergeometric error").decode(
                "utf-8", errors="replace")}}
        return {BOOST: {"values": [result[0]]}}
