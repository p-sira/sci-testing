"""GSL-quadrature-derived Kelvin ber, bei, berp, beip scalar adapters.

No native GSL Kelvin function is claimed: the C++ peer requests the existing
GSL-derived eight-component `kelvin` ABI and extracts the matching tuple lane.
"""

import ctypes as ct
import math


NAMES = frozenset(("ber", "bei", "berp", "beip"))


class GSLKelvinBerComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.gsl_kelvin_ber_error.argtypes = ()
        self.library.gsl_kelvin_ber_error.restype = ct.c_char_p
        self.library.gsl_kelvin_ber_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.gsl_kelvin_ber_supported.restype = ct.c_int
        self.library.gsl_kelvin_ber_output_count.argtypes = (ct.c_char_p,)
        self.library.gsl_kelvin_ber_output_count.restype = ct.c_int
        self.library.gsl_kelvin_ber_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
        )
        self.library.gsl_kelvin_ber_evaluate.restype = ct.c_int

    def evaluate(self, case):
        api = case.get("api", "")
        if api not in {"special." + name for name in NAMES}:
            return {}
        arguments = case.get("args", {})
        if not isinstance(arguments, dict) or set(arguments) != {"p0"}:
            return {}
        try:
            x = float(arguments["p0"])
        except (TypeError, ValueError, OverflowError):
            return {}

        name = api.rsplit(".", 1)[-1].encode("ascii")
        if self.library.gsl_kelvin_ber_supported(1, name) != 1:
            return {}
        inputs = (ct.c_double * 1)(x)
        output = (ct.c_double * 1)()
        status = self.library.gsl_kelvin_ber_evaluate(1, name, inputs, 1,
                                                      output, 1)
        if status:
            reason = self.library.gsl_kelvin_ber_error()
            return {"GSL": {"error": (reason or b"GSL Kelvin evaluation failed").decode(
                "utf-8", errors="replace")}}
        if not math.isfinite(output[0]):
            return {"GSL": {"error": "nonfinite GSL Kelvin quadrature result"}}
        return {"GSL": {"values": [output[0]]}}
