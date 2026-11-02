"""GSL-DERIVED Kelvin comparisons through public real GSL quadrature.

No native Kelvin function or complex Bessel K is attributed to GSL.
The native translation unit evaluates the defining real/imaginary
integrals of K_0(x exp(i*pi/4)) and I_0(x exp(i*pi/4)), including
their x derivatives, using gsl_integration_qag.

The strictly validated real domain is 0.1 <= x <= 8. The scalar APIs
ker, kei, kerp and keip return one value. The kelvin API returns
eight interleaved real/imaginary values, corresponding to SciPy's
(Be, Ke, Bep, Kep) tuple of four complex values.
"""

import ctypes as ct
import math


NAMES = frozenset(("ker", "kei", "kerp", "keip", "kelvin"))
LIBRARY_LABEL = "GSL"


class BesselKelvinKComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.bessel_kelvin_k_error.argtypes = ()
        self.library.bessel_kelvin_k_error.restype = ct.c_char_p
        self.library.bessel_kelvin_k_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.bessel_kelvin_k_supported.restype = ct.c_int
        self.library.bessel_kelvin_k_output_count.argtypes = (ct.c_char_p,)
        self.library.bessel_kelvin_k_output_count.restype = ct.c_int
        self.library.bessel_kelvin_k_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int)
        self.library.bessel_kelvin_k_evaluate.restype = ct.c_int

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("special.") or api.startswith("special.complex."):
            return {}
        name = api.rsplit(".", 1)[-1]
        if name not in NAMES:
            return {}
        arguments = case.get("args", {})
        if not isinstance(arguments, dict) or set(arguments) != {"p0"}:
            return {}
        try:
            x = float(arguments["p0"])
        except (ValueError, TypeError, OverflowError):
            return {}
        # Let the C ABI report explicit errors for scalar nonfinite and
        # out-of-domain real values, rather than silently omitting them.
        encoded = name.encode("ascii")
        output_count = self.library.bessel_kelvin_k_supported(1, encoded)
        if not output_count:
            return {}
        inputs = (ct.c_double * 1)(x)
        outputs = (ct.c_double * output_count)()
        status = self.library.bessel_kelvin_k_evaluate(
            1, encoded, inputs, 1, outputs, output_count)
        if status:
            error = self.library.bessel_kelvin_k_error()
            return {LIBRARY_LABEL: {
                "error": (error or b"GSL-derived Kelvin evaluation error").decode(
                    "utf-8", errors="replace")}}
        if not all(math.isfinite(value) for value in outputs):
            return {LIBRARY_LABEL: {"error": "nonfinite Kelvin quadrature output"}}
        return {LIBRARY_LABEL: {"values": list(outputs)}}
