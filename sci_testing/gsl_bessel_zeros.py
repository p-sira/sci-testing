"""GSL-derived integer-order Bessel/Neumann zero comparisons (not GSL zero APIs).

The C++ peer independently uses public GSL integer Jn/Yn, J zero and Brent
functions; Python only transports finite real arguments and flattened results.
"""
import ctypes as ct

SUPPORTED = frozenset("special." + name for name in (
    "y0_zeros", "y1_zeros", "y1p_zeros", "yn_zeros", "ynp_zeros",
    "jnjnp_zeros", "jnp_zeros", "jnyn_zeros"))


class GSLBesselZerosComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.gsl_bessel_zeros_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.gsl_bessel_zeros_supported.restype = ct.c_int
        self.library.gsl_bessel_zeros_output_count.argtypes = (
            ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int)
        self.library.gsl_bessel_zeros_output_count.restype = ct.c_int
        self.library.gsl_bessel_zeros_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int)
        self.library.gsl_bessel_zeros_evaluate.restype = ct.c_int
        self.library.gsl_bessel_zeros_error.argtypes = ()
        self.library.gsl_bessel_zeros_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if api not in SUPPORTED:
            return {}
        name = api.rsplit(".", 1)[-1].encode("ascii")
        arity = self.library.gsl_bessel_zeros_supported(1, name)
        if not arity:
            return {}
        values = case.get("args", {})
        try:
            inputs = (ct.c_double * arity)(
                *(float(values[f"p{i}"]) for i in range(arity)))
        except (KeyError, ValueError, TypeError, OverflowError):
            return {"GSL": {"error": f"expected {arity} real Bessel zero arguments"}}
        capacity = self.library.gsl_bessel_zeros_output_count(name, inputs, arity)
        if capacity <= 0:
            return {"GSL": {"error": "invalid GSL Bessel zero count/order"}}
        output = (ct.c_double * capacity)()
        status = self.library.gsl_bessel_zeros_evaluate(
            1, name, inputs, arity, output, capacity)
        if status:
            reason = self.library.gsl_bessel_zeros_error()
            return {"GSL": {"error": (reason or b"GSL Bessel root failure").decode(
                "utf-8", errors="replace")}}
        return {"GSL": {"values": list(output)}}
