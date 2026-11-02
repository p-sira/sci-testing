"""Boost.Math 1.83 public-2F3 comparisons for three real Struve integrals.

Link native/bessel_struve_integrals.cpp alongside the existing native shared
object, or compile it separately. This adapter uses only Boost's public
``hypergeometric_pFq`` and returns a scalar [value] with explicit errors for
unvalidated regions. GSL has no matching public 2F3 implementation.
"""

import ctypes as ct


SUPPORTED_NAMES = frozenset(("itstruve0", "itmodstruve0", "it2struve0"))


class BesselStruveIntegralsComparisons:
    """Evaluate scalar real-input Struve integrals from the separate C ABI."""

    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.bessel_struve_integrals_supported.argtypes = (
            ct.c_int, ct.c_char_p)
        self.library.bessel_struve_integrals_supported.restype = ct.c_int
        self.library.bessel_struve_integrals_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int)
        self.library.bessel_struve_integrals_evaluate.restype = ct.c_int
        self.library.bessel_struve_integrals_error.argtypes = ()
        self.library.bessel_struve_integrals_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("special."):
            return {}
        name = api.rsplit(".", 1)[-1]
        if name not in SUPPORTED_NAMES or api not in (
                f"special.{name}", f"special.bessel.{name}"):
            return {}

        args = case.get("args", {})
        if "p0" not in args:
            return {}
        encoded = name.encode("ascii")
        if self.library.bessel_struve_integrals_supported(0, encoded) != 1:
            return {}
        try:
            # A complex or vector argument cannot enter this scalar ABI.
            inputs = (ct.c_double * 1)(args["p0"])
        except (TypeError, ValueError, OverflowError):
            return {}

        output = (ct.c_double * 1)()
        status = self.library.bessel_struve_integrals_evaluate(
            0, encoded, inputs, 1, output, 1)
        if status:
            message = self.library.bessel_struve_integrals_error()
            return {"Boost.Math": {"error": (message or b"unknown Struve integral error")
                                    .decode("utf-8", errors="replace")}}
        return {"Boost.Math": {"values": [output[0]]}}
