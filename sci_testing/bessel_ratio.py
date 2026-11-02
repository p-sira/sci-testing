"""Boost.Math public modified-Bessel I ratios and their inverse.

The C++ peer links as native/bessel_ratio.cpp (independently or alongside
comparisons.cpp) and exports bessel_ratio_* symbols. Input shapes match
the saved Sci 0.3 special03 cases: {v,x} or {v,r}.

The mathematical large-argument asymptotic continuation is disclosed by the
C++ implementation; no Boost private functions or Sci internals are called.
Unsupported domains return explicit Boost.Math error records rather than
inventing coverage or a numerically unvalidated value.
"""

import ctypes as ct


SUPPORTED_NAMES = frozenset(("iv_ratio", "iv_ratio_c", "iv_ratioinv"))


class BesselRatioComparisons:
    """Wrap the unique Boost Bessel-ratio C ABI in .evaluate(case)."""

    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.bessel_ratio_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.bessel_ratio_supported.restype = ct.c_int
        self.library.bessel_ratio_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
        )
        self.library.bessel_ratio_evaluate.restype = ct.c_int
        self.library.bessel_ratio_error.argtypes = ()
        self.library.bessel_ratio_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("special.bessel."):
            return {}
        name = api.rsplit(".", 1)[-1]
        if name not in SUPPORTED_NAMES:
            return {}
        args = case.get("args", {})
        key = "r" if name == "iv_ratioinv" else "x"
        if "v" not in args or key not in args:
            return {}
        encoded = name.encode("ascii")
        if self.library.bessel_ratio_supported(0, encoded) != 2:
            return {}
        try:
            values = (ct.c_double * 2)(float(args["v"]), float(args[key]))
        except (ValueError, TypeError, OverflowError):
            return {}
        output = ct.c_double()
        status = self.library.bessel_ratio_evaluate(
            0, encoded, values, 2, ct.byref(output), 1)
        if status:
            message = self.library.bessel_ratio_error()
            return {"Boost.Math": {"error": (message or b"unknown Bessel-ratio error").decode(
                "utf-8", errors="replace")}}
        return {"Boost.Math": {"values": [output.value]}}
