"""Direct Boost.Math/GSL spherical Bessel comparisons for Riccati arrays.

SciPy's riccati_jn and riccati_yn return two order-0-through-n vectors:
values followed by their first derivatives. The native bridge constructs
these vectors from public spherical Bessel APIs and their analytic derivative
recurrence. Unlike several scalar native bridges, output size varies with n.
"""

import ctypes as ct
import math


LIBRARIES = ("Boost.Math", "GSL")
NAMES = frozenset(("riccati_jn", "riccati_yn"))


class BesselRiccatiComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.bessel_riccati_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.bessel_riccati_supported.restype = ct.c_int
        self.library.bessel_riccati_output_count.argtypes = (
            ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int)
        self.library.bessel_riccati_output_count.restype = ct.c_int
        self.library.bessel_riccati_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int)
        self.library.bessel_riccati_evaluate.restype = ct.c_int
        self.library.bessel_riccati_error.argtypes = ()
        self.library.bessel_riccati_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("special."):
            return {}
        name = api.rsplit(".", 1)[-1]
        if name not in NAMES:
            return {}
        args = case.get("args", {})
        if set(args) != {"p0", "p1"}:
            return {}
        try:
            # No vectors, arrays, or complex arguments can enter this real ABI.
            n = float(args["p0"])
            x = float(args["p1"])
        except (TypeError, ValueError, OverflowError):
            return {}
        if not math.isfinite(n) or not math.isfinite(x):
            return {}
        if n < 0 or n > 512 or not n.is_integer():
            return {}
        encoded = name.encode("ascii")
        params = (ct.c_double * 2)(n, x)
        capacity = self.library.bessel_riccati_output_count(encoded, params, 2)
        if capacity <= 0:
            return {}
        results = {}
        for lib, label in enumerate(LIBRARIES):
            if not self.library.bessel_riccati_supported(lib, encoded):
                continue
            output = (ct.c_double * capacity)()
            status = self.library.bessel_riccati_evaluate(
                lib, encoded, params, 2, output, capacity)
            if status:
                results[label] = {
                    "error": self.library.bessel_riccati_error().decode(
                        "utf-8", errors="replace")}
            elif not all(math.isfinite(value) for value in output):
                results[label] = {"error": "nonfinite Riccati result"}
            else:
                results[label] = {"values": list(output)}
        return results
