"""Public Boost.Math and GSL adapters for Sci's special probability cases.

The separate C ABI in native/probability_extra.cpp takes exactly the scalar
arguments from cases/generated/special/*.csv or *.csv.gz (p0, p1, ...).
Unavailable functions/libraries are omitted; failed evaluation returns an
explicit error. The noncentral F inverse degree-of-freedom functions may have
multiple mathematical roots. Their Boost static solvers can reject such
ambiguous inputs; a successful result is verified against the original CDF.
The generated SciPy sweeps select the same branch, but arbitrary parameter
domains are not guaranteed to select the same root as SciPy's inverse.
"""

import ctypes as ct
import math


LIBRARIES = ("Boost.Math", "GSL")


class ProbabilityExtraComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.probability_extra_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.probability_extra_supported.restype = ct.c_int
        self.library.probability_extra_error.argtypes = ()
        self.library.probability_extra_error.restype = ct.c_char_p
        self.library.probability_extra_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double),
        )
        self.library.probability_extra_evaluate.restype = ct.c_int

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("special.") or api.startswith("special.complex."):
            return {}
        name = api.rsplit(".", 1)[-1]
        encoded = name.encode("ascii")
        args = case.get("args", {})
        results = {}
        for index, label in enumerate(LIBRARIES):
            count = self.library.probability_extra_supported(index, encoded)
            if not count:
                continue
            try:
                values = [float(args[f"p{i}"]) for i in range(count)]
            except (ValueError, TypeError, KeyError):
                results[label] = {"error": f"expected {count} numeric p0...p{count-1} arguments"}
                continue
            inputs = (ct.c_double * count)(*values)
            output = ct.c_double()
            status = self.library.probability_extra_evaluate(
                index, encoded, inputs, count, ct.byref(output))
            if status:
                results[label] = {"error": self.library.probability_extra_error().decode()}
            elif not math.isfinite(output.value):
                results[label] = {"error": "nonfinite probability result"}
            else:
                results[label] = {"values": [output.value]}
        return results
