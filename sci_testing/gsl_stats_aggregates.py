"""Public GSL-derived logarithmic, power and central sample statistics.

Uses the standalone C++ peer ``native/gsl_stats_aggregates.cpp``.
GSL provides mean/variance/total-sum-of-squares primitives; the four
corresponding SciPy operations are mathematical reductions, not separate
GSL function exports. Only real 1-D, unweighted saved signatures are mapped.
"""

import ctypes as ct
import math


SUPPORTED_NAMES = frozenset(("gstd", "gzscore", "pmean", "moment"))


class StatsAggregatesGSLComparisons:
    """Evaluate the four real sample aggregates with the native GSL ABI."""

    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.gsl_stats_aggregates_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.gsl_stats_aggregates_supported.restype = ct.c_int
        self.library.gsl_stats_aggregates_output_count.argtypes = (ct.c_char_p, ct.c_int)
        self.library.gsl_stats_aggregates_output_count.restype = ct.c_int
        self.library.gsl_stats_aggregates_evaluate.argtypes = (
            ct.c_int, ct.c_char_p,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
        )
        self.library.gsl_stats_aggregates_evaluate.restype = ct.c_int
        self.library.gsl_stats_aggregates_error.argtypes = ()
        self.library.gsl_stats_aggregates_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("stats.") or api.count(".") != 1:
            return {}
        name = api.rsplit(".", 1)[-1]
        if name not in SUPPORTED_NAMES:
            return {}
        arguments = case.get("args", {})
        sample = arguments.get("x", [])
        y = arguments.get("y", [])
        p = arguments.get("p", [])
        g = arguments.get("g", [])
        if (not isinstance(sample, (list, tuple)) or not sample or
                not isinstance(y, (list, tuple)) or y or
                not isinstance(p, (list, tuple)) or len(p) != 1 or
                not isinstance(g, (list, tuple)) or g):
            return {}
        try:
            numbers = [float(v) for v in sample]
            parameters = [float(p[0])]
        except (TypeError, ValueError, OverflowError):
            return {}
        if (len(numbers) > 100000 or not all(math.isfinite(v) for v in
                                            (*numbers, *parameters))):
            return {}
        encoded = name.encode("ascii")
        if not self.library.gsl_stats_aggregates_supported(1, encoded):
            return {}
        capacity = self.library.gsl_stats_aggregates_output_count(encoded,len(numbers))
        if capacity <= 0:
            return {}
        data = (ct.c_double * len(numbers))(*numbers)
        extra = (ct.c_double * 1)(parameters[0])
        output = (ct.c_double * capacity)()
        status = self.library.gsl_stats_aggregates_evaluate(
            1, encoded, data, len(numbers), None, 0,
            extra, 1, None, 0, output, capacity)
        if status:
            message = self.library.gsl_stats_aggregates_error()
            return {"GSL": {"error": (message or b"unknown GSL-derived aggregate error")
                                  .decode("utf-8", errors="replace")}}
        return {"GSL": {"values": list(output)}}
