"""Public GSL-derived comparisons for ten Sci 0.3 statistical APIs.

Native companion: ``native/gsl_stats_more.cpp``. GSL supplies public means,
variances, total sums of squares, and chi-square/t/F distribution tails. The
ten functions are mathematical reductions of those primitives, not claimed
to be separate GSL exported functions.
"""

import ctypes as ct
import math


SUPPORTED_NAMES = frozenset((
    "gmean", "hmean", "circmean", "circstd", "circvar", "jarque_bera",
    "chisquare", "ttest_ind_from_stats", "f_oneway", "bartlett",
))

GROUP_NAMES = frozenset(("f_oneway", "bartlett"))
X_NAMES = frozenset(("gmean", "hmean", "circmean", "circstd", "circvar",
                     "jarque_bera", "chisquare"))
EXPECTED_P = {
    "gmean": 0, "hmean": 0, "circmean": 2, "circstd": 3,
    "circvar": 2, "jarque_bera": 0, "chisquare": 1,
    "ttest_ind_from_stats": 8, "f_oneway": 0, "bartlett": 0,
}


class StatsMoreGSLComparisons:
    """Evaluate a supported saved stats recipe using the GSL C ABI."""

    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.gsl_stats_more_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.gsl_stats_more_supported.restype = ct.c_int
        self.library.gsl_stats_more_evaluate.argtypes = (
            ct.c_int, ct.c_char_p,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_int), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
        )
        self.library.gsl_stats_more_evaluate.restype = ct.c_int
        self.library.gsl_stats_more_error.argtypes = ()
        self.library.gsl_stats_more_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("stats."):
            return {}
        name = api.rsplit(".", 1)[-1]
        if name not in SUPPORTED_NAMES or api != f"stats.{name}":
            return {}
        a = case.get("args", {})
        x, y, p, groups = (a.get(key, [] if key != "g" else [])
                           for key in ("x", "y", "p", "g"))
        if not all(isinstance(arr, (list, tuple)) for arr in (x, y, p, groups)):
            return {}
        if len(p) != EXPECTED_P[name]:
            return {}
        if name in GROUP_NAMES:
            if x or y or not groups or not all(isinstance(v, (list, tuple))
                                               for v in groups):
                return {}
        elif groups:
            return {}
        if name in X_NAMES and not x:
            return {}
        if name not in X_NAMES and x:
            return {}
        if name != "chisquare" and y:
            return {}
        try:
            sample = [float(value) for value in x]
            second = [float(value) for value in y]
            parameters = [float(value) for value in p]
            sizes = [len(row) for row in groups]
            flattened = [float(value) for row in groups for value in row]
        except (TypeError, ValueError, OverflowError):
            return {}
        if not all(math.isfinite(v) for v in (*sample, *second, *parameters,
                                             *flattened)):
            return {}
        if any(len(values) > 100000 for values in (sample, second, parameters,
                                                   flattened)):
            return {}

        def native(values, ctype=ct.c_double):
            return (ctype * len(values))(*values) if values else None

        encoded = name.encode("ascii")
        capacity = self.library.gsl_stats_more_supported(1, encoded)
        if not capacity:
            return {}
        output = (ct.c_double * capacity)()
        status = self.library.gsl_stats_more_evaluate(
            1, encoded,
            native(sample), len(sample), native(second), len(second),
            native(parameters), len(parameters), native(flattened), len(flattened),
            native(sizes, ct.c_int), len(sizes), output, capacity,
        )
        if status:
            message = self.library.gsl_stats_more_error()
            return {"GSL": {"error": (message or b"unknown GSL-derived statistical error")
                                  .decode("utf-8", errors="replace")}}
        return {"GSL": {"values": list(output)}}
