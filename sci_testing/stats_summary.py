"""Public Boost.Math/GSL vector statistics comparison bridge.

The native side lives in native/stats_summary.cpp and is linkable alongside
native/comparisons.cpp. It also builds standalone for isolated testing.

The preferred input schema is the saved Sci statistics recipe:
{'api': 'stats.descriptive.skew', 'args': {'x': [...], 'y': [], 'p': [1]}}.
Historical stats.variance cases carry {'values': [...]} and use ddof=1.
Only exact matching public library primitives or standard mathematical
conversions of their statistics are reported.
"""
import ctypes as ct
import math


LIBRARIES = ("Boost.Math", "GSL")
SCALAR_NAMES = frozenset({
    "mean", "variance", "std", "sem", "variation", "skew", "kurtosis",
    "median", "median_abs_deviation", "quantile", "scoreatpercentile",
    "iqr", "describe", "zscore",
})
PAIR_NAMES = frozenset({"pearsonr", "spearmanr", "spearmanrho", "pointbiserialr"})
SUPPORTED_NAMES = SCALAR_NAMES | PAIR_NAMES


class StatsSummaryComparisons:
    """Expose the native statistics C ABI through a loaded shared object.

    `library` can be a ctypes.CDLL or an existing NativeComparisons instance.
    """

    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.stats_summary_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.stats_summary_supported.restype = ct.c_int
        self.library.stats_summary_output_count.argtypes = (ct.c_char_p, ct.c_int)
        self.library.stats_summary_output_count.restype = ct.c_int
        self.library.stats_summary_evaluate.argtypes = (
            ct.c_int, ct.c_char_p,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
        )
        self.library.stats_summary_evaluate.restype = ct.c_int
        self.library.stats_summary_error.argtypes = ()
        self.library.stats_summary_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("stats."):
            return {}
        name = api.rsplit(".", 1)[-1]
        if name not in SUPPORTED_NAMES:
            return {}
        args = case.get("args", {})
        if "values" in args and name == "variance":
            # Pinned pre-0.3 stats.variance is sample variance with ddof=1.
            x, y, p = args["values"], [], [1]
        else:
            x, y, p = args.get("x"), args.get("y", []), args.get("p", [])
        if not isinstance(x, (list, tuple)) or not x:
            return {}
        if not isinstance(y, (list, tuple)) or not isinstance(p, (list, tuple)):
            return {}
        if name in PAIR_NAMES and (not y or len(x) != len(y)):
            return {}
        if name in SCALAR_NAMES and y:
            return {}
        # Numeric conversion happens once, before ctypes. Reject nonnumeric
        # sentinel values and protect the native function from malformed input.
        try:
            x = [float(value) for value in x]
            y = [float(value) for value in y]
            p = [float(value) for value in p]
        except (ValueError, TypeError, OverflowError):
            return {}
        if not all(math.isfinite(value) for value in (*x, *y, *p)):
            return {}
        encoded = name.encode("ascii")
        output_count = self.library.stats_summary_output_count(encoded, len(x))
        if output_count <= 0:
            return {}
        x_native = (ct.c_double * len(x))(*x)
        y_native = (ct.c_double * len(y))(*y) if y else None
        p_native = (ct.c_double * len(p))(*p) if p else None
        result = {}
        for index, library_name in enumerate(LIBRARIES):
            if not self.library.stats_summary_supported(index, encoded):
                continue
            output = (ct.c_double * output_count)()
            status = self.library.stats_summary_evaluate(
                index, encoded, x_native, len(x), y_native, len(y),
                p_native, len(p), output, output_count,
            )
            if status == 0:
                result[library_name] = {"values": list(output)}
            else:
                result[library_name] = {
                    "error": self.library.stats_summary_error().decode("utf-8", errors="replace")
                }
        return result
