"""Additional auditable Boost.Math-derived statistics comparisons."""
import ctypes as ct
import math

NAMES = frozenset({
    "gmean", "hmean", "gstd", "jarque_bera", "chisquare",
    "power_divergence", "mode", "binomtest", "fisher_exact",
})


class BoostStatsMoreComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.boost_stats_more_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.boost_stats_more_supported.restype = ct.c_int
        self.library.boost_stats_more_output_count.argtypes = (ct.c_char_p,)
        self.library.boost_stats_more_output_count.restype = ct.c_int
        self.library.boost_stats_more_error.argtypes = ()
        self.library.boost_stats_more_error.restype = ct.c_char_p
        self.library.boost_stats_more_evaluate.argtypes = (
            ct.c_int, ct.c_char_p,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
        )
        self.library.boost_stats_more_evaluate.restype = ct.c_int

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("stats.") or api.rsplit(".", 1)[-1] not in NAMES:
            return {}
        args = case.get("args", {})
        try:
            x = [float(v) for v in args.get("x", [])]
            y = [float(v) for v in args.get("y", [])]
            p = [float(v) for v in args.get("p", [])]
        except (ValueError, TypeError, OverflowError):
            return {}
        if not all(math.isfinite(v) for v in (*x, *y, *p)):
            return {}
        name = api.rsplit(".", 1)[-1].encode("ascii")
        if not self.library.boost_stats_more_supported(0, name):
            return {}
        output_count = self.library.boost_stats_more_output_count(name)
        def native(values):
            return (ct.c_double * len(values))(*values) if values else None
        output = (ct.c_double * output_count)()
        status = self.library.boost_stats_more_evaluate(
            0, name, native(x), len(x), native(y), len(y), native(p), len(p),
            output, output_count)
        return {"Boost.Math": {"values": list(output)}} if not status else {
            "Boost.Math": {"error": self.library.boost_stats_more_error().decode()}
        }
