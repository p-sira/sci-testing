"""Additional public GSL-derived oracles for saved statistics/special cases.

GSL supplies sorting, mean/TSS, weighted mean and real special-function
primitives. Compound statistical functions are mathematical reductions of
those primitives. The exports use a distinct C ABI so the existing adapters
can coexist without changes to their entry points.
"""

import ctypes as ct
import math

STATS_NAMES = frozenset((
    "entropy", "rankdata", "tiecorrect", "mode",
    "wasserstein_distance", "energy_distance", "boxcox_llf",
    "zmap", "expectile",
))
SPECIAL_ARITIES = {
    "exp2": 1, "exp10": 1, "cosm1": 1,
    "sindg": 1, "cosdg": 1, "tandg": 1, "cotdg": 1,
    "entr": 1, "xlogy": 2, "xlog1py": 2,
    "rel_entr": 2, "kl_div": 2,
}
SUPPORTED = frozenset({f"stats.{name}" for name in STATS_NAMES} |
                      {f"special.{name}" for name in SPECIAL_ARITIES})


class ExtraNativeStatsComparisons:
    """Evaluate supported saved cases through the independent GSL ABI."""

    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.extra_native_stats_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.extra_native_stats_supported.restype = ct.c_int
        self.library.extra_native_stats_output_count.argtypes = (
            ct.c_int, ct.c_char_p, ct.c_int, ct.c_int, ct.c_int)
        self.library.extra_native_stats_output_count.restype = ct.c_int
        self.library.extra_native_stats_evaluate.argtypes = (
            ct.c_int, ct.c_char_p,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
        )
        self.library.extra_native_stats_evaluate.restype = ct.c_int
        self.library.extra_native_stats_error.argtypes = ()
        self.library.extra_native_stats_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if api not in SUPPORTED:
            return {}
        args = case.get("args", {})
        if not isinstance(args, dict):
            return {}
        if api.startswith("special."):
            arity = SPECIAL_ARITIES[api.removeprefix("special.")]
            if set(args) != {f"p{i}" for i in range(arity)}:
                return {}
            x, y, p = [args[f"p{i}"] for i in range(arity)], [], []
        else:
            # Sci stats saved inputs encode x/y samples and p parameters;
            # d/g are empty for every operation supported here.
            if any(args.get(key) for key in ("d", "g")) or not all(
                    isinstance(args.get(key), (list, tuple)) for key in ("x", "y", "p")):
                return {}
            x, y, p = args["x"], args["y"], args["p"]
        try:
            x = [float(v) for v in x]
            y = [float(v) for v in y]
            p = [float(v) for v in p]
        except (ValueError, TypeError, OverflowError):
            return {}
        if (max(len(x), len(y), len(p)) > 100000 or
                not all(math.isfinite(v) for v in (*x, *y, *p))):
            return {}

        encoded = api.encode("ascii")
        provider_index = 1
        if not self.library.extra_native_stats_supported(provider_index, encoded):
            return {}
        capacity = self.library.extra_native_stats_output_count(
            provider_index, encoded, len(x), len(y), len(p))
        if capacity <= 0:
            return {}

        def buffer(values):
            return (ct.c_double * len(values))(*values) if values else None

        output = (ct.c_double * capacity)()
        status = self.library.extra_native_stats_evaluate(
            provider_index, encoded, buffer(x), len(x), buffer(y), len(y),
            buffer(p), len(p), output, capacity)
        if status:
            reason = self.library.extra_native_stats_error()
            return {"GSL": {"error": (reason or b"GSL-derived stats/special failure"
                                     ).decode("utf-8", errors="replace")}}
        if not all(math.isfinite(value) for value in output):
            return {"GSL": {"error": "nonfinite GSL-derived stats/special output"}}
        return {"GSL": {"values": list(output)}}
