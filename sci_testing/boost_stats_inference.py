"""Boost.Math derived inference and circular statistics on pinned Sci cases.

Native operations use public Boost.Math moment and distribution functions.
Multi-group inputs preserve the saved g matrix's individual group boundaries.
"""
import ctypes as ct
import math


NAMES = frozenset({
    "ttest_ind_from_stats", "f_oneway", "bartlett",
    "circmean", "circvar", "circstd",
})
GROUPED = frozenset({"f_oneway", "bartlett"})


class BoostStatsInferenceComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.boost_stats_inference_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.boost_stats_inference_supported.restype = ct.c_int
        self.library.boost_stats_inference_output_count.argtypes = (ct.c_char_p,)
        self.library.boost_stats_inference_output_count.restype = ct.c_int
        self.library.boost_stats_inference_error.argtypes = ()
        self.library.boost_stats_inference_error.restype = ct.c_char_p
        self.library.boost_stats_inference_evaluate.argtypes = (
            ct.c_int, ct.c_char_p,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_int), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
        )
        self.library.boost_stats_inference_evaluate.restype = ct.c_int

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("stats."):
            return {}
        name = api.rsplit(".", 1)[-1]
        if name not in NAMES:
            return {}
        args = case.get("args", {})
        x, y, p = args.get("x", []), args.get("y", []), args.get("p", [])
        g = args.get("g", [])
        if not all(isinstance(v, (list, tuple)) for v in (x, y, p, g)):
            return {}
        if name in GROUPED and (
            len(g) < 2 or any(not isinstance(row, (list, tuple)) for row in g)
        ):
            return {}
        if name not in GROUPED and g:
            return {}
        try:
            x, y, p = [float(v) for v in x], [float(v) for v in y], [float(v) for v in p]
            groups = [[float(v) for v in row] for row in g]
        except (ValueError, TypeError, OverflowError):
            return {}
        flat = [v for row in groups for v in row]
        if not all(math.isfinite(v) for v in (*x, *y, *p, *flat)):
            return {}
        encoded = name.encode("ascii")
        n_output = self.library.boost_stats_inference_output_count(encoded)
        if not n_output or not self.library.boost_stats_inference_supported(0, encoded):
            return {}

        def doubles(values):
            return (ct.c_double * len(values))(*values) if values else None

        sizes = [len(row) for row in groups]
        native_sizes = (ct.c_int * len(sizes))(*sizes) if sizes else None
        output = (ct.c_double * n_output)()
        status = self.library.boost_stats_inference_evaluate(
            0, encoded, doubles(x), len(x), doubles(y), len(y),
            doubles(p), len(p), doubles(flat), len(flat),
            native_sizes, len(sizes), output, n_output,
        )
        if status == 0:
            return {"Boost.Math": {"values": list(output)}}
        return {"Boost.Math": {
            "error": self.library.boost_stats_inference_error().decode(errors="replace")
        }}
