"""Public Boost.Math 1.83+ interpolation value comparators.

Only pairwise-distinct scalar nodes and function values are accepted.
Krogh is supported as polynomial interpolation for distinct
nodes, NOT as confluent Krogh/Hermite data with repeated derivative nodes.
Boost pchip requires >=4 sorted nodes and does not extrapolate.
"""
import ctypes as ct
import math


SUPPORTED = frozenset((
    "pchip_interpolate",
    "barycentric_interpolate",
    "krogh_interpolate",
))


def _real_vector(value):
    if not isinstance(value, (list, tuple)):
        return None
    if any(type(v) not in (float, int) for v in value):
        return None
    if not all(math.isfinite(v) for v in value):
        return None
    return value


class BoostInterpolationComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        lib = self.library
        lib.boost_interpolation_supported.argtypes = (ct.c_char_p,)
        lib.boost_interpolation_supported.restype = ct.c_int
        lib.boost_interpolation_evaluate.argtypes = (
            ct.c_char_p, ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),
            ct.c_int, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
        )
        lib.boost_interpolation_evaluate.restype = ct.c_int
        lib.boost_interpolation_error.argtypes = ()
        lib.boost_interpolation_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if not isinstance(api, str) or not api.startswith("interpolate."):
            return {}
        name = api.rsplit(".", 1)[-1]
        if name not in SUPPORTED:
            return {}
        args = case.get("args")
        if not isinstance(args, dict):
            return {}
        if set(args) not in (
            {"nodes", "values", "query"},
            {"nodes", "values", "query", "explicit"},
        ):
            return {}
        if "explicit" in args and args["explicit"] is not True:
            return {}
        nodes = _real_vector(args["nodes"])
        values = _real_vector(args["values"])
        query = _real_vector(args["query"])
        if nodes is None or values is None or query is None:
            return {}
        n, nq = len(nodes), len(query)
        if n < 2 or n > 16 or n != len(values) or nq < 1 or nq > 2048:
            return {}
        encoded = name.encode("ascii")
        if not self.library.boost_interpolation_supported(encoded):
            return {}
        inp_x = (ct.c_double * n)(*nodes)
        inp_y = (ct.c_double * n)(*values)
        inp_q = (ct.c_double * nq)(*query)
        out = (ct.c_double * nq)()
        status = self.library.boost_interpolation_evaluate(
            encoded, inp_x, inp_y, n, inp_q, nq, out, nq,
        )
        if status:
            err = self.library.boost_interpolation_error()
            return {"Boost.Math": {"error": (err or b"Boost interpolation error").decode(
                "utf-8", errors="replace")}}
        return {"Boost.Math": {"values": list(out)}}
