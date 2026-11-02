"""Boost.Math 1.83 public vector-norm comparisons for Sci spatial cases.

Scalar Euclidean/squared Euclidean/L1/L-infinity/Hamming and integer-p
Minkowski use public Boost.Math tools::*_distance functions. Pairwise
and weighted outputs are explicitly derived compositions, not independent
native Boost.Geometry matrices or weighted statistical metrics.
"""

import ctypes as ct
import math


NAMES = frozenset((
    "euclidean", "sqeuclidean", "cityblock", "chebyshev",
    "distance_matrix", "pdist", "cdist", "minkowski",
    "minkowski_distance", "minkowski_distance_p", "hamming",
))
PAIRED = frozenset(("distance_matrix", "pdist", "cdist"))


def _matrix(text):
    if not text:
        return []
    return [[float(value) for value in row.split(";")]
            for row in text.split(":")]


def _flat(rows):
    width = len(rows[0]) if rows else 0
    if any(len(row) != width for row in rows):
        raise ValueError("ragged spatial matrix")
    return [value for row in rows for value in row], width


class BoostSpatialComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.boost_spatial_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.boost_spatial_supported.restype = ct.c_int
        self.library.boost_spatial_output_count.argtypes = (
            ct.c_char_p, ct.c_int, ct.c_int, ct.c_int)
        self.library.boost_spatial_output_count.restype = ct.c_int
        self.library.boost_spatial_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int, ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int, ct.c_char_p, ct.c_double,
            ct.POINTER(ct.c_double), ct.c_int)
        self.library.boost_spatial_evaluate.restype = ct.c_int
        self.library.boost_spatial_error.argtypes = ()
        self.library.boost_spatial_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("spatial."):
            return {}
        name = api.rsplit(".", 1)[-1]
        if name not in NAMES:
            return {}
        arguments = case.get("args", {})
        try:
            if not isinstance(arguments, dict):
                return {}
            if isinstance(arguments.get("raw"), dict):
                raw = arguments["raw"]
                u = _matrix(raw.get("u", ""))
                v = _matrix(raw.get("v", ""))
                w = [float(value) for value in raw.get("w", "").split(";")
                     ] if raw.get("w") else []
                metric = str(raw.get("param", "")) if name in ("pdist", "cdist") else ""
                if name in ("distance_matrix", "minkowski",
                            "minkowski_distance", "minkowski_distance_p"):
                    p = float(raw.get("param", "2"))
                elif metric == "minkowski":
                    p = float(raw.get("extra", "2"))
                else:
                    p = 2.
            elif set(arguments) == {"u", "v"}:
                u = [list(map(float, arguments["u"]))]
                v = [list(map(float, arguments["v"]))]
                w, metric, p = [], "", 2.
            else:
                return {}
            if name not in PAIRED:
                # An empty vector is valid for scalar Lp distances. It has
                # one logical row, but zero feature columns.
                if not u: u = [[]]
                if not v: v = [[]]
            u_data, dim = _flat(u)
            v_data, vdim = _flat(v)
            if name == "pdist":
                if v:
                    raise ValueError("pdist expects a single point matrix")
            elif vdim != dim:
                raise ValueError("spatial input dimensions differ")
            if not all(math.isfinite(float(value))
                       for value in (*u_data, *v_data, *w)):
                raise ValueError("nonfinite real spatial coordinates/weights")
            if any(weight < 0. for weight in w):
                raise ValueError("negative spatial weight")
        except (TypeError, ValueError, OverflowError, KeyError) as exc:
            return {"Boost.Math": {"error": f"invalid spatial inputs: {exc}"}}

        encoded = name.encode("ascii")
        if not self.library.boost_spatial_supported(0, encoded):
            return {}
        count = self.library.boost_spatial_output_count(
            encoded, len(u), len(v), dim)
        if count <= 0:
            return {"Boost.Math": {"error": "invalid spatial output shape"}}
        inputs_u = (ct.c_double * max(1, len(u_data)))(*u_data)
        inputs_v = (ct.c_double * max(1, len(v_data)))(*v_data)
        weights = (ct.c_double * len(w))(*w) if w else None
        output = (ct.c_double * count)()
        status = self.library.boost_spatial_evaluate(
            0, encoded, inputs_u, len(u), inputs_v, len(v), dim, weights,
            len(w), metric.encode("ascii"), p, output, count,
        )
        if status:
            reason = self.library.boost_spatial_error()
            return {"Boost.Math": {"error": (reason or b"Boost.Math spatial error").decode(
                "utf-8", errors="replace")}}
        if not all(math.isfinite(value) for value in output):
            return {"Boost.Math": {"error": "nonfinite spatial result"}}
        return {"Boost.Math": {"values": list(output)}}
