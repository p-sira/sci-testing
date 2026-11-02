"""Public Boost.Graph comparisons for seven saved sparse.csgraph APIs.

This is not a general SciPy sparse.csgraph replacement. It supports only
the COO-encodable, real, finite, bounded graph fixtures and their exact
saved output contracts: distance matrices, component labels, and MST total
weight. SciPy predecessor matrices are deliberately not represented because
the saved references deliberately discard them.

For the sparse-parity reference, every finite scalar is tagged as [0,value],
positive/negative infinities as [sign,0] and NaNs as [2,0]. This adapter
preserves that *comparison* encoding, rather than pretending Inf is finite.
"""

import ctypes as ct
import math

NAME = "Boost.Graph"
NAMES = frozenset((
    "connected_components", "dijkstra", "bellman_ford",
    "floyd_warshall", "johnson", "minimum_spanning_tree", "shortest_path",
))
EXTRA_NAMES = frozenset((
    "breadth_first_order", "breadth_first_tree",
    "depth_first_order", "depth_first_tree",
    "maximum_flow", "maximum_bipartite_matching",
    "structural_rank", "laplacian",
))


def _tagged(values):
    result = []
    for value in values:
        if math.isfinite(value):
            result.extend((0., value))
        elif math.isnan(value):
            result.extend((2., 0.))
        else:
            result.extend((math.copysign(1., value), 0.))
    return result


class BoostGraphComparisons:
    """Typed C ABI over Boost.Graph public shortest-path/graph algorithms."""

    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.boost_graph_error.argtypes = ()
        self.library.boost_graph_error.restype = ct.c_char_p
        self.library.boost_graph_supported.argtypes = (ct.c_char_p,)
        self.library.boost_graph_supported.restype = ct.c_int
        self.library.boost_graph_output_count.argtypes = (
            ct.c_char_p, ct.c_int, ct.POINTER(ct.c_double), ct.c_int, ct.c_int)
        self.library.boost_graph_output_count.restype = ct.c_int
        self.library.boost_graph_evaluate.argtypes = (
            ct.c_char_p, ct.c_int,
            ct.POINTER(ct.c_int), ct.POINTER(ct.c_int), ct.POINTER(ct.c_double),
            ct.c_int, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_int), ct.c_int, ct.POINTER(ct.c_double), ct.c_int)
        self.library.boost_graph_evaluate.restype = ct.c_int
        # Rectangular bipartite matrices and traversal orders with variable
        # lengths require a separate, fully capacity-checked C entry point.
        self.library.boost_graph_extra_supported.argtypes = (ct.c_char_p,)
        self.library.boost_graph_extra_supported.restype = ct.c_int
        self.library.boost_graph_extra_output_capacity.argtypes = (
            ct.c_char_p, ct.c_int, ct.c_int)
        self.library.boost_graph_extra_output_capacity.restype = ct.c_int
        self.library.boost_graph_extra_evaluate.argtypes = (
            ct.c_char_p, ct.c_int, ct.c_int, ct.POINTER(ct.c_double),
            ct.c_int, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int, ct.POINTER(ct.c_int))
        self.library.boost_graph_extra_evaluate.restype = ct.c_int

    def _extra_evaluate(self, name, raw):
        matrix = raw.get("a")
        opts = raw.get("p")
        selected = raw.get("x", [])
        if not isinstance(matrix, list) or not matrix or len(matrix) > 256:
            return {}
        if not isinstance(matrix[0], list):
            return {}
        nr, nc = len(matrix), len(matrix[0])
        if not 1 <= nc <= 256 or any(
            not isinstance(row, list) or len(row) != nc for row in matrix
        ):
            return {}
        if not isinstance(opts, list) or selected != []:
            return {}
        arity = (3 if name == "laplacian" else
                 2 if name in ("breadth_first_order", "breadth_first_tree",
                               "depth_first_order", "depth_first_tree",
                               "maximum_flow") else 1)
        if len(opts) != arity or any(
            type(v) not in (int, float) or not math.isfinite(v) for v in opts
        ):
            return {}
        if raw.get("zero_edge"):
            # Dense input alone cannot communicate the distinction between
            # an explicit stored 0-weight edge and an absent sparse edge.
            return {}
        flat = []
        for row in matrix:
            for value in row:
                if type(value) not in (int, float) or not math.isfinite(value):
                    return {}
                flat.append(float(value))
        key = name.encode("ascii")
        capacity = self.library.boost_graph_extra_output_capacity(key, nr, nc)
        if capacity <= 0:
            return {}
        inputs = (ct.c_double * len(flat))(*flat)
        params = (ct.c_double * arity)(*opts)
        output = (ct.c_double * capacity)()
        written = ct.c_int(-1)
        status = self.library.boost_graph_extra_evaluate(
            key, nr, nc, inputs, len(flat), params, arity, output,
            capacity, ct.byref(written))
        if status:
            err = self.library.boost_graph_error()
            return {NAME: {"error": (err or b"Boost.Graph extended error").decode(
                "utf-8", errors="replace")}}
        if not 0 <= written.value <= capacity:
            return {NAME: {"error": "invalid extended graph output length"}}
        return {NAME: {"values": _tagged(output[:written.value])}}

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("sparse.csgraph."):
            return {}
        name = api.rsplit(".", 1)[-1]
        if name not in NAMES and name not in EXTRA_NAMES:
            return {}
        saved = case.get("args", {})
        if not isinstance(saved, dict) or saved.get("suite_name") != "sparse_parity":
            return {}
        raw = saved.get("raw")
        if not isinstance(raw, dict) or raw.get("api") != name:
            return {}
        if name in EXTRA_NAMES:
            return self._extra_evaluate(name, raw)
        matrix = raw.get("a")
        opts = raw.get("p")
        selected = raw.get("x", [])
        if not isinstance(matrix, list) or not matrix or len(matrix) > 256:
            return {}
        n = len(matrix)
        if any(not isinstance(row, list) or len(row) != n for row in matrix):
            return {}
        if not isinstance(opts, list) or not isinstance(selected, list):
            return {}
        if len(opts) not in ((2, 4) if name == "dijkstra" else (2,)):
            return {}
        if name != "dijkstra" and len(opts) != 2:
            return {}
        if name in ("connected_components", "minimum_spanning_tree") and selected:
            return {}
        # Do not silently coerce float vertex indices, complex edges,
        # malformed flags, booleans masquerading as node IDs, or NaN weights.
        if any(type(v) is not int or v < 0 or v >= n for v in selected):
            return {}
        if any(type(v) not in (int, float) or not math.isfinite(v) for v in opts):
            return {}
        if opts[0] not in (0, 1) or opts[1] not in (0, 1):
            return {}
        if len(opts) == 4 and (opts[2] not in (0, 1) or opts[3] < 0):
            return {}
        r, c, weight = [], [], []
        for i, row in enumerate(matrix):
            for j, value in enumerate(row):
                if type(value) not in (int, float) or not math.isfinite(value):
                    return {}
                if value != 0.:
                    r.append(i)
                    c.append(j)
                    weight.append(float(value))
        if raw.get("zero_edge"):
            if n < 2:
                return {}
            r.append(0)
            c.append(1)
            weight.append(0.)
        encoded = name.encode("ascii")
        count = self.library.boost_graph_output_count(
            encoded, n, (ct.c_double * len(opts))(*opts), len(opts),
            len(selected))
        if count < 1:
            return {}
        rows = (ct.c_int * len(r))(*r)
        cols = (ct.c_int * len(c))(*c)
        values = (ct.c_double * len(weight))(*weight)
        parameters = (ct.c_double * len(opts))(*opts)
        indices = (ct.c_int * len(selected))(*selected)
        output = (ct.c_double * count)()
        status = self.library.boost_graph_evaluate(
            encoded, n, rows, cols, values, len(r),
            parameters, len(opts), indices, len(selected), output, count)
        if status:
            err = self.library.boost_graph_error()
            return {NAME: {"error": (err or b"Boost.Graph error").decode(
                "utf-8", errors="replace")}}
        return {NAME: {"values": _tagged(output)}}
