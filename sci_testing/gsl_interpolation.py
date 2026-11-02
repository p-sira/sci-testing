"""Safe, independent public-GSL interpolation comparisons.

The polynomial GSL kernel matches SciPy barycentric interpolation and
non-confluent (distinct-abscissa) Krogh polynomial interpolation; GSL's
polynomial eval callback is used for extrapolation beyond the spline interval.
Sci PCHIP is NOT the GSL Steffen formula: only exactly affine data or
exact-node queries use the GSL Steffen spline.
Sci interpn uses multidimensional tensor interpolation. Its linear/nearest
algorithms are composed with GSL's public gsl_interp_bsearch interval search,
not misrepresented as a direct GSL n-D spline API.
"""

import ctypes as ct
import math


LIBRARY = "GSL"
SIMPLE = frozenset((
    "barycentric_interpolate", "krogh_interpolate", "pchip_interpolate"))
NAMES = SIMPLE | {"interpn"}


def _real(x):
    if type(x) not in (float, int):
        raise ValueError("numeric scalar required")
    number = float(x)
    if not math.isfinite(number):
        raise ValueError("nonfinite scalar")
    return number


def _vector(values):
    if not isinstance(values, list):
        raise ValueError("expected numerical vector")
    return [_real(value) for value in values]


def _count(value, maximum):
    number = _real(value)
    if number < 0 or number > maximum or number != int(number):
        raise ValueError("invalid integer count")
    return int(number)


def _simple(name, args):
    if set(args) not in ({"nodes", "values", "query"},
                         {"nodes", "values", "query", "explicit"}):
        raise ValueError("unsupported simple interpolation signature")
    if "explicit" in args and args["explicit"] is not True:
        raise ValueError("invalid explicit flag")
    nodes, values, query = (_vector(args[key]) for key in ("nodes", "values", "query"))
    n, nq = len(nodes), len(query)
    if n < 2 or n > 64 or nq < 1 or nq > 4096 or len(values) != n:
        raise ValueError("invalid interpolation vector dimensions")
    if not all(nodes[i] < nodes[i+1] for i in range(n-1)):
        raise ValueError("unsupported non-strictly-increasing nodes")
    if name == "pchip_interpolate":
        # Confluent/Steffen adaptations are explicitly excluded. GSL Steffen
        # and Sci PCHIP share all knot values and affine interpolants only.
        if n < 3 or any(q < nodes[0] or q > nodes[-1] for q in query):
            raise ValueError("unsupported PCHIP knot count/domain")
        slope = (values[1]-values[0])/(nodes[1]-nodes[0])
        affine = all((values[i]-values[i-1])/(nodes[i]-nodes[i-1]) == slope
                     for i in range(2,n))
        if not affine and not all(q in nodes for q in query):
            raise ValueError("GSL Steffen is not equivalent to Sci PCHIP")
    return [float(n), float(nq), *nodes, *values, *query]


def _interpn(args):
    if set(args) != {"fields", "raw", "suite", "suite_name"}:
        raise ValueError("unsupported interpn saved signature")
    if args["suite"] != "interpolation" or args["suite_name"] != "interpolation":
        raise ValueError("unsupported interpolation suite")
    raw = args["raw"]
    if not isinstance(raw, dict) or raw.get("api") != "interpn":
        raise ValueError("incorrect interpn raw payload")
    params = raw["args"]
    if not isinstance(params, list) or len(params) != 6:
        raise ValueError("unsupported interpn raw argument count")
    axes, values, query, method, bounds, fill = params
    if not isinstance(axes, list) or not 1 <= len(axes) <= 4:
        raise ValueError("invalid interpn axes")
    if method not in ("linear", "nearest") or type(bounds) is not bool:
        raise ValueError("invalid interpn method/bounds")
    if fill is not None:
        fill = _real(fill)
    axes = [_vector(axis) for axis in axes]
    lengths = [len(axis) for axis in axes]
    product = 1
    for axis in axes:
        if not 1 <= len(axis) <= 64:
            raise ValueError("invalid interpn grid length")
        if len(axis) > 1 and not (
            all(axis[i] < axis[i+1] for i in range(len(axis)-1))
            or all(axis[i] > axis[i+1] for i in range(len(axis)-1))
        ):
            raise ValueError("non-monotone regular-grid axis")
        product *= len(axis)
    if product > 32768:
        raise ValueError("regular grid too large")
    values = _vector(values)
    if len(values) != product:
        raise ValueError("regular grid shape mismatch")
    if not isinstance(query, list) or not 1 <= len(query) <= 4096:
        raise ValueError("invalid query count")
    points = [_vector(point) for point in query]
    if any(len(point) != len(axes) for point in points):
        raise ValueError("invalid regular grid query dimension")
    # A bounds_error=True outside query is deliberately not assigned a GSL
    # numerical comparator; SciPy itself raises an error for that contract.
    if bounds and any(
        point[d] < min(axis) or point[d] > max(axis)
        for point in points for d,axis in enumerate(axes)
    ):
        raise ValueError("interpn bounds_error=True outside the interpolation domain")
    return [
        float(len(axes)), *[float(n) for n in lengths],
        *(x for axis in axes for x in axis),
        *values,
        float(len(points)),
        *(x for point in points for x in point),
        0. if method == "linear" else 1.,
        float(bounds),
        float(fill is not None),
        0. if fill is None else fill,
    ]


class GslInterpolationComparisons:
    """One GSL provider, with per-case conservative numerical-domain guards."""

    def __init__(self, library):
        self.library = getattr(library, "library", library)
        lib = self.library
        lib.gsl_interpolation_supported.argtypes = (ct.c_int, ct.c_char_p)
        lib.gsl_interpolation_supported.restype = ct.c_int
        lib.gsl_interpolation_output_count.argtypes = (
            ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int)
        lib.gsl_interpolation_output_count.restype = ct.c_int
        lib.gsl_interpolation_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int)
        lib.gsl_interpolation_evaluate.restype = ct.c_int
        lib.gsl_interpolation_error.argtypes = ()
        lib.gsl_interpolation_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if not isinstance(api, str) or not api.startswith("interpolate."):
            return {}
        name = api.removeprefix("interpolate.")
        if name not in NAMES:
            return {}
        try:
            args = case.get("args", {})
            if not isinstance(args, dict):
                return {}
            params = _interpn(args) if name == "interpn" else _simple(name, args)
            if len(params) > 200000:
                return {}
        except (ValueError, TypeError, KeyError, OverflowError):
            return {}
        name_bytes = api.encode("ascii")
        if self.library.gsl_interpolation_supported(1, name_bytes) != 2:
            return {}
        input_values = (ct.c_double * len(params))(*params)
        count = self.library.gsl_interpolation_output_count(
            name_bytes, input_values, len(params))
        if count <= 0 or count > 4096:
            return {}
        output = (ct.c_double * count)()
        status = self.library.gsl_interpolation_evaluate(
            1, name_bytes, input_values, len(params), output, count)
        if status:
            error = self.library.gsl_interpolation_error()
            return {LIBRARY: {"error": (error or b"GSL interpolation failed").decode(
                "utf-8", errors="replace")}}
        if not all(math.isfinite(v) for v in output):
            return {LIBRARY: {"error": "nonfinite GSL interpolation result"}}
        return {LIBRARY: {"values": list(output)}}
