"""Public Boost.Numeric.uBLAS dense and sparse array competitors.

Dense inputs use the saved explicit or original linalg dictionaries. Sparse
inputs use the saved args.raw COO matrices, including explicit zeros and
duplicate entries. The native result carries row-major elements and shape;
the Python bridge transforms sparse values into Sci's exact (class, value)
finite-number serialization without changing their mathematical meaning.
"""

import ctypes as ct
import math


LIBRARY = "Boost.Numeric.uBLAS"
DENSE_NAMES = frozenset(("det", "solve", "dot", "matvec", "matmul", "transpose"))
SPARSE_NAMES = frozenset((
    "add", "multiply", "matvec", "matmul", "transpose", "from_dense", "to_dense"))
SPARSE_BINARY = frozenset(("add", "multiply", "matmul"))


def _real(value):
    # Reject complex, strings and booleans before entering the real C ABI.
    if type(value) not in (int, float):
        raise ValueError("input must be a real number")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("input must be finite")
    return number


def _vector(values):
    if not isinstance(values, list):
        raise ValueError("expected a vector")
    return [_real(v) for v in values]


def _matrix(values):
    if not isinstance(values, list) or not values or not all(
        isinstance(row, list) for row in values
    ):
        raise ValueError("expected nonempty matrix rows")
    cols = len(values[0])
    if cols == 0 or any(len(row) != cols for row in values):
        raise ValueError("ragged matrix")
    return len(values), cols, [v for row in values for v in _vector(row)]


def _dimension(value):
    number = _real(value)
    if number != int(number) or number < 0 or number > 128:
        raise ValueError("invalid dimension")
    return int(number)


def _coo(raw, key, rows, cols):
    triples = raw.get(key, [])
    if not isinstance(triples, list) or len(triples) > 32768:
        raise ValueError("invalid COO entries")
    encoded = [float(rows), float(cols), float(len(triples))]
    for triple in triples:
        if not isinstance(triple, list) or len(triple) != 3:
            raise ValueError("invalid COO triple")
        i, j, value = map(_real, triple)
        if i != int(i) or j != int(j) or not (0 <= i < rows and 0 <= j < cols):
            raise ValueError("invalid COO indices")
        encoded.extend((i, j, value))
    return encoded


def _dense_input(name, args):
    if name == "dot":
        if set(args) != {"a", "b", "u", "v", "explicit"}:
            raise ValueError("invalid saved dot inputs")
        u, v = _vector(args["u"]), _vector(args["v"])
        if len(u) != len(v):
            raise ValueError("dot size mismatch")
        return [len(u), *u, *v]
    if name in ("matmul", "matvec", "transpose"):
        if set(args) != {"a", "b", "u", "v", "explicit"}:
            raise ValueError("invalid saved dense inputs")
    elif set(args) != ({"a", "b"} if name == "solve" else {"a"}):
        raise ValueError("invalid saved linalg inputs")
    rows, cols, a = _matrix(args["a"])
    if name == "matmul":
        brows, bcols, b = _matrix(args["b"])
        if cols != brows:
            raise ValueError("invalid dense matmul shapes")
        return [rows, cols, *a, bcols, *b]
    if name == "matvec":
        x = _vector(args["u"])
        if len(x) != cols:
            raise ValueError("invalid dense matvec shapes")
        return [rows, cols, *a, *x]
    result = [rows, cols, *a]
    if name == "solve":
        b = _vector(args["b"])
        if rows != cols or len(b) != rows:
            raise ValueError("invalid solve input shape")
        result.extend(b)
    if name == "det" and rows != cols:
        raise ValueError("invalid determinant shape")
    return result


def _sparse_input(name, args):
    if set(args) != {"fields", "raw", "suite", "suite_name"}:
        raise ValueError("invalid saved sparse input format")
    if args["suite"] != "sparse" or args["suite_name"] != "sparse":
        raise ValueError("unsupported sparse suite")
    raw = args["raw"]
    if not isinstance(raw, dict) or raw.get("api") != name:
        raise ValueError("incorrect sparse raw operation")
    rows, cols = _dimension(raw["rows"]), _dimension(raw["cols"])
    data = _coo(raw, "a", rows, cols)
    if name in SPARSE_BINARY:
        brows, bcols = _dimension(raw["brows"]), _dimension(raw["bcols"])
        data.extend(_coo(raw, "b", brows, bcols))
    elif name == "matvec":
        x = _vector(raw["x"])
        if len(x) != cols:
            raise ValueError("invalid sparse matvec input shape")
        data.extend(x)
    return data


class BoostUblasComparisons:
    """Compatible with the existing NativeComparisons adapter protocol."""

    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.boost_ublas_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.boost_ublas_supported.restype = ct.c_int
        self.library.boost_ublas_output_count.argtypes = (
            ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int)
        self.library.boost_ublas_output_count.restype = ct.c_int
        self.library.boost_ublas_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int)
        self.library.boost_ublas_evaluate.restype = ct.c_int
        self.library.boost_ublas_error.argtypes = ()
        self.library.boost_ublas_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if not isinstance(api, str) or "." not in api:
            return {}
        domain, name = api.split(".", 1)
        if (domain == "linalg" and name in DENSE_NAMES):
            prepare = _dense_input
        elif (domain == "sparse" and name in SPARSE_NAMES):
            prepare = _sparse_input
        else:
            return {}
        try:
            data = prepare(name, case.get("args", {}))
            if len(data) > 100_000:
                return {}
        except (TypeError, ValueError, OverflowError, KeyError):
            return {}
        encoded = api.encode("ascii")
        values = (ct.c_double * len(data))(*data)
        if self.library.boost_ublas_supported(0, encoded) != 2:
            return {}
        count = self.library.boost_ublas_output_count(encoded, values, len(data))
        if count <= 0 or count > 40_000:
            return {}
        output = (ct.c_double * count)()
        status = self.library.boost_ublas_evaluate(
            0, encoded, values, len(data), output, count)
        if status:
            error = self.library.boost_ublas_error()
            return {LIBRARY: {"error": (error or b"unknown uBLAS error").decode(
                "utf-8", errors="replace")}}
        result = list(output)
        if not all(math.isfinite(v) for v in result):
            return {LIBRARY: {"error": "nonfinite uBLAS output"}}
        if domain == "sparse":
            # Sci sparse references serialize [rows,cols,data...] and tag
            # each scalar as (0,value) when finite.
            result = [part for value in result for part in (0., value)]
        return {LIBRARY: {"values": result}}
