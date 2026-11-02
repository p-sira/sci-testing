"""Public GSL sparse matrix and sparse BLAS comparisons for saved Sci cases.

Encoded COO entries preserve duplicate coordinates until native GSL coalesces
them by addition (the SciPy COO->CSR rule). Returned dense values are tagged
as (finite-class=0, value), including shape, exactly as saved sparse fixtures.
"""
import ctypes as ct
import math


LIBRARY = "GSL"
NAMES = frozenset({
    "matvec", "matrix_transpose", "transpose", "from_dense",
    "to_dense", "add", "matmul",
})
PAIRS = frozenset({"add", "matmul"})


def _real(value):
    if type(value) not in (int, float):
        raise ValueError("sparse operands must be real scalars")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("sparse operand must be finite")
    return value


def _dimension(value):
    number = _real(value)
    if not 0 < number <= 128 or number != int(number):
        raise ValueError("invalid sparse matrix dimensions")
    return int(number)


def _coo(raw, key, rows, cols):
    triples = raw.get(key, [])
    if not isinstance(triples, list) or len(triples) > 32768:
        raise ValueError("malformed COO entries")
    wire = [float(rows), float(cols), float(len(triples))]
    for entry in triples:
        if not isinstance(entry, list) or len(entry) != 3:
            raise ValueError("invalid COO triple")
        i, j, v = map(_real, entry)
        if (i != int(i) or j != int(j) or i < 0 or i >= rows
                or j < 0 or j >= cols):
            raise ValueError("out-of-bounds COO coordinate")
        wire.extend((i, j, v))
    return wire


class GslSparseComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        lib = self.library
        lib.gsl_sparse_supported.argtypes = (ct.c_int, ct.c_char_p)
        lib.gsl_sparse_supported.restype = ct.c_int
        lib.gsl_sparse_output_count.argtypes = (
            ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
        )
        lib.gsl_sparse_output_count.restype = ct.c_int
        lib.gsl_sparse_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
        )
        lib.gsl_sparse_evaluate.restype = ct.c_int
        lib.gsl_sparse_error.argtypes = ()
        lib.gsl_sparse_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if not isinstance(api, str) or not api.startswith("sparse."):
            return {}
        name = api.removeprefix("sparse.")
        if name not in NAMES:
            return {}
        args = case.get("args", {})
        if not isinstance(args, dict) or set(args) != {
                "fields", "raw", "suite", "suite_name"}:
            return {}
        if args["suite"] != "sparse" or args["suite_name"] != "sparse":
            return {}
        raw = args["raw"]
        if not isinstance(raw, dict) or raw.get("api") != name:
            return {}
        try:
            rows, cols = _dimension(raw["rows"]), _dimension(raw["cols"])
            wire = _coo(raw, "a", rows, cols)
            if name in PAIRS:
                brows, bcols = _dimension(raw["brows"]), _dimension(raw["bcols"])
                if name == "add" and (rows != brows or cols != bcols):
                    raise ValueError("incompatible addition shapes")
                if name == "matmul" and cols != brows:
                    raise ValueError("incompatible sparse multiplication shapes")
                wire.extend(_coo(raw, "b", brows, bcols))
            elif name == "matvec":
                vector = raw["x"]
                if not isinstance(vector, list) or len(vector) != cols:
                    raise ValueError("invalid sparse vector length")
                wire.extend(_real(item) for item in vector)
            if len(wire) > 160000:
                raise ValueError("sparse wire input too large")
        except (KeyError, TypeError, ValueError, OverflowError):
            return {}
        encoded = api.encode("ascii")
        if not self.library.gsl_sparse_supported(1, encoded):
            return {}
        native = (ct.c_double * len(wire))(*wire)
        count = self.library.gsl_sparse_output_count(
            encoded, native, len(wire))
        if count < 3 or count > 32770:
            return {}
        output = (ct.c_double * count)()
        status = self.library.gsl_sparse_evaluate(
            1, encoded, native, len(wire), output, count)
        if status:
            return {LIBRARY: {"error": self.library.gsl_sparse_error().decode(
                "utf-8", errors="replace")}}
        values = list(output)
        if not all(math.isfinite(v) for v in values):
            return {LIBRARY: {"error": "nonfinite sparse output"}}
        return {LIBRARY: {"values": [
            part for value in values for part in (0., value)]}}
