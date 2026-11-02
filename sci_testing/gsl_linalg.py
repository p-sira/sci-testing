"""Public GSL 2.7-era LU/BLAS/Cholesky/eigen/SVD comparisons.

Consumes exact Sci saved cases, including historical (a,b,u,v) fields and
the newer suite's encoded raw.args. All output arrays are C-order flattened.
"""
import ctypes as ct
import math


NAMES = frozenset({
    "det", "solve", "inv", "dot", "matvec", "matmul",
    "cholesky", "cho_solve", "eigvalsh", "svdvals",
})


def _matrix(x):
    if not isinstance(x, (list, tuple)) or not x:
        raise ValueError("matrix requires nonempty rows")
    if not all(isinstance(row, (list, tuple)) for row in x):
        raise ValueError("matrix rows must be sequences")
    cols = len(x[0])
    if not cols or any(len(row) != cols for row in x):
        raise ValueError("matrix must be nonempty rectangular")
    return [float(v) for row in x for v in row], len(x), cols


def _vector(x):
    if not isinstance(x, (list, tuple)):
        raise ValueError("vector must be a sequence")
    return [float(v) for v in x]


def _doubles(x):
    return (ct.c_double * len(x))(*x) if x else None


class GslLinalgComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.gsl_linalg_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.gsl_linalg_supported.restype = ct.c_int
        self.library.gsl_linalg_output_count.argtypes = (
            ct.c_char_p, ct.c_int, ct.c_int, ct.c_int, ct.c_int,
            ct.c_int, ct.c_int,
        )
        self.library.gsl_linalg_output_count.restype = ct.c_int
        self.library.gsl_linalg_error.argtypes = ()
        self.library.gsl_linalg_error.restype = ct.c_char_p
        self.library.gsl_linalg_evaluate.argtypes = (
            ct.c_int, ct.c_char_p,
            ct.POINTER(ct.c_double), ct.c_int, ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int, ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.c_int, ct.POINTER(ct.c_double), ct.c_int,
        )
        self.library.gsl_linalg_evaluate.restype = ct.c_int

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("linalg."):
            return {}
        name = api.rsplit(".", 1)[-1]
        if name not in NAMES:
            return {}
        args = case.get("args", {})
        if not isinstance(args, dict):
            return {}
        try:
            a = b = u = v = []
            ar = ac = br = bc = 0
            lower = 1
            raw = args.get("raw")
            if isinstance(raw, dict):
                saved = raw.get("args", [])
                if not isinstance(saved, list):
                    return {}
                if name not in {"inv", "cholesky", "cho_solve",
                                "eigvalsh", "svdvals"}:
                    return {}
                a, ar, ac = _matrix(saved[0])
                if name == "cho_solve":
                    b = _vector(saved[1])
                    br, bc = len(b), 1
                if name in {"cholesky", "cho_solve", "eigvalsh"}:
                    lower_flag = saved[-1]
                    if not isinstance(lower_flag, bool):
                        raise ValueError("lower flag must be boolean")
                    lower = int(lower_flag)
            else:
                if name not in {"det", "solve", "dot", "matvec", "matmul"}:
                    return {}
                a, ar, ac = _matrix(args["a"])
                if name == "solve":
                    b = _vector(args["b"])
                    br, bc = len(b), 1
                if name == "matmul":
                    b, br, bc = _matrix(args["b"])
                if name in {"dot", "matvec"}:
                    u = _vector(args["u"])
                if name == "dot":
                    v = _vector(args["v"])
            if not all(math.isfinite(item) for item in (*a, *b, *u, *v)):
                return {}
        except (TypeError, ValueError, KeyError, IndexError, OverflowError):
            return {}
        encoded = name.encode("ascii")
        if not self.library.gsl_linalg_supported(1, encoded):
            return {}
        amount = self.library.gsl_linalg_output_count(
            encoded, ar, ac, br, bc, len(u), len(v))
        if amount <= 0 or amount > 1000000:
            return {}
        out = (ct.c_double * amount)()
        status = self.library.gsl_linalg_evaluate(
            1, encoded, _doubles(a), ar, ac, _doubles(b), br, bc,
            _doubles(u), len(u), _doubles(v), len(v), lower, out, amount,
        )
        if status == 0:
            return {"GSL": {"values": list(out)}}
        return {"GSL": {"error": self.library.gsl_linalg_error().decode(
            errors="replace")}}
