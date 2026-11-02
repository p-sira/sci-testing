"""Public GSL 2.7+ fixed-quadrature nodes/weights for Sci special.roots_*.

Each saved roots_* signature accepts n and optional family parameters and
returns two arrays of length n: ascending nodes, then the actual weights
of its SciPy orthogonality measure. This adapter DOES NOT support the
separate optional SciPy mu=True output; no saved case requests it.

GSL is used only for genuinely matching public Gaussian quadrature rules.
Shifted, Chebyshev C/S, and probabilists' Hermite rules use explicit
mathematically equivalent transformations of those GSL rules.
"""

import ctypes as ct


ROOT_ARITIES = {
    "roots_legendre": 1,
    "roots_jacobi": 3,
    "roots_gegenbauer": 2,
    "roots_laguerre": 1,
    "roots_genlaguerre": 2,
    "roots_hermite": 1,
    "roots_hermitenorm": 1,
    "roots_chebyt": 1,
    "roots_chebyu": 1,
    "roots_chebyc": 1,
    "roots_chebys": 1,
    "roots_sh_legendre": 1,
    "roots_sh_jacobi": 3,
    "roots_sh_chebyt": 1,
    "roots_sh_chebyu": 1,
}


class GSLQuadratureRootsComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.gsl_quadrature_roots_supported.argtypes = (
            ct.c_int, ct.c_char_p)
        self.library.gsl_quadrature_roots_supported.restype = ct.c_int
        self.library.gsl_quadrature_roots_output_count.argtypes = (
            ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int)
        self.library.gsl_quadrature_roots_output_count.restype = ct.c_int
        self.library.gsl_quadrature_roots_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int)
        self.library.gsl_quadrature_roots_evaluate.restype = ct.c_int
        self.library.gsl_quadrature_roots_error.argtypes = ()
        self.library.gsl_quadrature_roots_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("special.") or api.startswith("special.complex."):
            return {}
        name = api.rsplit(".", 1)[-1]
        required = ROOT_ARITIES.get(name)
        if required is None:
            return {}
        args = case.get("args", {})
        if not isinstance(args, dict) or set(args) != {
                f"p{i}" for i in range(required)}:
            return {}
        try:
            params = (ct.c_double * required)(
                *(float(args[f"p{i}"]) for i in range(required)))
        except (TypeError, ValueError, OverflowError):
            # Complex and vector inputs cannot be represented in this ABI.
            return {}
        encoded = name.encode("ascii")
        if self.library.gsl_quadrature_roots_supported(1, encoded) != required:
            return {}
        count = self.library.gsl_quadrature_roots_output_count(
            encoded, params, required)
        if count <= 0:
            return {"GSL": {"error": "invalid roots degree or quadrature parameters"}}
        output = (ct.c_double * count)()
        status = self.library.gsl_quadrature_roots_evaluate(
            1, encoded, params, required, output, count)
        if status:
            error = self.library.gsl_quadrature_roots_error()
            return {"GSL": {"error": (error or b"GSL roots evaluation failed").decode(
                "utf-8", errors="replace")}}
        return {"GSL": {"values": list(output)}}
