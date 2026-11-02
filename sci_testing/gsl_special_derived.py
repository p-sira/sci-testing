"""Public GSL 2.7/2.8 special-function comparisons, direct and derived.

The two Airy zero collections use public GSL zero APIs for both Ai/Bi
and their derivatives; the four arrays are concatenated exactly as
SciPy's saved Sci facade: zeros, derivative zeros, function values at
derivative zeros, derivative values at zeros.

Other functions use mathematically explicit PUBLIC GSL-derived formulas,
not independent implementations substituted under a GSL label:
normal quantiles (erfinv/erfcinv), Gamma quantiles (gdtria), Dawson (erfi),
log(erfc) (erfcx/log_ndtr), Carlson RF/RD (elliprg), and GSL log
factorials (ln_multinomial). Domain errors are returned explicitly.

Native ABI is independent and links alongside native/comparisons.cpp.
"""
import ctypes as ct
import math


SCALARS = frozenset({
    "special.erfinv", "special.erfcinv", "special.erfi",
    "special.erfcx", "special.gdtria", "special.elliprg",
    "special.log_ndtr",
})
VECTORS = frozenset({"special.ai_zeros", "special.bi_zeros"})
COMBINATORICS = frozenset({"special.ln_multinomial"})
SUPPORTED = SCALARS | VECTORS | COMBINATORICS


class GslSpecialDerivedComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.gsl_special_derived_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.gsl_special_derived_supported.restype = ct.c_int
        self.library.gsl_special_derived_output_count.argtypes = (
            ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int)
        self.library.gsl_special_derived_output_count.restype = ct.c_int
        self.library.gsl_special_derived_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int)
        self.library.gsl_special_derived_evaluate.restype = ct.c_int
        self.library.gsl_special_derived_error.argtypes = ()
        self.library.gsl_special_derived_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if api not in SUPPORTED:
            return {}
        name = api.rsplit(".", 1)[-1].encode("ascii")
        arity = self.library.gsl_special_derived_supported(1, name)
        if not arity:
            return {}
        args = case.get("args", {})
        try:
            if api in COMBINATORICS:
                counts = args["p1"]
                if not isinstance(counts, (tuple, list)):
                    raise ValueError("expected a category-count vector")
                if set(args) != {"p0", "p1"}:
                    raise ValueError("expected n and category counts")
                values = [float(args["p0"]), *(float(value) for value in counts)]
            else:
                if set(args) != {f"p{i}" for i in range(arity)}:
                    raise ValueError("missing or extra special arguments")
                values = [float(args[f"p{i}"]) for i in range(arity)]
        except (KeyError, TypeError, ValueError, OverflowError) as exc:
            return {"GSL": {"error": f"invalid GSL special arguments: {exc}"}}
        if len(values) > 100001 or not all(math.isfinite(v) for v in values):
            return {"GSL": {"error": "nonfinite or oversized special input"}}
        inputs = (ct.c_double * len(values))(*values)
        capacity = self.library.gsl_special_derived_output_count(
            name, inputs, len(values))
        if capacity <= 0:
            return {"GSL": {"error": "invalid GSL-derived output shape/domain"}}
        output = (ct.c_double * capacity)()
        status = self.library.gsl_special_derived_evaluate(
            1, name, inputs, len(values), output, capacity)
        if status:
            error = self.library.gsl_special_derived_error()
            return {"GSL": {"error": (error or b"GSL-derived computation failed"
                                      ).decode("utf-8", errors="replace")}}
        return {"GSL": {"values": list(output)}}
