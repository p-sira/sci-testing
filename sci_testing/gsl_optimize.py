"""GSL 2.7 simplex2/BFGS2 comparisons for the saved analytic 2-D objective.

The SAVED callback is f(x,y)=3+(x-1)^2+(y+2)^2. Only the explicit
two-component start-vector recipe is accepted; no arbitrary Python callback
or unsupported minimize method is claimed. GSL is an independent optimizer,
and the output is exactly [x_opt, y_opt, f_opt], not algorithm metadata.
"""
import ctypes as ct
import math


NAMES = frozenset((
    "minimize", "try_minimize",
    "minimize_with_gradient", "try_minimize_with_gradient",
))


class GSLOptimizeComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.gsl_optimize_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.gsl_optimize_supported.restype = ct.c_int
        self.library.gsl_optimize_output_count.argtypes = (
            ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int)
        self.library.gsl_optimize_output_count.restype = ct.c_int
        self.library.gsl_optimize_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int)
        self.library.gsl_optimize_evaluate.restype = ct.c_int
        self.library.gsl_optimize_error.argtypes = ()
        self.library.gsl_optimize_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("optimize."):
            return {}
        name = api.rsplit(".", 1)[-1]
        if name not in NAMES:
            return {}
        arguments = case.get("args", {})
        if not isinstance(arguments, dict) or set(arguments) != {"explicit", "start"}:
            return {}
        if arguments["explicit"] is not True:
            return {}
        start = arguments["start"]
        if not isinstance(start, (list, tuple)) or len(start) != 2:
            return {}
        try:
            values = [float(value) for value in start]
        except (TypeError, ValueError, OverflowError):
            return {}
        if not all(math.isfinite(value) for value in values):
            return {"GSL": {"error": "GSL optimizer requires finite start coordinates"}}
        inputs = (ct.c_double * 2)(*values)
        encoded = name.encode("ascii")
        if self.library.gsl_optimize_supported(1, encoded) != 2:
            return {}
        length = self.library.gsl_optimize_output_count(encoded, inputs, 2)
        if length != 3:
            return {"GSL": {"error": "GSL optimizer requires bounded 2D start"}}
        output = (ct.c_double * length)()
        status = self.library.gsl_optimize_evaluate(
            1, encoded, inputs, 2, output, length)
        if status:
            reason = self.library.gsl_optimize_error()
            return {"GSL": {"error": (reason or b"GSL optimizer failed").decode(
                "utf-8", errors="replace")}}
        if not all(math.isfinite(v) for v in output):
            return {"GSL": {"error": "nonfinite optimized point/cost"}}
        return {"GSL": {"values": list(output)}}
