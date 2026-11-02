"""Public Boost.Math and GSL distribution comparison adapters.

Accepts the same typed `d`, `p`, `x` arguments as Sci's stats fixtures.
Every result comes from an actual public library API exposed by the native
bridge. Missing distributions and missing operations are omitted explicitly.
"""
import ctypes as ct
import math


OPERATIONS = frozenset(("pdf", "logpdf", "cdf", "sf", "ppf"))
LIBRARIES = ("Boost.Math", "GSL")


class StatsDistributionComparisons:
    def __init__(self, library):
        # NativeComparisons exposes .library; ctypes.CDLL can also be used
        # directly for isolated native tests.
        self.library = getattr(library, "library", library)
        self.library.stats_distribution_error.argtypes = ()
        self.library.stats_distribution_error.restype = ct.c_char_p
        self.library.stats_distribution_supported.argtypes = (
            ct.c_int, ct.c_char_p, ct.c_char_p)
        self.library.stats_distribution_supported.restype = ct.c_int
        self.library.stats_distribution_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.c_char_p, ct.POINTER(ct.c_double),
            ct.c_int, ct.c_double, ct.POINTER(ct.c_double))
        self.library.stats_distribution_evaluate.restype = ct.c_int

    def evaluate(self, case):
        api = case.get("api", "")
        if not (api.startswith("stats.") or api.startswith("stats.distributions.")):
            return {}
        operation = api.rsplit(".", 1)[-1]
        if operation not in OPERATIONS:
            return {}
        args = case["args"]
        name = args.get("d")
        parameters = args.get("p")
        x = args.get("x")
        if not isinstance(name, str) or not isinstance(parameters, (list, tuple)):
            return {}
        if not isinstance(x, (list, tuple)) or not x:
            return {}
        name_bytes = name.encode("ascii")
        op_bytes = operation.encode("ascii")
        results = {}
        for lib_index, label in enumerate(LIBRARIES):
            arity = self.library.stats_distribution_supported(lib_index, name_bytes, op_bytes)
            if not arity:
                continue
            if len(parameters) != arity:
                results[label] = {"error": f"expected {arity} parameters, got {len(parameters)}"}
                continue
            native_params = (ct.c_double * arity)(*parameters)
            values = []
            for observation in x:
                output = ct.c_double()
                status = self.library.stats_distribution_evaluate(
                    lib_index, name_bytes, op_bytes, native_params, arity,
                    float(observation), ct.byref(output))
                if status:
                    results[label] = {"error": self.library.stats_distribution_error().decode()}
                    break
                if not math.isfinite(output.value):
                    # Accuracy metrics require finite values; Sci's existing
                    # native comparator reports nonfinite results as errors.
                    results[label] = {"error": "nonfinite output"}
                    break
                values.append(output.value)
            else:
                results[label] = {"values": values}
        return results
