"""Public-GSL-derived detrend, median and Wiener signal comparisons.

The native implementation uses GSL linear fitting, statistics, sorting and
sorted-data median primitives. These are *compositions*, not named GSL DSP
equivalents. Only the finite real, saved 1D/2D SciPy signal contracts are
accepted; output is flattened in row-major order.
"""

import ctypes as ct
import math


NAMES = frozenset({"detrend", "medfilt", "medfilt2d", "wiener"})


class GSLSignalStatsComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.gsl_signal_stats_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.gsl_signal_stats_supported.restype = ct.c_int
        self.library.gsl_signal_stats_output_count.argtypes = (
            ct.c_char_p, ct.c_int, ct.c_int)
        self.library.gsl_signal_stats_output_count.restype = ct.c_int
        self.library.gsl_signal_stats_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_double,
            ct.POINTER(ct.c_double), ct.c_int)
        self.library.gsl_signal_stats_evaluate.restype = ct.c_int
        self.library.gsl_signal_stats_error.argtypes = ()
        self.library.gsl_signal_stats_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("signal."):
            return {}
        name = api.rsplit(".", 1)[-1]
        if name not in NAMES:
            return {}
        arguments = case.get("args", {})
        if not isinstance(arguments, dict) or set(arguments) != {"arguments"}:
            return {}
        values = arguments["arguments"]
        if not isinstance(values, (list, tuple)):
            return {}
        try:
            if name == "medfilt2d":
                if len(values) != 2:
                    return {}
                matrix, kernel = values
                if not isinstance(matrix, (list, tuple)) or not matrix:
                    return {}
                if not all(isinstance(row, (list, tuple)) for row in matrix):
                    return {}
                rows = len(matrix)
                cols = len(matrix[0])
                if cols == 0 or any(len(row) != cols for row in matrix):
                    return {}
                data = [float(item) for row in matrix for item in row]
                if not isinstance(kernel, (list, tuple)) or len(kernel) != 2:
                    return {}
                kernel_rows, kernel_cols = kernel
                mode, noise = 0, 0.
            else:
                expected = {"detrend": 2, "medfilt": 2, "wiener": 3}[name]
                if len(values) != expected:
                    return {}
                data = values[0]
                if not isinstance(data, (list, tuple)):
                    return {}
                data = [float(item) for item in data]
                rows, cols = 1, len(data)
                kernel_rows, kernel_cols = 1, 1
                mode, noise = 0, 0.
                if name == "detrend":
                    if not isinstance(values[1], bool):
                        return {}
                    mode = int(values[1])
                else:
                    kernel_cols = values[1]
                    if name == "wiener":
                        if values[2] is not None:
                            mode, noise = 1, float(values[2])
            if not isinstance(kernel_rows, int) or isinstance(kernel_rows, bool):
                return {}
            if not isinstance(kernel_cols, int) or isinstance(kernel_cols, bool):
                return {}
            if not all(math.isfinite(item) for item in data):
                return {"GSL": {"error": "GSL signal requires finite real inputs"}}
            if not math.isfinite(noise):
                return {"GSL": {"error": "GSL Wiener requires finite noise"}}
        except (OverflowError, ValueError, TypeError):
            return {}
        encoded = name.encode("ascii")
        if self.library.gsl_signal_stats_supported(1, encoded) != 1:
            return {}
        count = self.library.gsl_signal_stats_output_count(encoded, rows, cols)
        if count <= 0 or count != len(data):
            return {"GSL": {"error": "invalid GSL signal shape"}}
        inputs = (ct.c_double * count)(*data)
        output = (ct.c_double * count)()
        status = self.library.gsl_signal_stats_evaluate(
            1, encoded, inputs, rows, cols, kernel_rows, kernel_cols,
            mode, noise, output, count)
        if status:
            message = self.library.gsl_signal_stats_error()
            return {"GSL": {"error": (message or b"GSL signal computation failed").decode(
                "utf-8", errors="replace")}}
        if not all(math.isfinite(item) for item in output):
            return {"GSL": {"error": "nonfinite GSL signal result"}}
        return {"GSL": {"values": list(output)}}
