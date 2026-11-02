"""GSL 2.7+ native complex forward/inverse FFTs on saved SciPy 1D cases.

The public gsl_fft_complex_forward and gsl_fft_complex_inverse functions
accept interleaved real/imaginary scalars and support arbitrary positive
transform lengths. The GSL inverse has SciPy's default 1/n normalization.
"""

import ctypes as ct
import math


NAMES = frozenset(("fft", "ifft"))


class GSLFFTSimpleComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        lib = self.library
        lib.gsl_fft_simple_supported.argtypes = (ct.c_int, ct.c_char_p)
        lib.gsl_fft_simple_supported.restype = ct.c_int
        lib.gsl_fft_simple_output_count.argtypes = (ct.c_char_p, ct.c_int)
        lib.gsl_fft_simple_output_count.restype = ct.c_int
        lib.gsl_fft_simple_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int)
        lib.gsl_fft_simple_evaluate.restype = ct.c_int
        lib.gsl_fft_simple_error.argtypes = ()
        lib.gsl_fft_simple_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api")
        if api not in ("fft.fft", "fft.ifft"):
            return {}
        arguments = case.get("args")
        if not isinstance(arguments, dict) or set(arguments) != {"z"}:
            return {}
        z = arguments["z"]
        if not isinstance(z, (list, tuple)) or not 1 <= len(z) <= 65536:
            return {"GSL": {"error": "invalid GSL FFT length"}}
        data = []
        try:
            for pair in z:
                if not isinstance(pair, (list, tuple)) or len(pair) != 2:
                    raise ValueError("expected [real, imaginary] pairs")
                data.extend((float(pair[0]), float(pair[1])))
        except (ValueError, TypeError, OverflowError) as error:
            return {"GSL": {"error": f"invalid GSL FFT input: {error}"}}
        if not all(math.isfinite(value) for value in data):
            return {"GSL": {"error": "FFT requires finite complex input"}}
        encoded = api.encode("ascii")
        if self.library.gsl_fft_simple_supported(1, encoded) != 1:
            return {}
        count = self.library.gsl_fft_simple_output_count(encoded, len(z))
        if count != len(data):
            return {"GSL": {"error": "invalid native GSL FFT output count"}}
        source = (ct.c_double * count)(*data)
        output = (ct.c_double * count)()
        status = self.library.gsl_fft_simple_evaluate(
            1, encoded, source, len(z), output, count)
        if status:
            message = self.library.gsl_fft_simple_error()
            return {"GSL": {"error": (message or b"GSL FFT failed").decode(
                "utf-8", errors="replace")}}
        if not all(math.isfinite(value) for value in output):
            return {"GSL": {"error": "nonfinite GSL FFT output"}}
        return {"GSL": {"values": list(output)}}
