"""Independent GSL comparisons for Sci's complex gamma family.

Compiled companion: native/complex_extra.cpp.  Link it into the existing
native shared object or pass a standalone ctypes.CDLL to this adapter.
No Sci implementation or private numerical backend is imported.

GSL exposes a complex log-|Gamma| and wrapped phase. The native adapter
constructs Gamma and reciprocal Gamma with exponentials, and unwraps the
analytic log-Gamma branch using the Gamma recurrence. Only shared finite
domains are reported; pole errors are not presented as numerical values.
"""

import ctypes as ct


SUPPORTED_APIS = frozenset((
    "special.complex.gamma",
    "special.complex.loggamma",
    "special.complex.rgamma",
))


class ComplexExtraComparisons:
    """Call the separate complex_extra_* C symbols.

    `library` is either a ctypes.CDLL or an object exposing `.library`
    (e.g. NativeComparisons after its shared-object build includes this TU).
    """

    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.complex_extra_supported.argtypes = (ct.c_char_p,)
        self.library.complex_extra_supported.restype = ct.c_int
        self.library.complex_extra_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
        )
        self.library.complex_extra_evaluate.restype = ct.c_int
        self.library.complex_extra_error.argtypes = ()
        self.library.complex_extra_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if api not in SUPPORTED_APIS:
            return {}
        name = api.rsplit(".", 1)[-1].encode("ascii")
        if self.library.complex_extra_supported(name) != 2:
            return {}
        z = case.get("args", {}).get("z")
        if not isinstance(z, (list, tuple)) or len(z) != 2:
            raise ValueError(f"{api}: expected a complex [real, imag] pair")
        try:
            arguments = (ct.c_double * 2)(float(z[0]), float(z[1]))
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError(f"{api}: invalid complex components") from exc
        output = (ct.c_double * 2)()
        status = self.library.complex_extra_evaluate(
            1, name, arguments, 2, output, 2)
        if status:
            message = self.library.complex_extra_error()
            return {"GSL": {"error": (message or b"unknown GSL error").decode(
                "utf-8", errors="replace")}}
        return {"GSL": {"values": list(output)}}
