"""Boost.Math vector-special native comparisons.

Accepts saved Sci cases with vector-valued `p0`, including
`special.logsumexp`, `special.softmax` and `special.log_softmax`.
The C++ implementation uses the installed public
`boost::math::logsumexp(container)` (and its standard shifted transforms).
No Sci private implementations or reference fixtures are imported.
"""

import ctypes as ct


SUPPORTED_NAMES = frozenset(("logsumexp", "softmax", "log_softmax"))


class VectorSpecialComparisons:
    """A standalone/linked ctypes wrapper with a conventional .evaluate(case).

    `library` is a ctypes.CDLL or an object exposing `.library`, such as
    NativeComparisons after the shared-object build includes vector_special.cpp.
    """

    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.vector_special_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.vector_special_supported.restype = ct.c_int
        self.library.vector_special_output_count.argtypes = (ct.c_char_p, ct.c_int)
        self.library.vector_special_output_count.restype = ct.c_int
        self.library.vector_special_error.argtypes = ()
        self.library.vector_special_error.restype = ct.c_char_p
        self.library.vector_special_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
        )
        self.library.vector_special_evaluate.restype = ct.c_int

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("special.") or api.startswith("special.complex."):
            return {}
        name = api.rsplit(".", 1)[-1]
        if name not in SUPPORTED_NAMES:
            return {}
        args = case.get("args", {})
        x = args.get("p0")
        if not isinstance(x, (tuple, list)):
            return {}
        if set(args) != {"p0"}:
            return {}
        encoded = name.encode("ascii")
        if not self.library.vector_special_supported(0, encoded):
            return {}
        try:
            values = [float(v) for v in x]
            arguments = (ct.c_double * len(values))(*values) if values else None
        except (ValueError, TypeError, OverflowError):
            return {}
        size = self.library.vector_special_output_count(encoded, len(values))
        if size <= 0:
            return {"Boost.Math": {"error": "invalid input shape"}}
        output = (ct.c_double * size)()
        status = self.library.vector_special_evaluate(
            0, encoded, arguments, len(values), output, size,
        )
        if status:
            error = self.library.vector_special_error()
            return {"Boost.Math": {"error": (error or b"unknown vector-special error").decode(
                "utf-8", errors="replace")}}
        return {"Boost.Math": {"values": list(output)}}
