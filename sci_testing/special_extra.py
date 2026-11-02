"""Native Boost.Math/GSL comparisons for additional Sci special functions.

The C ABI lives in native/special_extra.cpp. This adapter deliberately checks
library availability per API; a missing public equivalent produces no result.
"""

import ctypes as ct


VERSIONS = ("Boost.Math", "GSL")


class SpecialExtraComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.special_extra_error.argtypes = ()
        self.library.special_extra_error.restype = ct.c_char_p
        self.library.special_extra_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.special_extra_supported.restype = ct.c_int
        self.library.special_extra_output_count.argtypes = (
            ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int)
        self.library.special_extra_output_count.restype = ct.c_int
        self.library.special_extra_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int)
        self.library.special_extra_evaluate.restype = ct.c_int

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("special."):
            return {}
        name = api.rsplit(".", 1)[-1].encode("ascii")
        args = case.get("args", {})
        results = {}
        for library, version in enumerate(VERSIONS):
            arity = self.library.special_extra_supported(library, name)
            if not arity:
                continue
            keys = [f"p{i}" for i in range(arity)]
            if any(key not in args for key in keys):
                continue
            try:
                parameters = (ct.c_double * arity)(*(args[key] for key in keys))
            except (TypeError, ValueError, OverflowError):
                # Complex and vector inputs need dedicated representation and
                # cannot be passed to this real-scalar ABI.
                continue
            capacity = self.library.special_extra_output_count(name, parameters, arity)
            if capacity <= 0:
                results[version] = {"error": "invalid output shape"}
                continue
            output = (ct.c_double * capacity)()
            status = self.library.special_extra_evaluate(
                library, name, parameters, arity, output, capacity)
            if status:
                results[version] = {"error": self.library.special_extra_error().decode()}
            else:
                results[version] = {"values": list(output)}
        return results
