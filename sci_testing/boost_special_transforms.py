"""Independent, domain-bounded Boost.Math-derived special comparators.

Operations in this module are compositions of public Boost.Math primitives,
NOT claims that Boost publishes identical named functions. Only Boost.Math
provider index 0 is registered by ``native/boost_special_transforms.cpp``.
The original saved Sci real-scalar signatures use ``p0``, ``p1``, ``p2``.
"""

import ctypes as ct
import math


NAMES = frozenset({
    "multigammaln", "boxcox", "boxcox1p", "softplus", "log_expit",
    "nrdtrimn", "nrdtrisd", "ndtri_exp", "log_ndtr", "factorialk", "agm",
})


class BoostSpecialTransformsComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.boost_special_transforms_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.boost_special_transforms_supported.restype = ct.c_int
        self.library.boost_special_transforms_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
        )
        self.library.boost_special_transforms_evaluate.restype = ct.c_int
        self.library.boost_special_transforms_error.argtypes = ()
        self.library.boost_special_transforms_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("special.") or api.startswith("special.complex."):
            return {}
        name = api.rsplit(".", 1)[-1]
        if name not in NAMES:
            return {}
        encoded = name.encode("ascii")
        count = self.library.boost_special_transforms_supported(0, encoded)
        if count <= 0:
            return {}
        params = case.get("args", {})
        try:
            arguments = [float(params[f"p{i}"]) for i in range(count)]
        except (ValueError, TypeError, OverflowError, KeyError):
            return {}
        if not all(math.isfinite(x) for x in arguments):
            return {"Boost.Math": {"error": "nonfinite real arguments are unsupported"}}
        native_args = (ct.c_double * count)(*arguments)
        output = (ct.c_double * 1)()
        status = self.library.boost_special_transforms_evaluate(
            0, encoded, native_args, count, output, 1,
        )
        if status:
            error = self.library.boost_special_transforms_error()
            return {"Boost.Math": {"error": (error or b"native transform error").decode(
                "utf-8", errors="replace")}}
        return {"Boost.Math": {"values": [output[0]]}}
