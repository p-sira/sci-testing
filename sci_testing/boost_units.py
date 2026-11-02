"""Independent unit-safe numeric comparisons backed by PUBLIC Boost.Units.

Companion C++: native/boost_units.cpp. The provider is deliberately labeled
Boost.Units, not Boost.Math: the former supplies the absolute temperature
conversion, SI dimensional quantities, and its exact speed-of-light constant.
Two dynamic text/key APIs (find, unit) are intentionally not claimed. The
value/precision APIs are supported ONLY for the exact speed-of-light CODATA
key where the earlier Boost and modern SciPy numeric values are identical.
"""

import ctypes as ct
import math


PROVIDER = "Boost.Units"
NUMERIC_APIS = frozenset({
    "constants.convert_temperature",
    "constants.lambda2nu",
    "constants.nu2lambda",
})
KEYED_APIS = frozenset({"constants.value", "constants.precision"})
SUPPORTED_APIS = NUMERIC_APIS | KEYED_APIS


class BoostUnitsComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.boost_units_supported.argtypes = (ct.c_char_p,)
        self.library.boost_units_supported.restype = ct.c_int
        self.library.boost_units_evaluate.argtypes = (
            ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
        )
        self.library.boost_units_evaluate.restype = ct.c_int
        self.library.boost_units_error.argtypes = ()
        self.library.boost_units_error.restype = ct.c_char_p
        self.library.boost_units_key_supported.argtypes = (ct.c_char_p, ct.c_char_p)
        self.library.boost_units_key_supported.restype = ct.c_int
        self.library.boost_units_evaluate_key.argtypes = (
            ct.c_char_p, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
        )
        self.library.boost_units_evaluate_key.restype = ct.c_int

    def evaluate(self, case):
        api = case.get("api", "")
        if api not in SUPPORTED_APIS:
            return {}
        name = api.rsplit(".", 1)[1].encode("ascii")
        if api in KEYED_APIS:
            key = case.get("args", {}).get("key")
            if not isinstance(key, str):
                return {}
            encoded = key.encode("utf-8")
            if self.library.boost_units_key_supported(name, encoded) != 1:
                return {}
            result = (ct.c_double * 1)()
            status = self.library.boost_units_evaluate_key(name, encoded, result, 1)
            if status == 1:
                return {}
            if status != 0:
                error = self.library.boost_units_error()
                return {PROVIDER: {"error": (error or b"unknown CODATA error").decode(
                    "utf-8", errors="replace")}}
            return {PROVIDER: {"values": [result[0]]}}
        if self.library.boost_units_supported(name) != 1:
            return {}
        args = case.get("args", {})
        try:
            x = float(args["x"])
        except (TypeError, ValueError, KeyError, OverflowError):
            return {}
        if not math.isfinite(x):
            return {}
        numeric = (ct.c_double * 1)(x)
        result = (ct.c_double * 1)()
        status = self.library.boost_units_evaluate(name, numeric, 1, result, 1)
        if status == 1:
            return {}
        if status != 0:
            error = self.library.boost_units_error()
            return {PROVIDER: {"error": (error or b"unknown Boost.Units error").decode(
                "utf-8", errors="replace")}}
        return {PROVIDER: {"values": [result[0]]}}
