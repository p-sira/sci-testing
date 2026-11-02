"""Numerically integrated Bessel functions using public Boost.Math/GSL APIs.

The C ABI is defined in native/bessel_integrals.cpp, compiled as part of the
shared object by NativeComparisons after prime integration.
"""
import ctypes as ct

VERSIONS = ("Boost.Math", "GSL")


class BesselIntegralComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.bessel_integrals_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.bessel_integrals_supported.restype = ct.c_int
        self.library.bessel_integrals_output_count.argtypes = (ct.c_char_p,)
        self.library.bessel_integrals_output_count.restype = ct.c_int
        self.library.bessel_integrals_error.argtypes = ()
        self.library.bessel_integrals_error.restype = ct.c_char_p
        self.library.bessel_integrals_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int)
        self.library.bessel_integrals_evaluate.restype = ct.c_int

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("special."):
            return {}
        name = api.rsplit(".", 1)[-1].encode("ascii")
        count = self.library.bessel_integrals_output_count(name)
        if not count:
            return {}
        args = case.get("args", {})
        results = {}
        for library, version in enumerate(VERSIONS):
            arity = self.library.bessel_integrals_supported(library, name)
            keys = [f"p{i}" for i in range(arity)]
            if not arity or any(key not in args for key in keys):
                continue
            try:
                params = (ct.c_double * arity)(*(args[key] for key in keys))
            except (TypeError, ValueError, OverflowError):
                continue
            output = (ct.c_double * count)()
            status = self.library.bessel_integrals_evaluate(
                library, name, params, arity, output, count)
            results[version] = (
                {"values": list(output)} if not status else
                {"error": self.library.bessel_integrals_error().decode()})
        return results
