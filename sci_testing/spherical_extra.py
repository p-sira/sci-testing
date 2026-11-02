"""Spherical special-function comparisons from public Boost.Math/GSL APIs.

Both libraries provide complex ``sph_harm_y`` as two floating-point values
``[real, imaginary]``; ``sph_legendre_p`` returns one real value. SciPy takes
the polar angle ``theta`` first, followed by the azimuth ``phi``. These
functions accept scalar fixture arguments ``p0..p3`` without transforming
angle order or dropping the Condon-Shortley sign for odd orders.

The native implementation lives in ``native/spherical_extra.cpp`` and exposes
separate ``spherical_extra_*`` symbols for linking with native/comparisons.cpp.
"""

import ctypes as ct


LIBRARIES = ("Boost.Math", "GSL")
APIS = frozenset(("special.sph_harm_y", "special.sph_legendre_p"))


class SphericalExtraComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.spherical_extra_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.spherical_extra_supported.restype = ct.c_int
        self.library.spherical_extra_output_count.argtypes = (ct.c_char_p,)
        self.library.spherical_extra_output_count.restype = ct.c_int
        self.library.spherical_extra_error.argtypes = ()
        self.library.spherical_extra_error.restype = ct.c_char_p
        self.library.spherical_extra_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int)
        self.library.spherical_extra_evaluate.restype = ct.c_int

    def evaluate(self, case):
        api = case.get("api", "")
        if api not in APIS:
            return {}
        name = api.rsplit(".", 1)[-1].encode("ascii")
        args = case.get("args", {})
        output_size = self.library.spherical_extra_output_count(name)
        results = {}
        for index, library in enumerate(LIBRARIES):
            arity = self.library.spherical_extra_supported(index, name)
            if not arity:
                continue
            try:
                values = [float(args[f"p{i}"]) for i in range(arity)]
            except (ValueError, TypeError, OverflowError, KeyError):
                results[library] = {
                    "error": f"expected {arity} numeric arguments p0..p{arity-1}"}
                continue
            native_args = (ct.c_double * arity)(*values)
            output = (ct.c_double * output_size)()
            status = self.library.spherical_extra_evaluate(
                index, name, native_args, arity, output, output_size)
            if status:
                err = self.library.spherical_extra_error()
                results[library] = {"error": (err or b"spherical evaluation failed").decode()}
            else:
                results[library] = {"values": list(output)}
        return results
