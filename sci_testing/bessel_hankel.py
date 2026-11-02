"""Independent Boost.Math/GSL real-argument complex Hankel comparisons.

Build ``native/bessel_hankel.cpp`` into the native shared object (or a
standalone .so) and use :class:`BesselHankelComparisons` with its CDLL.
Only the six scalar APIs in ``SUPPORTED_NAMES`` are exposed. SciPy's Hankel
functions also accept complex arguments; those require a separate oracle.

The native functions evaluate public real Bessel J/Y and assemble the
corresponding complex value. Successful output always contains two scalars,
``[real, imaginary]``, matching the saved Sci special test format.
"""

import ctypes as ct


VERSIONS = ("Boost.Math", "GSL")
SUPPORTED_NAMES = frozenset((
    "hankel1", "hankel1e", "hankel2", "hankel2e", "h1vp", "h2vp",
))


class BesselHankelComparisons:
    """Call the dedicated bessel_hankel_* ABI for scalar real inputs."""

    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.bessel_hankel_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.bessel_hankel_supported.restype = ct.c_int
        self.library.bessel_hankel_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
        )
        self.library.bessel_hankel_evaluate.restype = ct.c_int
        self.library.bessel_hankel_error.argtypes = ()
        self.library.bessel_hankel_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("special.") or api.count(".") != 1:
            return {}
        name = api.rsplit(".", 1)[-1]
        if name not in SUPPORTED_NAMES:
            return {}
        encoded = name.encode("ascii")
        arguments = case.get("args", {})
        results = {}
        for index, version in enumerate(VERSIONS):
            arity = self.library.bessel_hankel_supported(index, encoded)
            if not arity or any(f"p{i}" not in arguments for i in range(arity)):
                continue
            try:
                # Only scalar real-argument cases fit this C ABI. Complex x
                # (represented by [real,imag]) is deliberately unsupported.
                values = (ct.c_double * arity)(
                    *(arguments[f"p{i}"] for i in range(arity)))
            except (TypeError, ValueError, OverflowError):
                continue
            output = (ct.c_double * 2)()
            status = self.library.bessel_hankel_evaluate(
                index, encoded, values, arity, output, 2)
            if status:
                message = self.library.bessel_hankel_error()
                results[version] = {"error": (message or b"Hankel evaluation failed")
                                    .decode("utf-8", errors="replace")}
            else:
                results[version] = {"values": list(output)}
        return results
