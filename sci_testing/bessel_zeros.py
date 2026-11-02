"""Independent real Bessel zero collections via public Boost.Math routines.

The companion native/bessel_zeros.cpp exports a separate C ABI.
Supported outputs use Sci's standard flattened real vector convention:

* jnp_zeros, ynp_zeros: nt positive zeros.
* y0_zeros, y1_zeros, y1p_zeros: 2*nt, first roots,
  then companion real values (SciPy returns these in complex dtype with
  zero imaginary part for the default real-root mode).
* jnyn_zeros: 4*nt, Jn, Jn', Yn, Yn' zero sequences.
* jnjnp_zeros: 4*nt, the globally sorted merged root sequence
  followed by orders, serial numbers, derivative flags (NOT four
  independent root sequences); first entry is (0,0,0,1).
* ai_zeros, bi_zeros: 4*nt, negative real zeros, derivative zeros,
  function values at derivative zeros, derivative values at zeros.

Existing public Boost Jn/Yn root adapters in special_extra.cpp are not
duplicated. The optional complex-zero mode of Y0/Y1 is not supported.
"""

import ctypes as ct


SUPPORTED = frozenset({
    "special.jnp_zeros", "special.ynp_zeros", "special.y0_zeros",
    "special.y1_zeros", "special.y1p_zeros", "special.jnjnp_zeros",
    "special.jnyn_zeros", "special.ai_zeros", "special.bi_zeros",
})


class BesselZerosComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.bessel_zeros_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.bessel_zeros_supported.restype = ct.c_int
        self.library.bessel_zeros_output_count.argtypes = (
            ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int)
        self.library.bessel_zeros_output_count.restype = ct.c_int
        self.library.bessel_zeros_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int)
        self.library.bessel_zeros_evaluate.restype = ct.c_int
        self.library.bessel_zeros_error.argtypes = ()
        self.library.bessel_zeros_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if api not in SUPPORTED:
            return {}
        name = api.rsplit(".", 1)[-1].encode("ascii")
        count = self.library.bessel_zeros_supported(0, name)
        if not count:
            return {}
        args = case.get("args", {})
        try:
            inputs = (ct.c_double * count)(
                *(float(args[f"p{i}"]) for i in range(count)))
        except (KeyError, TypeError, ValueError, OverflowError):
            return {"Boost.Math": {"error": f"expected {count} numeric Bessel zero arguments"}}
        capacity = self.library.bessel_zeros_output_count(name, inputs, count)
        if capacity <= 0:
            return {"Boost.Math": {"error": "invalid Bessel zero order or count"}}
        output = (ct.c_double * capacity)()
        status = self.library.bessel_zeros_evaluate(
            0, name, inputs, count, output, capacity)
        if status:
            reason = self.library.bessel_zeros_error()
            return {"Boost.Math": {
                "error": (reason or b"Bessel zero evaluation failed").decode(
                    "utf-8", errors="replace")}}
        return {"Boost.Math": {"values": list(output)}}
