"""Public Boost.Math quadrature, finite differences and autodiff oracles.

These are independent mathematical comparisons for the saved Sci callbacks,
not assertions that Boost shares Sci's fixed-step/adaptive derivative method.
The native code uses only public Boost headers and returns provider Boost.Math.
"""

import ctypes as ct
import math


ARITIES = {
    "integrate.fixed_quad": 2,
    "differentiate.derivative": 2,
    "differentiate.gradient": 3,
    "differentiate.jacobian": 2,
    "differentiate.hessian": 2,
}
OUTPUT_COUNTS = {
    "integrate.fixed_quad": 1,
    "differentiate.derivative": 1,
    "differentiate.gradient": 2,
    "differentiate.jacobian": 4,
    "differentiate.hessian": 4,
}


class BoostNumericalCoreComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.boost_numerical_core_supported.argtypes = (
            ct.c_int, ct.c_char_p)
        self.library.boost_numerical_core_supported.restype = ct.c_int
        self.library.boost_numerical_core_output_count.argtypes = (
            ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int)
        self.library.boost_numerical_core_output_count.restype = ct.c_int
        self.library.boost_numerical_core_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int)
        self.library.boost_numerical_core_evaluate.restype = ct.c_int
        self.library.boost_numerical_core_error.argtypes = ()
        self.library.boost_numerical_core_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if api not in ARITIES:
            return {}
        name = api.rsplit(".", 1)[-1].encode("ascii")
        arity = self.library.boost_numerical_core_supported(0, name)
        if arity != ARITIES[api]:
            return {}
        args = case.get("args", {})
        try:
            if api == "integrate.fixed_quad":
                if set(args) != {"a", "b"}:
                    raise ValueError("expected bounds a and b")
                values = [float(args["a"]), float(args["b"])]
            elif api == "differentiate.derivative":
                if set(args) != {"x", "step"}:
                    raise ValueError("expected x and step")
                values = [float(args["x"]), float(args["step"])]
            elif api == "differentiate.gradient":
                if set(args) != {"point", "step"}:
                    raise ValueError("expected point and step")
                point = args["point"]
                if not isinstance(point, (list, tuple)) or len(point) != 2:
                    raise ValueError("gradient requires a two-element point")
                values = [float(v) for v in point] + [float(args["step"])]
            else:
                if set(args) != {"point", "explicit"} or args["explicit"] is not True:
                    raise ValueError("Jacobian/Hessian require saved explicit callback")
                point = args["point"]
                if not isinstance(point, (list, tuple)) or len(point) != 2:
                    raise ValueError("Jacobian/Hessian require a two-element point")
                values = [float(v) for v in point]
            if not all(math.isfinite(v) for v in values):
                raise ValueError("finite arguments required")
        except (KeyError, TypeError, ValueError, OverflowError) as error:
            return {"Boost.Math": {"error": f"Boost numerical input: {error}"}}
        inputs = (ct.c_double * arity)(*values)
        capacity = self.library.boost_numerical_core_output_count(
            name, inputs, arity)
        if capacity != OUTPUT_COUNTS[api]:
            return {"Boost.Math": {"error": "invalid Boost numerical input domain"}}
        output = (ct.c_double * capacity)()
        result = self.library.boost_numerical_core_evaluate(
            0, name, inputs, arity, output, capacity)
        if result:
            error = self.library.boost_numerical_core_error()
            return {"Boost.Math": {"error": (
                error or b"Boost numerical evaluation failed").decode(
                    "utf-8", errors="replace")}}
        return {"Boost.Math": {"values": list(output)}}
