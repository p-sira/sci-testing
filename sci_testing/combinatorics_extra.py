"""Boost.Math/GSL comparison adapter for Bernoulli and harmonic-number cases."""

import ctypes as ct


LIBRARIES = ("Boost.Math", "GSL")


class CombinatoricsExtraComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.combinatorics_extra_error.argtypes = ()
        self.library.combinatorics_extra_error.restype = ct.c_char_p
        self.library.combinatorics_extra_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.combinatorics_extra_supported.restype = ct.c_int
        self.library.combinatorics_extra_output_count.argtypes = (
            ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int)
        self.library.combinatorics_extra_output_count.restype = ct.c_int
        self.library.combinatorics_extra_evaluate.argtypes = (
            ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int)
        self.library.combinatorics_extra_evaluate.restype = ct.c_int

    def evaluate(self, case):
        name = case.get("api", "").rsplit(".", 1)[-1]
        if not case.get("api", "").startswith("special.") or name not in (
            "bernoulli", "harmonic", "gen_harmonic"):
            return {}
        args = case.get("args", {})
        encoded = name.encode("ascii")
        results = {}
        for index, library in enumerate(LIBRARIES):
            arity = self.library.combinatorics_extra_supported(index, encoded)
            if not arity:
                continue
            keys = [f"p{i}" for i in range(arity)]
            if any(key not in args for key in keys):
                continue
            try:
                parameters = (ct.c_double * arity)(*(args[key] for key in keys))
            except (TypeError, ValueError, OverflowError):
                continue
            capacity = self.library.combinatorics_extra_output_count(
                encoded, parameters, arity)
            if capacity <= 0:
                results[library] = {"error": "invalid finite output shape or input domain"}
                continue
            output = (ct.c_double * capacity)()
            status = self.library.combinatorics_extra_evaluate(
                index, encoded, parameters, arity, output, capacity)
            results[library] = (
                {"error": self.library.combinatorics_extra_error().decode()}
                if status else {"values": list(output)})
        return results
