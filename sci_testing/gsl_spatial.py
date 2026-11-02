"""Public GSL 2.7 BLAS/statistics-derived scalar Sci spatial distances.

GSL does not expose SciPy-style distance metric names; instead, the native
module uses gsl_blas_dnrm2/ddot/dasum/idamax and gsl_stats_correlation,
with explicit SciPy-compatible weighting and zero-dimensional contracts.
"""
import ctypes as ct
import math


NAMES = frozenset((
    "euclidean", "sqeuclidean", "cityblock", "chebyshev",
    "cosine", "correlation",
))


def _vector(raw):
    if raw == "":
        return []
    return [float(v) for v in raw.split(";")]


class GSLSpatialComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.gsl_spatial_supported.argtypes = (ct.c_int,ct.c_char_p)
        self.library.gsl_spatial_supported.restype = ct.c_int
        self.library.gsl_spatial_output_count.argtypes = (ct.c_char_p,ct.c_int)
        self.library.gsl_spatial_output_count.restype = ct.c_int
        self.library.gsl_spatial_evaluate.argtypes = (
            ct.c_int,ct.c_char_p,ct.POINTER(ct.c_double),
            ct.POINTER(ct.c_double),ct.c_int,ct.POINTER(ct.c_double),
            ct.c_int,ct.POINTER(ct.c_double),ct.c_int)
        self.library.gsl_spatial_evaluate.restype = ct.c_int
        self.library.gsl_spatial_error.argtypes = ()
        self.library.gsl_spatial_error.restype = ct.c_char_p

    def evaluate(self, case):
        api=case.get("api", "")
        if not api.startswith("spatial."):
            return {}
        name=api.rsplit(".",1)[-1]
        if name not in NAMES:
            return {}
        args=case.get("args",{})
        try:
            if not isinstance(args,dict):
                raise ValueError("expected spatial argument mapping")
            if "raw" in args:
                raw=args["raw"]
                if not isinstance(raw,dict):
                    raise ValueError("expected saved spatial raw fields")
                u=_vector(raw["u"])
                v=_vector(raw["v"])
                w=_vector(raw["w"]) if raw.get("w") else []
            elif set(args)=={"u","v"}:
                u=[float(v) for v in args["u"]]
                v=[float(v) for v in args["v"]]
                w=[]
            else:
                raise ValueError("expected paired vectors or saved spatial raw fields")
            if len(u)!=len(v):
                raise ValueError("spatial vectors must have matching dimensions")
            if w and len(w)!=len(u):
                raise ValueError("spatial weight count must match vector dimension")
            if not all(math.isfinite(value) for value in (*u,*v,*w)):
                raise ValueError("spatial coordinates and weights must be finite")
            if any(v<0. for v in w):
                raise ValueError("spatial weights must be nonnegative")
        except (TypeError,ValueError,KeyError,OverflowError) as e:
            return {"GSL":{"error":f"invalid GSL spatial input: {e}"}}
        encoded=name.encode("ascii")
        if self.library.gsl_spatial_supported(1,encoded)!=1:
            return {}
        capacity=self.library.gsl_spatial_output_count(encoded,len(u))
        if capacity!=1:
            return {"GSL":{"error":"unsupported spatial dimension"}}
        array_u=(ct.c_double*max(1,len(u)))(*u)
        array_v=(ct.c_double*max(1,len(v)))(*v)
        array_w=(ct.c_double*len(w))(*w) if w else None
        output=(ct.c_double*capacity)()
        status=self.library.gsl_spatial_evaluate(
            1,encoded,array_u,array_v,len(u),array_w,len(w),output,capacity)
        if status:
            error=self.library.gsl_spatial_error()
            return {"GSL":{"error":(error or b"GSL spatial evaluation failed").decode(
                "utf-8",errors="replace")}}
        return {"GSL":{"values":list(output)}}
