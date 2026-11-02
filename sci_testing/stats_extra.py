"""Additional Boost.Math/GSL public statistics comparisons.

Works with the saved Sci cases: x/y are sample vectors and p contains each
operation's exact Sci parameter list. Only real, finite input vectors pass
through the native ABI. Missing library APIs are omitted, not failures.

Direct primitives: GSL trimmed mean, selected min/max, pooled two-sample t,
linear fit, Dirichlet, multinomial and multivariate Gaussian densities;
Boost one-sample/paired t and ordinary least squares;
both libraries' mean, sample variance, correlation and t CDF.
Threshold filtering, Sci ddof, Welch, and alternative-tail conversions are
documented mathematical adaptations around those public primitives.
"""
import ctypes as ct
import math


NAMES = frozenset({
    "trim_mean", "tmean", "tvar", "tstd", "tsem", "tmin", "tmax",
    "ttest_1samp", "ttest_rel", "ttest_ind", "linregress",
    "dirichlet_pdf", "dirichlet_logpdf", "multinomial_pmf", "multinomial_logpmf",
    "multivariate_normal_pdf", "multivariate_normal_logpdf",
})
LIBRARIES = ("Boost.Math", "GSL")
PAIRS = frozenset({"ttest_ind", "ttest_rel", "linregress",
                   "dirichlet_pdf", "dirichlet_logpdf",
                   "multinomial_pmf", "multinomial_logpmf",
                   "multivariate_normal_pdf", "multivariate_normal_logpdf"})
GAUSSIAN = frozenset({"multivariate_normal_pdf", "multivariate_normal_logpdf"})


class StatsExtraComparisons:
    def __init__(self, library):
        self.library = getattr(library, "library", library)
        self.library.stats_extra_supported.argtypes = (ct.c_int, ct.c_char_p)
        self.library.stats_extra_supported.restype = ct.c_int
        self.library.stats_extra_output_count.argtypes = (ct.c_char_p,)
        self.library.stats_extra_output_count.restype = ct.c_int
        self.library.stats_extra_evaluate.argtypes = (
            ct.c_int, ct.c_char_p,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int,
            ct.POINTER(ct.c_double), ct.c_int)
        self.library.stats_extra_evaluate.restype = ct.c_int
        self.library.stats_extra_error.argtypes = ()
        self.library.stats_extra_error.restype = ct.c_char_p

    def evaluate(self, case):
        api = case.get("api", "")
        if not api.startswith("stats."):
            return {}
        name = api.rsplit(".", 1)[-1]
        if name not in NAMES:
            return {}
        args = case.get("args", {})
        x, y, p = args.get("x"), args.get("y", []), args.get("p", [])
        if not isinstance(x, (tuple, list)) or not x:
            return {}
        if not isinstance(y, (tuple, list)) or not isinstance(p, (tuple, list)):
            return {}
        if name in PAIRS and not y:
            return {}
        if name not in PAIRS and y:
            return {}
        if name in GAUSSIAN:
            covariance = args.get("g")
            if (not isinstance(covariance, (tuple, list)) or len(covariance) != len(x)
                    or not all(isinstance(row, (tuple, list)) and len(row) == len(x)
                               for row in covariance)):
                return {}
            p = [value for row in covariance for value in row]
        try:
            a = [float(v) for v in x]
            b = [float(v) for v in y]
            parameters = [float(v) for v in p]
        except (TypeError, ValueError, OverflowError):
            return {}
        if not all(math.isfinite(v) for v in (*a, *b, *parameters)):
            return {}
        encoded = name.encode("ascii")
        output_count = self.library.stats_extra_output_count(encoded)
        if not output_count:
            return {}
        native_a = (ct.c_double * len(a))(*a)
        native_b = (ct.c_double * len(b))(*b) if b else None
        native_p = (ct.c_double * len(parameters))(*parameters) if parameters else None
        results = {}
        for lib, label in enumerate(LIBRARIES):
            if not self.library.stats_extra_supported(lib, encoded):
                continue
            output = (ct.c_double * output_count)()
            status = self.library.stats_extra_evaluate(
                lib, encoded, native_a, len(a), native_b, len(b),
                native_p, len(parameters), output, output_count)
            if status:
                results[label] = {
                    "error": self.library.stats_extra_error().decode(errors="replace")}
            else:
                results[label] = {"values": list(output)}
        return results
