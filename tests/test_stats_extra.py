"""SciPy validation of remaining public Boost.Math/GSL statistical adapters."""
import ctypes as ct
import math
from pathlib import Path
import subprocess
import tempfile
import unittest
import warnings

from scipy import stats

from sci_testing.stats_extra import StatsExtraComparisons


class StatsExtraTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.temporary.cleanup)
        root = Path(__file__).resolve().parent.parent
        target = Path(cls.temporary.name) / "stats-extra.so"
        subprocess.run(["/usr/bin/clang++", "-std=c++17", "-O2", "-shared", "-fPIC",
                        str(root / "native/stats_extra.cpp"), "-lgsl", "-lgslcblas",
                        "-o", str(target)], check=True)
        cls.bridge = StatsExtraComparisons(ct.CDLL(str(target)))

    @staticmethod
    def case(name, x, p, y=()):
        return {"api": f"stats.{name}",
                "args": {"x": list(x), "y": list(y), "p": list(p), "g": [], "d": ""}}

    def assert_results(self, case, expected, labels=("Boost.Math", "GSL"),
                       rtol=3e-10, atol=1e-11):
        observed = self.bridge.evaluate(case)
        self.assertEqual(set(observed), set(labels))
        for label, result in observed.items():
            with self.subTest(api=case["api"], library=label, args=case["args"]):
                self.assertNotIn("error", result, result)
                self.assertEqual(len(result["values"]), len(expected))
                for actual, reference in zip(result["values"], expected):
                    self.assertTrue(math.isclose(actual, reference,
                                                rel_tol=rtol, abs_tol=atol),
                                    f"{label}: {actual} != {reference}")

    def test_exact_coverage(self):
        for name in sorted(("trim_mean", "tmean", "tvar", "tstd", "tsem",
                            "tmin", "tmax", "ttest_1samp", "ttest_rel",
                            "ttest_ind", "linregress", "dirichlet_pdf",
                            "dirichlet_logpdf", "multinomial_pmf",
                            "multinomial_logpmf", "multivariate_normal_pdf",
                            "multivariate_normal_logpdf")):
            for lib in (0, 1):
                with self.subTest(name=name, library=lib):
                    expected = int(lib == 1 or name not in {
                        "trim_mean", "tmin", "tmax", "dirichlet_pdf",
                        "dirichlet_logpdf", "multinomial_pmf", "multinomial_logpmf",
                        "multivariate_normal_pdf", "multivariate_normal_logpdf"})
                    self.assertEqual(self.bridge.library.stats_extra_supported(
                        lib, name.encode()), expected)
        self.assertEqual(self.bridge.evaluate(self.case("gmean", [1, 2], [])), {})
        self.assertEqual(self.bridge.evaluate(self.case("median", [1, 2], [])), {})

    def test_gsl_public_trimmed_mean(self):
        for sample in ([1., 4., 10., 100., 101.], [4., 1., 9., 2., 3., 5.],
                       [2.], [-3., -1., 2., 5., 10., 100.]):
            for fraction in (0., .1, .25, .4, .49):
                self.assert_results(self.case("trim_mean", sample, [fraction]),
                                    [stats.trim_mean(sample, fraction)], ("GSL",))

    def test_threshold_and_ddof_conventions(self):
        samples = ([-3., -1., 0., 2., 4., 10.], [1., 2., 3., 4., 5., 100.])
        for x in samples:
            for low, high, flags in ((0., 4., (True, True)),
                                     (-1., 4., (False, True)),
                                     (-3., 10., (True, False))):
                for name in ("tmean", "tvar", "tstd", "tsem"):
                    for correction in (0, 1):
                        parameters = [low, high, *map(float, flags), correction]
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            expected = float(getattr(stats, name)(
                                x, limits=(low, high), inclusive=flags,
                                **({} if name == "tmean" else {"ddof": correction})))
                        if math.isfinite(expected):
                            self.assert_results(self.case(name, x, parameters), [expected])
            for name in ("tmin", "tmax"):
                for bound in (-1., 2., 4.):
                    for inclusive in (False, True):
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            try:
                                expected = float(getattr(stats, name)(
                                    x, bound, inclusive=inclusive))
                            except ValueError:
                                self.assertIn("error", self.bridge.evaluate(
                                    self.case(name, x, [bound, float(inclusive)]))["GSL"])
                                continue
                        if math.isfinite(expected):
                            self.assert_results(self.case(
                                name, x, [bound, float(inclusive)]),
                                [expected], ("GSL",))

    def test_public_t_tests_and_alternative_tails(self):
        first = [1.5, 2., 3., 2.3, 4.1, 1.3, 3.8]
        second = [1.1, 3.2, 4.3, 4.7, 2.8, 2.2, 5.1]
        for alt, word in enumerate(("two-sided", "less", "greater")):
            result = stats.ttest_1samp(first, 2., alternative=word)
            self.assert_results(self.case("ttest_1samp", first, [2., alt]),
                                [result.statistic, result.pvalue, result.df])
            result = stats.ttest_rel(first, second, alternative=word)
            self.assert_results(self.case("ttest_rel", first, [alt], second),
                                [result.statistic, result.pvalue, result.df])
            for equal in (False, True):
                result = stats.ttest_ind(first, second, equal_var=equal, alternative=word)
                self.assert_results(self.case("ttest_ind", first, [float(equal), alt], second),
                                    [result.statistic, result.pvalue, result.df])

    def test_native_errors_are_explicit(self):
        case = self.case("trim_mean", [1., 2.], [-.5])
        self.assertIn("error", self.bridge.evaluate(case)["GSL"])
        case = self.case("tvar", [1.], [0., 2., 1., 1., 1.])
        self.assertTrue(all("error" in v for v in self.bridge.evaluate(case).values()))

    def test_public_linear_regression(self):
        samples = (
            ([1.,2.,3.,4.,5.,6.], [4.,3.,8.,9.,7.,12.]),
            ([1.,2.,3.,4.], [4.,3.,2.,1.]),
            ([1.,1.,2.,3.], [2.,3.,5.,5.]),
        )
        for x,y in samples:
            for alt, word in enumerate(("two-sided", "less", "greater")):
                r = stats.linregress(x,y, alternative=word)
                self.assert_results(self.case("linregress",x,[alt],y),
                                    [r.slope,r.intercept,r.rvalue,r.pvalue,
                                     r.stderr,r.intercept_stderr])
        result = self.bridge.evaluate(self.case(
            "linregress", [1.,2.,3.,4.], [0], [3.,3.,3.,3.]))
        self.assertTrue(all("error" in value for value in result.values()))

    def test_gsl_public_multivariate_densities(self):
        for probabilities, alphas in (
            ([.3,.7],[2.,3.]),
            ([.2,.3,.5],[1.5,2.5,3.]),
            ([.5,.5],[1.,1.]),
        ):
            density = stats.dirichlet(alphas)
            for api, operation in (("dirichlet_pdf","pdf"),
                                   ("dirichlet_logpdf","logpdf")):
                self.assert_results(self.case(api,probabilities,[],alphas),
                                    [float(getattr(density,operation)(probabilities))],
                                    ("GSL",))
        for counts, probabilities in (
            ([1.,2.],[.4,.6]),
            ([2.,1.,3.],[.25,.25,.5]),
            ([500000.,500000.],[.5,.5]),
        ):
            distribution = stats.multinomial(n=int(sum(counts)), p=probabilities)
            for api, operation in (("multinomial_pmf","pmf"),
                                   ("multinomial_logpmf","logpmf")):
                self.assert_results(self.case(api,counts,[],probabilities),
                                    [float(getattr(distribution,operation)(counts))],
                                    ("GSL",),rtol=1e-8,atol=2e-9)

    def test_gsl_multivariate_gaussian_public_density(self):
        for x, center, covariance in (
            ([0.,1.],[1.,2.],[[2.,.3],[.3,1.]]),
            ([1.,-2.,3.],[0.,0.,0.],[[2.,0.,0.],[0.,1.,0.],[0.,0.,3.]]),
            ([1.,2.],[0.,0.],[[1.,.99],[.99,1.]]),
        ):
            distribution = stats.multivariate_normal(mean=center,cov=covariance)
            for api, operation in (("multivariate_normal_pdf","pdf"),
                                   ("multivariate_normal_logpdf","logpdf")):
                case = self.case(api,x,[],center)
                case["args"]["g"] = covariance
                self.assert_results(case,[float(getattr(distribution,operation)(x))],
                                    ("GSL",))


if __name__ == "__main__":
    unittest.main()
