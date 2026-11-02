"""Independent SciPy checks for the Boost.Math/GSL stats distribution bridge."""
import ctypes as ct
import math
from pathlib import Path
import subprocess
import tempfile
import unittest

from scipy import stats

from sci_testing.stats_distributions import StatsDistributionComparisons


class StatsNativeDistributionsTests(unittest.TestCase):
    ALL_DISTRIBUTIONS = {
        "Normal", "Beta", "Gamma", "Erlang", "ChiSquared", "Chi",
        "StudentT", "F", "Exponential", "Cauchy", "Laplace", "LogNormal",
        "Uniform", "Weibull", "Pareto", "InverseGamma", "Gumbel",
        "Triangular", "Bernoulli", "Binomial", "Poisson", "Geometric",
        "NegativeBinomial", "Hypergeometric", "DiscreteUniform", "Levy",
    }
    @classmethod
    def setUpClass(cls):
        cls._tempdir = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls._tempdir.cleanup)
        shared = Path(cls._tempdir.name) / "stats-distributions.so"
        root = Path(__file__).resolve().parent.parent
        subprocess.run(["c++", "-std=c++17", "-O2", "-shared", "-fPIC",
                        str(root / "native/stats_distributions.cpp"),
                        "-lgsl", "-lgslcblas", "-o", str(shared)], check=True)
        cls.native = StatsDistributionComparisons(ct.CDLL(str(shared)))

    def case(self, name, parameters, operation, values):
        return {"api": f"stats.{operation}",
                "args": {"d": name, "p": parameters, "x": values}}

    def test_distributions_against_scipy(self):
        # Explicitly match Sci's generator conventions, including rate/scale
        # changes and counts for geometric and negative-binomial families.
        profiles = {
            "Normal": ([1.5, 2.25], stats.norm(loc=1.5, scale=2.25), 1.0),
            "Beta": ([2.5, 3.5], stats.beta(2.5, 3.5), .42),
            "Gamma": ([2.5, 1.5], stats.gamma(2.5, scale=1.5), 3.1),
            "Erlang": ([3, 2], stats.erlang(3, scale=2), 4.1),
            "ChiSquared": ([3.5], stats.chi2(3.5), 1.2),
            "Chi": ([3], stats.chi(3), 1.2),
            "StudentT": ([4.5], stats.t(4.5), .6),
            "F": ([3, 7], stats.f(3, 7), 1.2),
            "Exponential": ([2], stats.expon(scale=2), 1.2),
            "Cauchy": ([1, 2], stats.cauchy(loc=1, scale=2), .6),
            "Laplace": ([1, 2], stats.laplace(loc=1, scale=2), .6),
            "LogNormal": ([.3, .7], stats.lognorm(.7, scale=math.exp(.3)), 1.2),
            "Uniform": ([-2, 3], stats.uniform(loc=-2, scale=5), 1.2),
            "Weibull": ([1.5, 2], stats.weibull_min(1.5, scale=2), 1.2),
            "Pareto": ([1, 3], stats.pareto(3, scale=1), 1.2),
            "InverseGamma": ([3, 2], stats.invgamma(3, scale=2), 1.2),
            "Gumbel": ([1, 2], stats.gumbel_r(loc=1, scale=2), .6),
            "Triangular": ([-2, 3, 1], stats.triang(.6, loc=-2, scale=5), 1.2),
            "Bernoulli": ([.3], stats.bernoulli(.3), 0.),
            "Binomial": ([10, .3], stats.binom(10, .3), 3.),
            "Poisson": ([2], stats.poisson(2), 3.),
            "Geometric": ([.3], stats.geom(.3), 3.),
            "NegativeBinomial": ([3, .3], stats.nbinom(3, .3), 3.),
            "Hypergeometric": ([30, 10, 8], stats.hypergeom(30, 10, 8), 3.),
        }
        for name, (p, scipy_distribution, x) in profiles.items():
            for operation in ("pdf", "logpdf", "cdf", "sf", "ppf"):
                observation = .42 if operation == "ppf" else x
                scipy_operation = "pmf" if operation == "pdf" and name in {
                    "Bernoulli", "Binomial", "Poisson", "Geometric",
                    "NegativeBinomial", "Hypergeometric"} else (
                    "logpmf" if operation == "logpdf" and name in {
                        "Bernoulli", "Binomial", "Poisson", "Geometric",
                        "NegativeBinomial", "Hypergeometric"} else operation)
                expected = float(getattr(scipy_distribution, scipy_operation)(observation))
                observed = self.native.evaluate(self.case(name, p, operation, [observation]))
                for library, result in observed.items():
                    with self.subTest(name=name, operation=operation, library=library):
                        self.assertNotIn("error", result, result)
                        self.assertTrue(math.isclose(result["values"][0], expected,
                                          rel_tol=3e-9, abs_tol=1e-11),
                                        f"{name} {operation} {library}: {result} vs {expected}")

    def test_unavailable_distributions_are_omitted(self):
        for operation in ("pdf", "logpdf", "cdf", "sf", "ppf"):
            self.assertEqual(self.native.evaluate(self.case("Levy", [0, 1], operation, [1])), {})
        self.assertEqual(self.native.evaluate(
            self.case("Nonexistent", [1], "pdf", [1])), {})
        self.assertEqual(set(self.native.evaluate(
            self.case("Triangular", [-2, 3, 1], "pdf", [1]))), {"Boost.Math"})
        self.assertNotIn("GSL", self.native.evaluate(
            self.case("Poisson", [2], "ppf", [.5])))

    def test_exact_direct_adapter_inventory(self):
        bridge = self.native.library
        boost = self.ALL_DISTRIBUTIONS - {"Chi", "DiscreteUniform", "Levy"}
        gsl = self.ALL_DISTRIBUTIONS - {
            "InverseGamma", "Triangular", "DiscreteUniform", "Levy"}
        gsl_no_quantile = {
            "Bernoulli", "Binomial", "Poisson", "Geometric",
            "NegativeBinomial", "Hypergeometric"}
        for name in sorted(self.ALL_DISTRIBUTIONS):
            for operation in ("pdf", "logpdf", "cdf", "sf", "ppf"):
                for index, expected_set in ((0, boost), (1, gsl)):
                    expected = name in expected_set and not (
                        index == 1 and operation == "ppf" and name in gsl_no_quantile)
                    count = bridge.stats_distribution_supported(
                        index, name.encode(), operation.encode())
                    with self.subTest(name=name, operation=operation, library=index):
                        self.assertEqual(count > 0, expected)

    def test_vector_results_and_support_boundaries(self):
        case = self.case("Normal", [1.5, 2.25], "cdf", [-3., 0., 1.5, 4., math.inf])
        results = self.native.evaluate(case)
        for result in results.values():
            self.assertNotIn("error", result, result)
            for actual, expected in zip(result["values"],
                                        stats.norm(1.5, 2.25).cdf(case["args"]["x"])):
                self.assertAlmostEqual(actual, expected, places=11)
        for name, parameters in (("Binomial", [10, .3]), ("Geometric", [.3]),
                                 ("Hypergeometric", [30, 10, 8])):
            result = self.native.evaluate(self.case(name, parameters, "pdf", [-1., 1.5]))
            for library in result:
                self.assertEqual(result[library]["values"], [0., 0.])

    def test_invalid_parameters_return_native_errors(self):
        result = self.native.evaluate(self.case("Normal", [0, 0], "pdf", [1]))
        self.assertTrue(all("error" in item for item in result.values()))

    def test_observations_outside_support(self):
        cases = [
            ("Beta", [2.5, 3.5], -1., 0., 0., 1.),
            ("Beta", [2.5, 3.5], 2., 0., 1., 0.),
            ("Gamma", [2.5, 1.5], -1., 0., 0., 1.),
            ("ChiSquared", [3.5], -1., 0., 0., 1.),
            ("F", [3., 7.], -1., 0., 0., 1.),
            ("Weibull", [1.5, 2.], -1., 0., 0., 1.),
            ("LogNormal", [.3, .7], -1., 0., 0., 1.),
            ("Pareto", [1., 3.], .5, 0., 0., 1.),
            ("Bernoulli", [.3], 10., 0., 1., 0.),
        ]
        for name, p, x, density, cdf, sf in cases:
            for operation, expected in (("pdf", density), ("cdf", cdf), ("sf", sf)):
                results = self.native.evaluate(self.case(name, p, operation, [x]))
                for library, observed in results.items():
                    with self.subTest(name=name, operation=operation, library=library):
                        self.assertEqual(observed, {"values": [expected]})

    def test_degenerate_discrete_distributions(self):
        for name, params, point in (("Poisson", [0.], 0.),
                                    ("Geometric", [1.], 1.),
                                    ("NegativeBinomial", [3., 1.], 0.)):
            for operation, values, expected in (
                ("pdf", [point, point+1.], [1., 0.]),
                ("logpdf", [point, point+1.], [0., -math.inf]),
                ("cdf", [point-.1, point, point+1.], [0., 1., 1.]),
                ("sf", [point-.1, point, point+1.], [1., 0., 0.]),
                ("ppf", [.1, .5, .9], [point]*3),
            ):
                results = self.native.evaluate(self.case(name, params, operation, values))
                for library, observed in results.items():
                    with self.subTest(name=name, operation=operation, library=library):
                        if operation == "logpdf":
                            self.assertEqual(observed, {"error": "nonfinite output"})
                        else:
                            self.assertEqual(observed, {"values": expected})

    def test_boundary_density_and_inverse_endpoints(self):
        for name, p, x, expected in (
            ("Gamma", [1., 2.], 0., .5),
            ("Erlang", [1., 2.], 0., .5),
            ("ChiSquared", [2.], 0., .5),
        ):
            for library, result in self.native.evaluate(
                    self.case(name, p, "pdf", [x])).items():
                with self.subTest(name=name, library=library):
                    self.assertEqual(result, {"values": [expected]})
        for name, p, q, expected in (
            ("Binomial", [10, 1.], 0., -1.),
            ("Poisson", [2.], 0., -1.),
            ("Geometric", [.3], 0., 0.),
            ("InverseGamma", [3., 2.], 0., 0.),
        ):
            result = self.native.evaluate(self.case(name, p, "ppf", [q]))
            for library, observed in result.items():
                with self.subTest(name=name, library=library):
                    self.assertEqual(observed, {"values": [expected]})

    def test_hypergeometric_cdf_floors_noninteger_observations(self):
        dist = stats.hypergeom(30, 10, 8)
        for operation in ("cdf", "sf"):
            result = self.native.evaluate(
                self.case("Hypergeometric", [30, 10, 8], operation, [.5, 2.8, 8.5]))
            expected = [float(getattr(dist, operation)(math.floor(x)))
                        for x in (.5, 2.8, 8.5)]
            for library, observed in result.items():
                with self.subTest(library=library, operation=operation):
                    self.assertNotIn("error", observed)
                    for actual, reference in zip(observed["values"], expected):
                        self.assertAlmostEqual(actual, reference, places=12)


if __name__ == "__main__":
    unittest.main()
