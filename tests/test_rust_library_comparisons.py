"""RustFFT/nalgebra/statrs use public independent Rust numerical APIs."""

import csv
import gzip
import json
import math
from pathlib import Path
import unittest

import numpy as np
from scipy import stats

from sci_testing.case_csv import read_case_file
from sci_testing.rust_library_comparisons import VERSIONS, evaluate


ROOT = Path(__file__).resolve().parents[1]


def saved_cases(module, name):
    source = ROOT / "cases" / "generated" / module / f"{name}.csv.gz"
    if not source.is_file():
        source = source.with_suffix("")
    opener = gzip.open if source.suffix == ".gz" else open
    with opener(source, mode="rt", newline="") as stream:
        for row in csv.DictReader(stream):
            if module == "fft":
                args = {"z": [[float(x) for x in pair.split()]
                              for pair in row["z"].removeprefix("m:").split(";")]}
            elif module == "stats":
                args = {"x": float(row["x"].removeprefix("s:"))}
            elif "raw" in row:
                args = {"raw": json.loads(row["raw"].removeprefix("j:"))}
            else:
                def mat(encoded):
                    return [[float(x) for x in line.split()]
                            for line in encoded.removeprefix("m:").split(";")]

                def vec(encoded):
                    return [float(x) for x in encoded.removeprefix("v:").split()]

                args = {"a": mat(row["a"])}
                for key in ("u", "v", "b"):
                    if row.get(key):
                        args[key] = mat(row[key]) if row[key].startswith("m:") else vec(row[key])
            yield {"api": f"{module}.{name}", "args": args}


def samples(module, name, count=60):
    cases = list(saved_cases(module, name))
    # Include early boundary profiles and widely spaced saved inputs.
    indices = set(range(min(10, len(cases))))
    indices.update(np.linspace(0, len(cases) - 1, min(count, len(cases)), dtype=int))
    return [cases[int(i)] for i in sorted(indices)]


class RustLibraryComparisonTests(unittest.TestCase):
    def test_raw_fft_saved_cases_numpy_full_shapes(self):
        """Every public RustFFT raw fixture including ND and packed transforms."""
        names = ("rfft", "irfft", "fftn", "ifftn", "fft2", "ifft2",
                 "rfftn", "irfftn", "rfft2", "irfft2")
        total = 0
        for name in names:
            path = ROOT / "cases/generated/fft" / (name + ".csv.gz")
            cases = read_case_file(path, "fft/" + name + ".csv")
            for case in cases:
                raw = case["args"]["raw"]
                result = evaluate(case)
                self.assertEqual(list(result), ["RustFFT"], case["id"])
                actual = np.asarray(result["RustFFT"]["values"])
                re = np.asarray(raw["re"], float)
                im = np.asarray(raw.get("im", [0.] * len(re)), float)
                input_ = re if name in ("rfft", "rfftn", "rfft2") else re + 1j * im
                if name in ("rfft", "irfft"):
                    expected = (np.fft.rfft(input_) if name == "rfft" else
                                np.fft.irfft(input_, n=raw.get("n")))
                else:
                    shape = list(raw["shape"])
                    axes = raw.get("axes", [0, 1])
                    kwargs = {"axes": axes}
                    if name in ("irfftn", "irfft2"):
                        shape[axes[-1]] = shape[axes[-1]] // 2 + 1
                        kwargs["s"] = [raw["shape"][axis] for axis in axes]
                    expected = getattr(np.fft, name)(input_.reshape(shape), **kwargs)
                flat = np.asarray(expected).ravel()
                expected_values = np.column_stack([flat.real, flat.imag]).ravel()
                self.assertEqual(actual.shape, expected_values.shape, case["id"])
                np.testing.assert_allclose(actual, expected_values, rtol=2e-11, atol=2e-11,
                                           err_msg=case["id"])
                total += 1
        self.assertEqual(total, 20044)

    def test_statrs_distribution_public_apis_saved(self):
        """Check each reported public statrs comparison against SciPy distributions."""
        factories = {
            "Normal": lambda p: stats.norm(loc=p[0], scale=p[1]),
            "Beta": lambda p: stats.beta(*p),
            "Gamma": lambda p: stats.gamma(p[0], scale=p[1]),
            "ChiSquared": lambda p: stats.chi2(*p),
            "Chi": lambda p: stats.chi(*p),
            "StudentT": lambda p: stats.t(*p),
            "F": lambda p: stats.f(*p),
            "Exponential": lambda p: stats.expon(scale=p[0]),
            "Cauchy": lambda p: stats.cauchy(loc=p[0], scale=p[1]),
            "Laplace": lambda p: stats.laplace(loc=p[0], scale=p[1]),
            "LogNormal": lambda p: stats.lognorm(s=p[1], scale=math.exp(p[0])),
            "Uniform": lambda p: stats.uniform(loc=p[0], scale=p[1] - p[0]),
            "Weibull": lambda p: stats.weibull_min(p[0], scale=p[1]),
            "Pareto": lambda p: stats.pareto(p[1], scale=p[0]),
            "Gumbel": lambda p: stats.gumbel_r(loc=p[0], scale=p[1]),
            "InverseGamma": lambda p: stats.invgamma(p[0], scale=p[1]),
            "Erlang": lambda p: stats.erlang(p[0], scale=p[1]),
            "Bernoulli": lambda p: stats.bernoulli(*p),
            "Binomial": lambda p: stats.binom(*p),
            "Poisson": lambda p: stats.poisson(*p),
            "Geometric": lambda p: stats.geom(*p),
            "NegativeBinomial": lambda p: stats.nbinom(*p),
        }
        discrete = {"Bernoulli", "Binomial", "Poisson", "Geometric", "NegativeBinomial"}
        supported = 0
        for operation in ("pdf", "logpdf", "cdf", "sf", "ppf"):
            path = ROOT / "cases/generated/stats" / (operation + ".csv.gz")
            for case in read_case_file(path, "stats/" + operation + ".csv"):
                result = evaluate(case)
                if not result:
                    continue
                args = case["args"]
                self.assertEqual(list(result), ["statrs"], case["id"])
                actual = np.asarray(result["statrs"]["values"])
                reference_method = ("pmf" if operation == "pdf" else "logpmf"
                                    if operation == "logpdf" else operation
                                    ) if args["d"] in discrete else operation
                expected = np.atleast_1d(getattr(
                    factories[args["d"]](args["p"]), reference_method
                )(np.asarray(args["x"], float)))
                self.assertEqual(actual.shape, expected.shape, case["id"])
                np.testing.assert_allclose(actual, expected, rtol=2e-7, atol=2e-8,
                                           err_msg=case["id"])
                supported += 1
        self.assertGreater(supported, 8500)
    def test_fft_matches_independent_numpy_fft_saved_shapes(self):
        self.assertEqual(VERSIONS["RustFFT"], "6.4.0")
        for name, operation in (("fft", np.fft.fft), ("ifft", np.fft.ifft)):
            with self.subTest(name=name):
                total = 0
                for case in samples("fft", name, count=120):
                    result = evaluate(case)
                    self.assertEqual(list(result), ["RustFFT"])
                    z = np.array([complex(*pair) for pair in case["args"]["z"]])
                    expected = operation(z)
                    values = result["RustFFT"]["values"]
                    self.assertEqual(len(values), 2 * len(expected))
                    actual = np.array(values[::2]) + 1j * np.array(values[1::2])
                    np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=2e-12)
                    total += 1
                self.assertGreaterEqual(total, 100)

    def test_statrs_normal_distribution_saved_cases(self):
        self.assertEqual(VERSIONS["statrs"], "0.19.1")
        for case in samples("stats", "normal_cdf", count=120):
            result = evaluate(case)
            self.assertEqual(list(result), ["statrs"])
            (actual,) = result["statrs"]["values"]
            x = case["args"]["x"]
            expected = 0.5 * math.erfc(-x / math.sqrt(2.0))
            # statrs 0.19.1 CDF uses an erf approximation whose largest
            # observed saved-case absolute error is 2.44e-11.
            self.assertAlmostEqual(actual, expected, delta=3e-11)
        for case in samples("stats", "normal_pdf", count=120):
            result = evaluate(case)
            self.assertEqual(list(result), ["statrs"])
            (actual,) = result["statrs"]["values"]
            x = case["args"]["x"]
            expected = math.exp(-0.5 * x * x) / math.sqrt(2.0 * math.pi)
            self.assertAlmostEqual(actual, expected, delta=3e-15)
        # Independent quantile checkpoints (two-sided 95% interval and median).
        for p, expected in ((0.025, -1.959963984540054),
                            (0.5, 0.0), (0.975, 1.959963984540054)):
            case = {"api": "stats.normal_ppf", "args": {"x": p}}
            result = evaluate(case)
            self.assertAlmostEqual(result["statrs"]["values"][0], expected, delta=3e-12)
        for case in samples("stats", "normal_ppf", count=120):
            p = case["args"]["x"]
            if not 0 < p < 1:
                self.assertEqual(evaluate(case), {})
                continue
            result = evaluate(case)
            self.assertEqual(list(result), ["statrs"])
            (value,) = result["statrs"]["values"]
            self.assertAlmostEqual(0.5 * math.erfc(-value / math.sqrt(2)), p,
                                   delta=2e-12)

    def test_nalgebra_saved_dense_matrix_cases(self):
        self.assertEqual(VERSIONS["nalgebra"], "0.35.0")
        checked = 0
        for name in ("det", "solve", "inv", "dot", "matvec", "matmul",
                     "cholesky", "cho_solve", "eigvalsh", "svdvals"):
            for case in samples("linalg", name, count=100):
                args = case["args"]
                if "raw" in args:
                    a = np.asarray(args["raw"]["args"][0], dtype=float)
                else:
                    a = np.asarray(args["a"], dtype=float)
                if name == "det":
                    expected = np.array([np.linalg.det(a)])
                elif name == "solve":
                    if np.linalg.cond(a) > 1e8:
                        continue  # Ill-conditioned LU comparisons require MP arithmetic.
                    expected = np.linalg.solve(a, np.asarray(args["b"], dtype=float))
                elif name == "inv":
                    if np.linalg.cond(a) > 1e8:
                        continue
                    expected = np.linalg.inv(a).ravel()
                elif name == "cholesky":
                    lower = args["raw"]["args"][-1]
                    factor = np.linalg.cholesky(a)
                    expected = (factor if lower else factor.T).ravel()
                elif name == "cho_solve":
                    expected = np.linalg.solve(a, args["raw"]["args"][1])
                elif name == "eigvalsh":
                    lower = args["raw"]["args"][-1]
                    expected = np.linalg.eigvalsh(a, UPLO="L" if lower else "U")
                elif name == "svdvals":
                    expected = np.linalg.svd(a, compute_uv=False)
                elif name == "dot":
                    expected = np.array([np.dot(args["u"], args["v"])])
                elif name == "matvec":
                    expected = a @ np.asarray(args["u"], dtype=float)
                else:
                    expected = (a @ np.asarray(args["b"], dtype=float)).ravel()
                result = evaluate(case)
                self.assertEqual(list(result), ["nalgebra"], case)
                actual = np.asarray(result["nalgebra"]["values"])
                self.assertEqual(actual.shape, expected.shape, case)
                np.testing.assert_allclose(actual, expected, rtol=1e-8, atol=3e-11)
                checked += 1
        self.assertGreater(checked, 500)

    def test_unsupported_shapes_and_domains(self):
        invalid = (
            {"api": "fft.fft", "args": {"z": [[1, 2, 3]]}},
            {"api": "fft.fft", "args": {"z": [[float("nan"), 0]]}},
            {"api": "stats.normal_ppf", "args": {"x": 1}},
            {"api": "linalg.det", "args": {"a": [[1, 2], [3]]}},
            {"api": "linalg.solve", "args": {"a": [[1, 0], [0, 1]], "b": [1]}},
            {"api": "linalg.matmul", "args": {"a": [[1, 2]], "b": [[1, 2]]}},
            {"api": "linalg.inv", "args": {"raw": {"args": [[[1., 2.], [2., 4.]]]}}},
            {"api": "linalg.unsupported", "args": {}},
        )
        for case in invalid:
            with self.subTest(case=case):
                self.assertEqual(evaluate(case), {})


if __name__ == "__main__":
    unittest.main()
