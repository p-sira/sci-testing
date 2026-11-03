"""Public NumPy competitors must reproduce every saved output component."""

import math
import unittest
import warnings
from pathlib import Path

import numpy as np

from sci_testing.case_csv import case_paths, read_case_file
from sci_testing.extra_python_comparisons import evaluate
from sci_testing.registered_cases import reference
from sci_testing.report import competitors

ROOT = Path(__file__).resolve().parents[1]


class ExtraPythonComparisonsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.paths=case_paths(ROOT/"cases/generated")

    def test_public_api_examples_and_complete_shapes(self):
        names = (
            "fft/fft2", "fft/fftfreq", "fft/rfftn", "fft/irfft",
            "fft/fftshift", "signal/detrend", "signal/fftconvolve",
            "signal.windows/hann", "signal.windows/hamming",
            "spatial/cdist", "spatial/cosine", "spatial/distance_matrix",
            "spatial/minkowski", "spatial/is_valid_y", "spatial/is_valid_dm",
            "spatial/num_obs_y", "spatial/num_obs_dm", "spatial/squareform",
            "spatial/geometric_slerp", "linalg/pinv", "linalg/eigvalsh",
            "linalg/cholesky", "linalg/lstsq", "linalg/svdvals",
            "linalg/pinvh", "linalg/qr", "linalg/svd",
            "stats/mean", "stats/variance", "stats/std", "stats/quantile",
            "stats/median_abs_deviation", "stats/rankdata", "stats/iqr",
            "stats/percentileofscore", "stats/trimboth", "stats/gmean",
            "stats/gstd", "stats/hmean", "stats/mode", "stats/entropy",
            "stats/circmean", "stats/circvar", "stats/circstd",
            "stats/energy_distance", "stats/gzscore",
            "integrate/fixed_quad", "integrate/cumulative_trapezoid",
            "integrate/simpson", "differentiate/derivative",
            "differentiate/gradient", "differentiate/hessian",
            "differentiate/jacobian",
        )
        covered=set()
        for name in names:
            path=self.paths.get(name+".csv")
            if path is None:
                continue
            cases=read_case_file(path,name+".csv")
            for position in sorted({0,len(cases)//2,len(cases)-1}):
                case=cases[position]
                with self.subTest(case=case["id"]):
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        result=evaluate(case)
                        if not result:
                            # The adapter exposes only metrics with matching
                            # NumPy mathematical implementations.
                            continue
                        if (case["args"].get("raw") or case["args"].get("explicit")
                                or case["args"].get("arguments")):
                            expected=reference(case)
                        else:
                            expected=competitors(case)["SciPy"]["values"]
                    self.assertEqual(set(result),{"NumPy"})
                    values=np.asarray(result["NumPy"]["values"])
                    self.assertTrue(np.all(np.isfinite(values)))
                    # The same flattened length matters as much as the values:
                    # fft real/imag pairs, decomposition metadata, and ranks.
                    self.assertEqual(values.shape,np.asarray(expected).shape)
                    np.testing.assert_allclose(values,expected,atol=2e-10,rtol=1e-8)
                    covered.add(case["api"])
        self.assertGreaterEqual(len(covered),30)

    def test_manual_independent_algorithms(self):
        cases={
            "cluster.vq.vq":{"obs":[[0.,1.],[3.,4.]],
                             "codebook":[[0.,0.],[3.,5.]]},
            "cluster.vq.whiten":{"obs":[[0.,1.],[2.,3.]]},
            "differentiate.gradient":{"point":[1.,2.],"step":1e-6},
            "differentiate.jacobian":{"point":[1.,2.],"explicit":True},
            "differentiate.hessian":{"point":[1.,2.],"explicit":True},
            "integrate.fixed_quad":{"a":0.,"b":1.},
            "signal.convolve":{"arguments":[[1.,2.],[3.,4.],"full"]},
            "signal.correlate":{"arguments":[[1.,2.],[3.,4.],"full"]},
        }
        expected={
            "cluster.vq.vq":[0.,1.,1.,1.],
            "cluster.vq.whiten":[0.,1.,2.,3.],
            "differentiate.gradient":[8.,11.],
            "differentiate.jacobian":[2.,1.,1.,1.],
            "differentiate.hessian":[2.,3.,3.,4.],
            "integrate.fixed_quad":[math.e-1.],
            "signal.convolve":[3.,10.,8.],
            "signal.correlate":[4.,11.,6.],
        }
        for api,args in cases.items():
            with self.subTest(api=api):
                result=evaluate({"api":api,"args":args})
                self.assertEqual(set(result),{"NumPy"})
                np.testing.assert_allclose(result["NumPy"]["values"],expected[api],rtol=1e-14,atol=1e-14)

    def test_scipy_blas_supplies_explicit_linalg_alternatives(self):
        args={"explicit":True,"a":[[1.,2.],[3.,4.]],
              "b":[[2.,0.],[0.,3.]],"u":[5.,6.],"v":[7.,8.]}
        for api,expected in {
            "linalg.dot":[83.],
            "linalg.matvec":[17.,39.],
            "linalg.matmul":[2.,6.,6.,12.],
        }.items():
            with self.subTest(api=api):
                result=evaluate({"api":api,"args":args})
                self.assertEqual(set(result),{"SciPy"})
                np.testing.assert_allclose(result["SciPy"]["values"],expected)

    def test_only_supported_finite_complete_results(self):
        self.assertEqual(evaluate({"api":"special.gamma","args":{"p0":1.}}),{})
        # The report already provides a NumPy competitor for these base cases.
        self.assertEqual(evaluate({"api":"stats.variance","args":{"values":[1.,2.,3.]}}),{})
        self.assertEqual(evaluate({"api":"signal.convolve","args":{"x":[1.],"h":[2.]}}),{})
        self.assertEqual(evaluate({"api":"signal.find_peaks","args":{"arguments":[[1.,2.]]}}),{})
        self.assertEqual(evaluate({"api":"linalg.inv","args":{"raw":{"api":"inv","args":[[[0.,0.],[0.,0.]]]},"suite_name":"linalg"}}),{})
        self.assertEqual(evaluate({"api":"spatial.cosine","args":{"u":[0.,0.],"v":[1.,2.]}}),{})
        self.assertEqual(evaluate({"api":"fft.fht","args":{"raw":{"api":"fht","re":[1.]},"suite_name":"fft"}}),{})


if __name__=="__main__":
    unittest.main()
