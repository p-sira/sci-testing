"""Independent SciPy and analytic tests for GSL multimin simplex2 / BFGS2.

The saved explicit two-dimensional quadratic has a unique minimizer
([1,-2], cost=3). All 8,000 saved starts are checked against both that exact
solution and independently run SciPy Nelder-Mead/BFGS convergence.
No tests depend on Sci's Rust implementation or its optimizer internals.
"""

import ctypes as ct
from concurrent.futures import ThreadPoolExecutor
import math
import os
from pathlib import Path
import subprocess
import tempfile
import time
import unittest

import numpy as np
from scipy import optimize

from sci_testing.case_csv import read_case_file
from sci_testing.gsl_optimize import GSLOptimizeComparisons, NAMES


ROOT = Path(__file__).resolve().parents[1]
COUNT_PER_API = 2000
OBJECTIVE = lambda point: 3. + (point[0]-1.)**2 + (point[1]+2.)**2
GRADIENT = lambda point: np.array(
    [2.*(point[0]-1.), 2.*(point[1]+2.)], dtype=float)
OPTIMUM = np.asarray([1., -2., 3.])


def scipy_minimizer(name, start):
    """Same independent SciPy solver options as the saved fixture comparator."""
    derivative = "gradient" in name
    result = optimize.minimize(
        OBJECTIVE, np.asarray(start, dtype=float),
        method="BFGS" if derivative else "Nelder-Mead",
        jac=GRADIENT if derivative else None,
        tol=1e-8, options={"maxiter": 512})
    if not result.success:
        raise AssertionError(f"{name} SciPy did not converge: {result.message}")
    return np.asarray([*result.x, float(result.fun)])


class GSLOptimizeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        workspace = tempfile.TemporaryDirectory()
        cls.addClassCleanup(workspace.cleanup)
        filename = Path(workspace.name) / "gsl-optimize.so"
        environment = dict(os.environ)
        environment.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2",
            "-Wall", "-Wextra", "-Werror", "-Wpedantic", "-shared", "-fPIC",
            str(ROOT / "native/gsl_optimize.cpp"),
            "-lgsl", "-lgslcblas", "-o", str(filename),
        ], env=environment, check=True)
        cls.adapter = GSLOptimizeComparisons(ct.CDLL(str(filename)))

    @staticmethod
    def case(name, start):
        return {
            "api": "optimize." + name,
            "args": {"explicit": True, "start": list(start)},
        }

    def get(self, case):
        result = self.adapter.evaluate(case)
        self.assertEqual(set(result), {"GSL"}, (case, result))
        entry = result["GSL"]
        self.assertNotIn("error", entry, (case, entry))
        output = np.asarray(entry["values"], dtype=float)
        self.assertEqual(output.shape, (3,))
        self.assertTrue(np.all(np.isfinite(output)))
        return output

    def test_all_8000_saved_cases_independent_analytic_and_scipy(self):
        checked = 0
        error_by_api = {}
        start_time = time.perf_counter()
        for name in sorted(NAMES):
            path = ROOT / "cases/generated/optimize" / (name + ".csv")
            cases = read_case_file(path, "optimize/" + name + ".csv")
            self.assertEqual(len(cases), COUNT_PER_API, name)
            max_error = 0.
            max_disagreement = 0.
            for case in cases:
                args = case["args"]
                self.assertEqual(set(args), {"explicit", "start"})
                self.assertIs(args["explicit"], True)
                self.assertEqual(len(args["start"]), 2)
                result = self.get(case)
                independent = scipy_minimizer(name, args["start"])
                # Different converged derivative-free minimizers need not
                # coincide at the last bits. Preserve full output contract.
                with self.subTest(case=case["id"]):
                    np.testing.assert_allclose(result, OPTIMUM,
                                               rtol=0., atol=1e-7)
                    np.testing.assert_allclose(result, independent,
                                               rtol=0., atol=1.5e-7)
                    self.assertAlmostEqual(result[2], OBJECTIVE(result[:2]),
                                           delta=1e-12)
                max_error = max(max_error, float(np.max(abs(result-OPTIMUM))))
                max_disagreement = max(max_disagreement, float(
                    np.max(abs(result-independent))))
                checked += 1
            error_by_api[name] = (max_error, max_disagreement)
        self.assertEqual(checked, 8000)
        self.assertLess(time.perf_counter()-start_time, 70.,
                        "GSL/SciPy optimization case validation exceeded budget")
        for name,(error,disagreement) in error_by_api.items():
            self.assertLess(error, 1e-7, name)
            self.assertLess(disagreement, 1.5e-7, name)
        print(f"independent GSL optimization saved-case coverage: "
              f"{checked}; max(errors vs analytic, vs SciPy)={error_by_api}")

    def test_independent_quadratic_convergence_from_additional_starts(self):
        for name in NAMES:
            for start in (
                (-3.,-1.),(-1.5,2.75),(0.,1.),(3.,3.),
                (1.,-2.),(-.5,0.),(0.,0.),(.3,-3.5),
            ):
                with self.subTest(name=name,start=start):
                    result=self.get(self.case(name,start))
                    np.testing.assert_allclose(result,OPTIMUM,rtol=0.,atol=1e-7)
                    self.assertAlmostEqual(OBJECTIVE(result[:2]),result[2],
                                           delta=1e-13)

    def test_repeated_full_sweeps_use_same_immutable_standalone_library(self):
        """Regression for sporadic integrated objective-cache disagreement.

        The shared project .so is built elsewhere and must not be rebuilt or
        overwritten by this test; setUpClass keeps this unique temporary
        standalone library loaded for the entire suite.
        """
        cases = [
            case for name in sorted(NAMES)
            for case in read_case_file(
                ROOT / "cases/generated/optimize" / (name + ".csv"),
                "optimize/" + name + ".csv"
            )
        ]
        self.assertEqual(len(cases), 8000)
        for sweep in range(5):
            for case in cases:
                response=self.adapter.evaluate(case)
                entry=response.get("GSL",{})
                self.assertNotIn("error",entry,
                                 (sweep,case["id"],entry))
                output=entry["values"]
                self.assertEqual(len(output),3)
                x,y,cost=output
                self.assertLessEqual(abs(x-1.),1e-7)
                self.assertLessEqual(abs(y+2.),1e-7)
                self.assertEqual(cost, OBJECTIVE((x,y)))

    def test_threaded_reentrancy_and_repeated_first_saved_start(self):
        """Multiple ctypes threads must not corrupt GSL solver state."""
        names=sorted(NAMES)
        starts=((-3.,-1.),(1.,-2.),(0.,1.),(3.,3.),
                (-.5,0.),(.3,-3.5),(0.,0.))

        def workload(seed):
            for iteration in range(1500):
                # Every thread repeatedly exercises the reported failing
                # start, but also interleaves simplex and BFGS2 calls.
                name = names[(seed+iteration) % len(names)]
                start = starts[(seed+iteration*3) % len(starts)]
                if iteration % 7 == 0:
                    name, start="minimize",(-3.,-1.)
                response=self.adapter.evaluate(self.case(name,start))
                entry=response.get("GSL",{})
                if "error" in entry:
                    raise AssertionError((seed,iteration,name,start,entry))
                x,y,cost=entry["values"]
                if (abs(x-1.)>1e-7 or abs(y+2.)>1e-7 or
                        cost!=OBJECTIVE((x,y))):
                    raise AssertionError((seed,iteration,name,start,entry))
            return 1500

        with ThreadPoolExecutor(max_workers=12) as executor:
            self.assertEqual(sum(executor.map(workload,range(12))),18000)

    def test_no_fake_provider_or_arbitrary_callback_support(self):
        lib=self.adapter.library
        self.assertEqual(set(NAMES), {
            "minimize","try_minimize",
            "minimize_with_gradient","try_minimize_with_gradient",
        })
        for name in NAMES:
            encoded=name.encode("ascii")
            self.assertEqual(lib.gsl_optimize_supported(1,encoded),2)
            self.assertEqual(lib.gsl_optimize_supported(0,encoded),0)
        for case in (
            self.case("rosen", (2,3)),
            {"api":"optimize.minimize","args":{"start":[0.,0.]}},
            {"api":"optimize.minimize","args":{"explicit":False,"start":[0.,0.]}},
            {"api":"optimize.minimize","args":{"explicit":True,"start":[0.,0.,1.]}},
            {"api":"optimize.minimize","args":{"explicit":True,"start":0.}},
            {"api":"optimize.minimize","args":{"explicit":True,"start":["bad",3.]}},
            {"api":"optimize.minimize","args":{"explicit":True,"start":[[1.],3.]}},
        ):
            self.assertEqual(self.adapter.evaluate(case),{})

    def test_native_abi_invalid_inputs_leave_output_untouched(self):
        lib=self.adapter.library
        for name in NAMES:
            encoded=name.encode("ascii")
            valid=(ct.c_double*2)(-3.,-1.)
            out=(ct.c_double*4)(91.,92.,93.,94.)
            self.assertEqual(lib.gsl_optimize_output_count(encoded,valid,2),3)
            self.assertEqual(lib.gsl_optimize_output_count(encoded,valid,1),0)
            invalids=(
                (0,valid,2,3),(1,valid,1,3),
                (1,valid,2,2),(1,None,2,3),
                (1,(ct.c_double*2)(math.nan,0.),2,3),
                (1,(ct.c_double*2)(0.,math.inf),2,3),
                (1,(ct.c_double*2)(1001.,0.),2,3),
            )
            for provider,inputs,arity,capacity in invalids:
                with self.subTest(api=name,provider=provider,
                                  arity=arity,capacity=capacity):
                    status=lib.gsl_optimize_evaluate(
                        provider,encoded,inputs,arity,out,capacity)
                    self.assertNotEqual(status,0)
                    self.assertEqual(list(out),[91.,92.,93.,94.])
                    self.assertTrue(lib.gsl_optimize_error())
            status=lib.gsl_optimize_evaluate(1,encoded,valid,2,out,3)
            self.assertEqual(status,0)
            self.assertEqual(out[3],94.)
            self.assertEqual(lib.gsl_optimize_error(),b"")
            np.testing.assert_allclose(list(out[:3]),OPTIMUM,rtol=0.,atol=1e-7)
        self.assertEqual(lib.gsl_optimize_supported(1,b"rosen"),0)
        self.assertEqual(lib.gsl_optimize_output_count(b"rosen",
                          (ct.c_double*2)(0.,0.),2),0)

    def test_python_guards_on_nonfinite_or_out_of_domain_start(self):
        for name in NAMES:
            for start in ([math.nan,1.],[math.inf,1.],[1.,-math.inf],
                          [1001.,0.],[-1001.,2.]):
                result=self.adapter.evaluate(self.case(name,start))
                with self.subTest(name=name,start=start):
                    self.assertEqual(set(result),{"GSL"})
                    self.assertIn("error",result["GSL"])


if __name__=="__main__":
    unittest.main()
