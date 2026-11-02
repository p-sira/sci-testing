"""Saved-case and standalone ABI tests of public Boost.Math spatial norms."""

import ctypes as ct
import math
import os
from pathlib import Path
import subprocess
import tempfile
import time
import unittest

import numpy as np
from scipy import spatial
from scipy.spatial import distance

from sci_testing.boost_spatial import BoostSpatialComparisons, NAMES
from sci_testing.case_csv import read_case_file


ROOT = Path(__file__).resolve().parents[1]
COUNTS = {
    "euclidean": 2064,
    "sqeuclidean": 2004,
    "cityblock": 2004,
    "chebyshev": 2003,
    "distance_matrix": 2003,
    "pdist": 2019,
    "cdist": 2019,
    "minkowski": 2004,
    "minkowski_distance": 2004,
    "minkowski_distance_p": 2004,
    "hamming": 2003,
}
PAIRED_METRICS = frozenset((
    "euclidean", "sqeuclidean", "cityblock", "chebyshev",
    "minkowski", "hamming",
))


def parse_vector(value):
    return np.asarray([float(x) for x in value.split(";")], dtype=float) if value else np.array([])


def parse_matrix(value):
    return np.asarray([parse_vector(r) for r in value.split(":")],
                      dtype=float) if value else np.empty((0, 0))


def independent_reference(name, args):
    if "raw" not in args:
        return np.asarray([distance.euclidean(args["u"], args["v"])])
    raw = args["raw"]
    u, v = raw["u"], raw["v"]
    if name in ("pdist", "cdist"):
        metric = raw["param"]
        kwargs = {}
        if metric == "minkowski": kwargs["p"] = float(raw["extra"])
        first = parse_matrix(u)
        return distance.pdist(first, metric=metric, **kwargs) if name == "pdist" else (
            distance.cdist(first, parse_matrix(v), metric=metric, **kwargs))
    if name == "distance_matrix":
        return spatial.distance_matrix(parse_matrix(u), parse_matrix(v),
                                       p=float(raw["param"]))
    first, second = parse_vector(u), parse_vector(v)
    weights = parse_vector(raw["w"]) if raw["w"] else None
    if name == "minkowski_distance" or name == "minkowski_distance_p":
        return np.asarray([getattr(spatial, name)(
            first, second, p=float(raw["param"]))])
    if name == "minkowski":
        return np.asarray([distance.minkowski(
            first, second, p=float(raw["param"]), w=weights)])
    return np.asarray([getattr(distance, name)(first, second, w=weights)])


class BoostSpatialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        workspace = tempfile.TemporaryDirectory()
        cls.addClassCleanup(workspace.cleanup)
        target = Path(workspace.name) / "boost_spatial.so"
        environment = dict(os.environ)
        environment.pop("LD_LIBRARY_PATH", None)
        subprocess.run([
            os.environ.get("CXX", "c++"), "-std=c++17", "-O2",
            "-Wall", "-Wextra", "-Werror", "-Wpedantic", "-shared", "-fPIC",
            str(ROOT / "native/boost_spatial.cpp"), "-o", str(target),
        ], env=environment, check=True)
        cls.adapter = BoostSpatialComparisons(ct.CDLL(str(target)))

    @staticmethod
    def case(name, raw):
        return {"api": "spatial." + name, "args": {"raw": raw}}

    def test_all_8075_four_core_saved_cases(self):
        confirmed = 0
        worst_relative = 0.
        for name in ("euclidean", "sqeuclidean", "cityblock", "chebyshev"):
            directory = ROOT / "cases/generated/spatial"
            cases = read_case_file(directory / (name + ".csv.gz"),
                                   "spatial/" + name + ".csv")
            self.assertEqual(len(cases), COUNTS[name], name)
            for case in cases:
                response = self.adapter.evaluate(case)
                with self.subTest(case=case["id"]):
                    self.assertEqual(set(response), {"Boost.Math"})
                    self.assertNotIn("error", response["Boost.Math"], response)
                    value = response["Boost.Math"]["values"]
                    independent = np.asarray(
                        independent_reference(name, case["args"])).ravel()
                    self.assertEqual(len(value), len(independent))
                    np.testing.assert_allclose(value, independent,
                                               rtol=3e-12, atol=3e-13)
                    worst_relative = max(worst_relative, abs(
                        value[0]-independent[0]) / max(abs(independent[0]), 1e-300))
                confirmed += 1
        self.assertEqual(confirmed, 8075)
        self.assertLess(worst_relative, 5e-12)

    def test_all_additional_saved_cases_with_explicit_partial_coverage(self):
        counts = {}
        start = time.perf_counter()
        for name in (
            "distance_matrix", "pdist", "cdist", "minkowski",
            "minkowski_distance", "minkowski_distance_p", "hamming",
        ):
            directory = ROOT / "cases/generated/spatial"
            cases = read_case_file(directory / (name + ".csv.gz"),
                                   "spatial/" + name + ".csv")
            self.assertEqual(len(cases), COUNTS[name], name)
            ok, unsupported = 0, 0
            for case in cases:
                response = self.adapter.evaluate(case)
                self.assertEqual(set(response), {"Boost.Math"})
                entry = response["Boost.Math"]
                raw = case["args"]["raw"]
                is_unsupported = (
                    name in ("pdist", "cdist") and raw["param"] not in PAIRED_METRICS
                    or name in ("minkowski", "minkowski_distance",
                                "minkowski_distance_p") and
                    (not math.isinf(float(raw["param"])) and
                     not float(raw["param"]).is_integer())
                )
                with self.subTest(case=case["id"]):
                    if is_unsupported:
                        self.assertIn("error", entry, response)
                        unsupported += 1
                        continue
                    self.assertNotIn("error", entry, response)
                    actual = np.asarray(entry["values"])
                    expected = np.asarray(independent_reference(name,case["args"])).ravel()
                    self.assertEqual(actual.size,expected.size)
                    np.testing.assert_allclose(actual,expected,rtol=3e-12,atol=4e-13)
                    ok += 1
            counts[name] = (ok, unsupported)
        self.assertEqual(counts["distance_matrix"], (2003, 0))
        self.assertEqual(counts["hamming"], (2003, 0))
        for name in ("minkowski", "minkowski_distance", "minkowski_distance_p"):
            self.assertEqual(counts[name], (2003, 1))
        self.assertGreaterEqual(counts["pdist"][0], 2005)
        self.assertGreaterEqual(counts["cdist"][0], 2005)
        self.assertLess(time.perf_counter()-start, 45.)

    def test_pairwise_shape_and_metric_ordering(self):
        x = np.array([[0.,1.,2.],[2.,3.,1.],[4.,1.,3.]])
        y = np.array([[1.,2.,3.],[2.,4.,1.]])
        first = ":".join(";".join(map(str,row)) for row in x)
        second = ":".join(";".join(map(str,row)) for row in y)
        for metric in ("euclidean","sqeuclidean","cityblock",
                       "chebyshev","minkowski","hamming"):
            for name in ("pdist","cdist"):
                raw = dict(u=first, v=second if name=="cdist" else "",
                           param=metric, extra="3" if metric=="minkowski" else "",
                           w="")
                response = self.adapter.evaluate(self.case(name, raw))
                actual = np.asarray(response["Boost.Math"]["values"])
                kw={"p":3} if metric=="minkowski" else {}
                expected = (distance.cdist(x,y,metric=metric,**kw).ravel()
                            if name=="cdist" else distance.pdist(x,metric=metric,**kw))
                np.testing.assert_allclose(actual,expected,rtol=3e-13,atol=4e-14)

    def test_guard_unsupported_domains_and_wrong_shapes(self):
        self.assertEqual(set(NAMES), set(COUNTS))
        self.assertEqual(self.adapter.evaluate({"api":"spatial.procrustes","args":{}}),{})
        for name in NAMES:
            self.assertEqual(self.adapter.library.boost_spatial_supported(0,name.encode()),1)
            self.assertEqual(self.adapter.library.boost_spatial_supported(1,name.encode()),0)
        for api,raw in [
            ("minkowski", dict(u="1;2",v="3;4",w="",param="1.5")),
            ("minkowski_distance_p",dict(u="1;2",v="3;4",w="",param="0.5")),
            ("pdist",dict(u="1;2:3;4",v="",w="",param="correlation")),
            ("cdist",dict(u="1;2:3;4",v="1;2",w="",param="jensenshannon")),
            ("distance_matrix",dict(u="1;2",v="3;4",w="",param="1.5")),
            ("cityblock",dict(u="1;2",v="3;4",w="-1;2",param="2")),
            ("sqeuclidean",dict(u="1;2",v="3;4",w="1;2;3",param="2")),
            ("euclidean",dict(u="1;2",v="3",w="",param="2")),
        ]:
            with self.subTest(api=api,raw=raw):
                result=self.adapter.evaluate(self.case(api,raw))
                self.assertIn("error",result["Boost.Math"],result)

    def test_native_abi_capacity_error_does_not_mutate_output(self):
        lib=self.adapter.library
        u=(ct.c_double*3)(1.,2.,3.)
        v=(ct.c_double*3)(2.,4.,1.)
        out=(ct.c_double*2)(9876.,5432.)
        self.assertEqual(lib.boost_spatial_output_count(b"euclidean",1,1,3),1)
        self.assertEqual(lib.boost_spatial_output_count(b"cdist",3,2,3),6)
        self.assertEqual(lib.boost_spatial_output_count(b"pdist",3,0,3),3)
        self.assertEqual(lib.boost_spatial_output_count(b"unknown",1,1,3),0)
        for index, capacity in [(1,1),(0,0)]:
            status=lib.boost_spatial_evaluate(
                index,b"euclidean",u,1,v,1,3,None,0,b"",2.,out,capacity)
            self.assertNotEqual(status,0)
            self.assertEqual(list(out),[9876.,5432.])
            self.assertTrue(lib.boost_spatial_error())
        self.assertEqual(lib.boost_spatial_evaluate(
            0,b"euclidean",u,1,v,1,3,None,0,b"",2.,out,1),0)
        self.assertEqual(out[0],3.)
        self.assertEqual(out[1],5432.)
        self.assertEqual(lib.boost_spatial_error(),b"")


if __name__ == "__main__":
    unittest.main()
