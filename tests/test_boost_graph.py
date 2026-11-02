"""Public Boost.Graph vs independent SciPy sparse.csgraph saved oracles."""

import ctypes as ct
import math
import os
import random
from pathlib import Path
import subprocess
import tempfile
import time
import unittest

import numpy as np
from scipy import sparse
from scipy.sparse import csgraph

from sci_testing.boost_graph import BoostGraphComparisons, NAMES, EXTRA_NAMES
from sci_testing.case_csv import read_case_file
from sci_testing.sparse_cases import reference as scipy_reference


ROOT = Path(__file__).resolve().parents[1]


class BoostGraphTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tmp.cleanup)
        shared = Path(cls.tmp.name) / "boost_graph.so"
        subprocess.run([
            os.environ.get("CXX", "/usr/bin/clang++"), "-std=c++17", "-O2",
            "-shared", "-fPIC", str(ROOT/"native/boost_graph.cpp"),
            "-lstdc++", "-o", str(shared)], check=True,
            env={k:v for k,v in os.environ.items() if k!="LD_LIBRARY_PATH"})
        cls.adapter = BoostGraphComparisons(ct.CDLL(str(shared)))

    @classmethod
    def saved_cases(cls, name):
        return read_case_file(
            ROOT/"cases/generated/sparse.csgraph"/(name+".csv.gz"),
            "sparse.csgraph/"+name+".csv",
        )

    def check_scipy(self, case, *, abs_tolerance=3e-11):
        result = self.adapter.evaluate(case)
        self.assertIn("Boost.Graph", result, case["id"])
        self.assertNotIn("error", result["Boost.Graph"], (case["id"], result))
        values = result["Boost.Graph"]["values"]
        # This independent oracle uses actual SciPy special-path algorithms,
        # including its infinity-tagging and shape prefix encoding.
        expected = scipy_reference(case)
        self.assertEqual(len(values), len(expected), case["id"])
        for i, (got, ref) in enumerate(zip(values, expected)):
            self.assertTrue(
                math.isclose(got, ref, abs_tol=abs_tolerance, rel_tol=2e-11),
                f"{case['id']} output[{i}] Boost={got!r} SciPy={ref!r}")

    def test_all_saved_graph_cases(self):
        total = 0
        for name in sorted(NAMES):
            cases = self.saved_cases(name)
            self.assertIn(len(cases), (2002, 2003), name)
            for case in cases:
                self.check_scipy(case)
                total += 1
        self.assertEqual(total, 5*2003+2*2002)

    def test_all_saved_extended_graph_cases(self):
        total = 0
        for name in sorted(EXTRA_NAMES):
            cases = self.saved_cases(name)
            self.assertEqual(len(cases),
                             2003 if name == "laplacian" else 2002, name)
            for case in cases:
                self.check_scipy(case)
                total += 1
        self.assertEqual(total, 7*2002+2003)

    def test_extended_randomized_small_graphs_independent_scipy(self):
        # The saved 2,001-row sweeps largely reuse a single 4x4 DAG or 2x3
        # bipartite input. Check genuinely different graphs as well, including
        # reciprocal edges with different weights and unsymmetric Laplacians.
        rng = random.Random(12)
        for name in sorted(EXTRA_NAMES):
            template = self.saved_cases(name)[0]
            for index in range(40):
                if name in ("maximum_bipartite_matching", "structural_rank"):
                    nr, nc = 3+index%4, 4+(index//4)%3
                else:
                    nr = nc = 3+index%5
                matrix = [
                    [float(rng.choice((0,0,0,1,2,3,4,5)))
                     for _ in range(nc)]
                    for _ in range(nr)
                ]
                if name == "maximum_flow":
                    options = [0,nr-1]
                elif name == "laplacian":
                    options = [index%2,(index//2)%2,(index//4)%2]
                elif name in ("maximum_bipartite_matching","structural_rank"):
                    options = [index%2]
                else:
                    options = [index%nr,index%2]
                raw = dict(template["args"]["raw"], a=matrix,p=options,x=[])
                case = dict(template,args=dict(template["args"],raw=raw))
                with self.subTest(name=name,index=index):
                    self.check_scipy(case)

    def test_extended_capacity_and_invalid_inputs(self):
        library = self.adapter.library
        p2 = (ct.c_double*2)(0.,1.)
        p3 = (ct.c_double*3)(0.,0.,0.)
        for name in sorted(EXTRA_NAMES):
            key = name.encode()
            self.assertEqual(library.boost_graph_extra_supported(key),1)
            nr,nc=(2,3) if name in (
                "maximum_bipartite_matching","structural_rank") else (4,4)
            cap=library.boost_graph_extra_output_capacity(key,nr,nc)
            self.assertGreater(cap,0)
            out=(ct.c_double*cap)(*([-77.]*cap))
            mat=(ct.c_double*(nr*nc))()
            written=ct.c_int(-1)
            if name == "laplacian":
                params = p3
            elif name in ("breadth_first_order", "breadth_first_tree",
                          "depth_first_order", "depth_first_tree", "maximum_flow"):
                params = p2
            else:
                params = (ct.c_double*1)(0.)
            np = len(params)
            status=library.boost_graph_extra_evaluate(
                key,nr,nc,mat,nr*nc,params,np,out,cap-1,ct.byref(written))
            self.assertEqual(status,-1)
            self.assertIn(b"capacity",library.boost_graph_error())
            self.assertEqual(list(out),[-77.]*cap)
            self.assertEqual(written.value,-1)
            self.assertEqual(library.boost_graph_extra_evaluate(
                key,nr,nc,mat,nr*nc-1,params,np,out,cap,
                ct.byref(written)),-1)
            self.assertEqual(library.boost_graph_extra_evaluate(
                key,nr,nc,None,nr*nc,params,np,out,cap,
                ct.byref(written)),-1)
            self.assertEqual(library.boost_graph_extra_evaluate(
                key,nr,nc,mat,nr*nc,params,np,None,cap,
                ct.byref(written)),-1)
            self.assertEqual(list(out),[-77.]*cap)
        self.assertEqual(library.boost_graph_extra_supported(b"yen"),0)
        self.assertEqual(library.boost_graph_extra_output_capacity(
            b"breadth_first_order",257,257),0)
        # The dense ABI intentionally cannot represent explicit stored zero
        # edges, unlike the COO representation of the first seven functions.
        for name in sorted(EXTRA_NAMES):
            case=self.saved_cases(name)[0]
            raw=dict(case["args"]["raw"],zero_edge=True)
            self.assertEqual(self.adapter.evaluate(
                dict(case,args=dict(case["args"],raw=raw))),{})
        flow=(ct.c_double*4)(0.,1.5,0.,0.)
        output=(ct.c_double*1)(-77.)
        written=ct.c_int(-1)
        self.assertEqual(library.boost_graph_extra_evaluate(
            b"maximum_flow",2,2,flow,4,p2,2,output,1,ct.byref(written)),-1)
        self.assertIn(b"integral capacities",library.boost_graph_error())
        self.assertEqual(list(output),[-77.])

    def test_c_abi_and_output_capacity(self):
        native = self.adapter.library
        p = (ct.c_double * 2)(1, 0)
        indices = (ct.c_int * 2)(0, 2)
        for name in sorted(NAMES):
            encoded = name.encode()
            self.assertEqual(native.boost_graph_supported(encoded), 1)
            count = native.boost_graph_output_count(encoded, 4, p, 2, 2)
            self.assertEqual(count, 5 if name == "connected_components" else
                             1 if name == "minimum_spanning_tree" else 10)
            out = (ct.c_double * 10)(*([-77.]*10))
            self.assertEqual(native.boost_graph_evaluate(
                encoded, 4, None, None, None, 0, p, 2,
                indices, 2, out, count-1), -1)
            self.assertIn(b"capacity",native.boost_graph_error())
            self.assertEqual(list(out), [-77.]*10)
        self.assertEqual(native.boost_graph_supported(b"yen"), 0)
        self.assertEqual(native.boost_graph_output_count(
            b"dijkstra", 0, p, 2, 0), 0)
        self.assertEqual(native.boost_graph_output_count(
            b"dijkstra", 257, p, 2, 0), 0)

    def test_unsupported_and_invalid_graphs(self):
        case=self.saved_cases("dijkstra")[0]
        self.assertEqual(self.adapter.evaluate({"api":"sparse.csgraph.yen",
                                                 "args":case["args"]}),{})
        for replacement in (None, [], [1.,2.]):
            invalid=dict(case)
            invalid["args"]=dict(case["args"], raw=dict(case["args"]["raw"], a=replacement))
            self.assertEqual(self.adapter.evaluate(invalid),{})
        invalid=dict(case)
        invalid["args"]=dict(case["args"],raw=dict(case["args"]["raw"],p=[1,0,1,5]))
        # Still valid Dijkstra min_only+limit.
        self.assertIn("values",self.adapter.evaluate(invalid)["Boost.Graph"])
        invalid["args"]=dict(case["args"],raw=dict(case["args"]["raw"],p=[2,0]))
        self.assertEqual(self.adapter.evaluate(invalid),{})
        # Public BGL Dijkstra does not support negative-weight edges.
        negative=dict(case)
        graph=[[0.,-2.],[0.,0.]]
        negative["args"]=dict(case["args"],raw=dict(case["args"]["raw"],a=graph,x=[0]))
        self.assertIn("error",self.adapter.evaluate(negative)["Boost.Graph"])

    def test_explicit_zero_edge(self):
        case=self.saved_cases("dijkstra")[0]
        raw=dict(case["args"]["raw"])
        raw.update(a=[[0.,0.,0.],[0.,0.,1.],[0.,0.,0.]],
                   p=[1,0],x=[0],zero_edge=True)
        synthetic=dict(case,args=dict(case["args"],raw=raw))
        # CSR explicitly stores the zero-weight (0,1) edge.
        self.check_scipy(synthetic)

    def test_wider_graphs_and_negative_weights(self):
        case=self.saved_cases("shortest_path")[0]
        graphs=[
            [[0.,1.,0.],[0.,0.,2.],[0.,0.,0.]],
            [[0.,2.,8.],[0.,0.,-1.],[0.,0.,0.]],
            [[0.,4.,3.,0.],[1.,0.,0.,0.],[0.,0.,0.,2.],[0.,0.,0.,0.]],
        ]
        for name in sorted(NAMES):
            for i, matrix in enumerate(graphs):
                if name == "dijkstra" and i == 1:
                    continue
                for directed in (False, True):
                    if name == "minimum_spanning_tree" and not directed:
                        continue
                    raw=dict(case["args"]["raw"])
                    raw.update(api=name,a=matrix,p=[int(directed),0],x=[])
                    if name=="connected_components":raw["p"]=[int(directed),1]
                    if name=="minimum_spanning_tree":raw["p"]=[1,0]
                    mock=dict(case,api="sparse.csgraph."+name,
                              args=dict(case["args"],raw=raw))
                    with self.subTest(name=name,graph=i,directed=directed):
                        # A negative undirected edge is a two-edge negative
                        # cycle. BGL rejects it; SciPy raises NegativeCycleError.
                        if i == 1 and not directed and name in (
                                "bellman_ford", "floyd_warshall",
                                "johnson", "shortest_path"):
                            with self.assertRaises(csgraph.NegativeCycleError):
                                scipy_reference(mock)
                            out = self.adapter.evaluate(mock)
                            self.assertIn("error", out["Boost.Graph"])
                            self.assertIn("negative graph cycle",
                                          out["Boost.Graph"]["error"])
                        else:
                            self.check_scipy(mock)

    def test_2000_case_performance(self):
        cases=self.saved_cases("dijkstra")
        start=time.perf_counter()
        for case in cases[:2000]:
            self.assertIn("values",self.adapter.evaluate(case)["Boost.Graph"])
        self.assertLess(time.perf_counter()-start,15.)


if __name__ == "__main__":
    unittest.main()
