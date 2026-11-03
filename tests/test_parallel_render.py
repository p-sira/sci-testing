"""Parallel accuracy aggregation must publish the same authenticated report."""

import tempfile
import unittest
from pathlib import Path

from sci_testing import report
from sci_testing.case_csv import (
    MeasurementStore,
    directory_digest,
    measurement_files,
    write_cases,
    write_measurements,
)


class TinyCases:
    """Picklable case inventory for a real spawned-worker rendering test."""

    def __init__(self, cases):
        self.cases = cases
        self.apis = sorted({case["api"] for case in cases})

    def function(self, api):
        return [case for case in self.cases if case["api"] == api]

    def module_counts(self, module):
        apis = [api for api in self.apis if api.split(".")[0] == module]
        return len(apis), sum(len(self.function(api)) for api in apis)

    def __len__(self):
        return len(self.cases)


class ParallelRenderTests(unittest.TestCase):
    def test_parallel_and_serial_summaries_are_identical(self):
        cases = TinyCases([
            {"id": "special.erf/test", "api": "special.erf", "group": "sample",
             "args": {"x": 0.5}},
            {"id": "stats.variance/test", "api": "stats.variance", "group": "sample",
             "args": {"values": [1.0, 2.0]}},
            {"id": "stats.variance/test2", "api": "stats.variance", "group": "sample",
             "args": {"values": [2.0, 3.0]}},
        ])
        records = {
            "special.erf/test": {
                "id": "special.erf/test", "reference": ["0.52049987781304653768"],
                "libraries": {"Sci": {"values": [0.52]},
                              "mpmath": {"values": [0.5204998778130465]}},
            },
            "stats.variance/test": {
                "id": "stats.variance/test", "reference": ["0.25"],
                "libraries": {"Sci": {"values": [0.25]},
                              "SciPy": {"values": [0.25]}},
            },
            "stats.variance/test2": {
                "id": "stats.variance/test2", "reference": ["0.25"],
                "libraries": {"Sci": {"values": [0.25]},
                              "SciPy": {"error": "deliberate failure"}},
            },
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = root / "cases/generated"
            measurements = root / "results/measurements"
            write_cases(inputs, cases.cases, signature="parallel-render-test")
            write_measurements(measurements, measurement_files(cases.cases, records))
            originals = directory_digest(measurements)
            meta = {
                "sci_commit": "fixture", "testing_sources_sha256": "fixture",
                "cases_sha256": "fixture", "native_libraries": {},
                "module_libraries": {
                    "special": ["Sci", "mpmath", "GSL"],
                    "stats": ["Sci", "SciPy"],
                },
            }
            serial_meta = meta.copy()
            parallel_meta = meta.copy()
            exhaustive_meta = meta.copy()
            store = MeasurementStore(measurements)
            report.render_cached(cases, store, serial_meta, root / "serial",
                                 root=root, inputs=inputs)
            report.render_cached_parallel(cases, store, parallel_meta,
                                          root / "parallel", jobs=2,
                                          root=root, inputs=inputs)
            report.render_cached_parallel(cases, store, exhaustive_meta,
                                          root / "exhaustive", jobs=2,
                                          root=root, inputs=inputs,
                                          exhaustive=True)
            serial_files = {p.name: p.read_bytes()
                            for p in (root / "serial").iterdir()}
            parallel_files = {p.name: p.read_bytes()
                              for p in (root / "parallel").iterdir()}
            exhaustive_files = {p.name: p.read_bytes()
                                for p in (root / "exhaustive").iterdir()}
            self.assertEqual(serial_files, parallel_files)
            self.assertEqual(serial_files, exhaustive_files)
            self.assertEqual(serial_meta["accuracy_summaries_sha256"],
                             parallel_meta["accuracy_summaries_sha256"])
            self.assertEqual(serial_meta["accuracy_summaries_sha256"],
                             exhaustive_meta["accuracy_summaries_sha256"])
            self.assertIn("deliberate failure", serial_files["stats.md"].decode())
            self.assertEqual(originals, directory_digest(measurements))


if __name__ == "__main__":
    unittest.main()
