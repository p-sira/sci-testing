"""Reference provenance and selected generation preserve numerical evidence."""

import hashlib
import tempfile
import unittest
from collections.abc import Mapping, Sequence
from pathlib import Path

from sci_testing.reference_generation import (
    _regenerate_batch,  # pyright: ignore[reportPrivateUsage] — direct regression target
    _select,  # pyright: ignore[reportPrivateUsage] — direct regression target
    provenance_for_cases,
    recipe_source,
)
from sci_testing.reference_integrity import (
    ReferenceCase,
    ReferenceMeasurement,
    reference_hashes,
    serialize_provenance,
    verify_provenance,
)


def sample(api: str, **args: object) -> ReferenceCase:
    return {'api': api, 'id': f'{api}/test', 'group': 'test', 'args': args}


class FixtureCases:
    def __init__(self, batches: dict[str, list[ReferenceCase]]) -> None:
        self.apis: list[str] = sorted(batches)
        self.batches: dict[str, list[ReferenceCase]] = batches

    def function(self, api: str) -> list[ReferenceCase]:
        return self.batches[api]


def measurement(case: ReferenceCase, values: list[str],
                libraries: Mapping[str, object]) -> ReferenceMeasurement:
    return {'id': case['id'], 'reference': values, 'libraries': libraries}


class ReferenceGenerationTests(unittest.TestCase):
    def test_provenance_names_actual_recipe_code_paths(self):
        for case, backend, path in (
            (sample('complex.sqrt', z=[.5, 0.]), 'mpmath', 'complex_accuracy'),
            (sample('special.complex.hankel1', special_complex03_case=True),
             'mpmath', 'special_complex03_cases'),
            (sample('core.csr.from_coo', structural03_case=True),
             'SciPy/NumPy', 'structural03_cases'),
            (sample('io.wavfile.read', io_case=True), 'SciPy/NumPy', 'io_cases'),
            (sample('special.erf', p0=.5), 'SciPy/NumPy', 'registered_cases'),
            (sample('differentiate.derivative', x=.5), 'mpmath', 'report.high_precision'),
        ):
            provider, recipe, _ = recipe_source(case)
            self.assertIn(backend, provider)
            self.assertIn(path, recipe)

    def test_mixed_recipes_cannot_hide_historical_source(self):
        batch=[sample('special.erf', p0=.5),
               sample('special.erf', special03_case=True, x=.5)]
        batches=FixtureCases({'special.erf': batch})
        entries=provenance_for_cases(batches)
        self.assertEqual(len(entries['special.erf']),2)
        self.assertEqual(sum(row['cases'] for row in entries['special.erf']),2)
        self.assertTrue(all(row['status']=='historical_source_path_inferred'
                            for row in entries['special.erf']))
        packed=serialize_provenance(entries,inputs_sha256='inputs',
                                    measurements_sha256='results')
        self.assertIn(b'historical_source_path_inferred',packed)

    def test_reference_only_refresh_preserves_all_library_measurements(self):
        case=sample('special.erf', x=.5)
        libraries: dict[str, object] = {'Sci': {'values': [.5]},
                                        'Boost.Math': {'error': 'domain error'},
                                        'GSL': {'values': [.5000000000000001], 'warnings': ['warning']}}
        records = {case['id']: measurement(case, ['0.4'], libraries)}
        class Wolfram:
            def supports(self, _case: Mapping[str, object]) -> bool: return True
            def references(self, _batch: Sequence[Mapping[str, object]]) -> list[list[str]]:
                return [['0.5000000000000000000000000001']]
        refreshed=_regenerate_batch([case],records,'wolfram',Wolfram())
        self.assertEqual(refreshed[case['id']]['reference'],
                         ['0.5000000000000000000000000001'])
        self.assertEqual(refreshed[case['id']]['libraries'],libraries)
        self.assertEqual(records[case['id']]['reference'],['0.4'])

    def test_unsupported_provider_or_shape_fails(self):
        case=sample('special.erf', x=.5)
        record = {case['id']: measurement(case, ['1.0'], {'Sci': {'values': [1.]}})}
        class Unsupported:
            def supports(self, _case: Mapping[str, object]) -> bool: return False
            def references(self, _batch: Sequence[Mapping[str, object]]) -> list[list[str]]:
                raise AssertionError('called')
        class BadShape:
            def supports(self, _case: Mapping[str, object]) -> bool: return True
            def references(self, _batch: Sequence[Mapping[str, object]]) -> list[list[str]]:
                return [['1.0','2.0']]
        with self.assertRaisesRegex(ValueError,'No curated Wolfram'):
            _ = _regenerate_batch([case],record,'wolfram',Unsupported())
        with self.assertRaisesRegex(ValueError,'Invalid reference shape'):
            _ = _regenerate_batch([case],record,'wolfram',BadShape())

    def test_explicit_api_selection(self):
        cases=FixtureCases({api: [] for api in ('special.erf','special.erfc','stats.norm')})
        self.assertEqual(_select(cases,['special.erf']),['special.erf'])
        self.assertEqual(_select(cases,['special']),['special.erf','special.erfc'])
        with self.assertRaisesRegex(ValueError,'No saved API'):
            _ = _select(cases,['special.does_not_exist'])

    def test_provenance_digests_protect_reference_values(self):
        case=sample('special.erf', p0=.5)
        cases=FixtureCases({'special.erf': [case]})
        records = {case['id']: measurement(case, ['0.52049987781304653768'],
                                           {'Sci': {'values': [.5]}})}
        entries=provenance_for_cases(cases)
        hashes=reference_hashes(cases,records)
        content=serialize_provenance(entries,inputs_sha256='inputs',
                                     measurements_sha256='original',hashes=hashes)
        meta={'reference_provenance_sha256':hashlib.sha256(content).hexdigest(),
              'cases_sha256':'inputs','measurements_sha256':'original'}
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            (root/'results').mkdir()
            _ = (root/'results/reference-provenance.json').write_bytes(content)
            verify_provenance(cases,records,meta,root)
            meta['measurements_sha256']='replaced'
            with self.assertRaisesRegex(ValueError,'measurement digest mismatch'):
                verify_provenance(cases,records,meta,root)
            meta['measurements_sha256']='original'
            records[case['id']]['reference']=['0.5']
            with self.assertRaisesRegex(ValueError,'values differ'):
                verify_provenance(cases,records,meta,root)
            records[case['id']]['reference']=['0.52049987781304653768']
            _ = (root/'results/reference-provenance.json').write_bytes(content+b' ')
            with self.assertRaisesRegex(ValueError,'digest mismatch'):
                verify_provenance(cases,records,meta,root)


if __name__ == '__main__':
    _ = unittest.main()
