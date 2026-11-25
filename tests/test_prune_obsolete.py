"""Safe dry-run and reviewed staging for obsolete case snapshots."""

# These tests deliberately validate private snapshot integrity helpers.
# pyright: reportPrivateUsage=false

import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from typing import cast
from unittest.mock import patch

from sci_testing import prune_obsolete as pruning


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


class _Fixture:
    """Small fully authenticated snapshot with one unresolved measured alias."""

    def __init__(self, root: Path) -> None:
        self.root: Path = root
        self._write('api_inventory.csv',
                    '''api,canonical,source,line,status,evidence,kind
special.old,,legacy.rs,,unresolved_legacy,unmatched,legacy_alias
special.nomatch,,legacy.rs,,unresolved_legacy,unmatched,legacy_alias
special.renamed,special.renamed,modern.rs,,legacy_recipe_pending,retained,legacy_alias
special.other,special.other,modern.rs,,untested,modern,public_export
''')
        self._write('cases/recipes/special/old.csv',
                    'id,group,p0\nspecial.old/linspace-0,edge,s:0.5\n')
        self._write('cases/generated/special/old.csv',
                    '''id,group,p0
special.old/linspace-0,linspace,s:0.5
special.old/random-0,random,s:0.75
''')
        self._write('cases/generated/special/renamed.csv',
                    'id,group,p0\nspecial.renamed/linspace-0,linspace,s:0.5\n')
        self._write('results/measurements/special/old.csv',
                    '''id,library,reference,values,error
special.old/linspace-0,Sci,0.5,0.5,
special.old/random-0,Sci,0.75,0.75,
''')
        self._write('results/measurements/special/renamed.csv',
                    '''id,library,reference,values,error
special.renamed/linspace-0,Sci,0.5,0.5,
''')
        self._write('docs/accuracy/README.md', 'original report\n')
        self._write('results/accuracy-summaries.json', '{}\n')
        self._write('results/coverage.csv',
                    'api,canonical,cases\nspecial.old,,2\nspecial.renamed,special.renamed,1\n')
        inputs = pruning.case_paths(root / 'cases/generated')
        entries = {name: {'count': 2 if name.endswith('old.csv') else 1,
                          'sha256': pruning._sha256(path)} for name, path in inputs.items()}
        self._json('cases/generated/.cache.json',
                   {'format': 2, 'recipe_sha256': 'historical-fixture-signature',
                    'files': entries,
                    'provenance': {'kind': 'historical_snapshot_extension'}})
        api_names = ['special.old', 'special.renamed']
        provenance = {'schema': 2,
                      'cases_sha256': pruning.case_directory_digest(root / 'cases/generated'),
                      'measurements_sha256': pruning.directory_digest(root / 'results/measurements'),
                      'apis': {api: [{'provider': 'fixture', 'recipe': 'value',
                                      'precision': 'exact', 'status': 'historical',
                                      'cases': 2 if api == 'special.old' else 1}]
                               for api in api_names},
                      'reference_sha256': {api: 'a' * 64 for api in api_names}}
        self._json('results/reference-provenance.json', provenance)
        self._json('results/environment.json', {
            'cases_sha256': provenance['cases_sha256'],
            'measurements_sha256': provenance['measurements_sha256'],
            'reference_provenance_sha256': pruning._sha256(root / 'results/reference-provenance.json'),
            'coverage_sha256': pruning._sha256(root / 'results/coverage.csv'),
            'accuracy_summaries_sha256': pruning._sha256(root / 'results/accuracy-summaries.json'),
            'testing_sources_sha256': 'fixture-source-digest',
            'module_libraries': {'special': ['Sci']},
            'native_libraries': {'Boost.Math': '1_92', 'Boost.Graph': '1_92'},
            'sci_measurement_commits': {
                'special.old': 'revision-old',
                'special.renamed': 'revision-still-active',
                'special.nomatch': 'revision-unused',
                'special.other': 'revision-unrelated',
            },
        }, sort_keys=False)

    def _write(self, name: str, data: str) -> None:
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        _ = path.write_text(data, encoding='utf-8')

    def _json(self, name: str, value: object, *, sort_keys: bool = True) -> None:
        self._write(name, json.dumps(value, sort_keys=sort_keys) + '\n')


class PruneObsoleteTests(unittest.TestCase):
    def __init__(self, methodName: str = 'runTest') -> None:
        super().__init__(methodName)
        self.temp: tempfile.TemporaryDirectory[str] = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root: Path = Path(self.temp.name)
        _ = _Fixture(self.root)

    def review(self, plan: pruning.Plan) -> dict[str, object]:
        doc = pruning.review_template(plan)
        raw = doc['approvals']
        assert isinstance(raw, list)
        rows = cast(list[dict[str, object]], raw)
        target = rows[0]
        target['approved'] = True
        target['reason'] = 'No public export or usable signature; independently checked.'
        target['possible_exports_reviewed'] = True
        raw_inventory = doc['inventory_only_approvals']
        assert isinstance(raw_inventory, list)
        no_files = cast(list[dict[str, object]], raw_inventory)[0]
        no_files['approved'] = True
        no_files['reason'] = 'No current export and no active input data.'
        return doc

    def test_read_only_scan_and_hints(self) -> None:
        before = {name: pruning._sha256(self.root / name)
                  for name in ('cases/generated/.cache.json', 'api_inventory.csv',
                               'results/environment.json', 'results/coverage.csv')}
        plan = pruning.scan(self.root)
        self.assertEqual(len(plan.candidates), 1)
        self.assertEqual(plan.total_cases, 2)
        self.assertEqual(plan.unresolved_without_files, ('special.nomatch',))
        self.assertEqual(plan.candidates[0].api, 'special.old')
        self.assertEqual(plan.candidates[0].classification, 'unresolved_candidate')
        self.assertEqual(plan.candidates[0].authored_recipe_files,
                         ('cases/recipes/special/old.csv',))
        self.assertEqual(plan.as_json()['candidate_cases'], 2)
        self.assertEqual(before, {key: pruning._sha256(self.root / key) for key in before})

    def test_default_unapproved_and_reviewed_allowlist(self) -> None:
        plan = pruning.scan(self.root)
        template = pruning.review_template(plan)
        with self.assertRaisesRegex(ValueError, 'No candidates explicitly approved'):
            _ = pruning.approved_candidates(plan, template)
        review = self.review(plan)
        self.assertEqual([item.api for item in pruning.approved_candidates(plan, review)],
                         ['special.old'])
        self.assertEqual(pruning.approved_inventory_only(plan, review), ('special.nomatch',))
        raw = review['approvals']
        assert isinstance(raw, list)
        candidate = cast(list[dict[str, object]], raw)[0]
        candidate['input_sha256'] = '0' * 64
        with self.assertRaisesRegex(ValueError, 'Stale or altered'):
            _ = pruning.approved_candidates(plan, review)
        candidate['input_sha256'] = plan.candidates[0].input_sha256
        candidate['api'] = 'special.renamed'
        with self.assertRaisesRegex(ValueError, 'Unknown, mapped'):
            _ = pruning.approved_candidates(plan, review)

    def test_detects_cache_measurements_and_provenance_mismatch(self) -> None:
        path = self.root / 'results/measurements/special/old.csv'
        original = path.read_bytes()
        _ = path.write_bytes(original + b'garbage\n')
        with self.assertRaisesRegex(ValueError, 'measurement-directory digest'):
            _ = pruning.scan(self.root)
        _ = path.write_bytes(original)
        path = self.root / 'cases/generated/.cache.json'
        cache = pruning._json(path)
        files = cache['files']
        assert isinstance(files, dict)
        entry = cast(dict[str, dict[str, object]], files)['special/old.csv']
        entry['sha256'] = '0' * 64
        _ = path.write_text(json.dumps(cache))
        with self.assertRaisesRegex(ValueError, 'Stale case cache'):
            _ = pruning.scan(self.root)

    def test_refuses_input_measurement_identity_mismatch(self) -> None:
        path = self.root / 'results/measurements/special/old.csv'
        contents = path.read_text().replace('special.old/random-0', 'special.old/random-1')
        _ = path.write_text(contents)
        meta_path = self.root / 'results/environment.json'
        meta = pruning._json(meta_path)
        meta['measurements_sha256'] = pruning.directory_digest(self.root / 'results/measurements')
        provenance_path = self.root / 'results/reference-provenance.json'
        provenance = pruning._json(provenance_path)
        provenance['measurements_sha256'] = meta['measurements_sha256']
        _ = provenance_path.write_text(json.dumps(provenance) + '\n')
        meta['reference_provenance_sha256'] = pruning._sha256(provenance_path)
        _ = meta_path.write_text(json.dumps(meta))
        with self.assertRaisesRegex(ValueError, 'Mismatched measurement IDs'):
            _ = pruning.scan(self.root)

    def test_staging_never_modifies_original_and_cleans_inventory_only_when_approved(self) -> None:
        plan = pruning.scan(self.root)
        review = self.review(plan)
        selected = pruning.approved_candidates(plan, review)
        only_inventory = pruning.approved_inventory_only(plan, review)
        with self.assertRaisesRegex(ValueError, 'requires --drop-inventory'):
            _ = pruning.stage(plan, selected, self.root / 'badstage',
                              inventory_only=only_inventory)
        self.assertFalse((self.root / 'badstage').exists())
        self.assertTrue((self.root / 'cases/generated/special/old.csv').exists())

        def renderer(cases: object, records: object, meta: dict[str, object],
                     output: Path, **kwargs: object) -> None:
            del cases, records, kwargs
            output.mkdir(parents=True, exist_ok=True)
            _ = (output / 'README.md').write_text(
                f"- Native versions: {meta['native_libraries']!r}.\n")
            data = b'{"entries":{}}\n'
            _ = (output / 'accuracy-summaries.json').write_bytes(data)
            meta['accuracy_summaries_sha256'] = _sha(data)

        from sci_testing import coverage_audit, report

        with patch.object(pruning, '_staged_source_digest', return_value='fixture-source-digest'), \
             patch.object(coverage_audit, 'audit', return_value='api,canonical,cases\nspecial.renamed,special.renamed,1\n'), \
             patch.object(report, 'render_cached', side_effect=renderer):
            staged = pruning.stage(plan, selected, self.root / 'reviewed-stage',
                                   inventory_only=only_inventory, drop_inventory=True)
            self.assertTrue((self.root / 'cases/generated/special/old.csv').exists())
            self.assertFalse((staged / 'cases/generated/special/old.csv').exists())
            self.assertFalse((staged / 'results/measurements/special/old.csv').exists())
            self.assertEqual(pruning.case_paths(staged / 'cases/generated').keys(),
                             {'special/renamed.csv'})
            self.assertEqual(pruning._inventory(staged)[0]['api'], 'special.renamed')
            stage_meta = pruning._json(staged / 'results/environment.json')
            history = stage_meta['prune_history']
            assert isinstance(history, dict)
            self.assertEqual(history['removed_cases'], 2)
            self.assertEqual(history['removed_inventory_only_apis'],
                             ['special.nomatch'])
            self.assertEqual(stage_meta['sci_measurement_commits'], {
                'special.renamed': 'revision-still-active',
                'special.other': 'revision-unrelated',
            })
            # Staged metadata is read back by `report --check`. Its nested
            # dictionary order must reproduce the published README bytewise.
            stage_versions = stage_meta['native_libraries']
            self.assertEqual(list(cast(dict[str, str], stage_versions)),
                             ['Boost.Math', 'Boost.Graph'])
            self.assertEqual((staged / 'docs/accuracy/README.md').read_text(),
                             f'- Native versions: {stage_versions!r}.\n')
            self.assertEqual(pruning._check_stage(plan, staged)['approved'], ['special.old'])
            pruning.publish(plan, staged)
            self.assertFalse((self.root / 'cases/generated/special/old.csv').exists())
            self.assertFalse((self.root / 'results/measurements/special/old.csv').exists())
            self.assertTrue((self.root / 'cases/generated/special/renamed.csv').exists())
            self.assertEqual({row['api'] for row in pruning._inventory(self.root)},
                             {'special.other', 'special.renamed'})
            self.assertEqual(pruning.scan(self.root).total_cases, 0)

    def test_publish_rolls_back_every_previous_rename_on_failure(self) -> None:
        from sci_testing import coverage_audit, report

        plan = pruning.scan(self.root)
        review = self.review(plan)
        selected = pruning.approved_candidates(plan, review)
        originals = {
            name: pruning._sha256(self.root / name)
            for name in ('api_inventory.csv', 'cases/generated/.cache.json',
                         'results/environment.json', 'results/coverage.csv',
                         'results/reference-provenance.json',
                         'results/accuracy-summaries.json',
                         'docs/accuracy/README.md',
                         'results/measurements/special/old.csv')
        }

        def renderer(cases: object, records: object, meta: dict[str, object],
                     output: Path, **kwargs: object) -> None:
            del cases, records, kwargs
            output.mkdir(parents=True, exist_ok=True)
            _ = (output / 'README.md').write_text('new report\n')
            summary = b'{"entries":{}}\n'
            _ = (output / 'accuracy-summaries.json').write_bytes(summary)
            meta['accuracy_summaries_sha256'] = _sha(summary)

        with patch.object(pruning, '_staged_source_digest', return_value='fixture-source-digest'), \
             patch.object(coverage_audit, 'audit', return_value='api,canonical,cases\nspecial.renamed,special.renamed,1\n'), \
             patch.object(report, 'render_cached', side_effect=renderer):
            staged = pruning.stage(plan, selected, self.root / 'reviewed-stage',
                                   inventory_only=pruning.approved_inventory_only(plan, review),
                                   drop_inventory=True)

        replace = Path.replace

        def fail_after_several_moves(source: Path, target: Path) -> Path:
            if source == staged / 'results/accuracy-summaries.json':
                raise OSError('injected publish failure')
            return replace(source, target)

        with patch.object(pruning, '_staged_source_digest', return_value='fixture-source-digest'), \
             patch.object(Path, 'replace', new=fail_after_several_moves), \
             self.assertRaisesRegex(OSError, 'injected publish failure'):
            pruning.publish(plan, staged)

        self.assertEqual(originals, {name: pruning._sha256(self.root / name)
                                     for name in originals})
        self.assertTrue((self.root / 'cases/generated/special/old.csv').exists())
        self.assertTrue((staged / 'results/accuracy-summaries.json').exists())
        self.assertEqual(pruning.scan(self.root), plan)


if __name__ == '__main__':
    _ = unittest.main()
