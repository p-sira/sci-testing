"""Namespace migration inventories must never imply unmeasured coverage."""
import ast
import csv
import gzip
import json
import tempfile
import unittest
from pathlib import Path
from typing import cast

from sci_testing import ROOT, SCI
from sci_testing.api_registry import (
    Function,
    rust_registrations,
    sci_target,
    typed_manifest,
)
from sci_testing.inventory import (
    inventory,
    reconcile,
    refresh_inventory,
    resolve_legacy_api,
)

# These saved SciPy fixture IDs are deliberately NOT renamed to their Sci 0.3
# public paths: adapters preserve the original argument/observable contracts.
AUDITED_LEGACY_ALIASES = {
    'cluster.hierarchy.maxRstat': ('cluster.hierarchy.max_rstat', 2007),
    'interpolate.splprep': ('interpolate.fitting.make_splprep', 2007),
    'interpolate.splrep': ('interpolate.fitting.make_splrep', 2005),
    'linalg.solve_lyapunov': ('linalg.solvers.solve_continuous_lyapunov', 2002),
    'signal.check_COLA': ('signal.spectral.check_cola', 2009),
    'signal.check_NOLA': ('signal.spectral.check_nola', 2009),
    'signal.closest_STFT_dual_window': ('signal.utilities.closest_stft_dual_window', 2008),
    'signal.sosfreqz': ('signal.filtering.freqz_sos', 2002),
    'sparse.bmat': ('sparse.construct.block_array', 2002),
    'sparse.diags_array': ('sparse.construct.diags', 2002),
    'sparse.eye_array': ('sparse.construct.eye', 2005),
    'sparse.matrix_transpose': ('sparse.operations.transpose', 2002),
    'special.assoc_laguerre': ('special.orthogonal.eval_genlaguerre', 2001),
    'special.psi': ('special.gamma.digamma', 2001),
    'stats.kstest': ('stats.goodness.ks_2samp', 2012),
}


class InventoryTests(unittest.TestCase):
    def test_audited_legacy_aliases_target_current_public_exports_and_keep_case_ids(self):
        with (ROOT / 'api_inventory.csv').open(newline='') as stream:
            rows = {row['api']: row for row in csv.DictReader(stream)}
        total = 0
        for alias, (canonical, expected_count) in AUDITED_LEGACY_ALIASES.items():
            with self.subTest(alias=alias):
                old = rows[alias]
                target = rows[canonical]
                self.assertEqual(old['canonical'], canonical)
                self.assertEqual(old['status'], 'legacy_recipe_pending')
                self.assertEqual(old['kind'], 'legacy_alias')
                self.assertIn('Sci 0.3 alias disposition', old['evidence'])
                self.assertEqual(target['api'], target['canonical'])
                self.assertEqual(target['kind'], 'public_export')
                # Files/IDs belong to the historical fixtures, not the new paths.
                namespace, name = alias.rsplit('.', 1)
                plain = ROOT / 'cases/generated' / namespace / f'{name}.csv'
                compressed = plain.with_suffix('.csv.gz')
                self.assertTrue(plain.is_file() or compressed.is_file())
                stream = (gzip.open(compressed, 'rt', newline='') if compressed.is_file()
                          else plain.open(newline=''))
                with stream:
                    ids = [row['id'] for row in csv.DictReader(stream)]
                self.assertEqual(len(ids), expected_count)
                self.assertTrue(all(isinstance(case_id, str)
                                    and case_id.startswith(alias + '/') for case_id in ids))
                total += len(ids)
        self.assertEqual(total, 30074)

    def test_audited_legacy_aliases_agree_with_pinned_sci_disposition(self):
        source = ast.parse((SCI / 'tests/scipy/api_disposition.py').read_text())
        aliases = cast(dict[str, dict[str, str]],
                       next(ast.literal_eval(node.value) for node in source.body
                            if isinstance(node, ast.Assign)
                            and any(isinstance(target, ast.Name) and target.id == 'ALIASES'
                                    for target in node.targets)))
        for alias, (canonical, _) in AUDITED_LEGACY_ALIASES.items():
            with self.subTest(alias=alias):
                family, name = alias.split('.', 1)[0], alias.rsplit('.', 1)[-1]
                self.assertEqual(aliases[family][name], canonical.rsplit('.', 1)[-1])

    def test_kstest_is_two_sample_and_associated_laguerre_reorders_arguments(self):
        dispatch = (SCI / 'tests/statistics.rs').read_text()
        self.assertIn('"ks_2samp" | "kstest" => {', dispatch)
        self.assertIn('ks_2samp(&x, &y, alternative(p[0]))?', dispatch)
        self.assertNotIn('"ks_1samp" | "kstest"', dispatch)
        self.assertEqual(AUDITED_LEGACY_ALIASES['stats.kstest'][0],
                         'stats.goodness.ks_2samp')

        registrations = (ROOT / 'src/generated.rs').read_text().splitlines()
        entry = next(line for line in registrations
                     if '"special.assoc_laguerre" =>' in line)
        self.assertIn('sci::special::orthogonal::eval_genlaguerre;', entry)
        indices = [entry.index(f'a["p{i}"]') for i in (1, 2, 0)]
        self.assertEqual(indices, sorted(indices))

    def test_codegen_uses_pinned_submodule_with_no_sibling_checkout(self):
        self.assertEqual(sci_target().resolve(), SCI.resolve())
        self.assertTrue((sci_target()/'tests/fft_extended.rs').is_file())

    def test_rustdoc_resolves_public_reexports_through_private_modules(self):
        def item(name, visibility, inner, crate=0, span=None):
            return dict(name=name, visibility=visibility, crate_id=crate, inner=inner,
                        span=span)

        data = dict(root='0', index={
            '0': item('special', 'public', {'module': {'items': ['1', '2', '3', '4']}}),
            '1': item('implementation', 'private', {'module': {'items': ['5', '6']}}),
            '2': item('family', 'public', {'module': {'items': ['7', '8', '9']}}),
            '3': item('hidden', 'private', {'function': {}}),
            '4': item('external', 'public', {'use': {'id': '10', 'is_glob': False}}),
            '5': item('gamma', 'public', {'function': {}}, span={
                'filename': 'crates/sci-special/src/mod.rs', 'begin': [42, 0]}),
            '6': item('private_only', 'public', {'function': {}}),
            '7': item(None, 'public', {'use': {'id': '5', 'is_glob': False, 'name': 'gamma'}}),
            '8': item('all', 'public', {'use': {'id': '1', 'is_glob': True}}),
            '9': item('internal', 'private', {'function': {}}),
            '10': item('external', 'public', {'function': {}}, crate=1),
        })
        with tempfile.TemporaryDirectory() as directory:
            (Path(directory) / 'sci_special.json').write_text(json.dumps(data))
            rows = inventory(directory)
        self.assertEqual([r[0] for r in rows],
                         ['special.family.gamma', 'special.family.private_only'])
        self.assertEqual(rows[0][1:], ['special.family.gamma',
                                        'crates/sci-special/src/mod.rs', 42])

    def test_legacy_aliases_keep_pending_status_and_new_exports_are_untested(self):
        old = [dict(api='special.gamma', source='crates/sci-special/src/mod.rs'),
               dict(api='special.removed', source='crates/sci-special/src/old.rs')]
        new = [['special.gamma.gamma', 'special.gamma.gamma',
                'crates/sci-special/src/mod.rs', 42],
               ['special.wave.brand_new', 'special.wave.brand_new',
                'crates/sci-special/src/waves.rs', 73]]
        rows = {r[0]: r for r in reconcile(old, new)}
        self.assertEqual(rows['special.gamma'][1], 'special.gamma.gamma')
        self.assertEqual(rows['special.gamma'][4], 'legacy_recipe_pending')
        self.assertEqual(rows['special.gamma'][6], 'legacy_alias')
        self.assertEqual(rows['special.gamma.gamma'][4], 'legacy_recipe_pending')
        self.assertEqual(rows['special.gamma.gamma'][6], 'public_export')
        self.assertEqual(rows['special.wave.brand_new'][4], 'new_api_untested')
        self.assertEqual(rows['special.removed'][4], 'unresolved_legacy')
        self.assertEqual(rows['special.removed'][1], '')

    def test_refresh_moves_graph_exports_and_preserves_fixture_aliases(self):
        previous = [
            dict(api='sparse.csgraph.dijkstra', canonical='sparse.csgraph.dijkstra',
                 source='sci-sparse/src/csgraph.rs', line='10',
                 status='legacy_recipe_pending', evidence='historical graph fixture',
                 kind='public_export'),
            dict(api='special.erf', canonical='special.elementary.erf',
                 source='sci-special/src/mod.rs', line='2',
                 status='legacy_recipe_pending', evidence='historical scalar fixture',
                 kind='legacy_alias'),
            dict(api='special.elementary.erf', canonical='special.elementary.erf',
                 source='sci-special/src/mod.rs', line='20',
                 status='legacy_recipe_pending', evidence='historical scalar fixture',
                 kind='public_export'),
        ]
        current = [
            ['graph.dijkstra', 'graph.dijkstra', 'sci-graph/src/algorithms.rs', 12],
            ['special.elementary.erf', 'special.elementary.erf',
             'sci-special/src/mod.rs', 24],
            ['special.new_function', 'special.new_function',
             'sci-special/src/new.rs', 30],
        ]
        rows = {row[0]: row for row in refresh_inventory(previous, current)}
        self.assertEqual(rows['graph.dijkstra'][4], 'legacy_recipe_pending')
        self.assertEqual(rows['graph.dijkstra'][6], 'public_export')
        self.assertEqual(rows['sparse.csgraph.dijkstra'][1], 'graph.dijkstra')
        self.assertEqual(rows['sparse.csgraph.dijkstra'][6], 'legacy_alias')
        self.assertEqual(rows['special.erf'][1], 'special.elementary.erf')
        self.assertEqual(rows['special.elementary.erf'][3], 24)
        self.assertEqual(rows['special.new_function'][4], 'new_api_untested')
        with self.assertRaisesRegex(ValueError, 'disappeared'):
            refresh_inventory(previous, current[:1])

    def test_ambiguous_names_do_not_guess_public_callable(self):
        current = [['special.a.foo', 'special.a.foo', 'src/a.rs', 5],
                   ['special.b.foo', 'special.b.foo', 'src/b.rs', 6]]
        self.assertIsNone(resolve_legacy_api('special.foo', current))
        self.assertEqual(resolve_legacy_api('special.foo', current, 'src/b.rs'), current[1])
        self.assertIsNone(resolve_legacy_api('special.no_such_function', current))

    def test_duplicate_public_reexports_resolve_by_compiler_identity(self):
        source = 'crates/sci-linalg/src/qr.rs'
        current = [['linalg.decompositions.qr.qr_delete',
                    'linalg.decompositions.qr_delete', source, 335],
                   ['linalg.decompositions.qr_delete',
                    'linalg.decompositions.qr_delete', source, 335]]
        self.assertEqual(resolve_legacy_api('linalg.qr_delete', current), current[1])
        old = [{'api':'linalg.qr_delete','source':'crates/sci-linalg/src/old.rs'}]
        reconciled={row[0]:row for row in reconcile(old,current)}
        self.assertEqual(reconciled[current[0][0]][4],'legacy_recipe_pending')
        self.assertEqual(reconciled[current[1][0]][4],'legacy_recipe_pending')

    def test_dispatch_uses_new_path_but_retains_historical_case_id(self):
        f = Function('special.erf', (('x', 'f64'),), 'f64', 'src/old.rs')
        current = [['special.elementary.erf', 'special.elementary.erf', 'src/new.rs', 24]]
        manifest = '"special.erf" => sci::special::elementary::erf; <f64 as Argument>::read(&a["p0"])?;'
        self.assertEqual(typed_manifest(manifest)[0][0], f.api)
        source = rust_registrations(current,manifest=manifest)
        self.assertIn('"special.erf" => sci::special::elementary::erf;', source)
        self.assertNotIn('special::ErrorPolicy', source)

    def test_disallows_unresolved_manifest_without_losing_registrations(self):
        manifest = '"special.erf" => sci::special::no_longer_public; <f64 as Argument>::read(&a["p0"])?;'
        current = [['special.elementary.erf', 'special.elementary.erf', 'src/new.rs', 24]]
        migrated = rust_registrations(current,manifest=manifest)
        self.assertIn('"special.erf" => sci::special::elementary::erf;', migrated)
        with self.assertRaisesRegex(ValueError, 'no longer exported'):
            rust_registrations([['special.other.unrelated','special.other.unrelated','x',1]],manifest=manifest)

    def test_rejects_duplicate_dispatch_ids(self):
        manifest = '"special.erf" => sci::special::elementary::erf; 1;\n'*2
        with self.assertRaisesRegex(ValueError, 'duplicate API IDs'):
            typed_manifest(manifest)


if __name__ == '__main__':
    unittest.main()
