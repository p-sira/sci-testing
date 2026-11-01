"""Namespace migration inventories must never imply unmeasured coverage."""
import json
from pathlib import Path
import tempfile
import unittest

from sci_testing.api_registry import (Function, dispatch_paths,
                                      rust_registrations, typed_manifest)
from sci_testing.inventory import inventory, reconcile, resolve_legacy_api


class InventoryTests(unittest.TestCase):
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
