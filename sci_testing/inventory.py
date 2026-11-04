"""Discover compiler-resolved free functions from Rustdoc JSON.

An exported function is not necessarily supported by a comparison recipe.
``reconcile`` keeps historical case IDs visible while recording every current
export, including functions still awaiting independent reference evidence.
"""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path
from typing import TypedDict, cast

from . import ROOT
from .case_csv import table


class _RustdocItem(TypedDict):
    crate_id: int
    visibility: str
    name: str | None
    inner: dict[str, object]
    span: dict[str, str | list[int]] | None


class _MigrationSummary(TypedDict):
    status: list[tuple[str, int]]
    new_apis: list[str]
    namespace_moves: list[tuple[str, str]]
    unresolved_legacy: list[str]


HEADERS = ['api', 'canonical', 'source', 'line', 'status', 'evidence', 'kind']

# Several real public names resolve to one Rustdoc function item. Prefer the
# explicit mathematical names without rewriting historical case identifiers.
# These choices may precede the corresponding breaking name cleanup in Sci.
PREFERRED_CANONICAL_EXPORTS = frozenset({
    'linalg.basic.issymmetric',
    'special.elliptic.ellipkinc',
})

# Rustdoc 1.99 also traversed a private `qr` module behind the public `qr`
# *function* reexport. The six nested paths fail to compile (E0603); they are
# not historical fixtures and must not become synthetic legacy aliases.
REVIEWED_PRIVATE_MODULE_GHOSTS = {
    'linalg.decompositions.qr.khatri_rao': 'linalg.basic.khatri_rao',
    'linalg.decompositions.qr.qr_delete': 'linalg.decompositions.qr_delete',
    'linalg.decompositions.qr.qr_insert': 'linalg.decompositions.qr_insert',
    'linalg.decompositions.qr.qr_multiply': 'linalg.decompositions.qr_multiply',
    'linalg.decompositions.qr.qr_update': 'linalg.decompositions.qr_update',
    'linalg.decompositions.qr.rq': 'linalg.decompositions.rq',
}


def inventory(doc_root=None):
    """Return (export, canonical, source, line) from crate-local Rustdoc JSON.

    Publicly reexported functions (including those defined in private modules)
    are traversed, but a named reexport of a private module cannot be used to
    expose its members. External-crate items are excluded: each Sci crate is
    inspected separately, and the umbrella crate reexports those modules.
    """
    # /sci and /sci-testing share Cargo's parent target directory by default.
    doc_root = Path(doc_root) if doc_root is not None else ROOT.parent / 'target/doc'
    rows = []
    for path in sorted(doc_root.glob('sci_*.json')):
        data = cast(dict[str, object], json.loads(path.read_text()))
        module = path.stem.removeprefix('sci_')
        index = cast(dict[str, _RustdocItem], data['index'])
        found: defaultdict[str, set[str]] = defaultdict(set)

        # Bind per-crate maps, not the loop variables from a later Rustdoc file.
        def walk(item_id: str, prefix: str, visited: set[str],
                 index: dict[str, _RustdocItem] = index,
                 found: defaultdict[str, set[str]] = found) -> None:
            item_id = str(item_id)
            if item_id in visited:
                return
            item = index.get(item_id)
            if item is None or item['crate_id'] != 0:
                return
            inner = item['inner']
            if 'function' in inner:
                found[item_id].add(prefix)
            elif 'module' in inner:
                for child_id in cast(dict[str, list[str]], inner['module'])['items']:
                    child = index.get(str(child_id))
                    if child is None or child['visibility'] != 'public':
                        continue
                    body = child['inner']
                    if 'use' in body:
                        use = cast(dict[str, object], body['use'])
                        target = use.get('id')
                        if target is not None:
                            target_item = index.get(str(target))
                            # Rustdoc may treat a name collision between a
                            # private module and a public function reexport as
                            # a navigable module. Rust correctly rejects the
                            # resulting `public_path::private_module::f`.
                            if (not use['is_glob'] and target_item is not None
                                    and 'module' in target_item['inner']
                                    and target_item['visibility'] != 'public'):
                                continue
                            # Rustdoc 1.99 stores re-export names on the use
                            # body; the containing item has ``name: null``.
                            symbol = use.get('name') or child.get('name')
                            name = prefix if use['is_glob'] else f'{prefix}.{symbol}'
                            walk(str(target), name, visited | {item_id})
                    elif 'module' in body or 'function' in body:
                        walk(str(child_id), f'{prefix}.{child["name"]}', visited | {item_id})

        walk(str(data['root']), module, set())
        for item_id, names in found.items():
            names = sorted(names, key=lambda name: (name.count('.'), len(name), name))
            canonical = next((name for name in names
                              if name in PREFERRED_CANONICAL_EXPORTS), names[0])
            span = index[item_id].get('span') or {}
            for name in names:
                rows.append([name, canonical, span.get('filename', ''),
                             (span.get('begin') or [''])[0]])
    return sorted(rows)


def resolve_legacy_api(api, current_rows, old_source=''):
    """Find the current public callable for a historical recipe API.

    Require a unique path. For duplicate names, the Rust source span must
    distinguish the target; otherwise the recipe remains unresolved.
    Resolving a name does not establish that its signature or semantics agree.
    """
    # When reading a reconciled CSV, ignore synthetic historical aliases:
    # they are recipe identifiers and cannot be called as current Rust paths.
    current_rows = [row for row in current_rows if len(row) < 7 or row[6] == 'public_export']
    root, _, name = api.partition('.')
    candidates = [row for row in current_rows if row[0] == api]
    if len(candidates) == 1:
        return candidates[0]
    suffix = '.' + name
    candidates = [row for row in current_rows if row[0].startswith(root + '.')
                  and (row[0].endswith(suffix) or row[0].rsplit('.', 1)[-1] == name.rsplit('.', 1)[-1])]
    if len(candidates) == 1:
        return candidates[0]
    # Multiple *public* reexports may resolve to one Rustdoc function item.
    if candidates and len({row[1] for row in candidates}) == 1:
        return next((row for row in candidates if row[0] == row[1]), candidates[0])
    same_origin = [row for row in candidates if old_source and row[2] == old_source]
    if len(same_origin) == 1:
        return same_origin[0]
    return None


def reconcile(legacy_rows, current_rows):
    """Preserve old identifiers and classify new exports without asserting parity.

    ``legacy_recipe_pending`` means that a public callable was identified by
    source/path evidence, but the old fixture has not been run against 0.3.
    ``untested`` means that an independent recipe is not registered here.
    ``unresolved_legacy`` means that the old fixture cannot be dispatched
    without further review. This CSV is an audit, not a success report.
    """
    by_api = {row[0]: list(row) for row in current_rows}
    if len(by_api) != len(current_rows):
        raise ValueError('duplicate compiler-resolved API names')
    legacy: list[dict[str, str]] = [{'api': row['api'], 'source': row.get('source', '')} if isinstance(row, dict)
              else {'api': row[0], 'source': row[2]} for row in legacy_rows]
    old_names = {(row['api'].split('.')[0], row['api'].rsplit('.', 1)[-1]) for row in legacy}
    entries = {}
    for api, row in by_api.items():
        new_name = (api.split('.')[0], api.rsplit('.', 1)[-1]) not in old_names
        entries[api] = [api, row[1], row[2], row[3],
                        'new_api_untested' if new_name else 'untested',
                        ('new public function; independent 0.3 reference not validated' if new_name
                         else 'compiler export; independent 0.3 reference not validated'),
                        'public_export']
    for old in legacy:
        api, origin = old['api'], old['source']
        current = resolve_legacy_api(api, current_rows, origin)
        if current is None:
            entries[api] = [api, '', origin, '', 'unresolved_legacy',
                            'no unique current public callable; historical fixture retained',
                            'legacy_alias']
            continue
        target = current[0]
        entries[api] = [api, current[1], current[2], current[3],
                        'legacy_recipe_pending',
                        'name/path match; signature and independent 0.3 comparison pending',
                        'public_export' if api in by_api else 'legacy_alias']
        # A canonical export has a historical recipe, even when that recipe's
        # stored case identifier predates its module namespace.
        entries[target][4:6] = ['legacy_recipe_pending',
                                 'historical recipe mapped; 0.3 comparison pending']
    mapped_canonicals = {row[1] for row in entries.values()
                         if row[4] == 'legacy_recipe_pending'}
    for row in entries.values():
        if (row[6] == 'public_export' and row[4] in {'untested', 'new_api_untested'}
                and row[1] in mapped_canonicals):
            row[4:6] = ['legacy_recipe_pending',
                        'same compiler item as mapped historical recipe; 0.3 comparison pending']
    return [entries[key] for key in sorted(entries)]


def refresh_inventory(previous_rows, current_rows):
    """Refresh compiler exports while retaining audited historical recipe IDs.

    The Sci 0.3 graph extraction moves the former ``sparse.csgraph`` functions
    into ``sci::graph``. Their saved input IDs remain stable, but they are no
    longer public exports at their historical paths.
    """
    previous = {row['api']: row for row in previous_rows}
    if len(previous) != len(previous_rows):
        raise ValueError('duplicate previous inventory names')
    current = {row[0]: row for row in current_rows}
    if len(current) != len(current_rows):
        raise ValueError('duplicate compiler-resolved export names')

    previous_exports = {name for name, row in previous.items()
                        if row['kind'] == 'public_export'}
    moved = {name: f'graph.{name.rsplit(".", 1)[-1]}'
             for name in previous_exports - current.keys()
             if name.startswith('sparse.csgraph.')}
    if any(destination not in current for destination in moved.values()):
        raise ValueError('graph extraction has unresolved public exports')
    # Do not turn confirmed non-callable Rustdoc traversal artifacts into
    # historical aliases. Require their actual public compiler item to survive.
    ghosts = {name for name, canonical in REVIEWED_PRIVATE_MODULE_GHOSTS.items()
              if name in previous_exports - current.keys()
              and previous[name]['canonical'] == canonical
              and previous[name]['source'] == 'crates/sci-linalg/src/qr.rs'
              and canonical in current}
    missing: set[str] = previous_exports - current.keys() - moved.keys() - ghosts
    if missing:
        raise ValueError(f'public exports disappeared without reviewed mapping: {sorted(missing)}')
    reverse_moved = {new: old for old, new in moved.items()}
    if len(reverse_moved) != len(moved):
        raise ValueError('ambiguous moved graph exports')

    updated = {}
    for name, (api, canonical, source, line) in current.items():
        prior = previous.get(name) or previous.get(reverse_moved.get(name, ''))
        if prior is None:
            status = 'new_api_untested'
            evidence = 'new compiler export; independent reference not validated'
        else:
            status, evidence = prior['status'], prior['evidence']
        updated[name] = [api, canonical, source, line, status, evidence, 'public_export']

    for name, prior in previous.items():
        if name in current or name in ghosts:
            continue
        canonical = moved.get(prior['canonical'], prior['canonical'])
        if name in moved:
            canonical = current[moved[name]][1]
        if canonical and canonical not in current:
            raise ValueError(f'{name}: historical alias target is no longer exported: {canonical}')
        if canonical:
            # The preferred canonical for a shared compiler item can change
            # while its historical spelling is still exported. Rebind aliases
            # to the new canonical without touching their stored case IDs.
            canonical = cast(str, current[canonical][1])
        updated[name] = [name, canonical, prior['source'], prior['line'],
                         prior['status'], prior['evidence'], 'legacy_alias']
    return [updated[name] for name in sorted(updated)]


def migration_summary(rows: list[list[str]]) -> _MigrationSummary:
    """Separate namespace moves, truly new names and unverified comparisons."""
    from collections import Counter

    status = Counter(row[4] for row in rows)
    return {'status': sorted(status.items()),
            'new_apis': sorted(row[0] for row in rows if row[4] == 'new_api_untested'),
            'namespace_moves': sorted((row[0], row[1]) for row in rows
                                      if row[4] == 'legacy_recipe_pending' and row[0] != row[1]
                                      and row[0].count('.') < row[1].count('.')),
            'unresolved_legacy': sorted(row[0] for row in rows
                                        if row[4] == 'unresolved_legacy')}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--doc-root', type=Path, default=ROOT.parent / 'target/doc')
    parser.add_argument('--output', type=Path, default=ROOT / 'api_inventory.csv')
    parser.add_argument('--reconcile', action='store_true',
                        help='map old case IDs to the current exports and record untested APIs')
    parser.add_argument('--refresh', action='store_true',
                        help='refresh an existing reconciled inventory against compiler exports')
    parser.add_argument('--legacy', type=Path,
                        help='original pre-migration API registry, when reconciling')
    arguments = parser.parse_args()
    rows = inventory(arguments.doc_root)
    if not rows:
        parser.error(f'no sci_*.json Rustdoc inventories in {arguments.doc_root}')
    if arguments.reconcile and arguments.refresh:
        parser.error('--reconcile and --refresh are mutually exclusive')
    if arguments.refresh:
        with arguments.output.open(newline='') as stream:
            previous = list(csv.DictReader(stream))
        if not previous or 'kind' not in previous[0]:
            parser.error('--refresh requires an existing reconciled inventory')
        rows = refresh_inventory(previous, rows)
        headers = HEADERS
    elif arguments.reconcile:
        if arguments.legacy is None:
            parser.error('--reconcile requires --legacy pointing to the original pre-migration CSV')
        with arguments.legacy.open(newline='') as stream:
            old = list(csv.DictReader(stream))
        if any('status' in row or 'kind' in row for row in old):
            parser.error('--legacy must be the original pre-migration CSV; a prior '
                         'reconciliation is not a valid independent baseline')
        rows = reconcile(old, rows)
        headers = HEADERS
    else:
        headers = HEADERS[:4]
    arguments.output.write_text(table(headers, rows))
    print(f'{len(rows)} rows; {len({row[1] for row in rows if row[1]})} exported functions')
    if arguments.reconcile:
        summary = migration_summary(cast(list[list[str]], rows))
        print(f'{len(summary["new_apis"])} genuinely new names; '
              f'{len(summary["namespace_moves"])} matched namespace moves; '
              f'{len(summary["unresolved_legacy"])} unresolved legacy IDs')
