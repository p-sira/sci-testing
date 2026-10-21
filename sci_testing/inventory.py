"""Export compiler-resolved public free functions, including macro expansions.

Regenerate Rustdoc JSON with:
RUSTC_BOOTSTRAP=1 RUSTDOCFLAGS='-Z unstable-options --output-format json' \
    cargo doc --workspace --no-deps
The ordinary harness reads the checked-in CSV and needs only stable Rust.
"""
import csv
from collections import defaultdict
from .case_csv import table
from . import ROOT, SCI
import json


def inventory():
    rows = []
    for path in sorted((SCI/'target/doc').glob('sci_*.json')):
        data = json.loads(path.read_text())
        module = path.stem.removeprefix('sci_')
        index = data['index']
        found = defaultdict(set)
        def walk(item_id, prefix, visited):
            if item_id in visited: return
            item = index.get(str(item_id))
            if item is None or item['crate_id'] != 0: return
            kind, body = next(iter(item['inner'].items()))
            if kind == 'function':
                found[item_id].add(prefix)
            elif kind == 'module':
                for child_id in body['items']:
                    child = index[str(child_id)]
                    if child['visibility'] != 'public': continue
                    child_kind, child_body = next(iter(child['inner'].items()))
                    if child_kind == 'use':
                        target = child_body.get('id')
                        if target is not None:
                            walk(target, prefix if child_body['is_glob'] else f'{prefix}.{child_body["name"]}', visited|{item_id})
                    elif child_kind in ['module','function']:
                        walk(child_id,f'{prefix}.{child["name"]}',visited|{item_id})
        walk(data['root'],module,set())
        for item_id, names in found.items():
            names = sorted(names,key=lambda name:(name.count('.'),len(name),name))
            item = index[str(item_id)]
            span = item['span'] or {}
            for name in names:
                rows.append([name, names[0], span.get('filename',''), (span.get('begin') or [''])[0]])
    return sorted(rows)


if __name__ == '__main__':
    rows=inventory()
    (ROOT/'api_inventory.csv').write_text(table(['api','canonical','source','line'],rows))
    print(f'{len(rows)} exported names; {len({row[1] for row in rows})} distinct functions')
