"""Audit pinned recipes and generate typed public Rust registrations.

The 0.1 input profiles and their independent references remain pinned. For
Sci 0.3, a compiler-resolved inventory supplies the new callable paths. A
path match alone never establishes signature or numerical compatibility.
"""
from dataclasses import dataclass
import csv
from .rust_source import formatted
import re

from . import ROOT, SCI


def sci_target():
    """Use the pinned Sci 0.3 submodule in every checkout, including CI."""
    return SCI


@dataclass(frozen=True)
class Function:
    api: str
    arguments: tuple
    returns: str
    source: str


def split_arguments(text):
    start = depth = 0
    result = []
    for i, char in enumerate(text):
        if char in '(<[': depth += 1
        elif char in ')>]': depth -= 1
        elif char == ',' and depth == 0:
            result.append(text[start:i].strip()); start = i+1
    if text[start:].strip(): result.append(text[start:].strip())
    return result


def functions():
    for crate in sorted((SCI/'crates').glob('sci-*')):
        module = crate.name.removeprefix('sci-')
        for path in sorted((crate/'src').glob('*.rs')):
            text = path.read_text()
            for match in re.finditer(r'^pub fn (\w+)', text, re.M):
                begin = text.index('(', match.end())
                depth, end = 1, begin+1
                while depth:
                    if text[end] == '(': depth += 1
                    elif text[end] == ')': depth -= 1
                    end += 1
                tail = text[end:text.index('{', end)].strip()
                if not tail.startswith('->'): continue
                args = split_arguments(text[begin+1:end-1])
                if any(':' not in arg for arg in args): continue
                arguments = tuple(tuple(part.strip() for part in arg.split(':', 1)) for arg in args)
                namespace = f'{module}.{path.stem}' if path.stem in ['windows', 'hierarchy', 'vq', 'contingency', 'qmc'] else module
                yield Function(f'{namespace}.{match[1]}', arguments, tail[2:].strip(), path.relative_to(SCI).as_posix())


SCALAR_TYPES = {'f64', 'u32', 'i32', 'usize', 'u64', 'i64', 'bool', 'Complex', 'isize', 'ErrorPolicy', '&[f64]', '&[u64]'}


def supported_type(kind):
    if kind in SCALAR_TYPES: return True
    if kind.startswith('&[') and kind.endswith(']'): return supported_type(kind[2:-1])
    if kind.startswith('Vec<') and kind.endswith('>'): return supported_type(kind[4:-1])
    if kind.startswith('Option<') and kind.endswith('>'): return supported_type(kind[7:-1])
    return False


def special_functions():
    return [f for f in functions() if f.api.startswith('special.')
            and all(supported_type(t) for _,t in f.arguments)]


def registered_functions():
    return [f for f in functions() if f.api.split('.')[0] in ['special','complex']
            and all(supported_type(t) for _,t in f.arguments)]


def argument_expression(kind, index):
    value = f'&a["p{index}"]'
    if kind.startswith('&['):
        return f'&<Vec<{kind[2:-1]}> as Argument>::read({value})?'
    return f'<{kind} as Argument>::read({value})?'


def dispatch_paths(current_rows):
    """Map old fixture IDs to uniquely resolved 0.3 public paths.

    The returned omissions require review before generating comparisons;
    source/arity/behavior changes must still be checked in compiled Rust.
    """
    from .inventory import resolve_legacy_api

    entries, missing = [], []
    for function in registered_functions():
        mapped = resolve_legacy_api(function.api, current_rows, function.source)
        if mapped is None:
            missing.append(function.api)
        else:
            entries.append((function, mapped[0]))
    return entries, missing


REGISTRATION = re.compile(r'^\s*"([\w.]+)" => sci::([\w:]+);(.*);$', re.M)


def typed_manifest(source=None):
    """Keep the existing 390+ manually reviewed argument adapters intact.

    The historical source-signature scan sees only 221 of these registrations.
    The Rust dispatcher is therefore the current typed adapter manifest; no
    function or argument expression may disappear during namespace updates.
    """
    if source is None:
        source = (ROOT / 'src/generated.rs').read_text()
    entries = REGISTRATION.findall(source)
    if not entries or source.count(' => sci::') != len(entries):
        raise ValueError('typed Rust dispatch manifest has missing or malformed entries')
    names = [api for api, _, _ in entries]
    if len(set(names)) != len(names):
        raise ValueError('duplicate API IDs in typed Rust dispatch manifest')
    return entries


def typed_functions():
    """Typed signatures used by historical comparison sweeps.

    Preserve every manifest entry even when its original Rust source contains
    re-exports or generic methods the limited source regexp cannot see.
    Argument expressions are the already compiled converters in generated.rs.
    """
    old = {function.api: function for function in registered_functions()}
    result = []
    for api, path, arguments in typed_manifest():
        kinds = re.findall(r'(&?)<([\w<>]+) as Argument>::read\(&a\["p\d+"\]\)\?', arguments)
        if len(kinds) != arguments.count('Argument>::read'):
            raise ValueError(f'{api}: cannot recover every typed argument')
        inferred = tuple((f'p{i}', f'&[{kind[4:-1]}]' if borrow and kind.startswith('Vec<')
                          else kind) for i, (borrow, kind) in enumerate(kinds))
        previous = old.get(api)
        result.append(Function(api, inferred, previous.returns if previous else '',
                               previous.source if previous else path.replace('::', '/')))
    return result


def current_exports():
    """Read the reconciled compiler manifest without loading Rustdoc at runtime."""
    path = ROOT / 'api_inventory.csv'
    with path.open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    if not rows or 'kind' not in rows[0]:
        return None
    return [[row['api'], row['canonical'], row['source'], row['line']]
            for row in rows if row['kind'] == 'public_export']


def rust_registrations(current_rows=None, *, manifest=None):
    """Reconstruct typed Rust dispatch with verified public 0.3 callable paths.

    Pre-existing argument conversions are preserved literally. Every accepted
    callable must appear among compiler-resolved exports. Ambiguous renames
    fail for manual review instead of silently dropping tested historical IDs.
    """
    from .inventory import resolve_legacy_api

    if current_rows is None:
        current_rows = current_exports()
    if current_rows is None:
        raise ValueError('regenerate the compiler-resolved API inventory first')
    exports = {row[0] for row in current_rows}
    entries = typed_manifest(manifest)
    lines = ['// Generated by api_registry.py from the pinned public signatures.',
             'use super::Result;', 'use sci::complex::Complex;',
             'use serde_json::Value;', 'use crate::numeric::{Argument, Numbers};',
             'pub fn evaluate(api: &str, a: &Value) -> Result<Vec<f64>> {',
             '    macro_rules! registered {',
             '        ($($name:literal => $function:path; $($argument:expr),*;)*) => {',
             '            match api { $($name => $function($($argument),*).numbers(),)*',
             '                _ => Err(format!("unregistered API: {api}").into()), }',
             '        };', '    }', '    registered! {']
    for api, prior_path, args in entries:
        callable_api = prior_path.replace('::', '.')
        if callable_api not in exports:
            matched = resolve_legacy_api(api, current_rows)
            if matched is None:
                raise ValueError(f'{api}: prior callable {callable_api} is no longer exported; '
                                 'manual signature/path review is required')
            callable_api = matched[0]
        lines.append(f'        "{api}" => sci::{callable_api.replace(".", "::")};{args};')
    lines += ['    }', '}']
    return '\n'.join(lines)+'\n'


if __name__ == '__main__':
    root = ROOT
    source = formatted(rust_registrations())
    (root/'src/generated.rs').write_text(source)
    print(f'Registered {len(typed_manifest(source))} compiler-resolved typed adapters')
