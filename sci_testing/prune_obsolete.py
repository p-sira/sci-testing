"""Audit and stage explicitly reviewed obsolete Sci recipe snapshots.

Default execution is read-only.  ``--write-review-template FILE`` produces a
deliberately *unapproved* allowlist for human review.  ``--stage DIRECTORY
--approved FILE`` then creates a complete, isolated proposed snapshot using
hardlinks for retained raw data.  ``--publish DIRECTORY`` is a separately
guarded, reversible multi-path replacement; it is never implicit in staging.

An ``unresolved_legacy`` inventory status does not establish obsolescence:
several such names are plausible Sci 0.3 renames.  Every deletion therefore
requires an exact reviewed API/count/input-hash/measurement-hash approval.
Mapped historical aliases and current public exports cannot be approved.

The source generators still contain historical recipes.  The staged case cache
retains its authentic recipe signature, with generator equivalence explicitly
marked false.  ``report --full`` can regenerate those historical fixtures; do
not use it to refresh a pruned snapshot until the generators are reconciled.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
from collections.abc import Callable, Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import cast

from . import ROOT, case_csv
from .reference_integrity import ProvenanceEntry, ReferenceCase, ReferenceMeasurement

case_directory_digest = cast(Callable[[Path], str], case_csv.case_directory_digest)
case_paths = cast(Callable[[Path], dict[str, Path]], case_csv.case_paths)
csv_paths = cast(Callable[[Path], list[Path]], case_csv.csv_paths)
csv_text = cast(Callable[[Path], str], case_csv.csv_text)
directory_digest = cast(Callable[[Path], str], case_csv.directory_digest)
read_case_file = cast(Callable[[Path, str], list[ReferenceCase]], case_csv.read_case_file)
read_measurement_content = cast(Callable[[str], dict[str, ReferenceMeasurement]],
                                case_csv.read_measurement_content)


class _Arguments(argparse.Namespace):
    def __init__(self) -> None:
        super().__init__()
        self.root: Path = ROOT
        self.write_review_template: Path | None = None
        self.approved: Path | None = None
        self.stage: Path | None = None
        self.publish: Path | None = None
        self.drop_inventory: bool = False
        self.jobs: int = 1


@dataclass(frozen=True)
class Candidate:
    api: str
    cases: int
    input_path: str
    measurement_path: str
    input_sha256: str
    measurement_sha256: str
    input_bytes: int
    measurement_bytes: int
    possible_exports: tuple[str, ...]
    authored_recipe_files: tuple[str, ...]
    classification: str = 'unresolved_candidate'


@dataclass(frozen=True)
class Plan:
    source: str
    environment_sha256: str
    inventory_sha256: str
    cache_sha256: str
    cases_sha256: str
    measurements_sha256: str
    reference_provenance_sha256: str
    coverage_sha256: str
    candidates: tuple[Candidate, ...]
    unresolved_without_files: tuple[str, ...]

    @property
    def total_cases(self) -> int:
        return sum(item.cases for item in self.candidates)

    def as_json(self) -> dict[str, object]:
        return {**asdict(self), 'candidate_apis': len(self.candidates),
                'candidate_cases': self.total_cases,
                'policy': 'review-required; dry-run never deletes anything',
                'generator_sources_to_review': [
                    'sci_testing/registered_cases.py', 'sci_testing/explicit_cases.py',
                    'sci_testing/generate_cases.py', 'sci_testing/case_set.py',
                    'sci_testing/api_registry.py'],
                'generator_warning': 'report --full can regenerate old recipes until their sources are migrated'}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _json(path: Path) -> dict[str, object]:
    decoded = cast(object, json.loads(path.read_bytes()))
    if not isinstance(decoded, dict):
        raise TypeError(f'Expected a JSON object: {path}')
    return cast(dict[str, object], decoded)


def _mapping(value: object, name: str) -> dict[str, object]:
    if not isinstance(value, dict):
        raise TypeError(f'Invalid {name} mapping')
    entries = cast(dict[object, object], value)
    if any(not isinstance(k, str) for k in entries):
        raise ValueError(f'Invalid {name} mapping keys')
    return cast(dict[str, object], entries)


def _inventory(root: Path) -> list[dict[str, str]]:
    with (root / 'api_inventory.csv').open(newline='', encoding='utf-8') as stream:
        rows = list(csv.DictReader(stream))
    expected = {'api', 'canonical', 'status', 'kind'}
    if not rows or not expected.issubset(rows[0]):
        raise ValueError('Unrecognized reconciled API inventory')
    if len({row['api'] for row in rows}) != len(rows):
        raise ValueError('Duplicate inventory API names')
    return rows


def _key(api: str) -> str:
    module, name = api.rsplit('.', 1)
    return f'{module}/{name}.csv'


def _authored_recipe_files(root: Path, api: str) -> tuple[str, ...]:
    logical = _key(api)
    return tuple(str(file.relative_to(root)) for kind in ('base', 'recipes')
                 for suffix in ('', '.gz')
                 if (file := root / 'cases' / kind / (logical + suffix)).is_file())


def _possible_exports(api: str, rows: Sequence[dict[str, str]]) -> tuple[str, ...]:
    module = api.split('.', 1)[0]
    name = re.sub('[^a-z0-9]', '', api.rsplit('.', 1)[-1].casefold())
    return tuple(sorted(row['api'] for row in rows if row['kind'] == 'public_export'
                        and row['api'].split('.', 1)[0] == module
                        and re.sub('[^a-z0-9]', '', row['api'].rsplit('.', 1)[-1].casefold()) == name))


def _measurement_paths(root: Path) -> dict[str, Path]:
    result: dict[str, Path] = {}
    directory = root / 'results/measurements'
    for path in csv_paths(directory):
        logical = path.relative_to(directory).as_posix().removesuffix('.gz')
        if logical in result:
            raise ValueError(f'Duplicate compressed/plain measurement file: {logical}')
        result[logical] = path
    return result


def scan(root: Path = ROOT) -> Plan:
    """Read-only authenticated inventory of potentially obsolete active pairs."""
    root = root.resolve(strict=True)
    meta_path = root / 'results/environment.json'
    inventory_path = root / 'api_inventory.csv'
    cache_path = root / 'cases/generated/.cache.json'
    meta = _json(meta_path)
    cache = _json(cache_path)
    inventory = _inventory(root)
    if cache.get('format') not in (1, 2):
        raise ValueError('Unrecognized case cache format')
    entries = _mapping(cache.get('files'), 'case cache files')
    inputs = case_paths(root / 'cases/generated')
    measured = _measurement_paths(root)
    if set(entries) != set(inputs) or set(inputs) != set(measured):
        raise ValueError('Input manifest, generated input and measurement files differ')
    if meta.get('cases_sha256') != case_directory_digest(root / 'cases/generated'):
        raise ValueError('Saved case-directory digest mismatch')
    if meta.get('measurements_sha256') != directory_digest(root / 'results/measurements'):
        raise ValueError('Saved measurement-directory digest mismatch')
    for logical, path in inputs.items():
        entry = _mapping(entries[logical], f'case cache entry for {logical}')
        if entry.get('sha256') != _sha256(path) or type(entry.get('count')) is not int or cast(int, entry['count']) < 1:
            raise ValueError(f'Stale case cache: {logical}')
    coverage = root / 'results/coverage.csv'
    if _sha256(coverage) != meta.get('coverage_sha256'):
        raise ValueError('Saved coverage digest mismatch')
    provenance = root / 'results/reference-provenance.json'
    if not isinstance(meta.get('reference_provenance_sha256'), str):
        raise TypeError('Reference provenance must be present before pruning')
    if _sha256(provenance) != meta['reference_provenance_sha256']:
        raise ValueError('Saved reference provenance digest mismatch')
    prior = _json(provenance)
    if prior.get('schema') != 2 or prior.get('cases_sha256') != meta['cases_sha256'] or prior.get('measurements_sha256') != meta['measurements_sha256']:
        raise ValueError('Reference provenance does not match the saved snapshot')
    provenance_apis = _mapping(prior.get('apis'), 'reference provenance entries')
    reference_hashes = _mapping(prior.get('reference_sha256'), 'reference digests')
    apis = {_key_to_api(name) for name in entries}
    if set(provenance_apis) != apis or set(reference_hashes) != apis:
        raise ValueError('Reference provenance API coverage does not match saved inputs')
    candidates: list[Candidate] = []
    without: list[str] = []
    for row in inventory:
        if not (row['status'] == 'unresolved_legacy' and row['kind'] == 'legacy_alias'
                and row['canonical'] == ''):
            continue
        api = row['api']
        logical = _key(api)
        if logical not in inputs:
            without.append(api)
            continue
        path = inputs[logical]
        measurement = measured[logical]
        entry = _mapping(entries[logical], f'case cache entry {logical}')
        batch = read_case_file(path, logical)
        if len(batch) != entry['count'] or len({c['id'] for c in batch}) != len(batch):
            raise ValueError(f'Missing or duplicate input cases: {api}')
        if any(c['api'] != api or not c['id'].startswith(api + '/') for c in batch):
            raise ValueError(f'Invalid case identities: {api}')
        records = read_measurement_content(csv_text(measurement))
        if set(records) != {c['id'] for c in batch}:
            raise ValueError(f'Mismatched measurement IDs: {api}')
        if any(not record['reference'] or 'Sci' not in record['libraries'] for record in records.values()):
            raise ValueError(f'Incomplete measurements: {api}')
        candidates.append(Candidate(api, len(batch), path.relative_to(root).as_posix(),
                                    measurement.relative_to(root).as_posix(), _sha256(path),
                                    _sha256(measurement), path.stat().st_size,
                                    measurement.stat().st_size,
                                    _possible_exports(api, inventory),
                                    _authored_recipe_files(root, api)))
    return Plan(str(root), _sha256(meta_path), _sha256(inventory_path),
                _sha256(cache_path), cast(str, meta['cases_sha256']),
                cast(str, meta['measurements_sha256']),
                cast(str, meta['reference_provenance_sha256']),
                cast(str, meta['coverage_sha256']),
                tuple(sorted(candidates, key=lambda item: item.api)), tuple(sorted(without)))


def _key_to_api(logical: str) -> str:
    return logical.removesuffix('.csv').replace('/', '.')


def review_template(plan: Plan) -> dict[str, object]:
    """A reviewer must deliberately edit each candidate's decision and reason."""
    return {
        'schema': 1,
        'snapshot': {key: getattr(plan, key) for key in
                     ('environment_sha256', 'inventory_sha256', 'cache_sha256')},
        'approvals': [{'api': item.api, 'cases': item.cases,
                       'input_sha256': item.input_sha256,
                       'measurement_sha256': item.measurement_sha256,
                       'approved': False, 'reason': '',
                       'possible_exports_reviewed': False}
                      for item in plan.candidates],
        'inventory_only_approvals': [{'api': api, 'approved': False, 'reason': ''}
                                     for api in plan.unresolved_without_files],
    }


def approved_candidates(plan: Plan, review: Mapping[str, object]) -> tuple[Candidate, ...]:
    if review.get('schema') != 1:
        raise ValueError('Unsupported prune review schema')
    expected = review_template(plan)['snapshot']
    if review.get('snapshot') != expected:
        raise ValueError('Prune approval was prepared for a different snapshot')
    raw_rows = review.get('approvals')
    if not isinstance(raw_rows, list):
        raise TypeError('Prune review is missing the approvals list')
    rows = cast(list[object], raw_rows)
    known = {item.api: item for item in plan.candidates}
    used: set[str] = set()
    selected: list[Candidate] = []
    for raw in rows:
        row = _mapping(raw, 'prune approval')
        api = row.get('api')
        if not isinstance(api, str) or api not in known or api in used:
            raise ValueError(f'Unknown, mapped or duplicate prune approval API: {api}')
        used.add(api)
        item = known[api]
        if (type(row.get('cases')) is not int or row['cases'] != item.cases
                or row.get('input_sha256') != item.input_sha256
                or row.get('measurement_sha256') != item.measurement_sha256):
            raise ValueError(f'Stale or altered prune approval hashes/count: {api}')
        if row.get('approved') is True:
            if (row.get('possible_exports_reviewed') is not True
                    or not isinstance(row.get('reason'), str)
                    or len(cast(str, row['reason']).strip()) < 12):
                raise ValueError(f'Prune approval needs explicit rename review and explanation: {api}')
            selected.append(item)
        elif row.get('approved') is not False:
            raise ValueError(f'Invalid prune approval flag: {api}')
    if not selected:
        raise ValueError('No candidates explicitly approved for pruning')
    return tuple(sorted(selected, key=lambda item: item.api))


def approved_inventory_only(plan: Plan, review: Mapping[str, object]) -> tuple[str, ...]:
    """Optional explicit inventory cleanup for legacy aliases with zero files."""
    if review.get('snapshot') != review_template(plan)['snapshot']:
        raise ValueError('Prune approval was prepared for a different snapshot')
    raw_rows = review.get('inventory_only_approvals')
    if not isinstance(raw_rows, list):
        raise TypeError('Missing inventory-only review approvals')
    candidates = set(plan.unresolved_without_files)
    selected: list[str] = []
    seen: set[str] = set()
    for raw in cast(list[object], raw_rows):
        row = _mapping(raw, 'inventory-only approval')
        api = row.get('api')
        if not isinstance(api, str) or api not in candidates or api in seen:
            raise ValueError(f'Unknown, mapped or duplicate inventory-only API: {api}')
        seen.add(api)
        if row.get('approved') is True:
            if not isinstance(row.get('reason'), str) or len(cast(str, row['reason']).strip()) < 12:
                raise ValueError(f'Inventory-only removal needs a reviewed explanation: {api}')
            selected.append(api)
        elif row.get('approved') is not False:
            raise ValueError(f'Invalid inventory-only review decision: {api}')
    return tuple(sorted(selected))


def _write_json(path: Path, value: object, *, sort_keys: bool = True) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    _ = path.write_text(json.dumps(value, indent=2, sort_keys=sort_keys) + '\n',
                        encoding='utf-8')


def _hardlink_survivors(src: Path, dest: Path, excluded: set[str]) -> None:
    for path in csv_paths(src):
        logical = path.relative_to(src).as_posix().removesuffix('.gz')
        if logical in excluded:
            continue
        target = dest / path.relative_to(src)
        target.parent.mkdir(parents=True, exist_ok=True)
        os.link(path, target, follow_symlinks=False)


def _staged_source_digest(root: Path, staged_inventory: Path) -> str:
    """Mirror report.digest_sources with exactly one staged inventory override.

    This keeps staged ``testing_sources_sha256`` identical to the report's
    regular source digest *after* the inventory replacement is published.
    """
    names = {'.gitmodules', 'Cargo.toml', 'Cargo.lock', 'pyproject.toml',
             'uv.lock', 'sci-revision.txt', 'api_inventory.csv', 'README.md'}
    patterns = [
        '*.py', 'sci_testing/**/*.py', 'tests/**/*.py', '*.rs', 'src/*.rs',
        'native/*.cpp', 'rust_oracles/Cargo.toml', 'rust_oracles/Cargo.lock',
        'rust_oracles/src/*.rs', '.github/workflows/*.yml',
        'cases/base/**/*.csv', 'cases/base/**/*.csv.gz',
        'cases/recipes/**/*.csv', 'cases/recipes/**/*.csv.gz',
    ]
    names.update(path.relative_to(root).as_posix()
                 for pattern in patterns for path in root.glob(pattern))
    digest = hashlib.sha256()
    for name in sorted(names):
        path = staged_inventory if name == 'api_inventory.csv' else root / name
        digest.update(name.encode() + b'\0' + path.read_bytes() + b'\0')
    return digest.hexdigest()


def stage(plan: Plan, approvals: Sequence[Candidate], destination: Path, *,
          jobs: int = 1, inventory_only: Sequence[str] = (),
          drop_inventory: bool = False) -> Path:
    """Stage and verify a complete proposed snapshot; never mutate ``plan.source``.

    Uses report's authenticated summaries and renderer against retained inputs.
    Stage execution may be expensive because the coverage audit inspects cases.
    """
    if jobs < 1:
        raise ValueError('jobs must be positive')
    root = Path(plan.source)
    if scan(root) != plan:
        raise ValueError('The source snapshot changed after the review')
    selected = {item.api for item in approvals}
    if not selected or selected - {item.api for item in plan.candidates}:
        raise ValueError('Stage selection contains an unreviewed API')
    extra = set(inventory_only)
    if extra - set(plan.unresolved_without_files):
        raise ValueError('Stage contains an unreviewed inventory-only API')
    if not drop_inventory and extra:
        raise ValueError('Inventory-only cleanup requires --drop-inventory')
    destination = destination.absolute()
    if destination.exists() or destination.parent != root:
        raise ValueError('Stage must be a new direct subdirectory of the source repository')
    destination.mkdir(parents=True, exist_ok=False)
    staged_inventory = destination / 'api_inventory.csv'
    if drop_inventory:
        with (root / 'api_inventory.csv').open(newline='', encoding='utf-8') as stream:
            reader = csv.DictReader(stream)
            fieldnames = reader.fieldnames
            if fieldnames is None:
                raise ValueError('Missing inventory header')
            filtered = [row for row in reader if row['api'] not in selected | extra]
        with staged_inventory.open('w', newline='', encoding='utf-8') as stream:
            writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator='\n')
            _ = cast(object, writer.writeheader())
            writer.writerows(filtered)
    else:
        os.link(root / 'api_inventory.csv', staged_inventory)
    inputs = root / 'cases/generated'
    measurements = root / 'results/measurements'
    staged_inputs = destination / 'cases/generated'
    staged_measurements = destination / 'results/measurements'
    excluded = {_key(api) for api in selected}
    _hardlink_survivors(inputs, staged_inputs, excluded)
    _hardlink_survivors(measurements, staged_measurements, excluded)
    original_cache = _json(inputs / '.cache.json')
    entries = _mapping(original_cache['files'], 'case cache files')
    original_cache['files'] = {key: value for key, value in entries.items() if key not in excluded}
    original_cache['format'] = 2
    original_cache['provenance'] = {
        'kind': 'reviewed_legacy_prune',
        'historical_cases_sha256': plan.cases_sha256,
        'historical_measurements_sha256': plan.measurements_sha256,
        'historical_file_count': len(entries),
        'removed_file_count': len(selected),
        'generator_equivalence_verified': False,
    }
    _write_json(staged_inputs / '.cache.json', original_cache)
    summary = root / 'results/accuracy-summaries.json'
    staged_summary = destination / 'results/accuracy-summaries.json'
    staged_summary.parent.mkdir(parents=True, exist_ok=True)
    os.link(summary, staged_summary)
    new_case_digest = case_directory_digest(staged_inputs)
    new_measurement_digest = directory_digest(staged_measurements)
    prior_provenance = _json(root / 'results/reference-provenance.json')
    prior_entries = _mapping(prior_provenance['apis'], 'reference provenance entries')
    prior_hashes = _mapping(prior_provenance['reference_sha256'], 'reference digests')
    from .reference_integrity import serialize_provenance

    updated_provenance = serialize_provenance(
        cast(dict[str, Sequence[ProvenanceEntry]], {key: value for key, value in prior_entries.items()
                                                  if key not in selected}),
        inputs_sha256=new_case_digest, measurements_sha256=new_measurement_digest,
        hashes={key: cast(str, value) for key, value in prior_hashes.items()
                if key not in selected},
    )
    provenance_path = destination / 'results/reference-provenance.json'
    _ = provenance_path.write_bytes(updated_provenance)
    from unittest.mock import patch

    from . import case_set, coverage_audit, report

    load_cases = cast(Callable[..., case_set.CachedCaseSet], case_set.historical_cases)
    cases = load_cases(staged_inputs, expected_cases_sha256=new_case_digest,
                       measurements_directory=staged_measurements,
                       expected_measurements_sha256=new_measurement_digest)
    records = case_csv.MeasurementStore(staged_measurements)

    # Reference hashes for the retained API names are immutable: surviving
    # measurement files are hardlinked, and the source manifest was verified.
    audit_cases = cast(Callable[..., str], coverage_audit.audit)
    with patch.object(coverage_audit, 'inventory', lambda: _inventory(destination)):
        coverage = audit_cases(cases, records)
    _ = (destination / 'results/coverage.csv').write_text(coverage, encoding='utf-8')
    meta = _json(root / 'results/environment.json')
    meta['cases_sha256'] = new_case_digest
    meta['measurements_sha256'] = new_measurement_digest
    meta['coverage_sha256'] = _digest(coverage.encode())
    meta['reference_provenance_sha256'] = _digest(updated_provenance)
    meta['testing_sources_sha256'] = _staged_source_digest(root, staged_inventory)
    if 'sci_measurement_commits' in meta:
        revisions = meta['sci_measurement_commits']
        if not isinstance(revisions, dict):
            raise ValueError('Invalid Sci measurement revisions in source environment')
        typed_revisions = cast(dict[str, object], revisions)
        # An API with no retained measurement must not retain a live
        # per-function revision. The original revisions remain in Git.
        meta['sci_measurement_commits'] = {
            api: revision for api, revision in typed_revisions.items()
            if api not in selected | extra
        }
    meta['prune_history'] = {'kind': 'reviewed_legacy_prune',
                             'prior_cases_sha256': plan.cases_sha256,
                             'prior_measurements_sha256': plan.measurements_sha256,
                             'removed_apis': sorted(selected),
                             'removed_inventory_only_apis': sorted(extra),
                             'removed_cases': sum(item.cases for item in approvals),
                             'reconciliation': 'reviewed explicit API allowlist; old data in Git'}
    output = destination / 'docs/accuracy'
    with patch.object(coverage_audit, 'inventory', lambda: _inventory(destination)):
        if jobs > 1:
            renderer_parallel = cast(Callable[..., None], report.render_cached_parallel)
            renderer_parallel(cases, records, meta, output, jobs=jobs,
                              root=destination, inputs=staged_inputs)
        else:
            renderer = cast(Callable[..., None], report.render_cached)
            renderer(cases, records, meta, output,
                     root=destination, inputs=staged_inputs)
    publish_summary = cast(Callable[..., None], report.publish_summary)
    publish_summary(output, root=destination)
    # The accuracy README renders version dictionaries with their native
    # insertion order. Sorting nested environment keys here would change that
    # order after reload and make `report --check` reject the rendered README.
    _write_json(destination / 'results/environment.json', meta, sort_keys=False)
    # A cached file snapshot can validate the staged output without invoking
    # current generators, Sci or any reference engine.
    load_snapshot = cast(Callable[..., object], report.load_snapshot)
    _ = load_snapshot(cases, root=destination, cached=True)
    if scan(root) != plan:
        raise ValueError('The source snapshot changed while staging; refuse publishing')
    staged_record = {
        'schema': 1, 'source': str(root), 'approved': sorted(selected),
        'plan': {key: getattr(plan, key) for key in
                 ('environment_sha256', 'inventory_sha256', 'cache_sha256')},
        'staged': {'cases_sha256': new_case_digest,
                   'inventory_sha256': _sha256(staged_inventory),
                   'measurements_sha256': new_measurement_digest,
                   'environment_sha256': _sha256(destination / 'results/environment.json'),
                   'cache_sha256': _sha256(staged_inputs / '.cache.json'),
                   'coverage_sha256': meta['coverage_sha256'],
                   'reference_provenance_sha256': meta['reference_provenance_sha256'],
                   'accuracy_summaries_sha256': cast(str, meta['accuracy_summaries_sha256'])},
        'source_input_files': len(case_paths(inputs)),
        'surviving_input_files': len(case_paths(staged_inputs)),
        'removed_inventory_only': sorted(extra),
        'drop_inventory': drop_inventory,
        'report_files': {page.name: _sha256(page) for page in sorted(output.glob('*.md'))},
    }
    _write_json(destination / 'prune-stage.json', staged_record)
    return destination


def _check_stage(plan: Plan, staged: Path) -> dict[str, object]:
    manifest = _json(staged / 'prune-stage.json')
    if manifest.get('schema') != 1 or manifest.get('source') != plan.source:
        raise ValueError('Unknown or foreign prune stage')
    expected = {key: getattr(plan, key) for key in
                ('environment_sha256', 'inventory_sha256', 'cache_sha256')}
    if manifest.get('plan') != expected:
        raise ValueError('Stage no longer matches the active source snapshot')
    digests = _mapping(manifest.get('staged'), 'staged digests')
    inputs = staged / 'cases/generated'
    measured = staged / 'results/measurements'
    for key, path in [('cache_sha256', inputs / '.cache.json'),
                      ('inventory_sha256', staged / 'api_inventory.csv'),
                      ('environment_sha256', staged / 'results/environment.json'),
                      ('coverage_sha256', staged / 'results/coverage.csv'),
                      ('reference_provenance_sha256', staged / 'results/reference-provenance.json'),
                      ('accuracy_summaries_sha256', staged / 'results/accuracy-summaries.json')]:
        if _sha256(path) != digests.get(key):
            raise ValueError(f'Tampered or incomplete stage file: {path}')
    if case_directory_digest(inputs) != digests.get('cases_sha256') or directory_digest(measured) != digests.get('measurements_sha256'):
        raise ValueError('Staged raw data no longer match reviewed digests')
    approved_value = manifest.get('approved')
    extra_value = manifest.get('removed_inventory_only')
    if not isinstance(approved_value, list) or not isinstance(extra_value, list):
        raise TypeError('Stage is missing its reviewed API selection')
    approved_list = cast(list[object], approved_value)
    extra_list = cast(list[object], extra_value)
    if (any(not isinstance(x, str) for x in approved_list + extra_list)
            or len(set(map(str, approved_list))) != len(approved_list)
            or len(set(map(str, extra_list))) != len(extra_list)):
        raise ValueError('Invalid or duplicate staged API removals')
    approved = cast(set[str], set(approved_list))
    extra = cast(set[str], set(extra_list))
    if approved - {item.api for item in plan.candidates} or extra - set(plan.unresolved_without_files):
        raise ValueError('Stage tries to remove an unreviewed or mapped API')
    active_inputs = set(case_paths(Path(plan.source) / 'cases/generated'))
    if set(case_paths(inputs)) != active_inputs - {_key(api) for api in approved}:
        raise ValueError('Staged case files differ from reviewed removals')
    if set(_measurement_paths(staged)) != active_inputs - {_key(api) for api in approved}:
        raise ValueError('Staged measurements differ from reviewed removals')
    unchanged = {row['api']: row for row in _inventory(Path(plan.source))
                 if row['api'] not in (approved | extra if manifest.get('drop_inventory') else set())}
    staged_rows = {row['api']: row for row in _inventory(staged)}
    if unchanged != staged_rows:
        raise ValueError('Staged inventory differs from reviewed removals')
    reports = _mapping(manifest.get('report_files'), 'rendered report files')
    actual_reports = {page.name: _sha256(page) for page in (staged / 'docs/accuracy').glob('*.md')}
    if actual_reports != reports:
        raise ValueError('Staged report pages changed after rendering')
    staged_meta = _json(staged / 'results/environment.json')

    if staged_meta.get('testing_sources_sha256') != _staged_source_digest(Path(plan.source),
                                                                          staged / 'api_inventory.csv'):
        raise ValueError('Source code or inventory changed since report staging')
    return manifest


def publish(plan: Plan, staged: Path) -> None:
    """Publish a verified stage with old paths preserved for rollback on errors.

    A multi-directory filesystem update cannot be a single atomic operation.
    Each rename is atomic; a synchronous exception rolls back all completed
    renames.  As with other multi-file Git workflows, an abrupt power loss
    between renames requires restoring the recorded backups before retrying.
    """
    root = Path(plan.source)
    if scan(root) != plan:
        raise ValueError('Active snapshot changed since dry run; publish refused')
    staged = staged.resolve(strict=True)
    _ = _check_stage(plan, staged)
    if root not in staged.parents:
        raise ValueError('Publish stages must be stored within the same repository filesystem')
    backup = staged / 'prune-backup'
    if backup.exists():
        raise ValueError('Stage already has a backup; restore or review its previous transaction')
    backup.mkdir()
    moves = [
        ('cases/generated', True),
        ('results/measurements', True),
        ('api_inventory.csv', False),
        ('results/reference-provenance.json', False),
        ('results/accuracy-summaries.json', False),
        ('results/coverage.csv', False),
        ('docs/accuracy', True),
        ('results/environment.json', False),
    ]
    done: list[tuple[Path, Path, Path]] = []
    try:
        for name, _ in moves:
            current = root / name
            proposed = staged / name
            previous = backup / name
            previous.parent.mkdir(parents=True, exist_ok=True)
            if not proposed.exists() or not current.exists():
                raise ValueError(f'Missing old/new publish artifact: {name}')
            _ = current.replace(previous)
            try:
                _ = proposed.replace(current)
            except BaseException:
                _ = previous.replace(current)
                raise
            done.append((current, proposed, previous))
    except BaseException:
        for current, proposed, previous in reversed(done):
            _ = current.replace(proposed)
            _ = previous.replace(current)
        raise
    # Retain original files in a clearly named backup for manual verification;
    # their bytes may include local uncommitted measurements.  Git can remove
    # those historical files only after the operator verifies this transaction.


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    _ = parser.add_argument('--root', type=Path, default=ROOT)
    _ = parser.add_argument('--write-review-template', type=Path)
    _ = parser.add_argument('--approved', type=Path, help='reviewed hash-pinned allowlist JSON')
    _ = parser.add_argument('--stage', type=Path, help='new isolated stage directory')
    _ = parser.add_argument('--publish', type=Path, help='publish an already rendered and checked stage')
    _ = parser.add_argument('--jobs', type=int, default=1)
    _ = parser.add_argument('--drop-inventory', action='store_true',
                            help='remove explicitly reviewed inventory aliases from the staged snapshot')
    args = parser.parse_args(argv, namespace=_Arguments())
    if args.stage and args.publish or args.write_review_template and (args.stage or args.publish):
        parser.error('Choose only one of review-template, stage or publish')
    plan = scan(args.root)
    if args.write_review_template:
        if args.write_review_template.exists():
            parser.error('Refusing to overwrite a review template')
        _write_json(args.write_review_template, review_template(plan))
        print(f'Wrote UNAPPROVED hash-pinned review template: {args.write_review_template}')
    elif args.stage:
        if args.approved is None:
            parser.error('--stage requires --approved review.json')
        approval = _json(args.approved)
        selected = approved_candidates(plan, approval)
        no_files = approved_inventory_only(plan, approval)
        target = stage(plan, selected, args.stage, jobs=args.jobs,
                       inventory_only=no_files, drop_inventory=args.drop_inventory)
        print(f'Staged {len(selected)} reviewed APIs and {sum(c.cases for c in selected)} cases at {target}')
        print('Active cases, measurements and reports were not modified.')
    elif args.publish:
        if args.approved is not None:
            parser.error('--approved is used with --stage, not --publish')
        publish(plan, args.publish)
        print(f'Published reviewed prune stage {args.publish}; original paths retained in prune-backup/')
    elif args.approved is not None:
        parser.error('--approved requires --stage')
    else:
        print(json.dumps(plan.as_json(), indent=2, sort_keys=True))
        print('Review required: no unresolved legacy pair is automatically considered obsolete.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
