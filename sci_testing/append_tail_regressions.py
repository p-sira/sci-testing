"""Append independent tail fixtures to an authenticated, current Sci snapshot.

No historical input, reference or library result is regenerated. All changes
are staged and validated before replacing saved files.
"""
import argparse
from collections import defaultdict
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import mpmath as mp

from . import ROOT
from . import report, tail_regressions
from .accuracy_cache import sha256_file
from .case_csv import (case_directory_digest, case_files, case_paths, csv_text,
                       directory_digest, MeasurementStore, write_measurements)
from .case_set import cached_cases, historical_cases
from .coverage_audit import audit
from .generate_cases import case_signature
from .reference_generation import provenance_for_batch
from .reference_integrity import (MANIFEST, load_prior_provenance,
                                  reference_hashes, serialize_provenance)


def append_csv(old, new, *, measurements=False):
    """Retain exact old field strings; allow extra measurement warning columns."""
    old_reader = csv.DictReader(io.StringIO(old))
    new_reader = csv.DictReader(io.StringIO(new))
    old_fields = old_reader.fieldnames
    new_fields = new_reader.fieldnames
    if not old_fields or not new_fields:
        raise ValueError('Missing CSV headers')
    extra = set(new_fields) ^ set(old_fields)
    if extra and (not measurements or any(not k.startswith('warning_') for k in extra)):
        raise ValueError('Regression CSV schema differs from historical schema')
    old_rows, new_rows = list(old_reader), list(new_reader)
    if {r['id'] for r in old_rows} & {r['id'] for r in new_rows}:
        raise ValueError('Duplicate historical regression case ID')
    if not new_rows:
        raise ValueError('No regression rows')
    keys = ('id', 'library') if measurements else ('id',)
    identifiers = [tuple(r[k] for k in keys) for r in new_rows]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError('Duplicate regression rows')
    if old_fields == new_fields:
        # Retain all historical bytes including quoting and decimal spelling.
        return old + new.partition('\n')[2]
    fields = old_fields + [k for k in new_fields if k not in old_fields]
    output = io.StringIO(newline='')
    writer = csv.DictWriter(output, fields, lineterminator='\n')
    writer.writeheader()
    writer.writerows(old_rows)
    writer.writerows(new_rows)
    return output.getvalue()


def append_files(source, destination, additions, *, measurements=False):
    """Change only staged function files; never truncate historical hard links."""
    paths = case_paths(source)
    for name, content in additions.items():
        old = paths[name]
        combined = append_csv(csv_text(old), content, measurements=measurements)
        for path in (destination/name, destination/(name+'.gz')):
            path.unlink(missing_ok=True)
        write_measurements(destination, {name: combined})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    previous, records, _ = report.load_stored_snapshot()
    results = ROOT/'results'
    inputs = ROOT/'cases/generated'
    historic = historical_cases(inputs,
        expected_cases_sha256=previous['cases_sha256'],
        measurements_directory=records.directory,
        expected_measurements_sha256=previous['measurements_sha256'])
    pinned = (ROOT/'sci-revision.txt').read_text().strip()
    current = report.environment(report.SCI)
    if current['sci_commit'] != pinned or previous['sci_commit'] != pinned or current['sci_source_dirty']:
        raise ValueError('Complete a clean Sci-only refresh at the pinned revision first')
    if any(revision != pinned for revision in previous.get('sci_measurement_commits', {}).values()):
        raise ValueError('Historical Sci measurements contain mixed revisions')
    for key in ('numpy', 'scipy', 'mpmath', 'filterpy', 'python', 'rust', 'cargo', 'profile'):
        if current[key] != previous[key]:
            raise ValueError(f'Regression environment differs from Sci refresh: {key}')
    if any(k.startswith(('RUSTFLAGS', 'CARGO_ENCODED_RUSTFLAGS', 'CARGO_PROFILE_')) or k == 'CARGO_BUILD_TARGET' for k in os.environ):
        raise ValueError('Unset custom compiler/profile/target overrides')
    prior = load_prior_provenance(previous, results/MANIFEST)
    if prior is None or prior['measurements_sha256'] != previous['measurements_sha256'] or set(prior['apis']) != set(historic.apis):
        raise ValueError('Authenticated historical reference provenance is required')
    added = defaultdict(list)
    for case in tail_regressions.cases():
        api = case['api']
        if api not in historic.apis:
            raise ValueError(f'Regression API not already measured: {api}')
        added[api].append(case)
    if not added:
        raise ValueError('No new tail regressions')
    if reference_hashes(historic, records, apis=sorted(added)) != {api: prior['reference_sha256'][api] for api in added}:
        raise ValueError('Historical reference values changed')
    with tempfile.TemporaryDirectory(dir=results) as directory:
        stage = Path(directory)
        staged_inputs, staged = stage/'cases', stage/'measurements'
        shutil.copytree(inputs, staged_inputs, copy_function=os.link)
        shutil.copytree(records.directory, staged, copy_function=os.link)
        additions = case_files([c for batch in added.values() for c in batch])
        append_files(inputs, staged_inputs, additions)
        paths = case_paths(staged_inputs)
        index = json.loads((inputs/'.cache.json').read_text())
        for api, batch in added.items():
            name = api.rsplit('.', 1)[0]+'/'+api.rsplit('.', 1)[1]+'.csv'
            index['files'][name] = {'count': historic.counts[name]+len(batch),
                                    'sha256': sha256_file(paths[name])}
        index.update(format=2, recipe_sha256=case_signature(), provenance={
            'kind': 'historical_snapshot_case_append',
            'historical_cases_sha256': previous['cases_sha256'],
            'historical_measurements_sha256': previous['measurements_sha256'],
            'generator_equivalence_verified': False,
            'added_cases': sum(map(len, added.values()))})
        manifest = staged_inputs/'.cache.json'
        manifest.unlink()
        manifest.write_text(json.dumps(index, sort_keys=True, indent=2)+'\n')
        cases = cached_cases(staged_inputs, case_signature())
        for api, batch in added.items():
            if cases.function(api) != historic.function(api)+batch:
                raise ValueError(f'Historical inputs or regression ordering changed: {api}')
        subprocess.run(['cargo', 'build', '--locked'], cwd=ROOT, check=True)
        target = Path(json.loads(report.command('cargo', 'metadata', '--format-version', '1', '--no-deps'))['target_directory'])
        new_measurements = stage/'new-measurements'
        native = report.NativeComparisons()
        new_libraries = report.extend_references(added, target/'debug/sci-testing', new_measurements, native=native)
        new_records = MeasurementStore(new_measurements)
        for batch in added.values():
            for case in batch:
                saved = new_records[case['id']]
                if 'error' in saved['libraries']['Sci'] or not saved['reference']:
                    raise ValueError(f'Failed Sci or reference regression: {case["id"]}: {saved}')
        append_files(records.directory, staged,
                     {name: csv_text(path) for name, path in case_paths(new_measurements).items()}, measurements=True)
        staged_records = MeasurementStore(staged)
        if reference_hashes(historic, staged_records, apis=sorted(added)) != {api: prior['reference_sha256'][api] for api in added}:
            raise ValueError('Regression extension changed historical primary references')
        meta = previous.copy()
        meta['testing_sources_sha256'] = report.digest_sources()
        meta['cases_sha256'] = case_directory_digest(staged_inputs)
        meta['measurements_sha256'] = directory_digest(staged)
        entries = prior['apis'].copy()
        for api, batch in added.items():
            entries[api] = entries[api]+provenance_for_batch(batch, historical=False)
        hashes = {**prior['reference_sha256'], **reference_hashes(cases, staged_records, apis=sorted(added))}
        provenance = serialize_provenance(entries, inputs_sha256=meta['cases_sha256'],
            measurements_sha256=meta['measurements_sha256'], hashes=hashes)
        meta['reference_provenance_sha256'] = hashlib.sha256(provenance).hexdigest()
        meta['module_libraries'] = {module: sorted(set(names)|new_libraries.get(module, set()))
                                    for module, names in previous['module_libraries'].items()}
        meta['tail_regression_extensions'] = [*previous.get('tail_regression_extensions', []), {
            'cases': sum(map(len, added.values())), 'apis': sorted(added),
            'sci_commit': pinned, 'testing_sources_sha256': meta['testing_sources_sha256'],
            'historical_cases_sha256': previous['cases_sha256'],
            'historical_measurements_sha256': previous['measurements_sha256'],
            'environment': {key: current[key] for key in ('numpy', 'scipy', 'mpmath', 'filterpy', 'python', 'rust', 'cargo', 'profile', 'rust_oracle_libraries', 'platform', 'library_config')},
            'native_libraries': native.versions,
            'cxx': report.command(os.environ.get('CXX', 'c++'), '--version')}]
        coverage = audit(cases, staged_records, required_apis=added)
        meta['coverage_sha256'] = hashlib.sha256(coverage.encode()).hexdigest()
        pages = stage/'pages'
        with mp.workdps(report.METRIC_PRECISION_DPS):
            report.render_cached(cases, staged_records, meta, pages, inputs=staged_inputs)
        if (case_directory_digest(inputs) != previous['cases_sha256'] or
                directory_digest(records.directory) != previous['measurements_sha256']):
            raise ValueError('Historical snapshot changed during regression extension')
        # Only touched CSVs and metadata are replaced; the historical input and
        # measurement files of every other API remain byte-identical.
        for source, destination in ((staged_inputs, inputs), (staged, records.directory)):
            files = case_paths(source)
            for name in additions:
                path = files[name]
                target = destination/path.relative_to(source)
                path.replace(target)
                alternate = target.with_suffix('') if target.suffix == '.gz' else target.with_suffix('.csv.gz')
                alternate.unlink(missing_ok=True)
        manifest.replace(inputs/'.cache.json')
        (results/MANIFEST).write_bytes(provenance)
        (results/'coverage.csv').write_text(coverage)
        report.publish_summary(pages)
        (results/'environment.json').write_text(json.dumps(meta, indent=2)+'\n')
        for page in pages.glob('*.md'):
            page.replace(ROOT/'docs/accuracy'/page.name)
    print(f'Appended {sum(map(len, added.values()))} cases to {len(added)} measured APIs')


if __name__ == '__main__':
    main()
