"""Regenerate curated primary reference values from pinned mathematical recipes.

This is an independent, opt-in reference-only operation. Sci and competitor
measurements are preserved verbatim. WolframEngine is an optional primary source
for explicitly supported functions, and neither its results nor historical
recipes are represented as certified or independently corroborated.

The provenance manifest distinguishes verified regeneration from retrospective
identification of historical recipe code paths. All data changes are staged
until the complete snapshot and report have been assembled.
"""

import argparse
import hashlib
import json
import os
import shutil
import tempfile
from collections import Counter
from collections.abc import Callable, Mapping, Sequence
from contextlib import AbstractContextManager
from decimal import Decimal, InvalidOperation
from importlib import import_module
from pathlib import Path
from typing import Protocol, TypedDict, cast

from . import ROOT, case_csv, case_set
from .case_csv import MeasurementStore
from .reference_integrity import (
    MANIFEST,
    ProvenanceEntry,
    ReferenceCase,
    ReferenceCases,
    ReferenceMeasurement,
    ReferenceRecords,
    load_prior_provenance,
    reference_hashes,
    serialize_provenance,
)


class _Mpmath(Protocol):
    def nstr(self, value: object, n: int) -> str: ...

    def workdps(self, digits: int) -> AbstractContextManager[object]: ...


# mpmath has no distributed type stubs; spell out the two methods used here.
mp = cast(_Mpmath, cast(object, import_module('mpmath')))


class _WolframBatch(Protocol):
    def supports(self, case: Mapping[str, object], /) -> bool: ...

    def references(self, batch: Sequence[Mapping[str, object]], /) -> list[list[str]]: ...


class _WolframProvider(_WolframBatch, Protocol):
    def reference(self, case: Mapping[str, object]) -> list[str]: ...

    def provenance(self, *, digits: int) -> dict[str, object]: ...

    def definition(self, case: Mapping[str, object]) -> str: ...


class _Selectors(Protocol):
    @property
    def apis(self) -> Sequence[str]: ...


class _Report(Protocol):
    METRIC_PRECISION_DPS: int

    def stable_reference(self, case: ReferenceCase) -> Sequence[object]: ...

    def load_stored_snapshot(self) -> tuple[dict[str, object], MeasurementStore, object]: ...

    def digest_sources(self) -> str: ...

    def preserved_coverage(self, previous: Mapping[str, object], *,
                           cases: ReferenceCases, records: ReferenceRecords) -> str: ...

    def render_cached(self, cases: ReferenceCases, records: ReferenceRecords,
                      meta: dict[str, object], output: Path) -> None: ...

    def render_cached_parallel(self, cases: ReferenceCases, records: ReferenceRecords,
                               meta: dict[str, object], output: Path, *, jobs: int) -> None: ...

    def publish_summary(self, pages: Path) -> None: ...


class _RecipeSource(TypedDict):
    provider: str
    recipe: str
    precision: str
    status: str


def _metadata_str(metadata: Mapping[str, object], key: str) -> str:
    value = metadata.get(key)
    if not isinstance(value, str):
        raise TypeError(f'Missing or invalid snapshot metadata: {key}')
    return value


def _finite_decimals(values: Sequence[str]) -> bool:
    """Check decimal references without binary64 conversion or rounding."""
    try:
        return all(Decimal(x).is_finite() for x in values)
    except InvalidOperation:
        return False


def _verified_reference_values(case: ReferenceCase, values: object) -> list[str]:
    if not isinstance(values, (tuple, list)) or not values:
        raise ValueError(f'Empty primary reference for {case["id"]}')
    result = [str(x) for x in cast(Sequence[object], values)]
    if not _finite_decimals(result):
        raise ValueError(f'Nonfinite primary reference for {case["id"]}')
    return result


def recipe_source(case: ReferenceCase) -> tuple[str, str, str]:
    """Document the actual recipe dispatcher, without inferring its accuracy."""
    api, args = case['api'], case['args']
    if api.startswith('complex.'):
        return ('mpmath', 'sci_testing.complex_accuracy.reference',
                '200/320 decimal digits; escalated near difficult inputs')
    if args.get('special_complex03_case'):
        return ('mpmath', 'sci_testing.special_complex03_cases.reference',
                '100/160 decimal digits')
    if args.get('special03_case'):
        return ('mpmath', 'sci_testing.special03_cases.reference',
                '100/160 decimal digits')
    if args.get('structural03_case'):
        return ('SciPy/NumPy plus explicit identities',
                'sci_testing.structural03_cases.reference', 'binary64')
    if args.get('io_case'):
        return ('SciPy/NumPy and explicit binary fixtures',
                'sci_testing.io_cases.reference', 'binary64')
    if args.get('explicit'):
        backend = ('NumPy/Scipy plus exact identities'
                   if api.startswith(('differentiate.', 'optimize.'))
                   else 'SciPy/NumPy plus explicit identities')
        return (backend, 'sci_testing.explicit_cases.reference', 'binary64 or exact')
    if any(key in args for key in ('p0', 'p', 'arguments', 'raw')):
        return ('SciPy/NumPy plus independent fixture recipes',
                'sci_testing.registered_cases.reference', 'binary64')
    return ('mpmath / analytic mathematical recipe',
            'sci_testing.report.high_precision', '100/160 decimal digits')


def source_entry(case: ReferenceCase, *, historical: bool = False) -> _RecipeSource:
    backend, method, precision = recipe_source(case)
    return {'provider': backend, 'recipe': method, 'precision': precision,
            'status': 'historical_source_path_inferred' if historical else 'regenerated'}


def provenance_for_batch(batch: Sequence[ReferenceCase], *, historical: bool = True) -> list[ProvenanceEntry]:
    kinds: Counter[tuple[str, str, str, str]] = Counter()
    for case in batch:
        entry = source_entry(case, historical=historical)
        kinds[(entry['provider'], entry['recipe'], entry['precision'],
               entry['status'])] += 1
    return [{'provider': provider, 'recipe': recipe,
                        'precision': precision, 'status': status, 'cases': count}
                       for (provider, recipe, precision, status), count
                       in sorted(kinds.items())]


def provenance_for_cases(cases: ReferenceCases, *, historical: bool = True) -> dict[str, list[ProvenanceEntry]]:
    """Summarize provenance by public API; mixed recipes stay distinguishable."""
    return {api: provenance_for_batch(cases.function(api), historical=historical)
            for api in cases.apis}


def _select(cases: _Selectors, selectors: Sequence[str]) -> list[str]:
    result: set[str] = set()
    for selector in selectors:
        matches = {api for api in cases.apis
                   if api == selector or api.startswith(selector + '.')}
        if not matches:
            raise ValueError(f'No saved API matches --module {selector!r}')
        result.update(matches)
    return sorted(result)


def _reference_strings(case: ReferenceCase, source: str,
                       wolfram: _WolframProvider | None) -> list[str]:
    if source == 'wolfram':
        if wolfram is None:
            raise ValueError('Wolfram provider is unavailable')
        if not wolfram.supports(case):
            raise ValueError(f'No curated Wolfram expression for {case["id"]}')
        values = wolfram.reference(case)
    else:
        report = cast(_Report, cast(object, import_module('sci_testing.report')))
        values = [mp.nstr(x, 165) for x in report.stable_reference(case)]
    return _verified_reference_values(case, values)


def _regenerate_batch(batch: Sequence[ReferenceCase], records: ReferenceRecords,
                      source: str, wolfram: _WolframBatch | None) -> dict[str, ReferenceMeasurement]:
    """Never allow partial function datasets or changed library measurements."""
    new: dict[str, ReferenceMeasurement] = {}
    if source == 'wolfram':
        if wolfram is None:
            raise ValueError('Wolfram provider is unavailable')
        for case in batch:
            if not wolfram.supports(case):
                raise ValueError(f'No curated Wolfram expression for {case["id"]}')
        references = wolfram.references(batch)
        if len(references) != len(batch):
            raise ValueError('WolframEngine returned the wrong number of cases')
    else:
        references = [_reference_strings(case, source, None) for case in batch]
    for case, values in zip(batch, references):
        original = records[case['id']]
        values = [str(x) for x in values]
        if (len(values) != len(original['reference']) or not values
                or not _finite_decimals(values)):
            raise ValueError(f'Invalid reference shape or value: {case["id"]}')
        new[case['id']] = {'id': case['id'], 'reference': values,
                           'libraries': original['libraries']}
    return new


def regenerate(selectors: Sequence[str], *, source: str = 'recipe', jobs: int = 1,
               dry_run: bool = False) -> dict[str, list[str] | str | int]:
    """Stage and publish a reproducible reference-only snapshot with provenance."""
    if source not in ('recipe', 'wolfram'):
        raise ValueError('Unknown reference source')
    if jobs < 1:
        raise ValueError('jobs must be positive')
    report = cast(_Report, cast(object, import_module('sci_testing.report')))
    previous, records, _ = report.load_stored_snapshot()
    cases = cast(Callable[..., ReferenceCases], case_set.historical_cases)(
        ROOT/'cases/generated',
        expected_cases_sha256=_metadata_str(previous, 'cases_sha256'),
        measurements_directory=ROOT/'results/measurements',
        expected_measurements_sha256=_metadata_str(previous, 'measurements_sha256'))
    selected = _select(cases, selectors)
    wolfram: _WolframProvider | None = (
        cast(_WolframProvider, cast(object, import_module('sci_testing.wolfram_reference')))
        if source == 'wolfram' else None)
    if source == 'wolfram':
        assert wolfram is not None
        for api in selected:
            for case in cases.function(api):
                if not wolfram.supports(case):
                    raise ValueError(f'No curated Wolfram expression for {case["id"]}')
    if dry_run:
        return {'selected': selected,
                'cases': sum(len(cases.function(api)) for api in selected),
                'source': source}

    details: dict[str, object] = {'provider': 'existing mathematical recipes',
                                  'recipe': 'sci_testing.report.stable_reference',
                                  'precision': 'per-case recipe and stability check'}
    if source == 'wolfram':
        assert wolfram is not None
        details = wolfram.provenance(digits=160)
    prior_provenance = load_prior_provenance(previous, ROOT/'results'/MANIFEST)
    if prior_provenance is None:
        entries = provenance_for_cases(cases)
    else:
        entries = prior_provenance['apis'].copy()
    if set(entries) != set(cases.apis):
        raise ValueError('Reference provenance API inventory differs from inputs')
    record_lookup = cast(ReferenceRecords, records)
    if prior_provenance is not None and (reference_hashes(cases, record_lookup)
                                        != prior_provenance['reference_sha256']):
        raise ValueError('Previously recorded primary references have changed')
    with tempfile.TemporaryDirectory(dir=ROOT/'results') as tmp:
        stage = Path(tmp)
        staged = stage/'measurements'
        _ = shutil.copytree(ROOT/'results/measurements', staged, copy_function=os.link)
        changed: set[str] = set()
        for api in selected:
            batch = cases.function(api)
            outputs = _regenerate_batch(batch, record_lookup, source, wolfram)
            file_map = cast(
                Callable[[Sequence[ReferenceCase], Mapping[str, ReferenceMeasurement]], dict[str, str]],
                case_csv.measurement_files)(batch, outputs)
            for logical in file_map:
                (staged/logical).unlink(missing_ok=True)
                (staged/(logical+'.gz')).unlink(missing_ok=True)
            changed.update(cast(Callable[[Path, Mapping[str, str]], set[str]],
                                case_csv.write_measurements)(staged, file_map))
            if source == 'wolfram':
                assert wolfram is not None
                digits = details.get('digits')
                if not isinstance(digits, int):
                    raise ValueError('Wolfram provider supplied no reference precision')
                new_entry: ProvenanceEntry = {'provider': 'WolframEngine',
                             'recipe': 'sci_testing.wolfram_reference.references',
                             'expression_template': wolfram.definition(batch[0]),
                             'documentation': 'https://reference.wolfram.com/language/guide/EllipticIntegrals.html'
                             if api.startswith('special.ellip') else
                             'https://reference.wolfram.com/language/guide/MathematicalFunctions.html',
                             'precision': f'{digits} decimal digits',
                             'status': 'regenerated', 'cases': len(batch)}
                entries[api] = [new_entry]
            else:
                entries[api] = provenance_for_batch(batch, historical=False)
            print(f'Regenerated reference: {api} ({len(batch)} cases, {source})', flush=True)

        meta = previous.copy()
        meta['testing_sources_sha256'] = report.digest_sources()
        meta['reference_environment'] = previous.get('reference_environment', previous)
        meta['measurements_sha256'] = cast(Callable[[Path], str],case_csv.directory_digest)(staged)
        refreshes = meta.setdefault('reference_refreshes', [])
        if not isinstance(refreshes, list):
            raise TypeError('Reference refresh history is invalid')
        cast(list[object], refreshes).append({
            'apis': selected, 'source': source, 'details': details,
            'cases_sha256': _metadata_str(previous, 'cases_sha256'),
            'measurements_sha256': meta['measurements_sha256'],
            'testing_sources_sha256': meta['testing_sources_sha256']})
        provenance = serialize_provenance(entries,
            inputs_sha256=_metadata_str(meta, 'cases_sha256'),
            measurements_sha256=_metadata_str(meta, 'measurements_sha256'),
            hashes=reference_hashes(cases, cast(ReferenceRecords, MeasurementStore(staged))))
        meta['reference_provenance_sha256'] = hashlib.sha256(provenance).hexdigest()
        _ = (stage/MANIFEST).write_bytes(provenance)

        staged_records = MeasurementStore(staged)
        coverage = report.preserved_coverage(previous, cases=cases,
                                              records=cast(ReferenceRecords, staged_records))
        meta['coverage_sha256'] = hashlib.sha256(coverage.encode()).hexdigest()
        pages = stage/'pages'
        if jobs > 1:
            report.render_cached_parallel(cases, staged_records, meta, pages,
                                          jobs=jobs)
        else:
            with mp.workdps(report.METRIC_PRECISION_DPS):
                report.render_cached(cases, staged_records, meta, pages)

        # Do not publish until all selected references, provenance, and summaries
        # have been successfully calculated from the staged measurements.
        for name in sorted(changed):
            src = staged/name
            dst = ROOT/'results/measurements'/name
            dst.parent.mkdir(parents=True, exist_ok=True)
            _ = src.replace(dst)
            other = (dst.with_suffix('') if dst.suffix == '.gz'
                     else dst.with_suffix('.csv.gz'))
            other.unlink(missing_ok=True)
        _ = (ROOT/'results/coverage.csv').write_text(coverage)
        _ = (stage/MANIFEST).replace(ROOT/'results'/MANIFEST)
        report.publish_summary(pages)
        _ = (ROOT/'results/environment.json').write_text(json.dumps(meta, indent=2)+'\n')
        destination = ROOT/'docs/accuracy'
        destination.mkdir(parents=True, exist_ok=True)
        for page in pages.glob('*.md'):
            _ = page.replace(destination/page.name)
    return {'selected': selected, 'source': source,
            'cases': sum(len(cases.function(api)) for api in selected)}


class _Arguments(argparse.Namespace):
    def __init__(self) -> None:
        super().__init__()
        self.module: list[str] = []
        self.source: str = 'recipe'
        self.jobs: int = 1
        self.dry_run: bool = False


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    _ = parser.add_argument('--module', action='append', required=True,
                        help='saved API or namespace to regenerate; may repeat')
    _ = parser.add_argument('--source', choices=('recipe','wolfram'),default='recipe',
                        help='primary reference generator (no automatic fallback)')
    _ = parser.add_argument('--jobs', type=int, default=1,
                        help='parallel processes for full comparison report rendering')
    _ = parser.add_argument('--dry-run',action='store_true',
                        help='inspect selected saved APIs without changing data')
    args = parser.parse_args(namespace=_Arguments())
    print(json.dumps(regenerate(args.module, source=args.source,
                                jobs=args.jobs, dry_run=args.dry_run),indent=2))


if __name__ == '__main__':
    _ = main()
