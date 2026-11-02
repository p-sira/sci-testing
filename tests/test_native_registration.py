"""Regression guard for *previously measured* native comparison adapters.

The checked-in measurements, not the current C++ registry or a hand-maintained
list of supposed unavailable functions, supply the expectations. This detects
accidental loss of a working adapter during refactors. It does not discover
unimplemented Boost.Math/GSL operations or certify exhaustive coverage;
independent numerical tests live in each adapter's own test module.
"""

import csv
import gzip
import io
import json
import math
import unittest

from sci_testing import ROOT
from sci_testing.case_csv import case_paths, csv_text, decode
from sci_testing.native_compare import NativeComparisons


# These are the distinct *measured* native providers, not an unsupported-API
# taxonomy. The saved measurement snapshot independently determines the
# positive regression baselines for each provider/API pair.
PROVIDERS = frozenset({'Boost.Math', 'GSL', 'Boost.Graph',
                       'Boost.Numeric.uBLAS', 'Boost.Multiprecision',
                       'Boost.Units'})


def _saved_successful_samples():
    """One known-successful fixture per API/provider from the pinned snapshot."""
    entries = json.loads((ROOT / 'results/accuracy-summaries.json').read_text())['entries']
    inputs = case_paths(ROOT / 'cases/generated')
    measurements = case_paths(ROOT / 'results/measurements')
    for api, entry in sorted(entries.items()):
        providers = set(entry['key']['libraries']) & PROVIDERS
        if not providers:
            continue
        # Function names may themselves contain dots (e.g. special.bessel.iv_ratio).
        module, name = api.rsplit('.', 1)
        path = module + '/' + name + '.csv'
        if path not in inputs or path not in measurements:
            raise AssertionError(f'{api}: measured native adapter has no saved fixtures')
        selected = {}
        for row in csv.DictReader(io.StringIO(csv_text(measurements[path]))):
            provider = row['library']
            if (provider in providers and provider not in selected
                    and not row['error'] and row['values']):
                selected[provider] = (row['id'], len(row['values'].split()))
                if len(selected) == len(providers):
                    break
        # A column with no successful measurement cannot establish a working
        # adapter; retain its failures in the ordinary report, not this guard.
        if not selected:
            continue
        ids = {identifier for identifier, _ in selected.values()}
        found = {}
        source = inputs[path]
        opener = gzip.open if source.suffix == '.gz' else open
        with opener(source, 'rt', encoding='utf-8', newline='') as stream:
            for row in csv.DictReader(stream):
                if row['id'] in ids:
                    found[row['id']] = {
                        'api': api,
                        'args': {name: decode(value) for name, value in row.items()
                                 if name not in ('id', 'group')},
                    }
                    if len(found) == len(ids):
                        break
        if ids != set(found):
            raise AssertionError(f'{api}: saved successful input IDs not found')
        for provider, (identifier, width) in selected.items():
            yield api, provider, width, found[identifier]


class NativeRegistrationTests(unittest.TestCase):
    def test_previously_successful_comparators_remain_callable(self):
        native = NativeComparisons()
        checked = 0
        failures = []
        for api, provider, width, case in _saved_successful_samples():
            try:
                actual = (native.evaluate_complex(case)
                          if api.startswith('special.complex.') else native.evaluate(case))
                result = actual.get(provider)
                if (not result or result.get('error') or
                        len(result.get('values', ())) != width or
                        not all(math.isfinite(value) for value in result['values'])):
                    failures.append(f'{api} / {provider}: {result!r} (expected {width} values)')
            except (TypeError, ValueError, OSError) as exc:
                failures.append(f'{api} / {provider}: {exc!r}')
            checked += 1
        self.assertGreater(checked, 300, 'native snapshot unexpectedly has too few comparisons')
        self.assertFalse(failures, '\n' + '\n'.join(failures))


if __name__ == '__main__':
    unittest.main()
