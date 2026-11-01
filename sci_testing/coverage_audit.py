"""Audit measured successes against the compiler-resolved public function list."""
from collections import Counter, defaultdict
import csv
from .case_csv import table
import hashlib
import json

from . import ROOT


def inventory():
    with (ROOT/'api_inventory.csv').open(newline='') as stream:return list(csv.DictReader(stream))


def audit(cases,records,*,required_apis=(),require_complete=False):
    """Report incomplete 0.3 migration coverage explicitly, with optional strictness.

    Pending public APIs appear as rows with no Sci measurements. Historical
    aliases with unresolved mappings never contribute canonical coverage.
    """
    exports=inventory()
    aliases={row['api']:row['canonical'] for row in exports}
    if len(aliases)!=len(exports):raise ValueError('Duplicate public API in inventory')
    statuses={row['api']:row.get('status','') for row in exports}
    rows=[];successes=Counter()
    measured=set()
    measured_apis=set()
    io_rejects={'io.matlab.loadmat','io.matlab.whosmat',
                'io.matrix_market.mminfo','io.matrix_market.mmread',
                'io.wavfile.read','io.wavfile.signal'}
    for api in cases.apis:
        batch=cases.function(api)
        if api not in aliases:raise ValueError(f'Recipe is not a public Sci API: {api}')
        if not aliases[api] and statuses[api]!='unresolved_legacy':
            raise ValueError(f'Unresolved public API mapping: {api}')
        methods=Counter(c['id'].split('/')[-1].split('-')[0] for c in batch)
        io_fixture=api.startswith('io.')
        structural_fixture=bool(batch) and all(case['args'].get('structural03_case') for case in batch)
        complex_fixture=bool(batch) and all(case['args'].get('special_complex03_case') for case in batch)
        if io_fixture:
            if any(not case['args'].get('io_case') or not isinstance(case['args'].get('payload_hex'),str)
                   for case in batch):
                raise ValueError(f'{api}: missing independent I/O input payloads')
            streams={case['args']['payload_hex'] for case in batch}
            reject=sum(case['args'].get('operation')=='reject' for case in batch)
            if len(batch)<8 or len(streams)<8 or (api in io_rejects and not reject):
                raise ValueError(f'{api}: insufficient I/O fixture diversity or rejection coverage')
        elif structural_fixture:
            streams=set();reject=sum(case['args'].get('operation')=='reject' for case in batch)
            if len(batch)<1000 or len({json.dumps(case['args'],sort_keys=True) for case in batch})<100:
                raise ValueError(f'{api}: insufficient Sci 0.3 structural sample diversity')
        elif complex_fixture:
            streams=set();reject=0
            if len(batch)<8 or len({json.dumps(case['args'],sort_keys=True) for case in batch})<8:
                raise ValueError(f'{api}: insufficient complex special-function samples')
        else:
            streams=set();reject=0
            if methods['linspace']!=1000 or methods['random']!=1000:
                raise ValueError(f'{api}: incomplete dense samples')
        valid=sum('error' not in records[c['id']]['libraries']['Sci'] for c in batch)
        if aliases[api]:successes[aliases[api]]+=valid
        if valid:measured.add(api)
        measured_apis.add(api)
        unique=len({hashlib.sha256(json.dumps(c['args'],sort_keys=True,separators=(',',':')).encode()).digest() for c in batch})
        state=('unresolved_legacy' if not aliases[api] else
               'legacy_recipe_pending' if statuses[api]=='legacy_recipe_pending' else
               'measured' if valid else 'failed')
        rows.append([api,aliases[api],len(batch),methods['linspace'],methods['random'],unique,
                     len(streams),reject,valid,len(batch)-valid,state])
    for row in exports:
        if row['api'] in measured_apis:continue
        if row.get('status') in {'new_api_untested','untested','unresolved_legacy'}:
            rows.append([row['api'],row['canonical'],0,0,0,0,0,0,0,0,
                         'pending_'+row['status']])
    required=set(required_apis)
    missing_required=sorted(required-measured)
    if missing_required:raise ValueError(f'New public APIs without successful comparisons: {missing_required}')
    failed=sorted(api for api in measured_apis if statuses[api]!='unresolved_legacy'
                  and aliases[api] and not successes[aliases[api]])
    if failed:raise ValueError(f'Public functions without a successful comparison: {failed}')
    if require_complete:
        pending=sorted(row['api'] for row in exports if statuses[row['api']] in
                       {'new_api_untested','untested','unresolved_legacy'} and row['api'] not in measured
                       or statuses[row['api']]=='legacy_recipe_pending')
        if pending:raise ValueError(f'Unverified public Sci APIs: {pending}')
    return table(['api','canonical','cases','linspace','random','unique_inputs',
                  'io_payloads','io_rejects','sci_successes','sci_failures','coverage_state'],rows)
