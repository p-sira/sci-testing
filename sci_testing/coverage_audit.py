"""Audit measured successes against the compiler-resolved public function list."""
from collections import Counter, defaultdict
import csv
from .case_csv import table
import hashlib
import json

from . import ROOT


def inventory():
    with (ROOT/'api_inventory.csv').open(newline='') as stream:return list(csv.DictReader(stream))


def audit(cases,records):
    aliases={row['api']:row['canonical'] for row in inventory()}
    rows=[];successes=Counter()
    for api in cases.apis:
        batch=cases.function(api)
        if api not in aliases:raise ValueError(f'Recipe is not a public Sci API: {api}')
        methods=Counter(c['id'].split('/')[-1].split('-')[0] for c in batch)
        if methods['linspace']!=1000 or methods['random']!=1000:raise ValueError(f'{api}: incomplete dense samples')
        valid=sum('error' not in records[c['id']]['libraries']['Sci'] for c in batch)
        successes[aliases[api]]+=valid
        unique=len({hashlib.sha256(json.dumps(c['args'],sort_keys=True,separators=(',',':')).encode()).digest() for c in batch})
        rows.append([api,aliases[api],len(batch),methods['linspace'],methods['random'],unique,valid,len(batch)-valid])
    missing=sorted(set(aliases.values())-{a for a,n in successes.items() if n})
    if missing:raise ValueError(f'Public functions without a successful comparison: {missing}')
    return table(['api','canonical','cases','linspace','random','unique_inputs','sci_successes','sci_failures'],rows)
