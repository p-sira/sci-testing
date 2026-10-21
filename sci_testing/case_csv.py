"""Compact, lossless CSV storage, partitioned by public module and function.

Numerical cells use s:<number>, v:<space-separated numbers>, or
m:<semicolon-separated rows>. Decimal repr preserves binary64 signed zero and
round trips. CSV quoting handles textual IDs, group names, errors and warnings.
"""
from collections import defaultdict
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path

# Flattened multidimensional outputs and high-precision references exceed the
# csv module's 128 KiB default. Keep an explicit bound on individual cells.
csv.field_size_limit(16 * 1024 * 1024)


def table(columns, rows):
    stream = io.StringIO(newline='')
    writer = csv.writer(stream, lineterminator='\n')
    writer.writerow(columns)
    writer.writerows(rows)
    return stream.getvalue()


def encode(value):
    if value is None: return 'n:'
    if isinstance(value,str): return 't:' + value
    if isinstance(value,dict): return 'j:' + json.dumps(value,separators=(',',':'),allow_nan=False)
    if isinstance(value,list) and any(not isinstance(x,(int,float)) for x in value):
        if not all(isinstance(row,list) and all(isinstance(x,(int,float)) for x in row) for row in value):
            return 'j:' + json.dumps(value,separators=(',',':'),allow_nan=False)
    if isinstance(value, bool): return 'b:' + str(int(value))
    if isinstance(value, int): return 'i:' + str(value)
    if isinstance(value, list):
        if value and isinstance(value[0], list):
            return 'm:' + ';'.join(' '.join(map(repr, row)) for row in value)
        return 'v:' + ' '.join(map(repr, value))
    if not math.isfinite(value):
        raise ValueError('Inputs must be finite')
    return 's:' + repr(value)


def parse_number(text):
    if text in ['True','False']: return text == 'True'
    if text.lstrip('-').isdigit(): return int(text)
    return float(text)


def decode(cell):
    kind, text = cell.split(':', 1)
    if kind == 'n': return None
    if kind == 't': return text
    if kind == 'j': return json.loads(text)
    if kind == 'b':
        if text not in ['0', '1']: raise ValueError('Invalid boolean')
        return text == '1'
    if kind == 'i': return int(text)
    if kind == 's': value = float(text)
    elif kind == 'v': value = list(map(parse_number, text.split()))
    elif kind == 'm': value = [list(map(parse_number, row.split())) for row in text.split(';')]
    else: raise ValueError(f'Unknown numerical cell kind: {kind}')
    numbers = [value] if kind == 's' else value if kind == 'v' else [x for row in value for x in row]
    if not all(math.isfinite(x) for x in numbers): raise ValueError('Inputs must be finite')
    return value


def iter_case_files(cases):
    from itertools import groupby
    ordered=cases if hasattr(cases,'function') else sorted(cases,key=lambda c:(c['api'],c['id']))
    for api, iterator in groupby(ordered,key=lambda c:c['api']):
        batch=list(iterator)
        module, function = api.rsplit('.', 1)
        arguments = sorted(batch[0]['args'])
        if any(sorted(c['args']) != arguments for c in batch):
            raise ValueError(f'Inconsistent arguments for {api}')
        yield f'{module}/{function}.csv', table(['id', 'group', *arguments],
            ([c['id'], c['group'], *(encode(c['args'][a]) for a in arguments)] for c in batch))


def case_files(cases):return dict(iter_case_files(cases))


def write_cases(directory,cases):
    directory=Path(directory);names=set()
    for name,content in iter_case_files(cases):
        path=directory/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text(content);names.add(name)
    for path in directory.glob('*/*.csv'):
        if path.relative_to(directory).as_posix() not in names:path.unlink()


def directory_digest(directory):
    directory=Path(directory);digest=hashlib.sha256()
    for path in sorted(csv_paths(directory),key=lambda p:p.relative_to(directory).as_posix()):
        digest.update(path.relative_to(directory).as_posix().encode()+b'\0'+path.read_bytes()+b'\0')
    return digest.hexdigest()


def read_cases(directory):
    cases = []
    for path in sorted(Path(directory).glob('*/*.csv'),key=lambda p:f'{p.parent.name}.{p.stem}'):
        api = f'{path.parent.name}.{path.stem}'
        with path.open(newline='') as stream:
            for row in csv.DictReader(stream):
                cases.append(dict(id=row.pop('id'), api=api, group=row.pop('group'),
                                  args={k: decode(v) for k, v in row.items()}))
    return cases


def write_files(directory, files):
    directory = Path(directory)
    for name, content in files.items():
        path = directory / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
    for path in directory.glob('*/*.csv'):
        if path.relative_to(directory).as_posix() not in files: path.unlink()


def files_digest(files):
    digest = hashlib.sha256()
    for name, content in sorted(files.items()):
        digest.update(name.encode() + b'\0' + content.encode() + b'\0')
    return digest.hexdigest()


def csv_paths(directory):
    return [p for p in Path(directory).glob('*/*') if p.is_file() and
            (p.name.endswith('.csv') or p.name.endswith('.csv.gz'))]


def csv_text(path):
    data=path.read_bytes()
    return (gzip.decompress(data) if path.suffix=='.gz' else data).decode()


def disk_files(directory):
    directory=Path(directory)
    return {p.relative_to(directory).as_posix():csv_text(p) for p in csv_paths(directory)}


def write_measurements(directory, files, *, compress_above=16*1024*1024):
    """Keep large per-function CSVs compact, with reproducible gzip headers."""
    directory=Path(directory);names=set()
    for name,content in files.items():
        raw=content.encode();compressed=len(raw)>compress_above
        actual=name+'.gz' if compressed else name
        path=directory/actual;path.parent.mkdir(parents=True,exist_ok=True)
        path.write_bytes(gzip.compress(raw,compresslevel=6,mtime=0) if compressed else raw)
        (directory/(name if compressed else name+'.gz')).unlink(missing_ok=True)
        names.add(actual)
    return names


def measurement_files(cases, records):
    groups = defaultdict(list)
    for case in cases: groups[case['api']].append(records[case['id']])
    files = {}
    for api, batch in sorted(groups.items()):
        module, function = api.rsplit('.', 1)
        count = max((len(r.get('warnings', [])) for c in batch for r in c['libraries'].values()), default=0)
        rows = []
        for case in batch:
            for ordinal,(library,result) in enumerate(case['libraries'].items()):
                rows.append([case['id'], library, ' '.join(case['reference']) if ordinal==0 else '',
                             ' '.join(map(repr, result.get('values', []))), result.get('error', ''),
                             *(result.get('warnings', []) + [''] * (count-len(result.get('warnings', []))))])
        files[f'{module}/{function}.csv'] = table(
            ['id', 'library', 'reference', 'values', 'error', *(f'warning_{i}' for i in range(count))], rows)
    return files


def read_measurements(directory):
    records = {}
    for content in disk_files(directory).values():
        batch = read_measurement_content(content)
        if set(records) & set(batch): raise ValueError('Duplicate measurement IDs')
        records.update(batch)
    return records


class MeasurementStore:
    """Read one function at a time so dense suites do not retain all outputs."""
    def __init__(self, directory):
        self.directory = Path(directory)
        self.api = None
        self.records = {}

    def __getitem__(self, identifier):
        api = identifier.split('/')[0]
        if api != self.api:
            module, function = api.rsplit('.',1)
            path=self.directory/module/f'{function}.csv'
            compressed=path.with_suffix('.csv.gz')
            if path.exists() and compressed.exists():raise ValueError(f'Duplicate measurement files: {api}')
            content=csv_text(path if path.exists() else compressed)
            self.records = read_measurement_content(content)
            self.api = api
        return self.records[identifier]


def read_measurement_content(content):
    records = {}
    for row in csv.DictReader(io.StringIO(content)):
        record = records.setdefault(row['id'], dict(id=row['id'], reference=row['reference'].split(), libraries={}))
        if row['reference'] and record['reference'] != row['reference'].split(): raise ValueError('Inconsistent reference')
        if row['library'] in record['libraries']: raise ValueError('Duplicate measurement')
        result = {'error': row['error']} if row['error'] else {'values': list(map(float, row['values'].split()))}
        if row['error'] and row['values']:result['values']=list(map(float,row['values'].split()))
        captured = [v for k,v in row.items() if k.startswith('warning_') and v]
        if captured: result['warnings'] = captured
        record['libraries'][row['library']] = result
    return records
