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
import shutil

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


def _case_batches(cases):
    from itertools import groupby
    ordered=cases if hasattr(cases,'function') else sorted(cases,key=lambda c:(c['api'],c['id']))
    for api, iterator in groupby(ordered,key=lambda c:c['api']):
        batch=list(iterator)
        module, function = api.rsplit('.', 1)
        arguments = sorted(batch[0]['args'])
        if any(sorted(c['args']) != arguments for c in batch):
            raise ValueError(f'Inconsistent arguments for {api}')
        yield f'{module}/{function}.csv', table(['id', 'group', *arguments],
            ([c['id'], c['group'], *(encode(c['args'][a]) for a in arguments)] for c in batch)), len(batch)


def iter_case_files(cases):
    for name, content, _ in _case_batches(cases):
        yield name, content


def case_files(cases):return dict(iter_case_files(cases))


def _write_if_changed(path, data):
    if not path.exists() or path.read_bytes() != data:
        path.write_bytes(data)


def write_cases(directory, cases, *, compress_above=256 * 1024, signature=None):
    """Store large input CSVs in deterministic gzip; optionally index a reusable snapshot.

    The canonical CSV content and numerical repr do not change. A cache index is
    written only after every input file has been stored successfully.
    """
    directory=Path(directory);names=set();index={}
    manifest=directory/'.cache.json'
    manifest.unlink(missing_ok=True)
    for name,content,count in _case_batches(cases):
        raw=content.encode('utf-8')
        compressed=len(raw)>compress_above
        actual=name+'.gz' if compressed else name
        path=directory/actual;path.parent.mkdir(parents=True,exist_ok=True)
        data=gzip.compress(raw,compresslevel=6,mtime=0) if compressed else raw
        _write_if_changed(path,data)
        (directory/(name if compressed else name+'.gz')).unlink(missing_ok=True)
        names.add(actual)
        index[name]={'count':count,'sha256':hashlib.sha256(data).hexdigest()}
    for path in csv_paths(directory):
        if path.relative_to(directory).as_posix() not in names:path.unlink()
    if signature is not None:
        manifest.write_text(json.dumps({'format':1,'recipe_sha256':signature,'files':index},
                                       sort_keys=True,indent=2)+'\n')


def case_paths(directory):
    """Map canonical input names to CSV or CSV.gz paths, refusing ambiguities."""
    directory=Path(directory);paths={}
    for path in csv_paths(directory):
        name=path.relative_to(directory).as_posix().removesuffix('.gz')
        if name in paths:raise ValueError(f'Duplicate case input files: {name}')
        paths[name]=path
    return paths


def case_directory_digest(directory):
    """Hash logical input CSVs exactly as the original uncompressed snapshots did.

    This preserves cases_sha256 in saved environment metadata when input files are
    compressed, while measurement digests continue hashing their physical bytes.
    """
    digest=hashlib.sha256()
    for name,path in sorted(case_paths(directory).items()):
        digest.update(name.encode()+b'\0')
        opener=gzip.open if path.suffix=='.gz' else open
        with opener(path,'rb') as stream:
            for chunk in iter(lambda:stream.read(1024*1024),b''):
                digest.update(chunk)
        digest.update(b'\0')
    return digest.hexdigest()


def compact_existing_cases(directory, *, signature, expected_digest, compress_above=256*1024):
    """Compress a verified saved input snapshot without generating numerical cases.

    The original content digest must match the independent snapshot metadata
    before any files are changed. CSV parsing counts the rows even when fields
    contain embedded newlines; source bytes and all numeric strings are retained.
    """
    directory=Path(directory)
    if case_directory_digest(directory)!=expected_digest:
        raise ValueError('Existing case inputs do not match the saved snapshot digest')
    manifest=directory/'.cache.json'
    manifest.unlink(missing_ok=True)
    entries={}
    for name,path in sorted(case_paths(directory).items()):
        raw=path.read_bytes()
        if path.suffix=='.gz':raw=gzip.decompress(raw)
        with io.StringIO(raw.decode('utf-8'),newline='') as stream:
            count=sum(1 for _ in csv.reader(stream))-1
        if count<1:raise ValueError(f'No cases in {name}')
        compressed=len(raw)>compress_above
        target=directory/(name+'.gz' if compressed else name)
        data=gzip.compress(raw,compresslevel=6,mtime=0) if compressed else raw
        if target!=path:
            _write_if_changed(target,data)
            path.unlink()
        else:
            _write_if_changed(target,data)
        entries[name]={'count':count,'sha256':hashlib.sha256(data).hexdigest()}
    manifest.write_text(json.dumps({'format':1,'recipe_sha256':signature,'files':entries},
                                   sort_keys=True,indent=2)+'\n')


def stage_case_extension(historical, additions, destination, *, compress_above=256*1024):
    """Copy historical input bytes verbatim, add files for genuinely new APIs.

    `historical` must be a verified CachedCaseSet. No numerical inputs from the
    old generator are evaluated; the original data files remain untouched.
    This stage has no cache manifest until new measurements are validated.
    """
    destination=Path(destination)
    if destination.exists():raise ValueError('Case extension stage already exists')
    destination.mkdir(parents=True)
    for name,path in sorted(historical.paths.items()):
        target=destination/path.relative_to(historical.directory)
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(path,target)
    expected={}
    for name,content,count in _case_batches(additions):
        if name in historical.paths:raise ValueError(f'Changed historical case function: {name}')
        raw=content.encode('utf-8')
        compressed=len(raw)>compress_above
        target=destination/(name+'.gz' if compressed else name)
        target.parent.mkdir(parents=True,exist_ok=True)
        _write_if_changed(target,gzip.compress(raw,compresslevel=6,mtime=0) if compressed else raw)
        expected[name]=count
    if not expected:raise ValueError('No new case files to extend')
    return expected


def finalize_extended_cache(destination, historical, additions, *, signature,
                            historical_cases_sha256, historical_measurements_directory,
                            historical_measurements_sha256, measurements_directory):
    """Index an extended input snapshot only after validating its measurement files.

    The manifest records historical lineage rather than asserting that new
    generator code reproduces the historical input bytes.
    """
    destination=Path(destination)
    historical_measurements_directory=Path(historical_measurements_directory)
    measurements_directory=Path(measurements_directory)
    if (destination/'.cache.json').exists():
        raise ValueError('Staged cache manifest already exists')
    if case_directory_digest(historical.directory)!=historical_cases_sha256:
        raise ValueError('Historical case digest changed during extension')
    if directory_digest(historical_measurements_directory)!=historical_measurements_sha256:
        raise ValueError('Historical measurement digest changed during extension')
    paths=case_paths(destination)
    added={}
    for name,content,count in _case_batches(additions):
        if name in added or name in historical.paths:raise ValueError(f'Duplicate extension API: {name}')
        if name not in paths or csv_text(paths[name])!=content:
            raise ValueError(f'Changed new input file: {name}')
        added[name]=count
    if not added or set(paths)!=set(historical.paths)|set(added):
        raise ValueError('Staged input file set differs from historical plus new inputs')
    for name,path in historical.paths.items():
        staged=paths[name]
        if staged.name!=path.name or _sha256_file(staged)!=_sha256_file(path):
            raise ValueError(f'Changed historical input bytes: {name}')
    old_measurements={p.relative_to(historical_measurements_directory).as_posix():p
                      for p in csv_paths(historical_measurements_directory)}
    new_measurements={p.relative_to(measurements_directory).as_posix():p
                      for p in csv_paths(measurements_directory)}
    if {name.removesuffix('.gz') for name in new_measurements}!=set(paths) or len(new_measurements)!=len(paths):
        raise ValueError('Extended measurements do not cover exactly the input files')
    for name,path in old_measurements.items():
        if name not in new_measurements or _sha256_file(path)!=_sha256_file(new_measurements[name]):
            raise ValueError(f'Changed historical measurements: {name}')
    for name in added:
        measured=new_measurements.get(name) or new_measurements.get(name+'.gz')
        records=read_measurement_content(csv_text(measured))
        original=read_case_file(paths[name],name)
        ids=[case['id'] for case in original]
        if len(ids)!=len(set(ids)) or set(ids)!=set(records) or len(ids)!=added[name]:
            raise ValueError(f'Missing, extra or duplicate measurement IDs: {name}')
        if any('Sci' not in record['libraries'] or not record['reference'] for record in records.values()):
            raise ValueError(f'Missing new Sci measurement or reference: {name}')
    index={name:{'count':historical.counts[name] if name in historical.counts else added[name],
                 'sha256':_sha256_file(path)} for name,path in sorted(paths.items())}
    provenance={'kind':'historical_snapshot_extension','historical_cases_sha256':historical_cases_sha256,
                'historical_measurements_sha256':historical_measurements_sha256,
                'historical_file_count':len(historical.paths),
                'new_file_count':len(added), 'generator_equivalence_verified':False}
    (destination/'.cache.json').write_text(json.dumps(
        {'format':2,'recipe_sha256':signature,'files':index,'provenance':provenance},
        sort_keys=True,indent=2)+'\n')
    return case_directory_digest(destination)


def _sha256_file(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):
            digest.update(chunk)
    return digest.hexdigest()


def directory_digest(directory):
    directory=Path(directory);digest=hashlib.sha256()
    for path in sorted(csv_paths(directory),key=lambda p:p.relative_to(directory).as_posix()):
        digest.update(path.relative_to(directory).as_posix().encode()+b'\0'+path.read_bytes()+b'\0')
    return digest.hexdigest()


def read_case_file(path, name):
    api=name.removesuffix('.csv').replace('/','.')
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',newline='',encoding='utf-8') as stream:
        return [dict(id=row.pop('id'), api=api, group=row.pop('group'),
                     args={k:decode(v) for k,v in row.items()}) for row in csv.DictReader(stream)]


def read_cases(directory):
    cases=[]
    for name,path in sorted(case_paths(directory).items(),key=lambda item:item[0].replace('/','.')):
        cases.extend(read_case_file(path,name))
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


def write_measurements(directory, files, *, compress_above=256*1024):
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
