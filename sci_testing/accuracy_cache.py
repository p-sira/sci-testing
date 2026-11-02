"""Content-addressed per-API accuracy table summaries.

The saved summary is derived data, never an oracle.  Each entry binds the
physical input and measurement file, the current metric implementation, and
the set of competitor columns.  The full snapshot digests and per-case ID
checks remain the responsibility of report.load_snapshot.
"""
import hashlib
import json
from pathlib import Path

from .case_csv import case_paths, csv_paths


SCHEMA = 1
FILENAME = 'accuracy-summaries.json'


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


class AccuracyCache:
    def __init__(self, root, measurements, metric_source, *, inputs=None, previous=None,
                 previous_sha256=None, reuse=True, require_hits=False):
        self.root = Path(root)
        self.measurements = Path(measurements)
        self.metric_sha256 = hashlib.sha256((sha256_file(metric_source) +
                                            sha256_file(__file__)).encode()).hexdigest()
        self.reuse = reuse
        self.require_hits = require_hits
        self.entries = {}
        self.input_paths = case_paths(inputs if inputs is not None else self.root / 'cases/generated')
        self.measurement_paths = {
            p.relative_to(self.measurements).as_posix().removesuffix('.gz'): p
            for p in csv_paths(self.measurements)
        }
        self._hashes = {}
        if previous is not None and previous_sha256 is not None:
            content = Path(previous).read_bytes()
            if hashlib.sha256(content).hexdigest() != previous_sha256:
                raise ValueError('Accuracy summary digest mismatch')
            saved = json.loads(content)
            if saved.get('schema') == SCHEMA and saved.get('metric_sha256') == self.metric_sha256:
                self.entries = saved['entries'].copy()

    def _key(self, api, libraries):
        module, function = api.rsplit('.', 1)
        name = f'{module}/{function}.csv'
        inputs = self.input_paths[name]
        measurements = self.measurement_paths[name]
        def file_hash(path):
            if path not in self._hashes:
                self._hashes[path] = sha256_file(path)
            return self._hashes[path]
        return {'input': file_hash(inputs), 'measurement': file_hash(measurements),
                'libraries': list(libraries)}

    def get(self, api, libraries):
        key = self._key(api, libraries)
        entry = self.entries.get(api)
        if self.reuse and entry is not None and entry.get('key') == key:
            return entry['summary']
        if self.require_hits:
            raise ValueError(f'Stale or missing accuracy summary: {api}; rerun --render')
        return None

    def put(self, api, libraries, summary):
        self.entries[api] = {'key': self._key(api, libraries), 'summary': summary}

    def serialize(self, expected_apis):
        if set(self.entries) != set(expected_apis):
            raise ValueError('Accuracy summary API set differs from the saved inputs; rerun --render')
        return (json.dumps({'schema': SCHEMA, 'metric_sha256': self.metric_sha256,
                            'entries': self.entries}, sort_keys=True,
                           separators=(',', ':'), ensure_ascii=False) + '\n').encode()
