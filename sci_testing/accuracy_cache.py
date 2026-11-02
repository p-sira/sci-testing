"""Content-addressed per-API accuracy table summaries.

The saved summary is derived data, never an oracle.  Each entry binds the
physical input and measurement file, the current metric implementation, and
the set of competitor columns.  The full snapshot digests and per-case ID
checks remain the responsibility of report.load_snapshot.
"""
import ast
import hashlib
import json
from pathlib import Path

from .case_csv import case_paths, csv_paths


SCHEMA = 1
FILENAME = 'accuracy-summaries.json'

# Per-API numerical rows depend on these definitions, not page presentation
# prose or the code that invokes the comparison libraries. Row aggregation and
# row formatting are both included; Markdown page formatting is not.
METRIC_FUNCTIONS = frozenset({'symmetric_error', 'component_error', 'number',
                              'summarize_api'})
METRIC_CONSTANTS = frozenset({'EPS', 'TINY', 'ZERO', 'METRIC_PRECISION_DPS'})


def metric_source_sha256(path):
    """Fingerprint executable accuracy calculations, excluding report prose.

    Extract full Python ASTs rather than looking for an edit-sensitive range of
    source lines. Missing/renamed definitions fail closed instead of silently
    preserving stale summaries. The separate testing_sources_sha256 continues
    to authenticate *all* source bytes, including report notes and adapters.
    """
    tree = ast.parse(Path(path).read_text(encoding='utf-8'))
    functions = {}
    constants = {}
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name in METRIC_FUNCTIONS:
            functions[node.name] = ast.dump(node, include_attributes=False)
        elif isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            for target in targets:
                if isinstance(target, ast.Name) and target.id in METRIC_CONSTANTS:
                    constants[target.id] = ast.dump(node.value, include_attributes=False)
    if set(functions) != METRIC_FUNCTIONS or set(constants) != METRIC_CONSTANTS:
        raise ValueError('Missing numerical metric definition; cannot trust cached summaries')
    payload = json.dumps({'functions': functions, 'constants': constants},
                         sort_keys=True, separators=(',', ':'))
    return hashlib.sha256(payload.encode()).hexdigest()


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
        self.metric_sha256 = hashlib.sha256((metric_source_sha256(metric_source) +
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
