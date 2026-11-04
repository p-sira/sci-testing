"""Authenticate sourced primary reference values independently of report code."""

import hashlib
import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import NotRequired, Protocol, TypedDict, cast


class ReferenceCase(TypedDict):
    """Fields needed from a stored case; arguments retain their native values."""

    api: str
    id: str
    args: dict[str, object]
    group: NotRequired[str]


class ReferenceCases(Protocol):
    @property
    def apis(self) -> Sequence[str]: ...

    def function(self, api: str) -> list[ReferenceCase]: ...


class ReferenceMeasurement(TypedDict):
    id: str
    reference: list[str]
    libraries: Mapping[str, object]


class ReferenceRecords(Protocol):
    def __getitem__(self, key: str, /) -> ReferenceMeasurement: ...


class ProvenanceEntry(TypedDict):
    provider: str
    recipe: str
    precision: str
    status: str
    cases: int
    expression_template: NotRequired[str]
    documentation: NotRequired[str]


class PriorProvenance(TypedDict):
    schema: int
    cases_sha256: str
    measurements_sha256: str
    reference_sha256: dict[str, str]
    apis: dict[str, list[ProvenanceEntry]]

MANIFEST = 'reference-provenance.json'
SCHEMA = 2


def reference_hashes(cases: ReferenceCases, records: ReferenceRecords,
                     *, apis: Sequence[str] | None = None) -> dict[str, str]:
    """Hash reference decimals and IDs, excluding comparator output bytes.

    A selected API subset suffices when all other measurement files remain
    unchanged in a hardlinked staging snapshot.
    """
    hashes: dict[str, str] = {}
    for api in cases.apis if apis is None else apis:
        digest = hashlib.sha256()
        for case in cases.function(api):
            values = records[case['id']]['reference']
            digest.update(case['id'].encode() + b'\0')
            digest.update(json.dumps(values, separators=(',', ':')).encode() + b'\0')
        hashes[api] = digest.hexdigest()
    return hashes


def serialize_provenance(entries: Mapping[str, Sequence[ProvenanceEntry]], *,
                         inputs_sha256: str, measurements_sha256: str,
                         hashes: Mapping[str, str] | None = None) -> bytes:
    return (json.dumps({'schema': SCHEMA, 'cases_sha256': inputs_sha256,
                        'measurements_sha256': measurements_sha256,
                        'reference_sha256': hashes or {}, 'apis': entries},
                       sort_keys=True, indent=2) + '\n').encode()


def load_prior_provenance(meta: Mapping[str, object], path: Path) -> PriorProvenance | None:
    """Accept historical reference provenance only with authenticated bytes."""
    checksum = meta.get('reference_provenance_sha256')
    if not checksum:
        return None
    if not path.is_file():
        raise ValueError('Saved reference provenance is missing')
    content = path.read_bytes()
    if hashlib.sha256(content).hexdigest() != checksum:
        raise ValueError('Saved reference provenance digest mismatch')
    decoded = cast(object, json.loads(content))
    if not isinstance(decoded, dict):
        raise TypeError('Saved reference provenance is not a JSON object')
    saved = cast(dict[str, object], decoded)
    if saved.get('schema') != SCHEMA or saved['cases_sha256'] != meta['cases_sha256']:
        raise ValueError('Saved reference provenance case identity mismatch')
    entries = saved.get('apis')
    hashes = saved.get('reference_sha256')
    if not isinstance(entries, dict) or not isinstance(hashes, dict):
        raise TypeError('Saved reference provenance has invalid entries or hashes')
    raw_entries = cast(dict[object, object], entries)
    raw_hashes = cast(dict[object, object], hashes)
    if (any(not isinstance(api, str) or not isinstance(rows, list)
            for api, rows in raw_entries.items())
            or any(not isinstance(api, str) or not isinstance(digest, str)
                   for api, digest in raw_hashes.items())):
        raise ValueError('Saved reference provenance has invalid entries or hashes')
    if set(raw_hashes) != set(raw_entries):
        raise ValueError('Saved reference provenance hashes are incomplete')
    for rows in raw_entries.values():
        if not isinstance(rows, list):
            raise TypeError('Saved reference provenance entries must be lists')
        for raw in cast(list[object], rows):
            if not isinstance(raw, dict):
                raise TypeError('Saved reference provenance has invalid entries')
            item = cast(dict[str, object], raw)
            if (any(not isinstance(item.get(key), str)
                    for key in ('provider', 'recipe', 'precision', 'status'))
                    or type(item.get('cases')) is not int):
                raise ValueError('Saved reference provenance has invalid entries')
    return cast(PriorProvenance, cast(object, saved))


def verify_provenance(cases: ReferenceCases, records: ReferenceRecords,
                      meta: Mapping[str, object], root: Path) -> None:
    """Validate provenance against every stored reference without re-evaluation."""
    prior = load_prior_provenance(meta, root/'results'/MANIFEST)
    if prior is not None:
        if prior['measurements_sha256'] != meta.get('measurements_sha256'):
            raise ValueError('Saved reference provenance measurement digest mismatch')
        if reference_hashes(cases, records) != prior['reference_sha256']:
            raise ValueError('Saved reference provenance values differ from raw measurements')
