# core: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-core)

Measured binary64 (`f64`) results compared against the **saved per-case reference**.
Each implementation, including Sci-core, is evaluated against the same
reference values stored in the raw CSVs. This page
covers 3053 cases and 3 public API names;
coverage is limited to the operations and inputs below.

## Differences from saved references

Every implementation column reports **maximum / mean symmetric relative
difference from the saved reference in ε**, over successful evaluations.
Reference source and precision vary by function. In particular, binary64
reference values measure agreement with their generating method; these
tables do not presume that any provider represents mathematical ground truth.
`Cases` counts input cases; `Sci failed N/T` identifies Sci failures and
`Sci unavailable N/T` identifies absent Sci measurements, while
`failed N/T` counts other libraries' failures. A Sci failure does not exclude
the same input from other libraries' measurements.
`unavailable N/T` means the library has no compatible
public operation for those inputs or distribution families; they are excluded
from the reported mean and are not counted as numerical failures. A dash means
no saved evaluation or no successful scalar comparisons.
Means include all compared scalar components, including exact zeros.
Full outputs, references, warnings and errors remain in the per-function CSVs.
Tables are grouped by Sci 0.3 public submodule using the compiler-resolved
canonical API inventory. Historical flat recipe aliases appear under their
current public namespace, with their measured recipe names preserved.

### `sci::core::csr`

3 functions · 3053 cases

| Function | Cases | Sci | SciPy |
| --- | ---: | --- | --- |
| `csr.from_coo` | 1014 | 0.0 / 0.0 | 0.0 / 0.0 |
| `csr.to_dense` | 1020 | 0.0 / 0.0 | 0.0 / 0.0 |
| `csr.validate` | 1019 | 0.0 / 0.0 | 0.0 / 0.0 |


## Conventions and limitations

CSR constructors, validators and conversion routines are exercised against independent SciPy sparse and explicit storage invariants. Canonical pointer arrays, index ordering, zeros and rejection behavior are retained in the per-case measurements.

Inputs are generated deterministically into `sci-testing/cases/generated/<module>/<function>.csv` (optionally `.gz`); raw outputs, reference
values and the environment are stored in that repository's `results/`.
Reference provenance is recorded in `results/reference-provenance.json` when
generated with the curated reference provider; historical sources are described
by `results/environment.json` and the documented recipe definitions.
NaN, infinity, signed-zero contracts, invalid domains and work-limit enforcement
are outside this finite-input accuracy sample. Existing Sci API tests cover
additional behavior. The [coverage audit](../../results/coverage.csv) checks every public function
against the compiler-resolved inventory. These samples do not establish
worst-case bounds or full domain coverage. SciPy aliases, compatibility namespaces and Cython
variants are not counted as independent implementations.

## Reproduction

The [shared methodology and environment](README.md) and
`results/environment.json` record source revisions, hashes, toolchain versions,
compiler flags and per-function measurement provenance.
