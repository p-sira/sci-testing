# core: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-core)

Measured binary64 (`f64`) results using **Sci-core as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 3053 cases and 3 public API names;
coverage is limited to the operations and inputs below.

## Agreement with Sci

`Sci vs reference` is **maximum / mean difference between Sci and the saved
independent mathematical reference in ε**. Other library cells are **maximum /
mean difference from Sci in ε** over cases
where both evaluations succeed. `Cases` counts input cases; `Sci failed N/T`
identifies missing baseline evaluations. `failed N/T` counts that library's own
failures; `Sci unavailable N/T` counts its successful evaluations excluded
because Sci failed. `unavailable N/T` means the library has no compatible
public operation for those inputs or distribution families; they are excluded
from the reported mean and are not counted as numerical failures. A dash means
no saved comparison or no paired finite results.
Means include all compared scalar components, including exact zeros.
Full outputs, references, warnings and errors remain in the per-function CSVs.
Tables are grouped by Sci 0.3 public submodule using the compiler-resolved
canonical API inventory. Historical flat recipe aliases appear under their
current public namespace, with their measured recipe names preserved.

### `sci::core::csr` (3 functions, 3053 cases)

| Function | Cases | Sci vs reference | SciPy |
| --- | ---: | --- | --- |
| `csr.from_coo` | 1014 | 0.0 / 0.0 | 0.0 / 0.0 |
| `csr.to_dense` | 1020 | 0.0 / 0.0 | 0.0 / 0.0 |
| `csr.validate` | 1019 | 0.0 / 0.0 | 0.0 / 0.0 |


## Conventions and limitations

CSR constructors, validators and conversion routines are exercised against independent SciPy sparse and explicit storage invariants. Canonical pointer arrays, index ordering, zeros and rejection behavior are retained in the per-case measurements.

Inputs are generated deterministically into `sci-testing/cases/generated/<module>/<function>.csv` (optionally `.gz`); raw outputs, high-precision
references and the environment are stored in that repository's `results/`.
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
