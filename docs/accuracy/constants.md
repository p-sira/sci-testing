# constants: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-constants)

Measured binary64 (`f64`) results compared against the **saved per-case reference**.
Each implementation, including Sci-constants, is evaluated against the same
reference values stored in the raw CSVs. This page
covers 14252 cases and 7 public API names;
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

### `sci::constants`

7 functions · 14252 cases

| Function | Cases | Sci | mpmath | SciPy | NumPy | Boost.Units |
| --- | ---: | --- | --- | --- | --- | --- |
| `convert_temperature` | 2000 | 0.9076 / 0.031167 | 13.191 / 3.8327 | 0.0 / 0.0 | —; unavailable 2000/2000 | 13.309 / 3.8302 |
| `find` | 2000 | 0.0 / 0.0 | —; unavailable 2000/2000 | 0.0 / 0.0 | —; unavailable 2000/2000 | —; unavailable 2000/2000 |
| `lambda2nu` | 2126 | 0.48442 / 0.18108 | 0.48442 / 0.18108 | 0.48442 / 0.18108 | 0.48442 / 0.18108 | 0.48442 / 0.18108 |
| `nu2lambda` | 2126 | 0.48156 / 0.18018 | 0.48156 / 0.18018 | 0.48156 / 0.18018 | 0.48156 / 0.18018 | 0.48156 / 0.18018 |
| `precision` | 2000 | 0.0 / 0.0 | —; unavailable 2000/2000 | 0.0 / 0.0 | —; unavailable 2000/2000 | 0.0 / 0.0; unavailable 1998/2000 |
| `unit` | 2000 | 0.0 / 0.0 | —; unavailable 2000/2000 | 0.0 / 0.0 | —; unavailable 2000/2000 | —; unavailable 2000/2000 |
| `value` | 2000 | 0.0 / 0.0 | —; unavailable 2000/2000 | 0.0 / 0.0 | —; unavailable 2000/2000 | 0.0 / 0.0; unavailable 1996/2000 |


## Conventions and limitations

Vacuum conversions use metres and hertz and the exact SI value c = 299792458 m/s. Temperature recipes convert Celsius to Fahrenheit. CODATA values, uncertainties, unit strings and searches use the pinned SciPy constants table. Text metadata is encoded as UTF-8 lengths and bytes; Sci dtype f64 corresponds to NumPy float64. Boost.Units adds independently typed affine temperature and SI wavelength/frequency conversions; it is not part of Boost.Math. The currently sampled CODATA lookup/metadata tables are not assumed to match Boost.Units constants from other release years.

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
