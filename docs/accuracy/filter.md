# filter: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-filter)

Measured binary64 (`f64`) results compared against the **saved per-case reference**.
Each implementation, including Sci-filter, is evaluated against the same
reference values stored in the raw CSVs. This page
covers 20226 cases and 10 public API names;
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

### `sci::filter::kalman`

5 functions · 10156 cases

| Function | Cases | Sci | mpmath | SciPy | FilterPy |
| --- | ---: | --- | --- | --- | --- |
| `extended_kalman_predict` | 2003 | 0.70621 / 1.7611e-4 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 |
| `extended_kalman_update` | 2003 | 20.883 / 0.090124 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 |
| `kalman_predict` | 2067 | 62.462 / 0.34381 | 0.49622 / 0.1509 | —; unavailable 2067/2067 | 33.822 / 0.20742 |
| `kalman_update` | 2067 | 592.4 / 1.4033 | 0.4943 / 0.17981 | —; unavailable 2067/2067 | 592.4 / 0.9084 |
| `unscented_kalman_step` | 2016 | 160.0 / 0.18685 | —; unavailable 2016/2016 | 0.0 / 0.0 | —; unavailable 2016/2016 |

### `sci::filter::particle`

5 functions · 10070 cases

| Function | Cases | Sci | mpmath | SciPy | FilterPy |
| --- | ---: | --- | --- | --- | --- |
| `effective_sample_size` | 2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 |
| `particle_estimate` | 2006 | 0.61728 / 1.5371e-4 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 |
| `particle_predict` | 2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 |
| `particle_update` | 2028 | 1.8851 / 5.6209e-3 | —; unavailable 2028/2028 | 0.0 / 0.0 | —; unavailable 2028/2028 |
| `systematic_resample` | 2018 | 0.0 / 0.0 | —; unavailable 2018/2018 | 0.0 / 0.0 | —; unavailable 2018/2018 |


## Conventions and limitations

Gaussian outputs place means before row-major covariance; update and particle diagnostics follow the pinned filter fixture ordering. References use independent matrix identities, supplied deterministic particle draws and nonlinear callbacks. Original linear Kalman comparisons additionally call FilterPy.

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
