# filter: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-filter)

Measured binary64 (`f64`) results using **Sci-filter as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 20226 cases and 10 public API names;
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

### `sci::filter::kalman` (5 functions, 10156 cases)

| Function | Cases | Sci vs reference | SciPy | FilterPy |
| --- | ---: | --- | --- | --- |
| `extended_kalman_predict` | 2003 | 0.70621 / 1.7611e-4 | 0.70621 / 1.7611e-4 | — |
| `extended_kalman_update` | 2003 | 20.883 / 0.090124 | 20.883 / 0.090124 | — |
| `kalman_predict` | 2067 | 62.462 / 0.34381 | — | 62.501 / 0.18821 |
| `kalman_update` | 2067 | 592.4 / 1.4033 | — | 536.75 / 1.1847 |
| `unscented_kalman_step` | 2016 | 160.0 / 0.18685 | 160.0 / 0.18685 | — |

### `sci::filter::particle` (5 functions, 10070 cases)

| Function | Cases | Sci vs reference | SciPy | FilterPy |
| --- | ---: | --- | --- | --- |
| `effective_sample_size` | 2006 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `particle_estimate` | 2006 | 0.61728 / 1.5371e-4 | 0.61728 / 1.5371e-4 | — |
| `particle_predict` | 2012 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `particle_update` | 2028 | 1.8851 / 5.6209e-3 | 1.8851 / 5.6209e-3 | — |
| `systematic_resample` | 2018 | 0.0 / 0.0 | 0.0 / 0.0 | — |


## Conventions and limitations

Gaussian outputs place means before row-major covariance; update and particle diagnostics follow the pinned filter fixture ordering. References use independent matrix identities, supplied deterministic particle draws and nonlinear callbacks. Original linear Kalman comparisons additionally call FilterPy.

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
