# integrate: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-integrate)

Measured binary64 (`f64`) results using **Sci-integrate as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 14135 cases and 7 public API names;
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

### `sci::integrate::quadrature`

2 functions · 4068 cases

| Function | Cases | Sci vs reference | mpmath | SciPy | Boost.Math | GSL | NumPy |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `fixed_quad` | 2068 | 21.656 / 2.9595 | 21.46 / 2.9638 | 41.586 / 5.142 | 18.905 / 3.0513 | 18.598 / 2.976 | 24.526 / 3.0109 |
| `newton_cotes` | 2000 | 1.1988e+7 / 8920.3 | — | 1.1988e+7 / 8920.3 | — | — | — |

### `sci::integrate::sampled`

5 functions · 10067 cases

| Function | Cases | Sci vs reference | mpmath | SciPy | Boost.Math | GSL | NumPy |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `cumulative_simpson` | 2000 | 499.5 / 0.48561 | 386.38 / 0.36181 | 499.5 / 0.48561 | — | — | — |
| `cumulative_trapezoid` | 2000 | 666.0 / 0.18511 | 553.27 / 0.16959 | 666.0 / 0.18511 | — | — | 666.0 / 0.18511 |
| `romb` | 2000 | 1.984 / 0.1853 | 1.0035 / 0.20654 | 1.984 / 0.1853 | — | — | — |
| `simpson` | 2000 | 1.9521 / 0.361 | 0.9994 / 0.15247 | 1.9521 / 0.361 | — | — | 1.5799 / 0.32696 |
| `trapezoid` | 2067 | 130.12 / 0.43205 | 129.74 / 0.33265 | 102.39 / 0.34064 | — | — | 102.39 / 0.34064 |


## Conventions and limitations

Discrete-rule references integrate the supplied binary64 samples, excluding continuous sampling error. The dense added rules use dx = 0.25 and zero initial cumulative integrals. Newton-Cotes records weights followed by its error coefficient. Fixed Gaussian quadrature uses order 16 on exp(x), with analytic references. For the saved exponential callback only, Boost.Math and GSL independently provide fixed 16-point Gauss-Legendre quadrature. These are distinct from adaptive quadrature and the sampled-data trapezoid/Simpson operations.

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
