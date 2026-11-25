# integrate: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-integrate)

Measured binary64 (`f64`) results compared against the **saved per-case reference**.
Each implementation, including Sci-integrate, is evaluated against the same
reference values stored in the raw CSVs. This page
covers 14135 cases and 7 public API names;
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

### `sci::integrate::quadrature`

2 functions · 4068 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | GSL | NumPy |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `fixed_quad` | 2068 | 21.656 / 2.9595 | 0.49047 / 0.18336 | 42.606 / 5.0552 | 4.2838 / 0.6596 | 4.2838 / 0.65351 | 8.0788 / 1.4399 |
| `newton_cotes` | 2000 | 1.1988e+7 / 8920.3 | —; unavailable 2000/2000 | 0.0 / 0.0 | —; unavailable 2000/2000 | —; unavailable 2000/2000 | —; unavailable 2000/2000 |

### `sci::integrate::sampled`

5 functions · 10067 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | GSL | NumPy |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `cumulative_simpson` | 2000 | 499.5 / 0.48561 | 582.75 / 0.43248 | 0.0 / 0.0 | —; unavailable 2000/2000 | —; unavailable 2000/2000 | —; unavailable 2000/2000 |
| `cumulative_trapezoid` | 2000 | 666.0 / 0.18511 | 711.07 / 0.30292 | 0.0 / 0.0 | —; unavailable 2000/2000 | —; unavailable 2000/2000 | 0.0 / 0.0 |
| `romb` | 2000 | 1.984 / 0.1853 | 1.0035 / 0.28093 | 0.0 / 0.0 | —; unavailable 2000/2000 | —; unavailable 2000/2000 | —; unavailable 2000/2000 |
| `simpson` | 2000 | 1.9521 / 0.361 | 1.4503 / 0.29902 | 0.0 / 0.0 | —; unavailable 2000/2000 | —; unavailable 2000/2000 | 1.4919 / 0.25101 |
| `trapezoid` | 2067 | 130.12 / 0.43205 | 0.47777 / 0.17771 | 142.8 / 0.57778 | —; unavailable 2067/2067 | —; unavailable 2067/2067 | 142.8 / 0.57778 |


## Conventions and limitations

Discrete-rule references integrate the supplied binary64 samples, excluding continuous sampling error. The dense added rules use dx = 0.25 and zero initial cumulative integrals. Newton-Cotes records weights followed by its error coefficient. Fixed Gaussian quadrature uses order 16 on exp(x), with analytic references. For the saved exponential callback only, Boost.Math and GSL independently provide fixed 16-point Gauss-Legendre quadrature. These are distinct from adaptive quadrature and the sampled-data trapezoid/Simpson operations.

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
