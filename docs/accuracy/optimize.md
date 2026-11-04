# optimize: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-optimize)

Measured binary64 (`f64`) results compared against the **saved per-case reference**.
Each implementation, including Sci-optimize, is evaluated against the same
reference values stored in the raw CSVs. This page
covers 20152 cases and 10 public API names;
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

### `sci::optimize::derivatives`

6 functions · 12152 cases

| Function | Cases | Sci | mpmath | SciPy | GSL |
| --- | ---: | --- | --- | --- | --- |
| `approx_fprime` | 2012 | 32000.0 / 2851.2 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 |
| `check_grad` | 2008 | 6.0541e+6 / 1.3246e+6 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 |
| `rosen` | 2060 | 3.7821e+7 / 47584.0 | 0.49473 / 0.18029 | 3.7821e+7 / 47584.0 | —; unavailable 2060/2060 |
| `rosen_der` | 2060 | 4.5036e+15 / 8.5513e+11 | 0.49023 / 0.18134 | 4.5036e+15 / 8.5513e+11 | —; unavailable 2060/2060 |
| `rosen_hess` | 2006 | 1.28 / 0.098842 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 |
| `rosen_hess_prod` | 2006 | 1.2806 / 0.097823 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 |

### `sci::optimize::minimize`

4 functions · 8000 cases

| Function | Cases | Sci | mpmath | SciPy | GSL |
| --- | ---: | --- | --- | --- | --- |
| `minimize` | 2000 | 6.7039e+7 / 1.8008e+7 | 0.0 / 0.0 | 6.7039e+7 / 1.8008e+7 | 2.2317e+8 / 1.1964e+7 |
| `minimize_with_gradient` | 2000 | 1.6777e+7 / 4.403e+6 | 0.0 / 0.0 | 18.0 / 0.91633 | 1701.0 / 91.508 |
| `try_minimize` | 2000 | 6.7098e+7 / 1.7947e+7 | 0.0 / 0.0 | 6.7098e+7 / 1.7947e+7 | 2.7083e+8 / 1.181e+7 |
| `try_minimize_with_gradient` | 2000 | 1.6777e+7 / 4.4e+6 | 0.0 / 0.0 | 18.0 / 0.94533 | 1702.0 / 95.613 |


## Conventions and limitations

Objective helpers use their pinned Rosenbrock and finite-difference conventions. Minimizers use 3+(x-1)^2+(y+2)^2, a 512-iteration budget and tolerance 1e-8; only converged points and costs are compared. Iteration trajectories are not compared. Analytic minimizers supply the reference; SciPy runs Nelder-Mead or BFGS separately. GSL native multidimensional simplex and BFGS optimizers provide independent solutions of the saved convex two-dimensional quadratic objective, with guarded convergence, iteration budgets and function offset handling. Convergence trajectories and algorithms need not match Sci optimizers; only the saved final position and cost are compared.

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
