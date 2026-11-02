# optimize: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-optimize)

Measured binary64 (`f64`) results using **Sci-optimize as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 20152 cases and 10 public API names;
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

### `sci::optimize::derivatives` (6 functions, 12152 cases)

| Function | Cases | Sci vs reference | SciPy |
| --- | ---: | --- | --- |
| `approx_fprime` | 2012 | 32000.0 / 2851.2 | 32000.0 / 2851.2 |
| `check_grad` | 2008 | 6.0541e+6 / 1.3246e+6 | 6.0541e+6 / 1.3246e+6 |
| `rosen` | 2060 | 3.7821e+7 / 47584.0 | 2.5376 / 0.23879 |
| `rosen_der` | 2060 | 4.5036e+15 / 8.5513e+11 | 45.844 / 0.2412 |
| `rosen_hess` | 2006 | 1.28 / 0.098842 | 1.28 / 0.098842 |
| `rosen_hess_prod` | 2006 | 1.2806 / 0.097823 | 1.2806 / 0.097823 |

### `sci::optimize::minimize` (4 functions, 8000 cases)

| Function | Cases | Sci vs reference | SciPy |
| --- | ---: | --- | --- |
| `minimize` | 2000 | 6.7039e+7 / 1.8008e+7 | 91936.0 / 95.782 |
| `minimize_with_gradient` | 2000 | 1.6777e+7 / 4.403e+6 | 1.6777e+7 / 4.403e+6 |
| `try_minimize` | 2000 | 6.7098e+7 / 1.7947e+7 | 95031.0 / 108.9 |
| `try_minimize_with_gradient` | 2000 | 1.6777e+7 / 4.4e+6 | 1.6777e+7 / 4.4e+6 |


## Conventions and limitations

Objective helpers use their pinned Rosenbrock and finite-difference conventions. Minimizers use 3+(x-1)^2+(y+2)^2, a 512-iteration budget and tolerance 1e-8; only converged points and costs are compared. Iteration trajectories are not compared. Analytic minimizers supply the reference; SciPy runs Nelder-Mead or BFGS separately.

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
