# differentiate: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-differentiate)

Measured binary64 (`f64`) results using **Sci-differentiate as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 8277 cases and 4 public API names;
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

### `sci::differentiate` (4 functions, 8277 cases)

| Function | Cases | Sci vs reference | SciPy |
| --- | ---: | --- | --- |
| `derivative` | 2201 | 7.5061e+10 / 2.5002e+10 | 7.5061e+10 / 2.5002e+10 |
| `gradient` | 2076 | 2.0684e+8 / 1.0616e+6 | 2.0684e+8 / 1.0616e+6 |
| `hessian` | 2000 | 1.2274e+9 / 4.1672e+7 | 1.2274e+9 / 4.1672e+7 |
| `jacobian` | 2000 | 10000.0 / 662.56 | 10018.0 / 665.08 |


## Conventions and limitations

Callbacks are supplied explicitly: original sin and quadratic recipes, a quadratic Hessian and a two-output Jacobian. Fixed-step Sci stencils and adaptive SciPy differences have different work budgets. Analytic identities supply the callback references; raw inputs retain the selected steps and points.

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
