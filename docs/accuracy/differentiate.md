# differentiate: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-differentiate)

Measured binary64 (`f64`) results compared against the **saved per-case reference**.
Each implementation, including Sci-differentiate, is evaluated against the same
reference values stored in the raw CSVs. This page
covers 8277 cases and 4 public API names;
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

### `sci::differentiate`

4 functions · 8277 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | NumPy |
| --- | ---: | --- | --- | --- | --- | --- |
| `derivative` | 2201 | 7.5061e+10 / 2.5002e+10 | 0.48345 / 0.16759 | 14403.0 / 98.543 | 7.4001e+5 / 614.8 | 0.48345 / 0.16759 |
| `gradient` | 2076 | 2.0684e+8 / 1.0616e+6 | 0.5 / 0.10347 | 11880.0 / 18.355 | 49.318 / 0.18698 | 49.318 / 0.18698 |
| `hessian` | 2000 | 1.2274e+9 / 4.1672e+7 | 0.0 / 0.0 | 5918.0 / 399.64 | 0.0 / 0.0 | 0.0 / 0.0 |
| `jacobian` | 2000 | 10000.0 / 662.56 | 0.0 / 0.0 | 48.339 / 6.408 | 0.0 / 0.0 | 0.0 / 0.0 |


## Conventions and limitations

Callbacks are supplied explicitly: original sin and quadratic recipes, a quadratic Hessian and a two-output Jacobian. Fixed-step Sci stencils and adaptive SciPy differences have different work budgets. Analytic identities supply the callback references; raw inputs retain the selected steps and points. Boost.Math numerical finite differences and automatic differentiation evaluate the saved analytic callbacks independently. Their method and step policy differ from Sci fixed-step stencils, so agreement is numerical rather than a claim of identical algorithms.

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
