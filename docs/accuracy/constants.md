# constants: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-constants)

Measured binary64 (`f64`) results using **Sci-constants as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 14252 cases and 7 public API names;
coverage is limited to the operations and inputs below.

## Agreement with Sci

`Sci vs reference` is **maximum / mean difference between Sci and the saved
independent mathematical reference in ε**. Other library cells are **maximum /
mean difference from Sci in ε** over cases
where both evaluations succeed. `Cases` counts input cases; `Sci failed N/T`
identifies missing baseline evaluations. `failed N/T` counts that library's own
failures; `Sci unavailable N/T` counts its successful evaluations excluded
because Sci failed. A dash means no saved comparison or no paired finite results.
Means include all compared scalar components, including exact zeros.
Full outputs, references, warnings and errors remain in the per-function CSVs.
Tables are grouped by Sci 0.3 public submodule using the compiler-resolved
canonical API inventory. Historical flat recipe aliases appear under their
current public namespace, with their measured recipe names preserved.

### `sci::constants` (7 functions, 14252 cases)

| Function | Cases | Sci vs reference | SciPy | NumPy |
| --- | ---: | --- | --- | --- |
| `convert_temperature` | 2000 | 0.9076 / 0.031167 | 0.9076 / 0.031167 | — |
| `find` | 2000 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `lambda2nu` | 2126 | 0.48442 / 0.18108 | 0.0 / 0.0 | 0.0 / 0.0 |
| `nu2lambda` | 2126 | 0.48156 / 0.18018 | 0.0 / 0.0 | 0.0 / 0.0 |
| `precision` | 2000 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `unit` | 2000 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `value` | 2000 | 0.0 / 0.0 | 0.0 / 0.0 | — |


## Conventions and limitations

Vacuum conversions use metres and hertz and the exact SI value c = 299792458 m/s. Temperature recipes convert Celsius to Fahrenheit. CODATA values, uncertainties, unit strings and searches use the pinned SciPy constants table. Text metadata is encoded as UTF-8 lengths and bytes; Sci dtype f64 corresponds to NumPy float64.

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
