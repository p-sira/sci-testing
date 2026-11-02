# io: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-io)

Measured binary64 (`f64`) results using **Sci-io as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 413 cases and 10 public API names;
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

### `sci::io::matlab` (4 functions, 342 cases)

| Function | Cases | Sci vs reference | SciPy |
| --- | ---: | --- | --- |
| `matlab.loadmat` | 69 | 0.0 / 0.0 | 0.0 / 0.0 |
| `matlab.matfile_version` | 68 | 0.0 / 0.0 | 0.0 / 0.0 |
| `matlab.savemat` | 136 | 0.0 / 0.0 | 0.0 / 0.0 |
| `matlab.whosmat` | 69 | 0.0 / 0.0 | 0.0 / 0.0 |

### `sci::io::matrix_market` (3 functions, 33 cases)

| Function | Cases | Sci vs reference | SciPy |
| --- | ---: | --- | --- |
| `matrix_market.mminfo` | 12 | 0.0 / 0.0 | 0.0 / 0.0 |
| `matrix_market.mmread` | 11 | 0.0 / 0.0 | 0.0 / 0.0 |
| `matrix_market.mmwrite` | 10 | 0.0 / 0.0 | 0.0 / 0.0 |

### `sci::io::wavfile` (3 functions, 38 cases)

| Function | Cases | Sci vs reference | SciPy |
| --- | ---: | --- | --- |
| `wavfile.read` | 9 | 0.0 / 0.0 | 0.0 / 0.0 |
| `wavfile.signal` | 21 | 0.0 / 0.0 | 0.0 / 0.0 |
| `wavfile.write` | 8 | 0.0 / 0.0 | 0.0 / 0.0 |


## Conventions and limitations

Matrix Market, WAV and MATLAB comparisons use deterministic format fixtures and independent library references. Structured outputs follow the documented fixture ordering; raw input cases retain the shapes and format parameters.

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
