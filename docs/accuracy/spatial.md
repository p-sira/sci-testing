# spatial: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-spatial)

Measured binary64 (`f64`) results using **Sci-spatial as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 66272 cases and 33 public API names;
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

### `sci::spatial::distance` (30 functions, 60263 cases)

| Function | Cases | Sci vs reference | SciPy |
| --- | ---: | --- | --- |
| `braycurtis` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 |
| `canberra` | 2004 | 1.7148 / 0.19462 | 1.7148 / 0.19462 |
| `cdist` | 2019 | 4.5036e+15 / 3.7177e+11 | 4.5036e+15 / 3.7177e+11 |
| `chebyshev` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 |
| `cityblock` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 |
| `correlation` | 2003 | 3.9271 / 0.092261 | 3.9271 / 0.092261 |
| `cosine` | 2064 | 2.0455e+18 / 1.0604e+15 | 4.5036e+15 / 1.2108e+13 |
| `dice` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 |
| `directed_hausdorff` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 |
| `distance_matrix` | 2003 | 0.8165 / 1.1598e-4 | 0.8165 / 1.1598e-4 |
| `euclidean` | 2064 | 1.5729 / 0.2364 | 1.779 / 0.21919 |
| `hamming` | 2003 | 0.5 / 0.17374 | 0.5 / 0.17374 |
| `is_valid_dm` | 2008 | 0.0 / 0.0 | 0.0 / 0.0 |
| `is_valid_y` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 |
| `jaccard` | 2005 | 0.0 / 0.0 | 0.0 / 0.0 |
| `jensenshannon` | 2003 | 0.51954 / 2.5938e-4 | 0.51954 / 2.5938e-4 |
| `mahalanobis` | 2002 | 1.0006 / 0.21304 | 1.0006 / 0.21304 |
| `minkowski` | 2004 | 0.89267 / 0.18142 | 0.89267 / 0.18142 |
| `minkowski_distance` | 2004 | 0.66674 / 0.18563 | 0.66674 / 0.18563 |
| `minkowski_distance_p` | 2004 | 0.88883 / 4.4353e-4 | 0.88883 / 4.4353e-4 |
| `num_obs_dm` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 |
| `num_obs_y` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 |
| `pdist` | 2019 | 1.8783 / 0.13435 | 1.8783 / 0.13435 |
| `rogerstanimoto` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 |
| `russellrao` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 |
| `seuclidean` | 2002 | 1.2811 / 0.22172 | 1.2811 / 0.22172 |
| `sokalsneath` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 |
| `sqeuclidean` | 2004 | 1.7782 / 0.45407 | 1.7782 / 0.45407 |
| `squareform` | 2008 | 0.0 / 0.0 | 0.0 / 0.0 |
| `yule` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 |

### `sci::spatial::transform` (3 functions, 6009 cases)

| Function | Cases | Sci vs reference | SciPy |
| --- | ---: | --- | --- |
| `geometric_slerp` | 2004 | 4.5036e+15 / 3.7478e+11 | 4.5036e+15 / 3.7478e+11 |
| `procrustes` | 2003 | 8.9643e+16 / 1.1526e+15 | 8.9643e+16 / 1.1526e+15 |
| `tsearch` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 |


## Conventions and limitations

Metrics, pairwise arrays, geometry and validation helpers use explicit weights, variance/covariance data and shapes from the pinned spatial fixtures. Output ordering follows those fixtures. Near-zero distances and cancellation can amplify relative error; the CSV retains full outputs for absolute-error inspection.

Inputs are generated deterministically into `sci-testing/cases/generated/<module>/<function>.csv` (optionally `.gz`); raw outputs, high-precision
references and the environment are stored in that repository's `results/`.
NaN, infinity, signed-zero contracts, invalid domains and work-limit enforcement
are outside this finite-input accuracy sample. Existing Sci API tests cover
additional behavior. The [coverage audit](../../results/coverage.csv) checks every public function
against the compiler-resolved inventory. These samples do not establish
worst-case bounds or full domain coverage. SciPy aliases, compatibility namespaces and Cython
variants are not counted as independent implementations.

## Reproduction

- Sci source: `c781e8259b175ebeb3c7d27a522ef3f3a7e24998` (numerical sources clean: true).
- Testing source SHA-256: `3d1f538d12ef58c6738188ac9a526fb0459c4c5cb1023a3c7edaa3c23d5a3452`.
- NumPy 2.5.4; SciPy 1.18.1; FilterPy 1.4.5; mpmath 1.3.0.
- Native comparison versions: {'Boost.Math': '1_92', 'GSL': '2.8'}.
- rustc 1.99.1 (58cf8f9c3 2026-10-13); Python 3.13.3; dev; opt-level=0; debug=line-tables-only; no custom RUSTFLAGS.
- Linux-7.2.11-arch1-1-x86_64-with-glibc2.44; AMD Ryzen 5 4600H with Radeon Graphics.

See the [shared metric, provenance and reproduction instructions](README.md).

Sci-only refresh: independent references and other library outputs are retained from the saved snapshot. Their original environment is recorded in `results/environment.json` under `reference_environment`.

Sci results were refreshed by namespace. Per-function Sci revisions are recorded in `results/environment.json` under `sci_measurement_commits`.
