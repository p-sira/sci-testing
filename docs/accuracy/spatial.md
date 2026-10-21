# spatial: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-spatial)

Measured binary64 (`f64`) results using **Sci-spatial as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 66272 cases and 33 public API names;
coverage is limited to the operations and inputs below.

## Agreement with Sci

Each library cell is **maximum / mean difference from Sci in ε** over cases
where both evaluations succeed. `Cases` counts input cases; `Sci failed N/T`
identifies missing baseline evaluations. `failed N/T` counts that library's own
failures; `Sci unavailable N/T` counts its successful evaluations excluded
because Sci failed. A dash means unavailable or no paired finite results.
Means include all compared scalar components, including exact zeros.
Full outputs, references, warnings and errors remain in the per-function CSVs.

| Function | Cases | SciPy |
| --- | ---: | --- |
| `braycurtis` | 2003 | 0.0 / 0.0 |
| `canberra` | 2004 | 1.7148 / 0.19462 |
| `cdist` | 2019 | 4.5036e+15 / 3.7177e+11 |
| `chebyshev` | 2003 | 0.0 / 0.0 |
| `cityblock` | 2004 | 0.0 / 0.0 |
| `correlation` | 2003 | 3.9271 / 0.092261 |
| `cosine` | 2064 | 4.5036e+15 / 1.2108e+13 |
| `dice` | 2003 | 0.0 / 0.0 |
| `directed_hausdorff` | 2003 | 0.0 / 0.0 |
| `distance_matrix` | 2003 | 0.8165 / 1.1598e-4 |
| `euclidean` | 2064 | 1.779 / 0.21919 |
| `geometric_slerp` | 2004 | 4.5036e+15 / 3.7478e+11 |
| `hamming` | 2003 | 0.5 / 0.17374 |
| `is_valid_dm` | 2008 | 0.0 / 0.0 |
| `is_valid_y` | 2004 | 0.0 / 0.0 |
| `jaccard` | 2005 | 0.0 / 0.0 |
| `jensenshannon` | 2003 | 0.51954 / 2.5938e-4 |
| `mahalanobis` | 2002 | 1.0006 / 0.21304 |
| `minkowski` | 2004 | 0.89267 / 0.18142 |
| `minkowski_distance` | 2004 | 0.66674 / 0.18563 |
| `minkowski_distance_p` | 2004 | 0.88883 / 4.4353e-4 |
| `num_obs_dm` | 2004 | 0.0 / 0.0 |
| `num_obs_y` | 2003 | 0.0 / 0.0 |
| `pdist` | 2019 | 1.8783 / 0.13435 |
| `procrustes` | 2003 | 8.9643e+16 / 1.1526e+15 |
| `rogerstanimoto` | 2003 | 0.0 / 0.0 |
| `russellrao` | 2003 | 0.0 / 0.0 |
| `seuclidean` | 2002 | 1.2811 / 0.22172 |
| `sokalsneath` | 2003 | 0.0 / 0.0 |
| `sqeuclidean` | 2004 | 1.7782 / 0.45407 |
| `squareform` | 2008 | 0.0 / 0.0 |
| `tsearch` | 2002 | 0.0 / 0.0 |
| `yule` | 2004 | 0.0 / 0.0 |


## Conventions and limitations

Metrics, pairwise arrays, geometry and validation helpers use explicit weights, variance/covariance data and shapes from the pinned spatial fixtures. Output ordering follows those fixtures. Near-zero distances and cancellation can amplify relative error; the CSV retains full outputs for absolute-error inspection.

Inputs are generated deterministically into `sci-testing/cases/generated/<module>/<function>.csv`; raw outputs, high-precision
references and the environment are stored in that repository's `results/`.
NaN, infinity, signed-zero contracts, invalid domains and work-limit enforcement
are outside this finite-input accuracy sample. Existing Sci API tests cover
additional behavior. The [coverage audit](../../results/coverage.csv) checks every public function
against the compiler-resolved inventory. These samples do not establish
worst-case bounds or full domain coverage. SciPy aliases, compatibility namespaces and Cython
variants are not counted as independent implementations.

## Reproduction

- Sci source: `b52c598535920279dc0aadb14a9085b8545ba4c1` (numerical sources clean: true).
- Testing source SHA-256: `a1dbbae1205dbebed9bfb371b138d8f9a215cfec01568c6e3de356becdf58cd9`.
- NumPy 2.5.4; SciPy 1.18.1; FilterPy 1.4.5; mpmath 1.3.0.
- Native comparison versions: {'Boost.Math': '1_92', 'GSL': '2.8'}.
- rustc 1.98.1 (48a229cea 2026-09-01); Python 3.13.3; Cargo dev profile, opt-level 0.
- Linux-7.2.11-arch1-1-x86_64-with-glibc2.44; AMD Ryzen 5 4600H with Radeon Graphics.

See the [shared metric, provenance and reproduction instructions](README.md).
