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
because Sci failed. `unavailable N/T` means the library has no compatible
public operation for those inputs or distribution families; they are excluded
from the reported mean and are not counted as numerical failures. A dash means
no saved comparison or no paired finite results.
Means include all compared scalar components, including exact zeros.
Full outputs, references, warnings and errors remain in the per-function CSVs.
Tables are grouped by Sci 0.3 public submodule using the compiler-resolved
canonical API inventory. Historical flat recipe aliases appear under their
current public namespace, with their measured recipe names preserved.

### `sci::spatial::distance` (30 functions, 60263 cases)

| Function | Cases | Sci vs reference | SciPy | Boost.Math | GSL |
| --- | ---: | --- | --- | --- | --- |
| `braycurtis` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — |
| `canberra` | 2004 | 1.7148 / 0.19462 | 1.7148 / 0.19462 | — | — |
| `cdist` | 2019 | 4.5036e+15 / 3.7177e+11 | 4.5036e+15 / 3.7177e+11 | 1.3333 / 0.18284; failed 13/2019 | — |
| `chebyshev` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 |
| `cityblock` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 |
| `correlation` | 2003 | 3.9271 / 0.092261 | 3.9271 / 0.092261 | — | 3.0868e+16 / 1.5411e+13 |
| `cosine` | 2064 | 2.0455e+18 / 1.0604e+15 | 4.5036e+15 / 1.2108e+13 | — | 4.5036e+15 / 1.2108e+13 |
| `dice` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — |
| `directed_hausdorff` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — |
| `distance_matrix` | 2003 | 0.8165 / 1.1598e-4 | 0.8165 / 1.1598e-4 | 0.8165 / 1.1598e-4 | — |
| `euclidean` | 2064 | 1.5729 / 0.2364 | 1.779 / 0.21919 | 1.779 / 0.21318 | 1.779 / 0.21318 |
| `hamming` | 2003 | 0.5 / 0.17374 | 0.5 / 0.17374 | 0.0 / 0.0 | — |
| `is_valid_dm` | 2008 | 0.0 / 0.0 | 0.0 / 0.0 | — | — |
| `is_valid_y` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — | — |
| `jaccard` | 2005 | 0.0 / 0.0 | 0.0 / 0.0 | — | — |
| `jensenshannon` | 2003 | 0.51954 / 2.5938e-4 | 0.51954 / 2.5938e-4 | — | — |
| `mahalanobis` | 2002 | 1.0006 / 0.21304 | 1.0006 / 0.21304 | — | — |
| `minkowski` | 2004 | 0.89267 / 0.18142 | 0.89267 / 0.18142 | 0.66674 / 0.18106; failed 1/2004 | — |
| `minkowski_distance` | 2004 | 0.66674 / 0.18563 | 0.66674 / 0.18563 | 0.66674 / 0.18572; failed 1/2004 | — |
| `minkowski_distance_p` | 2004 | 0.88883 / 4.4353e-4 | 0.88883 / 4.4353e-4 | 0.88908 / 0.30133; failed 1/2004 | — |
| `num_obs_dm` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — | — |
| `num_obs_y` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — |
| `pdist` | 2019 | 1.8783 / 0.13435 | 1.8783 / 0.13435 | 1.3333 / 0.13399; failed 13/2019 | — |
| `rogerstanimoto` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — |
| `russellrao` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — |
| `seuclidean` | 2002 | 1.2811 / 0.22172 | 1.2811 / 0.22172 | — | — |
| `sokalsneath` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — |
| `sqeuclidean` | 2004 | 1.7782 / 0.45407 | 1.7782 / 0.45407 | 1.7782 / 0.39294 | 1.7782 / 0.45415 |
| `squareform` | 2008 | 0.0 / 0.0 | 0.0 / 0.0 | — | — |
| `yule` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — | — |

### `sci::spatial::transform` (3 functions, 6009 cases)

| Function | Cases | Sci vs reference | SciPy | Boost.Math | GSL |
| --- | ---: | --- | --- | --- | --- |
| `geometric_slerp` | 2004 | 4.5036e+15 / 3.7478e+11 | 4.5036e+15 / 3.7478e+11 | — | — |
| `procrustes` | 2003 | 8.9643e+16 / 1.1526e+15 | 8.9643e+16 / 1.1526e+15 | — | — |
| `tsearch` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — |

## Evaluation failures

- `spatial.cdist` (Boost.Math): 13 failed; metric has no public Boost.Math norm equivalent. Example: `spatial.cdist/recipe-0087`.
- `spatial.minkowski` (Boost.Math): 1 failed; fractional/invalid Minkowski p has no Boost lp_distance equivalent. Example: `spatial.minkowski/recipe-0028`.
- `spatial.minkowski_distance` (Boost.Math): 1 failed; fractional/invalid Minkowski p has no Boost lp_distance equivalent. Example: `spatial.minkowski_distance/recipe-0031`.
- `spatial.minkowski_distance_p` (Boost.Math): 1 failed; fractional/invalid Minkowski p has no Boost lp_distance equivalent. Example: `spatial.minkowski_distance_p/recipe-0034`.
- `spatial.pdist` (Boost.Math): 13 failed; metric has no public Boost.Math norm equivalent. Example: `spatial.pdist/recipe-0086`.


## Conventions and limitations

Metrics, pairwise arrays, geometry and validation helpers use explicit weights, variance/covariance data and shapes from the pinned spatial fixtures. Output ordering follows those fixtures. Near-zero distances and cancellation can amplify relative error; the CSV retains full outputs for absolute-error inspection. Boost.Math public vector norms independently implement bounded L1, L2, L-infinity and integer-p distance kernels; pairwise matrices and weighted distance conventions are documented derived compositions. Unsupported metric names and fractional powers remain unavailable. GSL-labelled distances are explicitly derived from public GSL BLAS norm/dot operations and sample covariance/mean statistics, with weighted Sci conventions applied by the adapter; GSL does not provide identical publicly named distance APIs.

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
