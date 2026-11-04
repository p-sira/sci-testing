# spatial: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-spatial)

Measured binary64 (`f64`) results compared against the **saved per-case reference**.
Each implementation, including Sci-spatial, is evaluated against the same
reference values stored in the raw CSVs. This page
covers 66272 cases and 33 public API names;
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

### `sci::spatial::distance`

30 functions · 60263 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | GSL | NumPy |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `braycurtis` | 2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | 0.0 / 0.0 |
| `canberra` | 2004 | 1.7148 / 0.19462 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | 0.0 / 0.0 |
| `cdist` | 2019 | 4.5036e+15 / 3.7177e+11 | —; unavailable 2019/2019 | 0.0 / 0.0 | 0.8 / 1.7725e-4; failed 13/2019 | —; unavailable 2019/2019 | 117.1 / 9.8548e-3; unavailable 15/2019 |
| `chebyshev` | 2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 |
| `cityblock` | 2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 |
| `correlation` | 2003 | 3.9271 / 0.092261 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 3.0868e+16 / 1.5411e+13 | 1.9635 / 1.3855e-3 |
| `cosine` | 2064 | 2.0455e+18 / 1.0604e+15 | 4.5036e+15 / 6.5459e+12 | 4.5036e+15 / 9.9237e+12 | —; unavailable 2064/2064 | 4.5036e+15 / 9.9237e+12 | 4.5036e+15 / 9.9237e+12 |
| `dice` | 2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `directed_hausdorff` | 2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `distance_matrix` | 2003 | 0.8165 / 1.1598e-4 | —; unavailable 2003/2003 | 0.0 / 0.0 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 |
| `euclidean` | 2064 | 1.5729 / 0.2364 | 0.47748 / 0.17168 | 1.0718 / 0.23593 | 1.2822 / 0.23568 | 1.2822 / 0.23568 | 1.0718 / 0.23593 |
| `hamming` | 2003 | 0.5 / 0.17374 | —; unavailable 2003/2003 | 0.0 / 0.0 | 0.5 / 0.17374 | —; unavailable 2003/2003 | 0.5 / 0.17374 |
| `is_valid_dm` | 2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | 0.0 / 0.0 |
| `is_valid_y` | 2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | 0.0 / 0.0 |
| `jaccard` | 2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | —; unavailable 2005/2005 | 0.0 / 0.0 |
| `jensenshannon` | 2003 | 0.51954 / 2.5938e-4 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `mahalanobis` | 2002 | 1.0006 / 0.21304 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `minkowski` | 2004 | 0.89267 / 0.18142 | —; unavailable 2004/2004 | 0.0 / 0.0 | 0.0 / 0.0; failed 1/2004 | —; unavailable 2004/2004 | 0.89267 / 4.4545e-4 |
| `minkowski_distance` | 2004 | 0.66674 / 0.18563 | —; unavailable 2004/2004 | 0.0 / 0.0 | 0.0 / 0.0; failed 1/2004 | —; unavailable 2004/2004 | 0.0 / 0.0 |
| `minkowski_distance_p` | 2004 | 0.88883 / 4.4353e-4 | —; unavailable 2004/2004 | 0.0 / 0.0 | 0.88908 / 0.30177; failed 1/2004 | —; unavailable 2004/2004 | 0.0 / 0.0; unavailable 1/2004 |
| `num_obs_dm` | 2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | 0.0 / 0.0 |
| `num_obs_y` | 2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | 0.0 / 0.0 |
| `pdist` | 2019 | 1.8783 / 1.4402e-3 | —; unavailable 2019/2019 | 0.0 / 0.0 | 0.66667 / 1.1078e-4; failed 13/2019 | —; unavailable 2019/2019 | 0.0 / 0.0; unavailable 15/2019 |
| `rogerstanimoto` | 2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `russellrao` | 2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `seuclidean` | 2002 | 1.2811 / 0.22172 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `sokalsneath` | 2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `sqeuclidean` | 2004 | 1.7782 / 0.45407 | —; unavailable 2004/2004 | 0.0 / 0.0 | 0.88908 / 0.30587 | 0.79012 / 7.0135e-4 | 0.0 / 0.0 |
| `squareform` | 2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | 0.0 / 0.0; unavailable 489/2008 |
| `yule` | 2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |

### `sci::spatial::transform`

3 functions · 6009 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | GSL | NumPy |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `geometric_slerp` | 2004 | 4.5036e+15 / 3.7478e+11 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | 4.5036e+15 / 5.6232e+11; unavailable 1/2004 |
| `procrustes` | 2003 | 8.9643e+16 / 1.1526e+15 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `tsearch` | 2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |

## Evaluation failures

- `spatial.cdist` (Boost.Math): 13 failed; metric has no public Boost.Math norm equivalent. Example: `spatial.cdist/recipe-0087`.
- `spatial.minkowski` (Boost.Math): 1 failed; fractional/invalid Minkowski p has no Boost lp_distance equivalent. Example: `spatial.minkowski/recipe-0028`.
- `spatial.minkowski_distance` (Boost.Math): 1 failed; fractional/invalid Minkowski p has no Boost lp_distance equivalent. Example: `spatial.minkowski_distance/recipe-0031`.
- `spatial.minkowski_distance_p` (Boost.Math): 1 failed; fractional/invalid Minkowski p has no Boost lp_distance equivalent. Example: `spatial.minkowski_distance_p/recipe-0034`.
- `spatial.pdist` (Boost.Math): 13 failed; metric has no public Boost.Math norm equivalent. Example: `spatial.pdist/recipe-0086`.


## Conventions and limitations

Metrics, pairwise arrays, geometry and validation helpers use explicit weights, variance/covariance data and shapes from the pinned spatial fixtures. Output ordering follows those fixtures. Near-zero distances and cancellation can amplify relative error; the CSV retains full outputs for absolute-error inspection. Boost.Math public vector norms independently implement bounded L1, L2, L-infinity and integer-p distance kernels; pairwise matrices and weighted distance conventions are documented derived compositions. Unsupported metric names and fractional powers remain unavailable. GSL-labelled distances are explicitly derived from public GSL BLAS norm/dot operations and sample covariance/mean statistics, with weighted Sci conventions applied by the adapter; GSL does not provide identical publicly named distance APIs.

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
