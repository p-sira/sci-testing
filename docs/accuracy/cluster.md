# cluster: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-cluster)

Measured binary64 (`f64`) results using **Sci-cluster as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 62217 cases and 31 public API names;
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

### `sci::cluster::hierarchy` (27 functions, 54114 cases)

| Function | Cases | Sci vs reference | SciPy |
| --- | ---: | --- | --- |
| `hierarchy.average` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 |
| `hierarchy.centroid` | 2004 | 1.254 / 0.010534 | 1.254 / 0.010534 |
| `hierarchy.complete` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 |
| `hierarchy.cophenet` | 2002 | 2.8446 / 0.099211 | 2.8446 / 0.099211 |
| `hierarchy.correspond` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 |
| `hierarchy.cut_tree` | 2006 | 0.0 / 0.0 | 0.0 / 0.0 |
| `hierarchy.fcluster` | 2015 | 0.0 / 0.0 | 0.0 / 0.0 |
| `hierarchy.fclusterdata` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 |
| `hierarchy.from_mlab_linkage` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 |
| `hierarchy.inconsistent` | 2007 | 4.189 / 0.20502 | 4.189 / 0.20502 |
| `hierarchy.is_isomorphic` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 |
| `hierarchy.is_monotonic` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 |
| `hierarchy.is_valid_im` | 2007 | 0.0 / 0.0 | 0.0 / 0.0 |
| `hierarchy.is_valid_linkage` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 |
| `hierarchy.leaders` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 |
| `hierarchy.leaves_list` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 |
| `hierarchy.linkage` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 |
| `hierarchy.maxRstat` | 2007 | 4.5826 / 0.31969 | 4.5826 / 0.31969 |
| `hierarchy.maxdists` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 |
| `hierarchy.maxinconsts` | 2003 | 3.6661 / 0.32916 | 3.6661 / 0.32916 |
| `hierarchy.median` | 2004 | 1.254 / 0.010586 | 1.254 / 0.010586 |
| `hierarchy.num_obs_linkage` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 |
| `hierarchy.single` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 |
| `hierarchy.to_mlab_linkage` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 |
| `hierarchy.to_tree` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 |
| `hierarchy.ward` | 2004 | 1.7727 / 0.024594 | 1.7727 / 0.024594 |
| `hierarchy.weighted` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 |

### `sci::cluster::vq` (4 functions, 8103 cases)

| Function | Cases | Sci vs reference | SciPy |
| --- | ---: | --- | --- |
| `vq.kmeans` | 2000 | 0.0 / 0.0 | 0.0 / 0.0 |
| `vq.kmeans2` | 2000 | 0.0 / 0.0 | 0.0 / 0.0 |
| `vq.vq` | 2052 | 0.87266 / 0.090455 | 0.99702 / 0.10237 |
| `vq.whiten` | 2051 | 1.5493 / 0.31883 | 1.9976 / 0.31042 |


## Conventions and limitations

Vector quantization uses Euclidean distances and first-minimum ties; whitening divides by population standard deviations. Hierarchy fixtures also check counts, leaf ordering and partitions. K-means receives an explicit initial codebook. Returned discrete indices and tree metadata are recorded alongside distances.

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
