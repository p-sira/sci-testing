# sparse: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-sparse)

Measured binary64 (`f64`) results using **Sci-sparse as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 198293 cases and 100 public API names;
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

### `sci::sparse::array` (7 functions, 12049 cases)

| Function | Cases | Sci vs reference | SciPy | GSL | Boost.Graph | Boost.Numeric.uBLAS |
| --- | ---: | --- | --- | --- | --- | --- |
| `array.csr_indices_i32` | 1019 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `array.csr_indices_i64` | 1019 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `as_array` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `expand_dims` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `permute_dims` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `random_array` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `swapaxes` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |

### `sci::sparse::construct` (16 functions, 32038 cases)

| Function | Cases | Sci vs reference | SciPy | GSL | Boost.Graph | Boost.Numeric.uBLAS |
| --- | ---: | --- | --- | --- | --- | --- |
| `block_array` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `block_diag` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `diags` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `eye` | 2005 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `from_coo` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `from_csc` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `from_csr` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `from_dense` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | — | 0.0 / 0.0 |
| `hstack` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `identity` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `kron` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `kronsum` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `random` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `spdiags` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `to_dense` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | — | 0.0 / 0.0 |
| `vstack` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |

### `sci::sparse::csgraph` (25 functions, 50056 cases)

| Function | Cases | Sci vs reference | SciPy | GSL | Boost.Graph | Boost.Numeric.uBLAS |
| --- | ---: | --- | --- | --- | --- | --- |
| `csgraph.bellman_ford` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | 0.0 / 0.0 | — |
| `csgraph.breadth_first_order` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | 0.0 / 0.0 | — |
| `csgraph.breadth_first_tree` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | 0.0 / 0.0 | — |
| `csgraph.connected_components` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | 0.0 / 0.0 | — |
| `csgraph.construct_dist_matrix` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `csgraph.csgraph_from_dense` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `csgraph.csgraph_from_masked` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `csgraph.csgraph_masked_from_dense` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `csgraph.csgraph_to_dense` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `csgraph.csgraph_to_masked` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `csgraph.depth_first_order` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | 0.0 / 0.0 | — |
| `csgraph.depth_first_tree` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | 0.0 / 0.0 | — |
| `csgraph.dijkstra` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | 0.0 / 0.0 | — |
| `csgraph.floyd_warshall` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | 0.0 / 0.0 | — |
| `csgraph.johnson` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | 0.0 / 0.0 | — |
| `csgraph.laplacian` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | 1.0 / 7.2728e-5 | — |
| `csgraph.maximum_bipartite_matching` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | 0.0 / 0.0 | — |
| `csgraph.maximum_flow` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | 0.0 / 0.0 | — |
| `csgraph.min_weight_full_bipartite_matching` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `csgraph.minimum_spanning_tree` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | 0.0 / 0.0 | — |
| `csgraph.reconstruct_path` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `csgraph.reverse_cuthill_mckee` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `csgraph.shortest_path` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | 0.0 / 0.0 | — |
| `csgraph.structural_rank` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | 0.0 / 0.0 | — |
| `csgraph.yen` | 2002 | 0.8008 / 0.067462 | 0.8008 / 0.067462 | — | — | — |

### `sci::sparse::io` (2 functions, 4004 cases)

| Function | Cases | Sci vs reference | SciPy | GSL | Boost.Graph | Boost.Numeric.uBLAS |
| --- | ---: | --- | --- | --- | --- | --- |
| `load_npz` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `save_npz` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |

### `sci::sparse::linalg` (32 functions, 64100 cases)

| Function | Cases | Sci vs reference | SciPy | GSL | Boost.Graph | Boost.Numeric.uBLAS |
| --- | ---: | --- | --- | --- | --- | --- |
| `linalg.aslinearoperator` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `linalg.bicg` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `linalg.bicgstab` | 2003 | 4.7972e+20 / 1.0563e+17 | 4.7972e+20 / 1.0563e+17 | — | — | — |
| `linalg.cg` | 2002 | 1.6855e+19 / 5.4627e+15 | 1.6855e+19 / 5.4627e+15 | — | — | — |
| `linalg.cgs` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `linalg.eigs` | 2003 | 2.5381 / 0.13993 | 2.5381 / 0.13993 | — | — | — |
| `linalg.eigsh` | 2002 | 3.3846 / 0.5829 | 3.3846 / 0.5829 | — | — | — |
| `linalg.expm` | 2002 | 40.587 / 2.8384 | 40.587 / 2.8384 | — | — | — |
| `linalg.expm_multiply` | 2002 | 36.57 / 4.502 | 36.57 / 4.502 | — | — | — |
| `linalg.factorized` | 2002 | 7.2572 / 0.41397 | 7.2572 / 0.41397 | — | — | — |
| `linalg.funm_multiply_krylov` | 2002 | 179.2 / 35.401 | 179.2 / 35.401 | — | — | — |
| `linalg.gcrotmk` | 2004 | 3.1525e+16 / 1.9457e+15 | 3.1525e+16 / 1.9457e+15 | — | — | — |
| `linalg.gmres` | 2006; Sci failed 1/2006 | 3.1525e+16 / 1.6664e+15 | 3.1525e+16 / 1.6664e+15; Sci unavailable 1/2006 | — | — | — |
| `linalg.inv` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `linalg.is_sptriangular` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `linalg.lgmres` | 2004 | 3.6029e+16 / 2.1864e+15 | 3.6029e+16 / 2.1864e+15 | — | — | — |
| `linalg.lobpcg` | 2002 | 3.3846 / 0.58839 | 3.3846 / 0.58839 | — | — | — |
| `linalg.lsmr` | 2004 | 6.5625 / 0.74906 | 6.5625 / 0.74906 | — | — | — |
| `linalg.lsqr` | 2004 | 9.375 / 0.64381 | 9.375 / 0.64381 | — | — | — |
| `linalg.matrix_power` | 2005 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `linalg.minres` | 2002 | 5.8547e+16 / 2.0425e+15 | 5.8547e+16 / 2.0425e+15 | — | — | — |
| `linalg.norm` | 2014 | 2.5383 / 0.34634 | 2.5383 / 0.34634 | — | — | — |
| `linalg.norm_axis` | 2002 | 1.5015 / 0.16459 | 1.5015 / 0.16459 | — | — | — |
| `linalg.onenormest` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `linalg.qmr` | 2003 | 3.1525e+16 / 7.8964e+14 | 3.1525e+16 / 7.8964e+14 | — | — | — |
| `linalg.spbandwidth` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `linalg.spilu` | 2002 | 4.5045 / 0.30706 | 4.5045 / 0.30706 | — | — | — |
| `linalg.splu` | 2002 | 7.2572 / 0.41883 | 7.2572 / 0.41883 | — | — | — |
| `linalg.spsolve` | 2002 | 7.2572 / 0.42219 | 7.2572 / 0.42219 | — | — | — |
| `linalg.spsolve_triangular` | 2004 | 0.75075 / 0.025823 | 0.75075 / 0.025823 | — | — | — |
| `linalg.svds` | 2003 | 2.538 / 0.34642 | 2.538 / 0.34642 | — | — | — |
| `linalg.tfqmr` | 2003 | 7.6561e+16 / 1.7852e+15 | 7.6561e+16 / 1.7852e+15 | — | — | — |

### `sci::sparse::operations` (11 functions, 22029 cases)

| Function | Cases | Sci vs reference | SciPy | GSL | Boost.Graph | Boost.Numeric.uBLAS |
| --- | ---: | --- | --- | --- | --- | --- |
| `add` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | — | 0.0 / 0.0 |
| `eliminate_zeros` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `find` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `matmul` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | — | 0.0 / 0.0 |
| `matvec` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | — | 0.0 / 0.0 |
| `multiply` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | 0.0 / 0.0 |
| `scale` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `transpose` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | — | 0.0 / 0.0 |
| `tril` | 2005 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `triu` | 2005 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `validate` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |

### `sci::sparse` (7 functions, 14017 cases)

| Function | Cases | Sci vs reference | SciPy | GSL | Boost.Graph | Boost.Numeric.uBLAS |
| --- | ---: | --- | --- | --- | --- | --- |
| `bmat` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `diags_array` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `eye_array` | 2005 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `get_index_dtype` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `issparse` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `matrix_transpose` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | — | — |
| `safely_cast_index_arrays` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |

## Evaluation failures

- `sparse.linalg.gmres` (Sci): 1 failed; gmres_scale_1e-20 residual 0.0000000000000000000013153341044116428. Example: `sparse.linalg.gmres/review-scale-1e-20`.


## Conventions and limitations

Canonical CSR and COO arrays retain explicit zeros. Matrix outputs include dimensions; graph distances retain infinities. Eigenpairs are compared by ordered values with independent residual checks; matching cardinality, permutation validity, flow conservation and random shape/count use invariants where outputs are nonunique. Sparse LU, iterative solves and matrix functions use pinned public SciPy references. Each scalar uses (class,value): class 0 for finite data, ±1 for signed infinity and 2 for NaN, with zero nonfinite payloads. This preserves unreachable vertices; exact tags contribute to component means. RNG streams differ. Dense samples preserve graph structure, integer capacities and dependent dimensions; fixed discrete profiles may repeat. Rust signature adaptations and bounded dense fallbacks are documented in Sci sparse coverage. Boost.Graph comparisons use public traversal, connected-component, shortest-path and minimum-spanning-tree algorithms on saved adjacency matrices; Boost.Numeric.uBLAS supplies sparse matrix arithmetic. Nonfinite distances use Sci's explicit (class,value) encoding, and matrix shape metadata is preserved. These adapters are not a general replacement for arbitrary Sci graph options or sparse formats.

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
