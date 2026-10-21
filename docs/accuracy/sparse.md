# sparse: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-sparse)

Measured binary64 (`f64`) results using **Sci-sparse as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 196255 cases and 98 public API names;
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
| `add` | 2002 | 0.0 / 0.0 |
| `as_array` | 2002 | 0.0 / 0.0 |
| `block_array` | 2002 | 0.0 / 0.0 |
| `block_diag` | 2002 | 0.0 / 0.0 |
| `bmat` | 2002 | 0.0 / 0.0 |
| `csgraph.bellman_ford` | 2003 | 0.0 / 0.0 |
| `csgraph.breadth_first_order` | 2002 | 0.0 / 0.0 |
| `csgraph.breadth_first_tree` | 2002 | 0.0 / 0.0 |
| `csgraph.connected_components` | 2002 | 0.0 / 0.0 |
| `csgraph.construct_dist_matrix` | 2002 | 0.0 / 0.0 |
| `csgraph.csgraph_from_dense` | 2002 | 0.0 / 0.0 |
| `csgraph.csgraph_from_masked` | 2002 | 0.0 / 0.0 |
| `csgraph.csgraph_masked_from_dense` | 2002 | 0.0 / 0.0 |
| `csgraph.csgraph_to_dense` | 2002 | 0.0 / 0.0 |
| `csgraph.csgraph_to_masked` | 2002 | 0.0 / 0.0 |
| `csgraph.depth_first_order` | 2002 | 0.0 / 0.0 |
| `csgraph.depth_first_tree` | 2002 | 0.0 / 0.0 |
| `csgraph.dijkstra` | 2003 | 0.0 / 0.0 |
| `csgraph.floyd_warshall` | 2003 | 0.0 / 0.0 |
| `csgraph.johnson` | 2003 | 0.0 / 0.0 |
| `csgraph.laplacian` | 2003 | 0.0 / 0.0 |
| `csgraph.maximum_bipartite_matching` | 2002 | 0.0 / 0.0 |
| `csgraph.maximum_flow` | 2002 | 0.0 / 0.0 |
| `csgraph.min_weight_full_bipartite_matching` | 2002 | 0.0 / 0.0 |
| `csgraph.minimum_spanning_tree` | 2002 | 0.0 / 0.0 |
| `csgraph.reconstruct_path` | 2002 | 0.0 / 0.0 |
| `csgraph.reverse_cuthill_mckee` | 2002 | 0.0 / 0.0 |
| `csgraph.shortest_path` | 2003 | 0.0 / 0.0 |
| `csgraph.structural_rank` | 2002 | 0.0 / 0.0 |
| `csgraph.yen` | 2002 | 0.8008 / 0.067462 |
| `diags` | 2002 | 0.0 / 0.0 |
| `diags_array` | 2002 | 0.0 / 0.0 |
| `eliminate_zeros` | 2002 | 0.0 / 0.0 |
| `expand_dims` | 2002 | 0.0 / 0.0 |
| `eye` | 2005 | 0.0 / 0.0 |
| `eye_array` | 2005 | 0.0 / 0.0 |
| `find` | 2002 | 0.0 / 0.0 |
| `from_coo` | 2002 | 0.0 / 0.0 |
| `from_csc` | 2002 | 0.0 / 0.0 |
| `from_csr` | 2002 | 0.0 / 0.0 |
| `from_dense` | 2002 | 0.0 / 0.0 |
| `get_index_dtype` | 2002 | 0.0 / 0.0 |
| `hstack` | 2002 | 0.0 / 0.0 |
| `identity` | 2003 | 0.0 / 0.0 |
| `issparse` | 2002 | 0.0 / 0.0 |
| `kron` | 2002 | 0.0 / 0.0 |
| `kronsum` | 2003 | 0.0 / 0.0 |
| `linalg.aslinearoperator` | 2002 | 0.0 / 0.0 |
| `linalg.bicg` | 2003 | 0.0 / 0.0 |
| `linalg.bicgstab` | 2003 | 4.7972e+20 / 1.0563e+17 |
| `linalg.cg` | 2002 | 1.6855e+19 / 5.4627e+15 |
| `linalg.cgs` | 2003 | 0.0 / 0.0 |
| `linalg.eigs` | 2003 | 2.5381 / 0.13993 |
| `linalg.eigsh` | 2002 | 3.3846 / 0.5829 |
| `linalg.expm` | 2002 | 40.587 / 2.8384 |
| `linalg.expm_multiply` | 2002 | 36.57 / 4.502 |
| `linalg.factorized` | 2002 | 7.2572 / 0.41397 |
| `linalg.funm_multiply_krylov` | 2002 | 179.2 / 35.401 |
| `linalg.gcrotmk` | 2004 | 3.1525e+16 / 1.9457e+15 |
| `linalg.gmres` | 2006; Sci failed 1/2006 | 3.1525e+16 / 1.6664e+15; Sci unavailable 1/2006 |
| `linalg.inv` | 2002 | 0.0 / 0.0 |
| `linalg.is_sptriangular` | 2002 | 0.0 / 0.0 |
| `linalg.lgmres` | 2004 | 3.6029e+16 / 2.1864e+15 |
| `linalg.lobpcg` | 2002 | 3.3846 / 0.58839 |
| `linalg.lsmr` | 2004 | 6.5625 / 0.74906 |
| `linalg.lsqr` | 2004 | 9.375 / 0.64381 |
| `linalg.matrix_power` | 2005 | 0.0 / 0.0 |
| `linalg.minres` | 2002 | 5.8547e+16 / 2.0425e+15 |
| `linalg.norm` | 2014 | 2.5383 / 0.34634 |
| `linalg.norm_axis` | 2002 | 1.5015 / 0.16459 |
| `linalg.onenormest` | 2002 | 0.0 / 0.0 |
| `linalg.qmr` | 2003 | 3.1525e+16 / 7.8964e+14 |
| `linalg.spbandwidth` | 2002 | 0.0 / 0.0 |
| `linalg.spilu` | 2002 | 4.5045 / 0.30706 |
| `linalg.splu` | 2002 | 7.2572 / 0.41883 |
| `linalg.spsolve` | 2002 | 7.2572 / 0.42219 |
| `linalg.spsolve_triangular` | 2004 | 0.75075 / 0.025823 |
| `linalg.svds` | 2003 | 2.538 / 0.34642 |
| `linalg.tfqmr` | 2003 | 7.6561e+16 / 1.7852e+15 |
| `load_npz` | 2002 | 0.0 / 0.0 |
| `matmul` | 2003 | 0.0 / 0.0 |
| `matrix_transpose` | 2002 | 0.0 / 0.0 |
| `matvec` | 2002 | 0.0 / 0.0 |
| `multiply` | 2002 | 0.0 / 0.0 |
| `permute_dims` | 2002 | 0.0 / 0.0 |
| `random` | 2003 | 0.0 / 0.0 |
| `random_array` | 2003 | 0.0 / 0.0 |
| `safely_cast_index_arrays` | 2002 | 0.0 / 0.0 |
| `save_npz` | 2002 | 0.0 / 0.0 |
| `scale` | 2002 | 0.0 / 0.0 |
| `spdiags` | 2002 | 0.0 / 0.0 |
| `swapaxes` | 2002 | 0.0 / 0.0 |
| `to_dense` | 2002 | 0.0 / 0.0 |
| `transpose` | 2002 | 0.0 / 0.0 |
| `tril` | 2005 | 0.0 / 0.0 |
| `triu` | 2005 | 0.0 / 0.0 |
| `validate` | 2002 | 0.0 / 0.0 |
| `vstack` | 2002 | 0.0 / 0.0 |

## Evaluation failures

- `sparse.linalg.gmres` (Sci): 1 failed; gmres_scale_1e-20 residual 0.0000000000000000000013153341044116428. Example: `sparse.linalg.gmres/review-scale-1e-20`.


## Conventions and limitations

Canonical CSR and COO arrays retain explicit zeros. Matrix outputs include dimensions; graph distances retain infinities. Eigenpairs are compared by ordered values with independent residual checks; matching cardinality, permutation validity, flow conservation and random shape/count use invariants where outputs are nonunique. Sparse LU, iterative solves and matrix functions use pinned public SciPy references. Each scalar uses (class,value): class 0 for finite data, ±1 for signed infinity and 2 for NaN, with zero nonfinite payloads. This preserves unreachable vertices; exact tags contribute to component means. RNG streams differ. Dense samples preserve graph structure, integer capacities and dependent dimensions; fixed discrete profiles may repeat. Rust signature adaptations and bounded dense fallbacks are documented in Sci sparse coverage.

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
