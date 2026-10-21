# linalg: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-linalg)

Measured binary64 (`f64`) results using **Sci-linalg as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 218680 cases and 109 public API names;
coverage is limited to the operations and inputs below.

## Agreement with Sci

Each library cell is **maximum / mean difference from Sci in ε** over cases
where both evaluations succeed. `Cases` counts input cases; `Sci failed N/T`
identifies missing baseline evaluations. `failed N/T` counts that library's own
failures; `Sci unavailable N/T` counts its successful evaluations excluded
because Sci failed. A dash means unavailable or no paired finite results.
Means include all compared scalar components, including exact zeros.
Full outputs, references, warnings and errors remain in the per-function CSVs.

| Function | Cases | SciPy | NumPy |
| --- | ---: | --- | --- |
| `bandwidth` | 2003 | 0.0 / 0.0 | — |
| `block_diag` | 2003 | 0.0 / 0.0 | — |
| `cdf2rdf` | 2002 | 0.0 / 0.0 | — |
| `cho_factor` | 2004 | 1.0005 / 0.38293 | — |
| `cho_solve` | 2004 | 10.688 / 1.0401 | — |
| `cho_solve_banded` | 2004 | 10.688 / 0.7667 | — |
| `cholesky` | 2004 | 1.0005 / 0.38417 | — |
| `cholesky_banded` | 2004 | 1.3329 / 0.3594 | — |
| `circulant` | 2003 | 0.0 / 0.0 | — |
| `clarkson_woodruff_transform` | 2003 | 0.0 / 0.0 | — |
| `companion` | 2004 | 0.0 / 0.0 | — |
| `convolution_matrix` | 2009 | 0.0 / 0.0 | — |
| `coshm` | 2004 | 11.805 / 0.83849 | — |
| `cosm` | 2004 | 7.4006 / 0.82619 | — |
| `cossin` | 2010 | 6.0957 / 0.074127 | — |
| `det` | 2058 | 4.3784e+8 / 2.1295e+5 | 4.3784e+8 / 2.1296e+5 |
| `dft` | 2009 | 9.2638e+17 / 3.085e+15 | — |
| `diagsvd` | 2003 | 0.0 / 0.0 | — |
| `dot` | 2000 | — | 0.0 / 0.0 |
| `eig` | 2004 | 2.0 / 0.18063 | — |
| `eig_banded` | 2004 | 6.75 / 0.041538 | — |
| `eigh` | 2006 | 37.787 / 8.3398 | — |
| `eigh_tridiagonal` | 2002 | 8.3971 / 0.011478 | — |
| `eigvals` | 2004 | 12.603 / 6.7129e-3 | — |
| `eigvals_banded` | 2004 | 1.0 / 0.11048 | — |
| `eigvalsh` | 2006 | 3.3844 / 0.93313 | — |
| `eigvalsh_tridiagonal` | 2002 | 2.0 / 2.6287e-3 | — |
| `expm` | 2005 | 4286.1 / 8.4365 | — |
| `expm_cond` | 2004 | 21.521 / 3.9973 | — |
| `expm_frechet` | 2003 | 5.5524 / 2.5624 | — |
| `fiedler` | 2003 | 0.0 / 0.0 | — |
| `fiedler_companion` | 2006 | 0.0 / 0.0 | — |
| `find_best_blas_type` | 2002 | 0.0 / 0.0 | — |
| `fractional_matrix_power` | 2027 | 3.4751e+16 / 6.8747e+12 | — |
| `funm` | 2009 | 4.5036e+15 / 1.6761e+12 | — |
| `get_blas_funcs` | 2005 | 0.0 / 0.0 | — |
| `get_lapack_funcs` | 2004 | 3.9997 / 0.57605 | — |
| `hadamard` | 2003 | 0.0 / 0.0 | — |
| `hankel` | 2003 | 0.0 / 0.0 | — |
| `helmert` | 2006 | 0.0 / 0.0 | — |
| `hessenberg` | 2005 | 4.6022e+17 / 1.9136e+15 | — |
| `hilbert` | 2004 | 0.0 / 0.0 | — |
| `interpolative.estimate_rank` | 2004 | 0.0 / 0.0 | — |
| `interpolative.estimate_spectral_norm` | 2004 | 4.7852 / 1.4 | — |
| `interpolative.estimate_spectral_norm_diff` | 2004 | 5.0648 / 1.0155 | — |
| `interpolative.id_to_svd` | 2003 | 7.0 / 0.94978 | — |
| `interpolative.interp_decomp` | 2008 | 4.7852 / 0.67558 | — |
| `interpolative.reconstruct_interp_matrix` | 2004 | 0.0 / 0.0 | — |
| `interpolative.reconstruct_matrix_from_id` | 2004 | 0.0 / 0.0 | — |
| `interpolative.reconstruct_skel_matrix` | 2003 | 0.0 / 0.0 | — |
| `interpolative.svd` | 2008 | 6.5 / 5.4349e-3 | — |
| `inv` | 2003 | 1.0 / 1.473e-3 | — |
| `invhilbert` | 2004 | 0.95272 / 1.2437e-3 | — |
| `invpascal` | 2009 | 0.0 / 0.0 | — |
| `ishermitian` | 2004 | 0.0 / 0.0 | — |
| `issymmetric` | 2004 | 0.0 / 0.0 | — |
| `khatri_rao` | 2002 | 0.0 / 0.0 | — |
| `ldl` | 2008 | 0.0 / 0.0 | — |
| `leslie` | 2002 | 0.0 / 0.0 | — |
| `logm` | 2009 | 4.5036e+15 / 4.9661e+12 | — |
| `lstsq` | 2012 | 4.5036e+15 / 1.1172e+12 | — |
| `lu` | 2014 | 0.75 / 1.1964e-4 | — |
| `lu_factor` | 2007 | 0.75 / 1.8437e-4 | — |
| `lu_solve` | 2006 | 1.6842 / 1.8257e-3 | — |
| `matmul` | 2000 | — | 0.0 / 0.0 |
| `matmul_toeplitz` | 2002 | 1.6015 / 0.37223 | — |
| `matrix_balance` | 2006 | 5.2355 / 0.17031 | — |
| `matvec` | 2000 | — | 1.0 / 0.077224 |
| `norm` | 2024 | 1.4368 / 7.099e-4 | — |
| `null_space` | 2017 | 7.5 / 1.2374 | — |
| `ordqz` | 2008 | 1.0273e+17 / 2.7242e+13 | — |
| `orth` | 2017 | 5.1857e+16 / 1.0322e+13 | — |
| `orthogonal_procrustes` | 2002 | 4.5036e+15 / 1.1116e+14 | — |
| `pascal` | 2009 | 0.0 / 0.0 | — |
| `pinv` | 2017 | 13.674 / 2.3156 | — |
| `pinvh` | 2018 | 53.341 / 2.1188 | — |
| `polar` | 2006 | 3.1481e+18 / 6.0773e+15 | — |
| `qr` | 2005 | 1.0845e+17 / 8.6269e+15 | — |
| `qr_delete` | 2003 | 1.5 / 2.7053e-4 | — |
| `qr_insert` | 2003 | 1.0852e+19 / 6.7709e+14 | — |
| `qr_multiply` | 2004 | 27.123 / 4.3737 | — |
| `qr_update` | 2002 | 10.019 / 1.1203 | — |
| `qz` | 2004 | 9.2181e+16 / 1.4036e+13 | — |
| `rq` | 2003 | 1.1699e+16 / 3.873e+12 | — |
| `rsf2csf` | 2002 | 4.5036e+15 / 6.7993e+14 | — |
| `schur` | 2005 | 9.694e+17 / 2.3762e+15 | — |
| `signm` | 2007 | 1.8495e+16 / 1.7097e+12 | — |
| `sinhm` | 2004 | 29.297 / 0.6562 | — |
| `sinm` | 2004 | 4.2546 / 0.68412 | — |
| `solve` | 2058 | 7.4052e+7 / 47552.0 | 7.835e+8 / 4.4136e+5 |
| `solve_banded` | 2003 | 0.66733 / 0.069514 | — |
| `solve_circulant` | 2002 | 2.625 / 1.3099e-3 | — |
| `solve_continuous_are` | 2003 | 5.6569 / 0.55421 | — |
| `solve_continuous_lyapunov` | 2002 | 0.0 / 0.0 | — |
| `solve_discrete_are` | 2003 | 4.5036e+15 / 4.4834e+12 | — |
| `solve_discrete_lyapunov` | 2002 | 1.8601 / 0.20185 | — |
| `solve_lyapunov` | 2002 | 0.0 / 0.0 | — |
| `solve_sylvester` | 2003 | 2.0 / 0.15424 | — |
| `solve_toeplitz` | 2002 | 3.4997 / 0.42147 | — |
| `solve_triangular` | 2008 | 1.3321 / 0.031869 | — |
| `solveh_banded` | 2004 | 8.4375 / 0.77894 | — |
| `sqrtm` | 2009 | 4.5036e+15 / 4.4684e+12 | — |
| `subspace_angles` | 2004 | 0.0 / 0.0 | — |
| `svd` | 2005 | 1.0525e+18 / 1.5578e+16 | — |
| `svdvals` | 2005 | 3.6493 / 0.8805 | — |
| `tanhm` | 2004 | 29.297 / 0.71794 | — |
| `tanm` | 2004 | 4.8006 / 0.63213 | — |
| `toeplitz` | 2003 | 0.0 / 0.0 | — |
| `transpose` | 2000 | — | 0.0 / 0.0 |


## Conventions and limitations

Nonunique decompositions are compared through ordered eigenvalues, projectors and reconstructions, retaining independent structural invariants from the Sci fixture dispatcher. Interpolative decomposition compares spectral residual and rank; its SVD compares reconstruction and singular values. Type metadata uses the Sci f64 spelling. Dense sample sizes and constructor orders stay bounded.

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
