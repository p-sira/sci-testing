# linalg: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-linalg)

Measured binary64 (`f64`) results using **Sci-linalg as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 218680 cases and 109 public API names;
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

### `sci::linalg::basic`

13 functions · 26133 cases

| Function | Cases | Sci vs reference | SciPy | GSL | NumPy | Boost.Numeric.uBLAS |
| --- | ---: | --- | --- | --- | --- | --- |
| `bandwidth` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `det` | 2058 | 2.5178e+8 / 1.2242e+5 | 4.3784e+8 / 2.1295e+5 | 2.2925e+8 / 1.1153e+5 | 4.3784e+8 / 2.1296e+5 | 2.518e+8 / 1.2243e+5 |
| `dot` | 2000 | 0.0 / 0.0 | — | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 |
| `inv` | 2003 | 1.0 / 1.473e-3 | 1.0 / 1.473e-3 | 1.0 / 1.473e-3 | — | — |
| `ishermitian` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `issymmetric` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `khatri_rao` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `matmul` | 2000 | 0.0 / 0.0 | — | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 |
| `matvec` | 2000 | 1.0 / 0.077224 | — | 0.0 / 0.0 | 1.0 / 0.077224 | 0.0 / 0.0 |
| `norm` | 2024 | 1.4368 / 7.099e-4 | 1.4368 / 7.099e-4 | — | — | — |
| `pinv` | 2017 | 13.674 / 2.3156 | 13.674 / 2.3156 | — | — | — |
| `pinvh` | 2018 | 53.341 / 2.1188 | 53.341 / 2.1188 | — | — | — |
| `transpose` | 2000 | 0.0 / 0.0 | — | — | 0.0 / 0.0 | 0.0 / 0.0 |

### `sci::linalg::decompositions`

20 functions · 40105 cases

| Function | Cases | Sci vs reference | SciPy | GSL | NumPy | Boost.Numeric.uBLAS |
| --- | ---: | --- | --- | --- | --- | --- |
| `cho_factor` | 2004 | 1.0005 / 0.38293 | 1.0005 / 0.38293 | — | — | — |
| `cholesky` | 2004 | 1.0005 / 0.38417 | 1.0005 / 0.38417 | 0.0 / 0.0 | — | — |
| `cholesky_banded` | 2004 | 1.3329 / 0.3594 | 1.3329 / 0.3594 | — | — | — |
| `cossin` | 2010 | 6.0957 / 0.074127 | 6.0957 / 0.074127 | — | — | — |
| `hessenberg` | 2005 | 4.6022e+17 / 1.9136e+15 | 4.6022e+17 / 1.9136e+15 | — | — | — |
| `ldl` | 2008 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `lu` | 2014 | 0.75 / 1.1964e-4 | 0.75 / 1.1964e-4 | — | — | — |
| `lu_factor` | 2007 | 0.75 / 1.8437e-4 | 0.75 / 1.8437e-4 | — | — | — |
| `matrix_balance` | 2006 | 5.2355 / 0.17031 | 5.2355 / 0.17031 | — | — | — |
| `polar` | 2006 | 3.1481e+18 / 6.0773e+15 | 3.1481e+18 / 6.0773e+15 | — | — | — |
| `qr` | 2005 | 1.0845e+17 / 8.6269e+15 | 1.0845e+17 / 8.6269e+15 | — | — | — |
| `qr_delete` | 2003 | 1.5 / 2.7053e-4 | 1.5 / 2.7053e-4 | — | — | — |
| `qr_insert` | 2003 | 1.0852e+19 / 6.7709e+14 | 1.0852e+19 / 6.7709e+14 | — | — | — |
| `qr_multiply` | 2004 | 27.123 / 4.3737 | 27.123 / 4.3737 | — | — | — |
| `qr_update` | 2002 | 10.019 / 1.1203 | 10.019 / 1.1203 | — | — | — |
| `rq` | 2003 | 1.1699e+16 / 3.873e+12 | 1.1699e+16 / 3.873e+12 | — | — | — |
| `rsf2csf` | 2002 | 4.5036e+15 / 6.7993e+14 | 4.5036e+15 / 6.7993e+14 | — | — | — |
| `schur` | 2005 | 9.694e+17 / 2.3762e+15 | 9.694e+17 / 2.3762e+15 | — | — | — |
| `svd` | 2005 | 1.0525e+18 / 1.5578e+16 | 1.0525e+18 / 1.5578e+16 | — | — | — |
| `svdvals` | 2005 | 3.6493 / 0.8805 | 3.6493 / 0.8805 | 3.6516 / 0.84989 | — | — |

### `sci::linalg::eigen`

11 functions · 22046 cases

| Function | Cases | Sci vs reference | SciPy | GSL | NumPy | Boost.Numeric.uBLAS |
| --- | ---: | --- | --- | --- | --- | --- |
| `cdf2rdf` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `eig` | 2004 | 2.0 / 0.18063 | 2.0 / 0.18063 | — | — | — |
| `eig_banded` | 2004 | 6.75 / 0.041538 | 6.75 / 0.041538 | — | — | — |
| `eigh` | 2006 | 37.787 / 8.3398 | 37.787 / 8.3398 | — | — | — |
| `eigh_tridiagonal` | 2002 | 8.3971 / 0.011478 | 8.3971 / 0.011478 | — | — | — |
| `eigvals` | 2004 | 12.603 / 6.7129e-3 | 12.603 / 6.7129e-3 | — | — | — |
| `eigvals_banded` | 2004 | 1.0 / 0.11048 | 1.0 / 0.11048 | — | — | — |
| `eigvalsh` | 2006 | 3.3844 / 0.93313 | 3.3844 / 0.93313 | 5.0704 / 0.9063 | — | — |
| `eigvalsh_tridiagonal` | 2002 | 2.0 / 2.6287e-3 | 2.0 / 2.6287e-3 | — | — | — |
| `ordqz` | 2008 | 1.0273e+17 / 2.7242e+13 | 1.0273e+17 / 2.7242e+13 | — | — | — |
| `qz` | 2004 | 9.2181e+16 / 1.4036e+13 | 9.2181e+16 / 1.4036e+13 | — | — | — |

### `sci::linalg::interpolative`

9 functions · 18042 cases

| Function | Cases | Sci vs reference | SciPy | GSL | NumPy | Boost.Numeric.uBLAS |
| --- | ---: | --- | --- | --- | --- | --- |
| `interpolative.estimate_rank` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `interpolative.estimate_spectral_norm` | 2004 | 4.7852 / 1.4 | 4.7852 / 1.4 | — | — | — |
| `interpolative.estimate_spectral_norm_diff` | 2004 | 5.0648 / 1.0155 | 5.0648 / 1.0155 | — | — | — |
| `interpolative.id_to_svd` | 2003 | 7.0 / 0.94978 | 7.0 / 0.94978 | — | — | — |
| `interpolative.interp_decomp` | 2008 | 4.7852 / 0.67558 | 4.7852 / 0.67558 | — | — | — |
| `interpolative.reconstruct_interp_matrix` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `interpolative.reconstruct_matrix_from_id` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `interpolative.reconstruct_skel_matrix` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `interpolative.svd` | 2008 | 6.5 / 5.4349e-3 | 6.5 / 5.4349e-3 | — | — | — |

### `sci::linalg::matrix_functions`

14 functions · 28097 cases

| Function | Cases | Sci vs reference | SciPy | GSL | NumPy | Boost.Numeric.uBLAS |
| --- | ---: | --- | --- | --- | --- | --- |
| `coshm` | 2004 | 11.805 / 0.83849 | 11.805 / 0.83849 | — | — | — |
| `cosm` | 2004 | 7.4006 / 0.82619 | 7.4006 / 0.82619 | — | — | — |
| `expm` | 2005 | 4286.1 / 8.4365 | 4286.1 / 8.4365 | — | — | — |
| `expm_cond` | 2004 | 21.521 / 3.9973 | 21.521 / 3.9973 | — | — | — |
| `expm_frechet` | 2003 | 5.5524 / 2.5624 | 5.5524 / 2.5624 | — | — | — |
| `fractional_matrix_power` | 2027 | 3.4751e+16 / 6.8747e+12 | 3.4751e+16 / 6.8747e+12 | — | — | — |
| `funm` | 2009 | 4.5036e+15 / 1.6761e+12 | 4.5036e+15 / 1.6761e+12 | — | — | — |
| `logm` | 2009 | 4.5036e+15 / 4.9661e+12 | 4.5036e+15 / 4.9661e+12 | — | — | — |
| `signm` | 2007 | 1.8495e+16 / 1.7097e+12 | 1.8495e+16 / 1.7097e+12 | — | — | — |
| `sinhm` | 2004 | 29.297 / 0.6562 | 29.297 / 0.6562 | — | — | — |
| `sinm` | 2004 | 4.2546 / 0.68412 | 4.2546 / 0.68412 | — | — | — |
| `sqrtm` | 2009 | 4.5036e+15 / 4.4684e+12 | 4.5036e+15 / 4.4684e+12 | — | — | — |
| `tanhm` | 2004 | 29.297 / 0.71794 | 29.297 / 0.71794 | — | — | — |
| `tanm` | 2004 | 4.8006 / 0.63213 | 4.8006 / 0.63213 | — | — | — |

### `sci::linalg::solvers`

15 functions · 30116 cases

| Function | Cases | Sci vs reference | SciPy | GSL | NumPy | Boost.Numeric.uBLAS |
| --- | ---: | --- | --- | --- | --- | --- |
| `cho_solve` | 2004 | 10.688 / 1.0401 | 10.688 / 1.0401 | 0.0 / 0.0 | — | — |
| `cho_solve_banded` | 2004 | 10.688 / 0.7667 | 10.688 / 0.7667 | — | — | — |
| `lstsq` | 2012 | 4.5036e+15 / 1.1172e+12 | 4.5036e+15 / 1.1172e+12 | — | — | — |
| `lu_solve` | 2006 | 1.6842 / 1.8257e-3 | 1.6842 / 1.8257e-3 | — | — | — |
| `solve` | 2058 | 4.6362e+8 / 2.6474e+5 | 7.4052e+7 / 47552.0 | 4.1439e+8 / 2.3323e+5 | 7.835e+8 / 4.4136e+5 | 4.6366e+8 / 2.6476e+5 |
| `solve_banded` | 2003 | 0.66733 / 0.069514 | 0.66733 / 0.069514 | — | — | — |
| `solve_circulant` | 2002 | 2.625 / 1.3099e-3 | 2.625 / 1.3099e-3 | — | — | — |
| `solve_continuous_are` | 2003 | 5.6569 / 0.55421 | 5.6569 / 0.55421 | — | — | — |
| `solve_continuous_lyapunov` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `solve_discrete_are` | 2003 | 4.5036e+15 / 4.4834e+12 | 4.5036e+15 / 4.4834e+12 | — | — | — |
| `solve_discrete_lyapunov` | 2002 | 1.8601 / 0.20185 | 1.8601 / 0.20185 | — | — | — |
| `solve_sylvester` | 2003 | 2.0 / 0.15424 | 2.0 / 0.15424 | — | — | — |
| `solve_toeplitz` | 2002 | 3.4997 / 0.42147 | 3.4997 / 0.42147 | — | — | — |
| `solve_triangular` | 2008 | 1.3321 / 0.031869 | 1.3321 / 0.031869 | — | — | — |
| `solveh_banded` | 2004 | 8.4375 / 0.77894 | 8.4375 / 0.77894 | — | — | — |

### `sci::linalg::structured`

18 functions · 36085 cases

| Function | Cases | Sci vs reference | SciPy | GSL | NumPy | Boost.Numeric.uBLAS |
| --- | ---: | --- | --- | --- | --- | --- |
| `block_diag` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `circulant` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `companion` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `convolution_matrix` | 2009 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `dft` | 2009 | 9.2638e+17 / 3.085e+15 | 9.2638e+17 / 3.085e+15 | — | — | — |
| `diagsvd` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `fiedler` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `fiedler_companion` | 2006 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `hadamard` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `hankel` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `helmert` | 2006 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `hilbert` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `invhilbert` | 2004 | 0.95272 / 1.2437e-3 | 0.95272 / 1.2437e-3 | — | — | — |
| `invpascal` | 2009 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `leslie` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `matmul_toeplitz` | 2002 | 1.6015 / 0.37223 | 1.6015 / 0.37223 | — | — | — |
| `pascal` | 2009 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `toeplitz` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |

### `sci::linalg::subspaces`

5 functions · 10043 cases

| Function | Cases | Sci vs reference | SciPy | GSL | NumPy | Boost.Numeric.uBLAS |
| --- | ---: | --- | --- | --- | --- | --- |
| `clarkson_woodruff_transform` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `null_space` | 2017 | 7.5 / 1.2374 | 7.5 / 1.2374 | — | — | — |
| `orth` | 2017 | 5.1857e+16 / 1.0322e+13 | 5.1857e+16 / 1.0322e+13 | — | — | — |
| `orthogonal_procrustes` | 2002 | 4.5036e+15 / 1.1116e+14 | 4.5036e+15 / 1.1116e+14 | — | — | — |
| `subspace_angles` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |

### `sci::linalg`

4 functions · 8013 cases

| Function | Cases | Sci vs reference | SciPy | GSL | NumPy | Boost.Numeric.uBLAS |
| --- | ---: | --- | --- | --- | --- | --- |
| `find_best_blas_type` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `get_blas_funcs` | 2005 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |
| `get_lapack_funcs` | 2004 | 3.9997 / 0.57605 | 3.9997 / 0.57605 | — | — | — |
| `solve_lyapunov` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — |


## Conventions and limitations

Nonunique decompositions are compared through ordered eigenvalues, projectors and reconstructions, retaining independent structural invariants from the Sci fixture dispatcher. Interpolative decomposition compares spectral residual and rank; its SVD compares reconstruction and singular values. Type metadata uses the Sci f64 spelling. Dense sample sizes and constructor orders stay bounded. Boost.Numeric.uBLAS provides public matrix/vector arithmetic and LU solves; GSL provides BLAS, LU, Cholesky, symmetric eigenvalue and SVD primitives. Output order is adapted to the saved Sci contracts. Ill-conditioned solves can disagree with SciPy rounding while still being closer to independent high-precision mathematical references.

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
