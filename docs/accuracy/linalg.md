# linalg: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-linalg)

Measured binary64 (`f64`) results compared against the **saved per-case reference**.
Each implementation, including Sci-linalg, is evaluated against the same
reference values stored in the raw CSVs. This page
covers 212669 cases and 106 public API names;
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

### `sci::linalg::basic`

13 functions · 26133 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy | Boost.Numeric.uBLAS | nalgebra |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `bandwidth` | 2003 | 0.0 / 0.0 | 0.0 / 0.0; unavailable 1/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `det` | 2058 | 2.5178e+8 / 1.2242e+5 | 0.4742 / 0.16594 | 1.8607e+8 / 90530.0 | 2.2527e+7 / 10999.0 | 1.8607e+8 / 90533.0 | 24148.0 / 11.908 | 2.2527e+7 / 10999.0 |
| `dot` | 2000 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 |
| `inv` | 2003 | 1.0 / 1.473e-3 | 2.5 / 4.4191e-3 | 0.0 / 0.0 | 0.0 / 0.0 | 1.0 / 4.9652e-4 | —; unavailable 2003/2003 | 2.5 / 4.4191e-3 |
| `ishermitian` | 2004 | 0.0 / 0.0 | 0.0 / 0.0; unavailable 1/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `issymmetric` | 2004 | 0.0 / 0.0 | 0.0 / 0.0; unavailable 1/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `khatri_rao` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `matmul` | 2000 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 |
| `matvec` | 2000 | 1.0 / 0.077224 | 0.0 / 0.0 | 1.0 / 0.077224 | 1.0 / 0.077224 | 0.0 / 0.0 | 1.0 / 0.077224 | 1.0 / 0.077224 |
| `norm` | 2024 | 1.4368 / 7.099e-4 | 0.0 / 0.0; unavailable 18/2024 | 0.0 / 0.0 | —; unavailable 2024/2024 | 0.8008 / 0.17594; unavailable 19/2024 | —; unavailable 2024/2024 | —; unavailable 2024/2024 |
| `pinv` | 2017 | 13.674 / 2.3156 | —; unavailable 2017/2017 | 0.0 / 0.0 | —; unavailable 2017/2017 | 6.5651 / 0.3707 | —; unavailable 2017/2017 | —; unavailable 2017/2017 |
| `pinvh` | 2018 | 53.341 / 2.1188 | —; unavailable 2018/2018 | 0.0 / 0.0 | —; unavailable 2018/2018 | 7.9023 / 0.44323 | —; unavailable 2018/2018 | —; unavailable 2018/2018 |
| `transpose` | 2000 | 0.0 / 0.0 | 0.0 / 0.0 | —; unavailable 2000/2000 | —; unavailable 2000/2000 | 0.0 / 0.0 | 0.0 / 0.0 | —; unavailable 2000/2000 |

### `sci::linalg::decompositions`

20 functions · 40105 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy | Boost.Numeric.uBLAS | nalgebra |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `cho_factor` | 2004 | 1.0005 / 0.38293 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `cholesky` | 2004 | 1.0005 / 0.38417 | 0.0 / 0.0 | 0.0 / 0.0 | 1.0005 / 0.38417 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 |
| `cholesky_banded` | 2004 | 1.3329 / 0.3594 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `cossin` | 2010 | 6.0957 / 0.074127 | —; unavailable 2010/2010 | 0.0 / 0.0 | —; unavailable 2010/2010 | —; unavailable 2010/2010 | —; unavailable 2010/2010 | —; unavailable 2010/2010 |
| `hessenberg` | 2005 | 4.6022e+17 / 1.9136e+15 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | —; unavailable 2005/2005 | —; unavailable 2005/2005 | —; unavailable 2005/2005 |
| `ldl` | 2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `lu` | 2014 | 0.75 / 1.1964e-4 | —; unavailable 2014/2014 | 0.0 / 0.0 | —; unavailable 2014/2014 | —; unavailable 2014/2014 | —; unavailable 2014/2014 | —; unavailable 2014/2014 |
| `lu_factor` | 2007 | 0.75 / 1.8437e-4 | —; unavailable 2007/2007 | 0.0 / 0.0 | —; unavailable 2007/2007 | —; unavailable 2007/2007 | —; unavailable 2007/2007 | —; unavailable 2007/2007 |
| `matrix_balance` | 2006 | 5.2355 / 0.17031 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `polar` | 2006 | 3.1481e+18 / 6.0773e+15 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `qr` | 2005 | 1.0845e+17 / 8.6269e+15 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | —; unavailable 2005/2005 |
| `qr_delete` | 2003 | 1.5 / 2.7053e-4 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `qr_insert` | 2003 | 1.0852e+19 / 6.7709e+14 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `qr_multiply` | 2004 | 27.123 / 4.3737 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `qr_update` | 2002 | 10.019 / 1.1203 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `rq` | 2003 | 1.1699e+16 / 3.873e+12 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `rsf2csf` | 2002 | 4.5036e+15 / 6.7993e+14 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `schur` | 2005 | 9.694e+17 / 2.3762e+15 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | —; unavailable 2005/2005 | —; unavailable 2005/2005 | —; unavailable 2005/2005 |
| `svd` | 2005 | 1.0525e+18 / 1.5578e+16 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | —; unavailable 2005/2005 |
| `svdvals` | 2005 | 3.6493 / 0.8805 | 2.7382 / 0.40666 | 0.0 / 0.0 | 2.7394 / 0.53977 | 0.0 / 0.0 | —; unavailable 2005/2005 | 3.6541 / 0.44829 |

### `sci::linalg::eigen`

11 functions · 22046 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy | Boost.Numeric.uBLAS | nalgebra |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `cdf2rdf` | 2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `eig` | 2004 | 2.0 / 0.18063 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `eig_banded` | 2004 | 6.75 / 0.041538 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `eigh` | 2006 | 37.787 / 8.3398 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `eigh_tridiagonal` | 2002 | 8.3971 / 0.011478 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `eigvals` | 2004 | 12.603 / 6.7129e-3 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `eigvals_banded` | 2004 | 0.8453 / 7.541e-4 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `eigvalsh` | 2006 | 3.3844 / 0.93313 | 2.0019 / 0.40171 | 0.0 / 0.0 | 5.5181 / 1.1351 | 0.0 / 0.0 | —; unavailable 2006/2006 | 3.1555 / 0.86814 |
| `eigvalsh_tridiagonal` | 2002 | 2.0 / 2.6287e-3 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `ordqz` | 2008 | 1.0273e+17 / 2.7242e+13 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `qz` | 2004 | 9.2181e+16 / 1.4036e+13 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |

### `sci::linalg::interpolative`

9 functions · 18042 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy | Boost.Numeric.uBLAS | nalgebra |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `interpolative.estimate_rank` | 2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `interpolative.estimate_spectral_norm` | 2004 | 4.7852 / 1.4 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `interpolative.estimate_spectral_norm_diff` | 2004 | 5.0648 / 1.0155 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `interpolative.id_to_svd` | 2003 | 7.0 / 0.94978 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `interpolative.interp_decomp` | 2008 | 4.7852 / 0.67558 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `interpolative.reconstruct_interp_matrix` | 2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `interpolative.reconstruct_matrix_from_id` | 2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `interpolative.reconstruct_skel_matrix` | 2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `interpolative.svd` | 2008 | 6.5 / 5.4349e-3 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |

### `sci::linalg::matrix_functions`

14 functions · 28097 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy | Boost.Numeric.uBLAS | nalgebra |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `coshm` | 2004 | 11.805 / 0.83849 | 11.253 / 0.80384 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `cosm` | 2004 | 7.4006 / 0.82619 | 5.8352 / 0.32745 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `expm` | 2005 | 4286.1 / 8.4365 | 4295.0 / 7.9752 | 0.0 / 0.0 | —; unavailable 2005/2005 | —; unavailable 2005/2005 | —; unavailable 2005/2005 | —; unavailable 2005/2005 |
| `expm_cond` | 2004 | 21.521 / 3.9973 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `expm_frechet` | 2003 | 5.5524 / 2.5624 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `fractional_matrix_power` | 2027 | 3.4751e+16 / 6.8747e+12 | 8.5718e+92 / 9.978e+88 | 0.0 / 0.0 | —; unavailable 2027/2027 | —; unavailable 2027/2027 | —; unavailable 2027/2027 | —; unavailable 2027/2027 |
| `funm` | 2009 | 4.5036e+15 / 1.6761e+12 | 2.3768 / 0.046284; unavailable 3/2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | —; unavailable 2009/2009 | —; unavailable 2009/2009 | —; unavailable 2009/2009 |
| `logm` | 2009 | 4.5036e+15 / 4.9661e+12 | 1.5388e+92 / 9.5425e+87 | 0.0 / 0.0 | —; unavailable 2009/2009 | —; unavailable 2009/2009 | —; unavailable 2009/2009 | —; unavailable 2009/2009 |
| `signm` | 2007 | 1.8495e+16 / 1.7097e+12 | —; unavailable 2007/2007 | 0.0 / 0.0 | —; unavailable 2007/2007 | —; unavailable 2007/2007 | —; unavailable 2007/2007 | —; unavailable 2007/2007 |
| `sinhm` | 2004 | 29.297 / 0.6562 | 28.32 / 0.46695 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `sinm` | 2004 | 4.2546 / 0.68412 | 2.5524 / 0.38377 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `sqrtm` | 2009 | 4.5036e+15 / 4.4684e+12 | 9.0072e+15 / 5.0295e+12; unavailable 1/2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | —; unavailable 2009/2009 | —; unavailable 2009/2009 | —; unavailable 2009/2009 |
| `tanhm` | 2004 | 29.297 / 0.71794 | 28.32 / 0.50757 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `tanm` | 2004 | 4.8006 / 0.63213 | 3.1977 / 0.47259 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |

### `sci::linalg::solvers`

16 functions · 32118 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy | Boost.Numeric.uBLAS | nalgebra |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `cho_solve` | 2004 | 10.688 / 1.0401 | 1.9999 / 0.51843 | 0.0 / 0.0 | 10.688 / 1.0401 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | 1.0 / 0.16631 |
| `cho_solve_banded` | 2004 | 10.688 / 0.7667 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `lstsq` | 2012 | 4.5036e+15 / 1.1172e+12 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | 3.0 / 7.6912e-3 | —; unavailable 2012/2012 | —; unavailable 2012/2012 |
| `lu_solve` | 2006 | 1.6842 / 1.8257e-3 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `solve` | 2058 | 4.6362e+8 / 2.6474e+5 | 0.5 / 0.17901 | 3.8957e+8 / 2.1752e+5 | 4.9228e+7 / 31689.0 | 3.1988e+8 / 1.7663e+5 | 37637.0 / 22.162 | 4.9228e+7 / 31689.0 |
| `solve_banded` | 2003 | 0.66733 / 0.069514 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `solve_circulant` | 2002 | 2.625 / 1.3099e-3 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `solve_continuous_are` | 2003 | 5.6569 / 0.55421 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `solve_continuous_lyapunov` | 2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `solve_discrete_are` | 2003 | 4.5036e+15 / 4.4834e+12 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `solve_discrete_lyapunov` | 2002 | 1.8601 / 0.20185 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `solve_lyapunov` | 2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `solve_sylvester` | 2003 | 2.0 / 0.15424 | 1.0 / 0.19832 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `solve_toeplitz` | 2002 | 3.4997 / 0.42147 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `solve_triangular` | 2008 | 1.3321 / 0.031869 | 1.0358e+17 / 6.2914e+13 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `solveh_banded` | 2004 | 8.4375 / 0.77894 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |

### `sci::linalg::structured`

18 functions · 36085 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy | Boost.Numeric.uBLAS | nalgebra |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `block_diag` | 2003 | 0.0 / 0.0 | 0.0 / 0.0; unavailable 656/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `circulant` | 2003 | 0.0 / 0.0 | 0.0 / 0.0; unavailable 1/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `companion` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `convolution_matrix` | 2009 | 0.0 / 0.0 | 0.0 / 0.0; failed 1/2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | —; unavailable 2009/2009 | —; unavailable 2009/2009 | —; unavailable 2009/2009 |
| `dft` | 2009 | 9.2638e+17 / 3.085e+15 | 1.1786e+93 / 1.1624e+91 | 0.0 / 0.0 | —; unavailable 2009/2009 | —; unavailable 2009/2009 | —; unavailable 2009/2009 | —; unavailable 2009/2009 |
| `diagsvd` | 2003 | 0.0 / 0.0 | 0.0 / 0.0; unavailable 1/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `fiedler` | 2003 | 0.0 / 0.0 | 0.0 / 0.0; unavailable 1/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `fiedler_companion` | 2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `hadamard` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `hankel` | 2003 | 0.0 / 0.0 | 0.0 / 0.0; unavailable 1/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `helmert` | 2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `hilbert` | 2004 | 0.0 / 0.0 | 0.0 / 0.0; unavailable 31/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `invhilbert` | 2004 | 0.95272 / 1.2437e-3 | 0.95272 / 6.5441e-3; unavailable 28/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `invpascal` | 2009 | 0.0 / 0.0 | 0.0 / 0.0; unavailable 45/2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | —; unavailable 2009/2009 | —; unavailable 2009/2009 | —; unavailable 2009/2009 |
| `leslie` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `matmul_toeplitz` | 2002 | 1.6015 / 0.37223 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `pascal` | 2009 | 0.0 / 0.0 | 0.0 / 0.0; unavailable 48/2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | —; unavailable 2009/2009 | —; unavailable 2009/2009 | —; unavailable 2009/2009 |
| `toeplitz` | 2003 | 0.0 / 0.0 | 0.0 / 0.0; unavailable 1/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |

### `sci::linalg::subspaces`

5 functions · 10043 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy | Boost.Numeric.uBLAS | nalgebra |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `clarkson_woodruff_transform` | 2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `null_space` | 2017 | 7.5 / 1.2374 | —; unavailable 2017/2017 | 0.0 / 0.0 | —; unavailable 2017/2017 | —; unavailable 2017/2017 | —; unavailable 2017/2017 | —; unavailable 2017/2017 |
| `orth` | 2017 | 5.1857e+16 / 1.0322e+13 | —; unavailable 2017/2017 | 0.0 / 0.0 | —; unavailable 2017/2017 | —; unavailable 2017/2017 | —; unavailable 2017/2017 | —; unavailable 2017/2017 |
| `orthogonal_procrustes` | 2002 | 4.5036e+15 / 1.1116e+14 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `subspace_angles` | 2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |

## Evaluation failures

- `linalg.convolution_matrix` (mpmath): 1 failed; shape mismatch: expected 6 scalar outputs, received 4. Example: `linalg.convolution_matrix/recipe-0066`.


## Conventions and limitations

Nonunique decompositions are compared through ordered eigenvalues, projectors and reconstructions, retaining independent structural invariants from the Sci fixture dispatcher. Interpolative decomposition compares spectral residual and rank; its SVD compares reconstruction and singular values. Type metadata uses the Sci f64 spelling. Dense sample sizes and constructor orders stay bounded. The nalgebra comparisons call a separate pinned Rust numerical crate and preserve matrix output shape and row-major flattening. High-precision mpmath computations offer additional bounded matrix comparisons with distinct rounding behavior from binary64 BLAS/LAPACK solvers. Boost.Numeric.uBLAS provides public matrix/vector arithmetic and LU solves; GSL provides BLAS, LU, Cholesky, symmetric eigenvalue and SVD primitives. Output order is adapted to the saved Sci contracts. Ill-conditioned solves can disagree with SciPy rounding while still being closer to independent high-precision mathematical references.

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
