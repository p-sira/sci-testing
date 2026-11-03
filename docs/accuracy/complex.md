# complex: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-complex)

Measured binary64 (`f64`) results using **Sci-complex as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 70919 cases and 34 public API names;
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

### `sci::complex`

34 functions · 70919 cases

| Function | Cases | Sci vs reference | mpmath | NumPy | Boost.Multiprecision |
| --- | ---: | --- | --- | --- | --- |
| `abs` | 2007 | 0.44813 / 0.2228 | 0.89471 / 8.9117e-4 | 0.89621 / 0.22006 | 0.89471 / 8.9117e-4 |
| `acos` | 2282 | 1.6697 / 0.34281 | 4.5036e+15 / 1.9735e+13 | 1.837 / 0.30072 | 4.5036e+15 / 2.5325e+13; unavailable 148/2282 |
| `acosh` | 2282 | 1.7985 / 0.34303 | 9.0072e+15 / 5.1312e+13 | 1.837 / 0.30628 | 4.5036e+15 / 1.6883e+13; unavailable 148/2282 |
| `add` | 2004 | 0.5 / 0.063365 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0; unavailable 1/2004 |
| `arg` | 2008 | 0.3091 / 0.22019 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0; unavailable 1/2008 |
| `asin` | 2282 | 1.6697 / 0.38443 | 4.5036e+15 / 1.9735e+13 | 1.9971 / 0.39866 | 4.5036e+15 / 1.6883e+13; unavailable 148/2282 |
| `asinh` | 2278 | 1.5756 / 0.35362 | 1.6737 / 0.32549 | 2.535 / 0.45588 | 1.6737 / 0.34373; unavailable 146/2278 |
| `atan` | 2274 | 1.2918 / 0.27359 | 9.2244e+219 / 3.7783e+217 | 1.2478 / 0.22307 | 4.5036e+15 / 2.1324e+13; unavailable 162/2274 |
| `atanh` | 2278 | 1.3523 / 0.29791 | 4.5036e+15 / 4.7448e+13 | 1.0466 / 0.098362 | 1.0012 / 0.2346; unavailable 144/2278 |
| `cbrt` | 2001 | 1.3824 / 0.34354 | 0.9864 / 0.29166 | 0.9864 / 0.18385 | 0.9864 / 0.29166 |
| `conjugate` | 2006 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 |
| `cos` | 2042 | 0.98771 / 0.29688 | 1.1053 / 0.27335 | 0.97888 / 6.3954e-3 | 1.1053 / 0.27335 |
| `cosh` | 2042 | 1.1471 / 0.32963 | 0.97888 / 0.25372 | 0.97888 / 6.3954e-3 | 0.97888 / 0.25372 |
| `divide` | 2007 | 2.5321 / 0.4649 | 2.5004 / 0.44196 | 2.501 / 0.54263 | 2.5004 / 0.44237; unavailable 3/2007 |
| `exp` | 2074 | 1.1771 / 0.29136 | 1.3584 / 0.23519 | 0.0 / 0.0 | 1.3584 / 0.23519 |
| `from_polar` | 2004 | 0.21543 / 0.061282 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 |
| `imag` | 2006 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 |
| `ln` | 2033 | 1.9907 / 0.43268 | 1.7253 / 0.31504 | 2.5873 / 0.33484 | 4.5036e+15 / 2.2273e+12; unavailable 11/2033 |
| `multiply` | 2047 | 0.36364 / 0.065703 | 0.0 / 0.0 | 4.5036e+15 / 1.5569e+13 | 0.0 / 0.0; unavailable 20/2047 |
| `negate` | 2006 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 |
| `norm_sqr` | 2054 | 0.79642 / 0.25099 | 0.8032 / 0.16621 | 0.8032 / 0.16149 | 0.8032 / 0.16719; unavailable 12/2054 |
| `powc` | 2094 | 4.5036e+15 / 4.3014e+12 | 4.5036e+15 / 4.3014e+12 | 4.5036e+15 / 7.5275e+12 | 9.7612e+64 / 2.4161e+61; unavailable 74/2094 |
| `powf` | 2094 | 4.5036e+15 / 4.3014e+12 | 4.5036e+15 / 4.3014e+12 | 4.5036e+15 / 9.6782e+12 | 9.7612e+64 / 2.4161e+61; unavailable 74/2094 |
| `powi` | 2013 | 0.875 / 0.21001 | 0.78125 / 0.081463 | 4.5036e+15 / 1.1186e+12 | 4.5036e+15 / 4.479e+12; unavailable 2/2013 |
| `powu` | 2001 | 3.0863 / 0.50822 | 2.9999 / 0.46339 | 3.9996 / 0.46508 | 2.9999 / 0.46339 |
| `real` | 2006 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 |
| `reciprocal` | 2004 | 1.0 / 0.25171 | 0.70833 / 0.24468 | 1.2517 / 0.27961 | 0.70833 / 0.24468 |
| `sin` | 2042 | 0.97471 / 0.28871 | 1.1292 / 0.26367 | 0.97888 / 6.3954e-3 | 1.1292 / 0.26367 |
| `sinh` | 2042 | 1.148 / 0.3403 | 0.99232 / 0.26579 | 0.97888 / 6.3954e-3 | 0.99232 / 0.26579 |
| `sqrt` | 2072 | 1.513 / 0.36094 | 1.9527 / 0.32091 | 1.9905 / 0.31184 | 1.9527 / 0.32169; unavailable 5/2072 |
| `subtract` | 2004 | 0.33336 / 0.040504 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0; unavailable 1/2004 |
| `tan` | 2262 | 2.5351 / 0.44255 | 2.4071 / 0.39983 | 2.4046 / 0.38774 | 2.4071 / 0.39978; unavailable 28/2262 |
| `tanh` | 2262 | 1.7828 / 0.37789 | 1.7757 / 0.35699 | 2.0589 / 0.37394 | 1.7214 / 0.35641; unavailable 28/2262 |
| `to_polar` | 2006 | 0.44807 / 0.22314 | 0.89471 / 4.4591e-4 | 0.89621 / 0.11276 | 0.89471 / 4.4591e-4 |


## Conventions and limitations

Complex outputs flatten real then imaginary components. Functions use principal branches; powers preserve their explicit real, complex or integer exponent conventions. All 34 complex APIs use high-precision mpmath mathematical references, with explicit signed-zero branch banks. NumPy provides the public binary64 complex operations; SciPy, Boost.Math and GSL do not expose equivalent general-purpose complex arithmetic. Signed-zero and nonfinite contracts remain in the Sci API tests. Boost.Multiprecision cpp_complex_50 is a separate public Boost complex arithmetic backend, not Boost.Math and not std::complex relabeled as Boost. Its transcendental and arithmetic comparisons are deliberately bounded around signed-zero banks, branch cuts and overflow; unavailable samples are not silently substituted with Sci or NumPy values.

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
