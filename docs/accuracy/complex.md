# complex: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-complex)

Measured binary64 (`f64`) results compared against the **saved per-case reference**.
Each implementation, including Sci-complex, is evaluated against the same
reference values stored in the raw CSVs. This page
covers 70919 cases and 34 public API names;
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

### `sci::complex`

34 functions · 70919 cases

| Function | Cases | Sci | mpmath | NumPy | Boost.Multiprecision |
| --- | ---: | --- | --- | --- | --- |
| `abs` | 2007 | 0.44813 / 0.2228 | 0.44735 / 0.2228 | 0.66566 / 0.27807 | 0.44735 / 0.2228 |
| `acos` | 2282 | 1.6697 / 0.34281 | 4.5036e+15 / 1.9735e+13 | 1.1686 / 0.29718 | 4.5036e+15 / 2.5325e+13; unavailable 148/2282 |
| `acosh` | 2282 | 1.7985 / 0.34303 | 9.0072e+15 / 5.1312e+13 | 1.1686 / 0.29555 | 4.5036e+15 / 1.6883e+13; unavailable 148/2282 |
| `add` | 2004 | 0.5 / 0.063365 | 0.5 / 0.063365 | 0.5 / 0.063365 | 0.5 / 0.063397; unavailable 1/2004 |
| `arg` | 2008 | 0.3091 / 0.22019 | 0.3091 / 0.22019 | 0.3091 / 0.22019 | 0.3091 / 0.22021; unavailable 1/2008 |
| `asin` | 2282 | 1.6697 / 0.38443 | 4.5036e+15 / 1.9735e+13 | 1.4661 / 0.3475 | 4.5036e+15 / 1.6883e+13; unavailable 148/2282 |
| `asinh` | 2278 | 1.5756 / 0.35362 | 0.46106 / 0.13699 | 1.8722 / 0.37143 | 0.36427 / 0.13624; unavailable 146/2278 |
| `atan` | 2274 | 1.2918 / 0.27359 | 9.2244e+219 / 3.7783e+217 | 1.2381 / 0.23881 | 4.5036e+15 / 2.1324e+13; unavailable 162/2274 |
| `atanh` | 2278 | 1.3523 / 0.29791 | 4.5036e+15 / 4.7448e+13 | 1.2025 / 0.28017 | 0.49977 / 0.2064; unavailable 144/2278 |
| `cbrt` | 2001 | 1.3824 / 0.34354 | 0.49256 / 0.20294 | 1.227 / 0.30502 | 0.49256 / 0.20294 |
| `conjugate` | 2006 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 |
| `cos` | 2042 | 0.98771 / 0.29688 | 0.44788 / 0.13292 | 0.98771 / 0.29522 | 0.44788 / 0.13292 |
| `cosh` | 2042 | 1.1471 / 0.32963 | 0.48446 / 0.23332 | 1.1471 / 0.32796 | 0.48446 / 0.23332 |
| `divide` | 2007 | 2.5321 / 0.4649 | 0.39729 / 0.15729 | 1.7133 / 0.3959 | 0.31252 / 0.15721; unavailable 3/2007 |
| `exp` | 2074 | 1.1771 / 0.29136 | 0.493 / 0.18277 | 1.1771 / 0.29136 | 0.493 / 0.18277 |
| `from_polar` | 2004 | 0.21543 / 0.061282 | 0.21543 / 0.061282 | 0.21543 / 0.061282 | 0.21543 / 0.061282 |
| `imag` | 2006 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 |
| `ln` | 2033 | 1.9907 / 0.43268 | 0.43016 / 0.21731 | 1.1688 / 0.32069 | 4.5036e+15 / 2.2273e+12; unavailable 11/2033 |
| `multiply` | 2047 | 0.36364 / 0.065703 | 0.36364 / 0.065703 | 4.5036e+15 / 1.5569e+13 | 0.36364 / 0.065569; unavailable 20/2047 |
| `negate` | 2006 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 |
| `norm_sqr` | 2054 | 0.79642 / 0.25099 | 0.40019 / 0.19405 | 0.4974 / 0.19807 | 0.40019 / 0.19353; unavailable 12/2054 |
| `powc` | 2094 | 4.5036e+15 / 4.3014e+12 | 0.49978 / 0.18067 | 4.5036e+15 / 3.2261e+12 | 9.7612e+64 / 2.4161e+61; unavailable 74/2094 |
| `powf` | 2094 | 4.5036e+15 / 4.3014e+12 | 0.49999 / 0.17376 | 4.5036e+15 / 5.3768e+12 | 9.7612e+64 / 2.4161e+61; unavailable 74/2094 |
| `powi` | 2013 | 0.875 / 0.21001 | 0.49906 / 0.1755 | 4.5036e+15 / 1.1186e+12 | 4.5036e+15 / 4.479e+12; unavailable 2/2013 |
| `powu` | 2001 | 3.0863 / 0.50822 | 0.4996 / 0.18713 | 2.9964 / 0.56174 | 0.4996 / 0.18713 |
| `real` | 2006 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 |
| `reciprocal` | 2004 | 1.0 / 0.25171 | 0.31279 / 0.15772 | 0.71136 / 0.24812 | 0.31279 / 0.15772 |
| `sin` | 2042 | 0.97471 / 0.28871 | 0.43772 / 0.13328 | 0.97471 / 0.28704 | 0.43772 / 0.13328 |
| `sinh` | 2042 | 1.148 / 0.3403 | 0.49523 / 0.23719 | 1.148 / 0.33864 | 0.49523 / 0.23719 |
| `sqrt` | 2072 | 1.513 / 0.36094 | 0.49297 / 0.18005 | 1.1326 / 0.24025 | 0.49297 / 0.18012; unavailable 5/2072 |
| `subtract` | 2004 | 0.33336 / 0.040504 | 0.33336 / 0.040504 | 0.33336 / 0.040504 | 0.33336 / 0.040524; unavailable 1/2004 |
| `tan` | 2262 | 2.5351 / 0.44255 | 0.49583 / 0.20709 | 2.121 / 0.41594 | 0.49583 / 0.20782; unavailable 28/2262 |
| `tanh` | 2262 | 1.7828 / 0.37789 | 0.4437 / 0.13551 | 2.0017 / 0.36596 | 0.4437 / 0.13535; unavailable 28/2262 |
| `to_polar` | 2006 | 0.44807 / 0.22314 | 0.44735 / 0.22314 | 0.66581 / 0.25156 | 0.44735 / 0.22314 |


## Conventions and limitations

Complex outputs flatten real then imaginary components. Functions use principal branches; powers preserve their explicit real, complex or integer exponent conventions. All 34 complex APIs use high-precision mpmath mathematical references, with explicit signed-zero branch banks. NumPy provides the public binary64 complex operations; SciPy, Boost.Math and GSL do not expose equivalent general-purpose complex arithmetic. Signed-zero and nonfinite contracts remain in the Sci API tests. Boost.Multiprecision cpp_complex_50 is a separate public Boost complex arithmetic backend, not Boost.Math and not std::complex relabeled as Boost. Its transcendental and arithmetic comparisons are deliberately bounded around signed-zero banks, branch cuts and overflow; unavailable samples are not silently substituted with Sci or NumPy values.

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
