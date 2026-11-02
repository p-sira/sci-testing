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

### `sci::complex` (34 functions, 70919 cases)

| Function | Cases | Sci vs reference | NumPy |
| --- | ---: | --- | --- |
| `abs` | 2007 | 0.44813 / 0.2228 | 0.89621 / 0.22006 |
| `acos` | 2282 | 1.6697 / 0.34281 | 1.837 / 0.30072 |
| `acosh` | 2282 | 1.7985 / 0.34303 | 1.837 / 0.30628 |
| `add` | 2004 | 0.5 / 0.063365 | 0.0 / 0.0 |
| `arg` | 2008 | 0.3091 / 0.22019 | 0.0 / 0.0 |
| `asin` | 2282 | 1.6697 / 0.38443 | 1.9971 / 0.39866 |
| `asinh` | 2278 | 1.5756 / 0.35362 | 2.535 / 0.45588 |
| `atan` | 2274 | 1.2918 / 0.27359 | 1.2478 / 0.22307 |
| `atanh` | 2278 | 1.3523 / 0.29791 | 1.0466 / 0.098362 |
| `cbrt` | 2001 | 1.3824 / 0.34354 | 0.9864 / 0.18385 |
| `conjugate` | 2006 | 0.0 / 0.0 | 0.0 / 0.0 |
| `cos` | 2042 | 0.98771 / 0.29688 | 0.97888 / 6.3954e-3 |
| `cosh` | 2042 | 1.1471 / 0.32963 | 0.97888 / 6.3954e-3 |
| `divide` | 2007 | 2.5321 / 0.4649 | 2.501 / 0.54263 |
| `exp` | 2074 | 1.1771 / 0.29136 | 0.0 / 0.0 |
| `from_polar` | 2004 | 0.21543 / 0.061282 | 0.0 / 0.0 |
| `imag` | 2006 | 0.0 / 0.0 | 0.0 / 0.0 |
| `ln` | 2033 | 1.9907 / 0.43268 | 2.5873 / 0.33484 |
| `multiply` | 2047 | 0.36364 / 0.065703 | 4.5036e+15 / 1.5569e+13 |
| `negate` | 2006 | 0.0 / 0.0 | 0.0 / 0.0 |
| `norm_sqr` | 2054 | 0.79642 / 0.25099 | 0.8032 / 0.16149 |
| `powc` | 2094 | 4.5036e+15 / 4.3014e+12 | 4.5036e+15 / 7.5275e+12 |
| `powf` | 2094 | 4.5036e+15 / 4.3014e+12 | 4.5036e+15 / 9.6782e+12 |
| `powi` | 2013 | 0.875 / 0.21001 | 4.5036e+15 / 1.1186e+12 |
| `powu` | 2001 | 3.0863 / 0.50822 | 3.9996 / 0.46508 |
| `real` | 2006 | 0.0 / 0.0 | 0.0 / 0.0 |
| `reciprocal` | 2004 | 1.0 / 0.25171 | 1.2517 / 0.27961 |
| `sin` | 2042 | 0.97471 / 0.28871 | 0.97888 / 6.3954e-3 |
| `sinh` | 2042 | 1.148 / 0.3403 | 0.97888 / 6.3954e-3 |
| `sqrt` | 2072 | 1.513 / 0.36094 | 1.9905 / 0.31184 |
| `subtract` | 2004 | 0.33336 / 0.040504 | 0.0 / 0.0 |
| `tan` | 2262 | 2.5351 / 0.44255 | 2.4046 / 0.38774 |
| `tanh` | 2262 | 1.7828 / 0.37789 | 2.0589 / 0.37394 |
| `to_polar` | 2006 | 0.44807 / 0.22314 | 0.89621 / 0.11276 |


## Conventions and limitations

Complex outputs flatten real then imaginary components. Functions use principal branches; powers preserve their explicit real, complex or integer exponent conventions. All 34 complex APIs use high-precision mpmath mathematical references, with explicit signed-zero branch banks. NumPy provides the public binary64 complex operations; SciPy, Boost.Math and GSL do not expose equivalent general-purpose complex arithmetic. Signed-zero and nonfinite contracts remain in the Sci API tests.

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
