# complex: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-complex)

Measured binary64 (`f64`) results using **Sci-complex as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 70919 cases and 34 public API names;
coverage is limited to the operations and inputs below.

## Agreement with Sci

Each library cell is **maximum / mean difference from Sci in ε** over cases
where both evaluations succeed. `Cases` counts input cases; `Sci failed N/T`
identifies missing baseline evaluations. `failed N/T` counts that library's own
failures; `Sci unavailable N/T` counts its successful evaluations excluded
because Sci failed. A dash means unavailable or no paired finite results.
Means include all compared scalar components, including exact zeros.
Full outputs, references, warnings and errors remain in the per-function CSVs.

| Function | Cases | NumPy | stdlib |
| --- | ---: | --- | --- |
| `abs` | 2007 | — | 0.0 / 0.0 |
| `acos` | 2282 | — | 2.6673 / 0.36619 |
| `acosh` | 2282 | — | 2.6673 / 0.36781 |
| `add` | 2004 | — | 0.0 / 0.0 |
| `arg` | 2008 | — | 0.0 / 0.0 |
| `asin` | 2282 | — | 2.6673 / 0.38283 |
| `asinh` | 2278 | — | 2.2315 / 0.40473 |
| `atan` | 2274 | — | 1.2478 / 0.13413 |
| `atanh` | 2278 | — | 1.0466 / 0.10341 |
| `cbrt` | 2001 | — | 0.9864 / 0.1801 |
| `conjugate` | 2006 | — | 0.0 / 0.0 |
| `cos` | 2042 | — | 0.97888 / 5.0407e-3 |
| `cosh` | 2042 | — | 0.97888 / 5.0407e-3 |
| `divide` | 2007 | — | 3.1237 / 0.5162 |
| `exp` | 2074 | 0.0 / 0.0 | 0.0 / 0.0 |
| `from_polar` | 2004 | — | 0.0 / 0.0 |
| `imag` | 2006 | — | 0.0 / 0.0 |
| `ln` | 2033 | — | 2.5873 / 0.33484 |
| `multiply` | 2047 | — | 3.5911e+14 / 3.5087e+11 |
| `negate` | 2006 | — | 0.0 / 0.0 |
| `norm_sqr` | 2054 | — | 0.8 / 0.07929 |
| `powc` | 2094 | — | 4.5036e+15 / 1.613e+13 |
| `powf` | 2094 | — | 4.5036e+15 / 1.8281e+13 |
| `powi` | 2013 | — | 4.5036e+15 / 1.1186e+12 |
| `powu` | 2001 | — | 3.9996 / 0.46508 |
| `real` | 2006 | — | 0.0 / 0.0 |
| `reciprocal` | 2004 | — | 1.2517 / 0.27961 |
| `sin` | 2042 | — | 0.97888 / 5.0407e-3 |
| `sinh` | 2042 | — | 0.97888 / 5.0407e-3 |
| `sqrt` | 2072 | 1.9905 / 0.31184 | 1.9905 / 0.31184 |
| `subtract` | 2004 | — | 0.0 / 0.0 |
| `tan` | 2262 | — | 3.2055 / 0.64917 |
| `tanh` | 2262 | — | 3.1531 / 0.5934 |
| `to_polar` | 2006 | — | 0.0 / 0.0 |


## Conventions and limitations

Complex outputs flatten real then imaginary components. Functions use principal branches; powers preserve their explicit real, complex or integer exponent conventions. The original exp/sqrt recipes use high-precision mpmath references; additional complex comparisons use the Python standard library as a binary64 reference. Signed-zero and nonfinite contracts remain in the Sci API tests.

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
