# integrate: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-integrate)

Measured binary64 (`f64`) results using **Sci-integrate as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 14135 cases and 7 public API names;
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
| `cumulative_simpson` | 2000 | 499.5 / 0.48561 | — |
| `cumulative_trapezoid` | 2000 | 666.0 / 0.18511 | — |
| `fixed_quad` | 2068 | 41.586 / 5.142 | — |
| `newton_cotes` | 2000 | 1.1988e+7 / 8920.3 | — |
| `romb` | 2000 | 1.984 / 0.1853 | — |
| `simpson` | 2000 | 1.9521 / 0.361 | — |
| `trapezoid` | 2067 | 102.39 / 0.34064 | 102.39 / 0.34064 |


## Conventions and limitations

Discrete-rule references integrate the supplied binary64 samples, excluding continuous sampling error. The dense added rules use dx = 0.25 and zero initial cumulative integrals. Newton-Cotes records weights followed by its error coefficient. Fixed Gaussian quadrature uses order 16 on exp(x), with analytic references.

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
