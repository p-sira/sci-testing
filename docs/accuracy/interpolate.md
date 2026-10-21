# interpolate: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-interpolate)

Measured binary64 (`f64`) results using **Sci-interpolate as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 50224 cases and 25 public API names;
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
| `approximate_taylor_polynomial` | 2004 | 2.7762e+18 / 4.9684e+14 |
| `barycentric_interpolate` | 2050 | 37692.0 / 28.937 |
| `bisplev` | 2003 | 3.125 / 0.049255 |
| `bisplrep` | 2002 | 7.6198e+16 / 4.876e+12 |
| `generate_knots` | 2004 | 0.0 / 0.0 |
| `griddata` | 2007 | 11.871 / 0.4387 |
| `insert` | 2007 | 43339.0 / 6.3195 |
| `interpn` | 2014 | 0.57199 / 0.024425 |
| `krogh_interpolate` | 2050 | 632.8 / 4.2342 |
| `lagrange` | 2000 | 1.3104e+6 / 338.44 |
| `make_interp_spline` | 2005 | 4.5036e+15 / 4.9901e+11 |
| `make_lsq_spline` | 2002 | 151.94 / 3.2348 |
| `make_smoothing_spline` | 2004 | 4.5036e+15 / 1.0215e+11 |
| `make_splprep` | 2007 | 2.674e+16 / 3.004e+14 |
| `make_splrep` | 2005 | 4.5036e+15 / 1.2476e+11 |
| `pade` | 2000 | 0.66574 / 0.051564 |
| `pchip_interpolate` | 2000 | 1.7699 / 0.33489 |
| `spalde` | 2003 | 3.1525e+16 / 9.2356e+11 |
| `splantider` | 2010 | 3.8791 / 0.022945 |
| `splder` | 2007 | 3.8947 / 0.021381 |
| `splev` | 2011 | 5.0 / 0.11639 |
| `splint` | 2013 | 1.0508e+16 / 5.2203e+12 |
| `splprep` | 2007 | 3.1662e+19 / 1.7246e+17 |
| `splrep` | 2005 | 1.993e+16 / 5.5211e+11 |
| `sproot` | 2004 | 12.5 / 2.601 |


## Conventions and limitations

Original barycentric/Krogh recipes sample an exactly represented quartic. Added fitting and spline recipes compare evaluation, knots, coefficients and degree in the established fixture order. Queries, boundary modes, repeated knots and work sizes remain explicit in the CSV inputs. Equivalent spline representations can differ at roundoff scale.

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
