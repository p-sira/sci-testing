# interpolate: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-interpolate)

Measured binary64 (`f64`) results compared against the **saved per-case reference**.
Each implementation, including Sci-interpolate, is evaluated against the same
reference values stored in the raw CSVs. This page
covers 50224 cases and 25 public API names;
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

### `sci::interpolate::fitting`

9 functions · 18041 cases

| Function | Cases | Sci | SciPy | Boost.Math | GSL |
| --- | ---: | --- | --- | --- | --- |
| `bisplrep` | 2002 | 7.6198e+16 / 4.876e+12 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `generate_knots` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `make_interp_spline` | 2005 | 10.495 / 0.12661 | 0.0 / 0.0 | —; unavailable 2005/2005 | —; unavailable 2005/2005 |
| `make_lsq_spline` | 2002 | 151.94 / 3.2348 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `make_smoothing_spline` | 2004 | 4.5036e+15 / 1.0215e+11 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `make_splprep` | 2007 | 2.674e+16 / 3.004e+14 | 0.0 / 0.0 | —; unavailable 2007/2007 | —; unavailable 2007/2007 |
| `make_splrep` | 2005 | 4.5036e+15 / 1.2476e+11 | 0.0 / 0.0 | —; unavailable 2005/2005 | —; unavailable 2005/2005 |
| `splprep` | 2007 | 3.1662e+19 / 1.7246e+17 | 0.0 / 0.0 | —; unavailable 2007/2007 | —; unavailable 2007/2007 |
| `splrep` | 2005 | 1.993e+16 / 5.5211e+11 | 0.0 / 0.0 | —; unavailable 2005/2005 | —; unavailable 2005/2005 |

### `sci::interpolate::multidimensional`

3 functions · 6024 cases

| Function | Cases | Sci | SciPy | Boost.Math | GSL |
| --- | ---: | --- | --- | --- | --- |
| `bisplev` | 2003 | 3.125 / 0.049255 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `griddata` | 2007 | 11.871 / 0.4387 | 0.0 / 0.0 | —; unavailable 2007/2007 | —; unavailable 2007/2007 |
| `interpn` | 2014 | 0.57199 / 0.024425 | 0.0 / 0.0 | —; unavailable 2014/2014 | 0.57199 / 0.024425 |

### `sci::interpolate::polynomial`

6 functions · 12104 cases

| Function | Cases | Sci | SciPy | Boost.Math | GSL |
| --- | ---: | --- | --- | --- | --- |
| `approximate_taylor_polynomial` | 2004 | 2.7762e+18 / 4.9684e+14 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `barycentric_interpolate` | 2050 | 23392.0 / 19.974 | 24623.0 / 15.962 | 15585.0 / 10.65 | 434.4 / 2.6421 |
| `krogh_interpolate` | 2050 | 414.66 / 2.7897 | 764.79 / 3.5273 | 2941.1 / 3.16 | 414.66 / 2.7897 |
| `lagrange` | 2000 | 1.3104e+6 / 338.44 | 0.0 / 0.0 | —; unavailable 2000/2000 | —; unavailable 2000/2000 |
| `pade` | 2000 | 0.66574 / 0.051564 | 0.0 / 0.0 | —; unavailable 2000/2000 | —; unavailable 2000/2000 |
| `pchip_interpolate` | 2000 | 1.7699 / 0.33489 | 0.0 / 0.0 | 2.059 / 0.3513 | —; unavailable 2000/2000 |

### `sci::interpolate::spline`

7 functions · 14055 cases

| Function | Cases | Sci | SciPy | Boost.Math | GSL |
| --- | ---: | --- | --- | --- | --- |
| `insert` | 2007 | 43339.0 / 6.3195 | 0.0 / 0.0 | —; unavailable 2007/2007 | —; unavailable 2007/2007 |
| `spalde` | 2003 | 3.1525e+16 / 9.2356e+11 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `splantider` | 2010 | 3.8791 / 0.022945 | 0.0 / 0.0 | —; unavailable 2010/2010 | —; unavailable 2010/2010 |
| `splder` | 2007 | 3.8947 / 0.021381 | 0.0 / 0.0 | —; unavailable 2007/2007 | —; unavailable 2007/2007 |
| `splev` | 2011 | 5.0 / 0.11639 | 0.0 / 0.0 | —; unavailable 2011/2011 | —; unavailable 2011/2011 |
| `splint` | 2013 | 1.0508e+16 / 5.2203e+12 | 0.0 / 0.0 | —; unavailable 2013/2013 | —; unavailable 2013/2013 |
| `sproot` | 2004 | 12.5 / 2.601 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |


## Conventions and limitations

Original barycentric/Krogh recipes sample an exactly represented quartic. Added fitting and spline recipes compare evaluation, knots, coefficients and degree in the established fixture order. Queries, boundary modes, repeated knots and work sizes remain explicit in the CSV inputs. Equivalent spline representations can differ at roundoff scale. Public Boost.Math PCHIP and barycentric polynomial interpolators provide independent bounded real-data comparisons. SciPy-compatible endpoint slopes must be supplied to Boost PCHIP; nontrivial extrapolation and repeated-node Krogh/Hermite derivatives are not automatically equivalent to the saved polynomial interpolation profile. GSL polynomial interpolation and public interval search provide independent saved Barycentric/Krogh and limited tensor-grid comparisons; general Sci PCHIP slope rules are not identical to GSL Steffen splines.

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
