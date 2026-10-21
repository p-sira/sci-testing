# fft: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-fft)

Measured binary64 (`f64`) results using **Sci-fft as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 90336 cases and 45 public API names;
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
| `cc_diff` | 2004 | 8.7137 / 0.015853 | — |
| `cs_diff` | 2004 | 7.3222 / 0.040276 | — |
| `dct` | 2008 | 1.3511e+16 / 3.9133e+12 | — |
| `dctn` | 2008 | 4.5036e+15 / 6.5528e+14 | — |
| `diff` | 2012 | 3.1833e+16 / 3.6827e+13 | — |
| `dst` | 2008 | 4.954e+16 / 7.2675e+12 | — |
| `dstn` | 2008 | 1.8014e+16 / 8.1393e+14 | — |
| `fft` | 2040 | 2115.0 / 1.8359 | 2115.0 / 1.8359 |
| `fft2` | 2004 | 4.0532e+16 / 3.708e+15 | — |
| `fftfreq` | 2005 | 0.0 / 0.0 | — |
| `fftn` | 2005 | 1.0085e+17 / 3.6571e+15 | — |
| `fftshift` | 2005 | 0.0 / 0.0 | — |
| `fht` | 2006 | 619.05 / 0.53128 | — |
| `fhtoffset` | 2006 | 0.0 / 0.0 | — |
| `hfft` | 2005 | 2.435 / 8.2451e-4 | — |
| `hfft2` | 2004 | 6.7214 / 0.50543 | — |
| `hfftn` | 2004 | 18.192 / 0.50645 | — |
| `hilbert` | 2004 | 28.745 / 0.04542 | — |
| `idct` | 2008 | 1.5763e+16 / 4.2627e+12 | — |
| `idctn` | 2008 | 1.8014e+16 / 6.5444e+14 | — |
| `idst` | 2008 | 4.954e+16 / 7.2675e+12 | — |
| `idstn` | 2008 | 6.5869e+30 / 2.7282e+26 | — |
| `ifft` | 2040 | 2115.0 / 1.9261 | 2115.0 / 1.9261 |
| `ifft2` | 2004 | 5.2326e+16 / 3.4286e+15 | — |
| `ifftn` | 2005 | 7.676e+16 / 3.4467e+15 | — |
| `ifftshift` | 2005 | 0.0 / 0.0 | — |
| `ifht` | 2006 | 612.77 / 0.27612 | — |
| `ihfft` | 2007 | 3.0 / 4.8078e-3 | — |
| `ihfft2` | 2004 | 2.8523e+16 / 2.0101e+15 | — |
| `ihfftn` | 2004 | 2.8523e+16 / 2.0549e+15 | — |
| `ihilbert` | 2004 | 28.745 / 0.04542 | — |
| `irfft` | 2005 | 2.435 / 9.1757e-4 | — |
| `irfft2` | 2004 | 7.5616 / 0.56249 | — |
| `irfftn` | 2004 | 17.908 / 0.56237 | — |
| `itilbert` | 2004 | 86.898 / 0.09966 | — |
| `next_fast_len` | 2016 | 0.0 / 0.0 | — |
| `prev_fast_len` | 2016 | 0.0 / 0.0 | — |
| `rfft` | 2005 | 3.1474 / 4.3585e-3 | — |
| `rfft2` | 2004 | 4.0532e+16 / 2.5891e+15 | — |
| `rfftfreq` | 2005 | 0.0 / 0.0 | — |
| `rfftn` | 2004 | 4.0532e+16 / 2.6188e+15 | — |
| `sc_diff` | 2005 | 407.7 / 0.2497 | — |
| `shift` | 2004 | 5.4813 / 0.016628 | — |
| `ss_diff` | 2005 | 18.388 / 0.03237 | — |
| `tilbert` | 2004 | 9.6962 / 0.029734 | — |


## Conventions and limitations

Transforms use unnormalized forward and inverse-divided-by-length conventions. Real transforms, multidimensional shapes, transform types, FFTLog and periodic operators follow the pinned fixture conventions. FFTLog uses SciPy fft; distinct periodic operators use its legacy fftpack implementation. Real and imaginary components are interleaved, including exact imaginary zeros.

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
