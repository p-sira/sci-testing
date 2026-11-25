# fft: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-fft)

Measured binary64 (`f64`) results compared against the **saved per-case reference**.
Each implementation, including Sci-fft, is evaluated against the same
reference values stored in the raw CSVs. This page
covers 90336 cases and 45 public API names;
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

### `sci::fft::fourier`

6 functions · 12098 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy | RustFFT |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `fft` | 2040 | 9.0072e+15 / 1.8055e+13 | 4.5036e+15 / 2.6019e+13 | 9.0072e+15 / 1.8055e+13 | 9.0072e+15 / 2.1958e+13 | 9.0072e+15 / 1.8055e+13 | 9.0072e+15 / 1.8055e+13 |
| `fft2` | 2004 | 4.0532e+16 / 3.708e+15 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | 4.0532e+16 / 3.708e+15 |
| `fftn` | 2005 | 1.0085e+17 / 3.6571e+15 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | 1.8711e+16 / 2.4416e+12 | 1.0085e+17 / 3.6571e+15 |
| `ifft` | 2040 | 9.0072e+15 / 1.8055e+13 | 4.5036e+15 / 2.8384e+13 | 9.0072e+15 / 1.8055e+13 | 9.0072e+15 / 2.2076e+13 | 9.0072e+15 / 1.8055e+13 | 9.0072e+15 / 1.8055e+13 |
| `ifft2` | 2004 | 5.2326e+16 / 3.4286e+15 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 2.0909e+16 / 9.6253e+11 | 5.2326e+16 / 3.4293e+15 |
| `ifftn` | 2005 | 7.676e+16 / 3.4467e+15 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | 3.0203e+16 / 2.5721e+12 | 7.676e+16 / 3.4467e+15 |

### `sci::fft::hankel`

3 functions · 6018 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy | RustFFT |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `fht` | 2006 | 619.05 / 0.53128 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `fhtoffset` | 2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `ifht` | 2006 | 612.77 / 0.27612 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |

### `sci::fft::real`

20 functions · 40118 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy | RustFFT |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `dct` | 2008 | 1.3511e+16 / 3.9133e+12 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `dctn` | 2008 | 4.5036e+15 / 6.5528e+14 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `dst` | 2008 | 4.954e+16 / 7.2675e+12 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `dstn` | 2008 | 1.8014e+16 / 8.1393e+14 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `hfft` | 2005 | 2.435 / 8.2451e-4 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 |
| `hfft2` | 2004 | 6.7214 / 0.50543 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `hfftn` | 2004 | 18.192 / 0.50645 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `idct` | 2008 | 1.5763e+16 / 4.2627e+12 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `idctn` | 2008 | 1.8014e+16 / 6.5444e+14 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `idst` | 2008 | 4.954e+16 / 7.2675e+12 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `idstn` | 2008 | 6.5869e+30 / 2.7282e+26 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `ihfft` | 2007 | 3.0 / 4.8078e-3 | —; unavailable 2007/2007 | 0.0 / 0.0 | —; unavailable 2007/2007 | 0.0 / 0.0 | —; unavailable 2007/2007 |
| `ihfft2` | 2004 | 2.8523e+16 / 2.0101e+15 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `ihfftn` | 2004 | 2.8523e+16 / 2.0549e+15 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `irfft` | 2005 | 2.435 / 9.1757e-4 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | 0.0 / 0.0 | 2.435 / 7.6906e-4 |
| `irfft2` | 2004 | 7.5616 / 0.56249 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 7.2188 / 3.0006e-4 | 4.5036e+15 / 3.744e+11 |
| `irfftn` | 2004 | 17.908 / 0.56237 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | 17.908 / 0.50638 |
| `rfft` | 2005 | 3.1474 / 4.3585e-3 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | 0.0 / 0.0 | 3.1474 / 4.1725e-3 |
| `rfft2` | 2004 | 4.0532e+16 / 2.5891e+15 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | 4.0532e+16 / 2.5891e+15 |
| `rfftn` | 2004 | 4.0532e+16 / 2.6188e+15 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | 4.0532e+16 / 2.6188e+15 |

### `sci::fft::utilities`

16 functions · 32102 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy | RustFFT |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `cc_diff` | 2004 | 8.7137 / 0.015853 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `cs_diff` | 2004 | 7.3222 / 0.040276 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `diff` | 2012 | 3.1833e+16 / 3.6827e+13 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | —; unavailable 2012/2012 | —; unavailable 2012/2012 |
| `fftfreq` | 2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 |
| `fftshift` | 2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 |
| `hilbert` | 2004 | 28.745 / 0.04542 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `ifftshift` | 2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 |
| `ihilbert` | 2004 | 28.745 / 0.04542 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `itilbert` | 2004 | 86.898 / 0.09966 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `next_fast_len` | 2016 | 0.0 / 0.0 | —; unavailable 2016/2016 | 0.0 / 0.0 | —; unavailable 2016/2016 | —; unavailable 2016/2016 | —; unavailable 2016/2016 |
| `prev_fast_len` | 2016 | 0.0 / 0.0 | —; unavailable 2016/2016 | 0.0 / 0.0 | —; unavailable 2016/2016 | —; unavailable 2016/2016 | —; unavailable 2016/2016 |
| `rfftfreq` | 2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 |
| `sc_diff` | 2005 | 407.7 / 0.2497 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | —; unavailable 2005/2005 | —; unavailable 2005/2005 |
| `shift` | 2004 | 5.4813 / 0.016628 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `ss_diff` | 2005 | 18.388 / 0.03237 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | —; unavailable 2005/2005 | —; unavailable 2005/2005 |
| `tilbert` | 2004 | 9.6962 / 0.029734 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |


## Conventions and limitations

Transforms use unnormalized forward and inverse-divided-by-length conventions. Real transforms, multidimensional shapes, transform types, FFTLog and periodic operators follow the pinned fixture conventions. FFTLog uses SciPy fft; distinct periodic operators use its legacy fftpack implementation. Real and imaginary components are interleaved, including exact imaginary zeros. RustFFT evaluations use a separate pinned Rust cdylib and the library's complex FFT kernels, preserving the Sci fixture's result ordering and forward/inverse normalization. NumPy FFT comparisons are independent public NumPy operations on matching input shapes. GSL public complex forward/inverse FFT routines independently evaluate the saved finite, one-dimensional fft/ifft profiles, including non-power-of-two lengths and inverse scaling. Real, multidimensional, cosine/sine and arbitrary-transform-option variants are not thereby covered.

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
