# fft: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-fft)

Measured binary64 (`f64`) results using **Sci-fft as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 90336 cases and 45 public API names;
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

### `sci::fft::fourier`

6 functions · 12098 cases

| Function | Cases | Sci vs reference | mpmath | SciPy | GSL | NumPy | RustFFT |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `fft` | 2040 | 9.0072e+15 / 1.8055e+13 | 1.0864e+91 / 1.1059e+87 | 2115.0 / 1.8359 | 4.5036e+15 / 3.9028e+12 | 2115.0 / 1.8359 | 0.0 / 0.0 |
| `fft2` | 2004 | 4.0532e+16 / 3.708e+15 | — | 4.0532e+16 / 3.708e+15 | — | 4.0532e+16 / 3.708e+15 | 0.0 / 0.0 |
| `fftn` | 2005 | 1.0085e+17 / 3.6571e+15 | — | 1.0085e+17 / 3.6571e+15 | — | 4.0532e+16 / 3.6541e+15 | 0.0 / 0.0 |
| `ifft` | 2040 | 9.0072e+15 / 1.8055e+13 | 1.6296e+91 / 1.3096e+87 | 2115.0 / 1.9261 | 4.5036e+15 / 4.0211e+12 | 2115.0 / 1.9261 | 0.99999 / 0.1004 |
| `ifft2` | 2004 | 5.2326e+16 / 3.4286e+15 | — | 5.2326e+16 / 3.4286e+15 | — | 5.2326e+16 / 3.4287e+15 | 4.5036e+15 / 3.5568e+11 |
| `ifftn` | 2005 | 7.676e+16 / 3.4467e+15 | — | 7.676e+16 / 3.4467e+15 | — | 5.2326e+16 / 3.4444e+15 | 0.99999 / 0.20041 |

### `sci::fft::hankel`

3 functions · 6018 cases

| Function | Cases | Sci vs reference | mpmath | SciPy | GSL | NumPy | RustFFT |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `fht` | 2006 | 619.05 / 0.53128 | — | 619.05 / 0.53128 | — | — | — |
| `fhtoffset` | 2006 | 0.0 / 0.0 | — | 0.0 / 0.0 | — | — | — |
| `ifht` | 2006 | 612.77 / 0.27612 | — | 612.77 / 0.27612 | — | — | — |

### `sci::fft::real`

20 functions · 40118 cases

| Function | Cases | Sci vs reference | mpmath | SciPy | GSL | NumPy | RustFFT |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `dct` | 2008 | 1.3511e+16 / 3.9133e+12 | — | 1.3511e+16 / 3.9133e+12 | — | — | — |
| `dctn` | 2008 | 4.5036e+15 / 6.5528e+14 | — | 4.5036e+15 / 6.5528e+14 | — | — | — |
| `dst` | 2008 | 4.954e+16 / 7.2675e+12 | — | 4.954e+16 / 7.2675e+12 | — | — | — |
| `dstn` | 2008 | 1.8014e+16 / 8.1393e+14 | — | 1.8014e+16 / 8.1393e+14 | — | — | — |
| `hfft` | 2005 | 2.435 / 8.2451e-4 | — | 2.435 / 8.2451e-4 | — | 2.435 / 8.2451e-4 | — |
| `hfft2` | 2004 | 6.7214 / 0.50543 | — | 6.7214 / 0.50543 | — | — | — |
| `hfftn` | 2004 | 18.192 / 0.50645 | — | 18.192 / 0.50645 | — | — | — |
| `idct` | 2008 | 1.5763e+16 / 4.2627e+12 | — | 1.5763e+16 / 4.2627e+12 | — | — | — |
| `idctn` | 2008 | 1.8014e+16 / 6.5444e+14 | — | 1.8014e+16 / 6.5444e+14 | — | — | — |
| `idst` | 2008 | 4.954e+16 / 7.2675e+12 | — | 4.954e+16 / 7.2675e+12 | — | — | — |
| `idstn` | 2008 | 6.5869e+30 / 2.7282e+26 | — | 6.5869e+30 / 2.7282e+26 | — | — | — |
| `ihfft` | 2007 | 3.0 / 4.8078e-3 | — | 3.0 / 4.8078e-3 | — | 3.0 / 4.8078e-3 | — |
| `ihfft2` | 2004 | 2.8523e+16 / 2.0101e+15 | — | 2.8523e+16 / 2.0101e+15 | — | — | — |
| `ihfftn` | 2004 | 2.8523e+16 / 2.0549e+15 | — | 2.8523e+16 / 2.0549e+15 | — | — | — |
| `irfft` | 2005 | 2.435 / 9.1757e-4 | — | 2.435 / 9.1757e-4 | — | 2.435 / 9.1757e-4 | 0.6 / 1.4851e-4 |
| `irfft2` | 2004 | 7.5616 / 0.56249 | — | 7.5616 / 0.56249 | — | 7.5616 / 0.56263 | 4.5036e+15 / 3.744e+11 |
| `irfftn` | 2004 | 17.908 / 0.56237 | — | 17.908 / 0.56237 | — | 17.908 / 0.56237 | 0.94737 / 0.056184 |
| `rfft` | 2005 | 3.1474 / 4.3585e-3 | — | 3.1474 / 4.3585e-3 | — | 3.1474 / 4.3585e-3 | 0.75 / 1.8601e-4 |
| `rfft2` | 2004 | 4.0532e+16 / 2.5891e+15 | — | 4.0532e+16 / 2.5891e+15 | — | 4.0532e+16 / 2.5891e+15 | 0.0 / 0.0 |
| `rfftn` | 2004 | 4.0532e+16 / 2.6188e+15 | — | 4.0532e+16 / 2.6188e+15 | — | 4.0532e+16 / 2.6188e+15 | 0.0 / 0.0 |

### `sci::fft::utilities`

16 functions · 32102 cases

| Function | Cases | Sci vs reference | mpmath | SciPy | GSL | NumPy | RustFFT |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `cc_diff` | 2004 | 8.7137 / 0.015853 | — | 8.7137 / 0.015853 | — | — | — |
| `cs_diff` | 2004 | 7.3222 / 0.040276 | — | 7.3222 / 0.040276 | — | — | — |
| `diff` | 2012 | 3.1833e+16 / 3.6827e+13 | — | 3.1833e+16 / 3.6827e+13 | — | — | — |
| `fftfreq` | 2005 | 0.0 / 0.0 | — | 0.0 / 0.0 | — | 0.0 / 0.0 | — |
| `fftshift` | 2005 | 0.0 / 0.0 | — | 0.0 / 0.0 | — | 0.0 / 0.0 | — |
| `hilbert` | 2004 | 28.745 / 0.04542 | — | 28.745 / 0.04542 | — | — | — |
| `ifftshift` | 2005 | 0.0 / 0.0 | — | 0.0 / 0.0 | — | 0.0 / 0.0 | — |
| `ihilbert` | 2004 | 28.745 / 0.04542 | — | 28.745 / 0.04542 | — | — | — |
| `itilbert` | 2004 | 86.898 / 0.09966 | — | 86.898 / 0.09966 | — | — | — |
| `next_fast_len` | 2016 | 0.0 / 0.0 | — | 0.0 / 0.0 | — | — | — |
| `prev_fast_len` | 2016 | 0.0 / 0.0 | — | 0.0 / 0.0 | — | — | — |
| `rfftfreq` | 2005 | 0.0 / 0.0 | — | 0.0 / 0.0 | — | 0.0 / 0.0 | — |
| `sc_diff` | 2005 | 407.7 / 0.2497 | — | 407.7 / 0.2497 | — | — | — |
| `shift` | 2004 | 5.4813 / 0.016628 | — | 5.4813 / 0.016628 | — | — | — |
| `ss_diff` | 2005 | 18.388 / 0.03237 | — | 18.388 / 0.03237 | — | — | — |
| `tilbert` | 2004 | 9.6962 / 0.029734 | — | 9.6962 / 0.029734 | — | — | — |


## Conventions and limitations

Transforms use unnormalized forward and inverse-divided-by-length conventions. Real transforms, multidimensional shapes, transform types, FFTLog and periodic operators follow the pinned fixture conventions. FFTLog uses SciPy fft; distinct periodic operators use its legacy fftpack implementation. Real and imaginary components are interleaved, including exact imaginary zeros. RustFFT evaluations use a separate pinned Rust cdylib and the library's complex FFT kernels, preserving the Sci fixture's result ordering and forward/inverse normalization. NumPy FFT comparisons are independent public NumPy operations on matching input shapes. GSL public complex forward/inverse FFT routines independently evaluate the saved finite, one-dimensional fft/ifft profiles, including non-power-of-two lengths and inverse scaling. Real, multidimensional, cosine/sine and arbitrary-transform-option variants are not thereby covered.

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
