# signal: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-signal)

Measured binary64 (`f64`) results compared against the **saved per-case reference**.
Each implementation, including Sci-signal, is evaluated against the same
reference values stored in the raw CSVs. This page
covers 352489 cases and 176 public API names;
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

### `sci::signal::advanced_design`

10 functions · 20151 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy |
| --- | ---: | --- | --- | --- | --- | --- |
| `band_stop_obj` | 2008 | 318.36 / 83.42 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `bessel` | 2024 | 5.5871 / 0.043311 | —; unavailable 2024/2024 | 0.0 / 0.0 | —; unavailable 2024/2024 | —; unavailable 2024/2024 |
| `besselap` | 2039 | 137.88 / 3.1479 | —; unavailable 2039/2039 | 0.0 / 0.0 | —; unavailable 2039/2039 | —; unavailable 2039/2039 |
| `ellip` | 2024 | 4.954e+16 / 8.6359e+12 | —; unavailable 2024/2024 | 0.0 / 0.0 | —; unavailable 2024/2024 | —; unavailable 2024/2024 |
| `ellipap` | 2012 | 29.688 / 2.4795 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | —; unavailable 2012/2012 |
| `ellipord` | 2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | —; unavailable 2005/2005 |
| `findfreqs` | 2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `gammatone` | 2008 | 43.952 / 2.3769 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `iircomb` | 2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `iirdesign` | 2020 | 10.876 / 3.6349 | —; unavailable 2020/2020 | 0.0 / 0.0 | —; unavailable 2020/2020 | —; unavailable 2020/2020 |

### `sci::signal::convolution`

8 functions · 16168 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy |
| --- | ---: | --- | --- | --- | --- | --- |
| `convolve` | 2057 | 3.3777e+15 / 2.3581e+11 | 0.49074 / 0.17831 | 3.3777e+15 / 2.3581e+11 | —; unavailable 2057/2057 | 3.3777e+15 / 2.3581e+11 |
| `correlate` | 2057 | 580.78 / 0.5823 | 0.49433 / 0.17742 | 580.78 / 0.5823 | —; unavailable 2057/2057 | 580.78 / 0.5823 |
| `correlation_lags` | 2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | —; unavailable 2012/2012 |
| `decimate` | 2012 | 89.245 / 1.3162 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | —; unavailable 2012/2012 |
| `deconvolve` | 2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `fftconvolve` | 2012 | 4.5036e+15 / 1.6821e+12 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | 1.5 / 0.1226; failed 1/2012 |
| `resample_poly` | 2010 | 34.154 / 4.0757e-3 | —; unavailable 2010/2010 | 0.0 / 0.0 | —; unavailable 2010/2010 | —; unavailable 2010/2010 |
| `upfirdn` | 2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | —; unavailable 2005/2005 |

### `sci::signal::design`

28 functions · 56165 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy |
| --- | ---: | --- | --- | --- | --- | --- |
| `bilinear` | 2002 | 0.65691 / 0.26587 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `bilinear_zpk` | 2002 | 1.5627 / 0.12371 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `bode` | 2002 | 28.571 / 1.9145 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `buttap` | 2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `butter` | 2032 | 5.4464 / 0.033037 | —; unavailable 2032/2032 | 0.0 / 0.0 | —; unavailable 2032/2032 | —; unavailable 2032/2032 |
| `buttord` | 2004 | 0.68007 / 3.3769e-4 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `cheb1ap` | 2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `cheb1ord` | 2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `cheb2ap` | 2004 | 1.0416 / 4.2034e-3 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `cheb2ord` | 2004 | 1.0251 / 2.5575e-4 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `cheby1` | 2032 | 15.347 / 0.040186 | —; unavailable 2032/2032 | 0.0 / 0.0 | —; unavailable 2032/2032 | —; unavailable 2032/2032 |
| `cheby2` | 2032 | 1.2981e+16 / 3.0529e+12 | —; unavailable 2032/2032 | 0.0 / 0.0 | —; unavailable 2032/2032 | —; unavailable 2032/2032 |
| `dbode` | 2003 | 14.983 / 0.69538 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `freqs_zpk` | 2002 | 1.0165 / 0.17286 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `freqz_zpk` | 2002 | 3.2497 / 0.48942 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `iirfilter` | 2006 | 1.2981e+16 / 8.0798e+11 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `lp2bp` | 2002 | 7.0 / 2.3208e-3 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `lp2bp_zpk` | 2002 | 0.53333 / 0.056134 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `lp2bs` | 2002 | 14.0 / 0.010471 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `lp2bs_zpk` | 2002 | 0.90212 / 5.8197e-5 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `lp2hp` | 2002 | 14.0 / 7.3015e-3 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `lp2hp_zpk` | 2002 | 0.8 / 6.1477e-5 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `lp2lp` | 2002 | 14.0 / 3.8276e-3 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `lp2lp_zpk` | 2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `normalize` | 2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `sos2zpk` | 2002 | 1.0 / 0.30555 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `tf2zpk` | 2003 | 11.0 / 4.1576e-3 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `zpk2tf` | 2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |

### `sci::signal::filtering`

25 functions · 50115 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy |
| --- | ---: | --- | --- | --- | --- | --- |
| `filtfilt` | 2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `firwin` | 2010 | 7.4555 / 1.7444 | —; unavailable 2010/2010 | 0.0 / 0.0 | —; unavailable 2010/2010 | —; unavailable 2010/2010 |
| `freqs` | 2002 | 1.9994 / 0.43274 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `freqz` | 2002 | 1.7141 / 0.37285 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `freqz_sos` | 2002 | 6.6952 / 0.42815 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `group_delay` | 2002 | 33.903 / 4.5991 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `iirnotch` | 2002 | 2.0423 / 1.6387e-3 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `iirpeak` | 2002 | 96.547 / 0.049031 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `kaiser_atten` | 2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `kaiser_beta` | 2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `kaiserord` | 2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `lfilter` | 2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `lfilter_with_state` | 2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `lfilter_zi` | 2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `lfiltic` | 2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `medfilt` | 2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 | 0.0 / 0.0 | —; unavailable 2003/2003 |
| `order_filter` | 2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | —; unavailable 2009/2009 |
| `savgol_coeffs` | 2016 | 1.6413e+16 / 8.0553e+12 | —; unavailable 2016/2016 | 0.0 / 0.0 | —; unavailable 2016/2016 | —; unavailable 2016/2016 |
| `savgol_filter` | 2012 | 1.8014e+16 / 7.4612e+11 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | —; unavailable 2012/2012 |
| `sos2tf` | 2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `sosfilt` | 2002 | 8.3429 / 0.26274 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `sosfilt_zi` | 2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `sosfiltfilt` | 2004 | 3.3461 / 0.33608 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `sosfreqz` | 2002 | 6.6952 / 0.42457 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `wiener` | 2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | 0.0 / 0.0 | 0.86957 / 7.9271e-5 | —; unavailable 2005/2005 |

### `sci::signal::fir`

5 functions · 10038 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy |
| --- | ---: | --- | --- | --- | --- | --- |
| `firls` | 2009 | 4313.1 / 258.56 | —; unavailable 2009/2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | —; unavailable 2009/2009 |
| `firwin2` | 2009 | 1.3324e+15 / 2.7324e+13 | —; unavailable 2009/2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | —; unavailable 2009/2009 |
| `firwin_2d` | 2002 | 4.8098 / 1.1811 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `minimum_phase` | 2012 | 229.3 / 67.865 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | —; unavailable 2012/2012 |
| `remez` | 2006 | 3.8424e+8 / 1.5651e+7 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |

### `sci::signal::multidimensional`

5 functions · 10045 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy |
| --- | ---: | --- | --- | --- | --- | --- |
| `convolve2d` | 2018 | 0.0 / 0.0 | —; unavailable 2018/2018 | 0.0 / 0.0 | —; unavailable 2018/2018 | —; unavailable 2018/2018 |
| `correlate2d` | 2018 | 0.0 / 0.0 | —; unavailable 2018/2018 | 0.0 / 0.0 | —; unavailable 2018/2018 | —; unavailable 2018/2018 |
| `hilbert2` | 2004 | 2.3685e+16 / 2.7049e+12 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `medfilt2d` | 2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 | 0.0 / 0.0 | —; unavailable 2003/2003 |
| `sepfir2d` | 2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |

### `sci::signal::peaks`

6 functions · 12082 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy |
| --- | ---: | --- | --- | --- | --- | --- |
| `argrelextrema` | 2016 | 0.0 / 0.0 | —; unavailable 2016/2016 | 0.0 / 0.0 | —; unavailable 2016/2016 | —; unavailable 2016/2016 |
| `argrelmax` | 2016 | 0.0 / 0.0 | —; unavailable 2016/2016 | 0.0 / 0.0 | —; unavailable 2016/2016 | —; unavailable 2016/2016 |
| `argrelmin` | 2016 | 0.0 / 0.0 | —; unavailable 2016/2016 | 0.0 / 0.0 | —; unavailable 2016/2016 | —; unavailable 2016/2016 |
| `find_peaks` | 2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `peak_prominences` | 2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `peak_widths` | 2024 | 0.0 / 0.0 | —; unavailable 2024/2024 | 0.0 / 0.0 | —; unavailable 2024/2024 | —; unavailable 2024/2024 |

### `sci::signal::representations`

7 functions · 14041 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy |
| --- | ---: | --- | --- | --- | --- | --- |
| `invres` | 2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `invresz` | 2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `residue` | 2008 | 19.0 / 8.2256e-3 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `residuez` | 2008 | 19.0 / 8.8041e-3 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `tf2sos` | 2004 | 94.688 / 0.31332 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `unique_roots` | 2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | —; unavailable 2009/2009 |
| `zpk2sos` | 2004 | 110.68 / 1.0752 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |

### `sci::signal::spectral`

17 functions · 34140 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy |
| --- | ---: | --- | --- | --- | --- | --- |
| `check_COLA` | 2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | —; unavailable 2009/2009 |
| `check_NOLA` | 2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | —; unavailable 2009/2009 |
| `check_cola` | 2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | —; unavailable 2009/2009 |
| `check_nola` | 2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | —; unavailable 2009/2009 |
| `coherence` | 2003 | 34.906 / 2.6874 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `csd` | 2012 | 50.358 / 1.7422 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | —; unavailable 2012/2012 |
| `czt` | 2004 | 22.931 / 1.9899 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `czt_points` | 2006 | 2.4271 / 4.3636e-3 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `hilbert` | 2003 | 4.5036e+15 / 1.1203e+12 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `istft` | 2008 | 1.0058e+17 / 2.1743e+15 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `lombscargle` | 2004 | 64.822 / 8.0263 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `periodogram` | 2012 | 6.7554e+15 / 2.5913e+14 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | —; unavailable 2012/2012 |
| `resample` | 2018 | 4.5036e+15 / 2.1767e+12 | —; unavailable 2018/2018 | 0.0 / 0.0 | —; unavailable 2018/2018 | —; unavailable 2018/2018 |
| `spectrogram` | 2012 | 5.5014 / 0.020721 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | —; unavailable 2012/2012 |
| `stft` | 2008 | 4.5036e+15 / 1.2082e+11 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `welch` | 2012 | 6.0447 / 0.39974 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | —; unavailable 2012/2012 |
| `zoom_fft` | 2002 | 220.68 / 9.9508 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |

### `sci::signal::splines`

10 functions · 20053 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy |
| --- | ---: | --- | --- | --- | --- | --- |
| `cspline1d` | 2008 | 5.0063 / 0.039147 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `cspline1d_eval` | 2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `cspline2d` | 2004 | 809.5 / 0.061668 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `qspline1d` | 2004 | 1.8378 / 5.9731e-3 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `qspline1d_eval` | 2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `qspline2d` | 2002 | 12.471 / 0.03535 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `spline_filter` | 2004 | 194.76 / 0.049308 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `symiirorder1` | 2003 | 1.9875 / 3.6143e-4 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `symiirorder2` | 2006 | 961.31 / 0.85295 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `whittaker_henderson` | 2018 | 103.27 / 0.051909 | —; unavailable 2018/2018 | 0.0 / 0.0 | —; unavailable 2018/2018 | —; unavailable 2018/2018 |

### `sci::signal::systems`

13 functions · 26058 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy |
| --- | ---: | --- | --- | --- | --- | --- |
| `abcd_normalize` | 2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `cont2discrete` | 2012 | 2.0687 / 0.22102 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | —; unavailable 2012/2012 |
| `dimpulse` | 2004 | 0.625 / 8.3112e-5 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `dlsim` | 2003 | 0.625 / 0.1041 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `dstep` | 2004 | 0.625 / 8.3112e-5 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `impulse` | 2002 | 3.3738 / 0.27516 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `lsim` | 2008 | 19.703 / 0.24737 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `place_poles` | 2002 | 6.5062 / 1.0856 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `ss2tf` | 2004 | 4.5036e+15 / 1.1216e+12 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `ss2zpk` | 2007 | 4.5036e+15 / 7.4452e+11 | —; unavailable 2007/2007 | 0.0 / 0.0 | —; unavailable 2007/2007 | —; unavailable 2007/2007 |
| `step` | 2004 | 2.8941 / 0.22663 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `tf2ss` | 2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `zpk2ss` | 2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |

### `sci::signal::utilities`

6 functions · 11128 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy |
| --- | ---: | --- | --- | --- | --- | --- |
| `closest_STFT_dual_window` | 2008 | 58.508 / 2.6203 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `envelope` | 2072 | 6.0 / 0.48479 | —; unavailable 2072/2072 | 0.0 / 0.0 | —; unavailable 2072/2072 | —; unavailable 2072/2072 |
| `find_peaks_cwt` | 2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `max_len_seq` | 2013 | 0.0 / 0.0 | —; unavailable 2013/2013 | 0.0 / 0.0 | —; unavailable 2013/2013 | —; unavailable 2013/2013 |
| `oaconvolve` | 2012 | 9.0072e+15 / 7.6696e+13 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | 4.5036e+15 / 1.1402e+14; failed 1/2012 |
| `utilities.closest_stft_dual_window` | 1015 | 4.5036e+15 / 2.6721e+11 | —; unavailable 1015/1015 | 0.0 / 0.0 | —; unavailable 1015/1015 | —; unavailable 1015/1015 |

### `sci::signal::waveforms`

9 functions · 18092 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy |
| --- | ---: | --- | --- | --- | --- | --- |
| `chirp` | 2016 | 1.9247 / 9.5472e-4 | —; unavailable 2016/2016 | 0.0 / 0.0 | —; unavailable 2016/2016 | —; unavailable 2016/2016 |
| `detrend` | 2006 | 4.5036e+15 / 4.3606e+15 | —; unavailable 2006/2006 | 0.0 / 0.0 | 4.5036e+15 / 4.3606e+15 | 4.5036e+15 / 4.3614e+15 |
| `gauss_spline` | 2009 | 0.60225 / 0.1313 | —; unavailable 2009/2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | —; unavailable 2009/2009 |
| `gausspulse` | 2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `sawtooth` | 2024 | 2.2254 / 0.13805 | —; unavailable 2024/2024 | 0.0 / 0.0 | —; unavailable 2024/2024 | —; unavailable 2024/2024 |
| `square` | 2024 | 0.0 / 0.0 | —; unavailable 2024/2024 | 0.0 / 0.0 | —; unavailable 2024/2024 | —; unavailable 2024/2024 |
| `sweep_poly` | 2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `unit_impulse` | 2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `vectorstrength` | 2003 | 9.5247e+15 / 2.3776e+12 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |

### `sci::signal::windows`

27 functions · 54213 cases

| Function | Cases | Sci | mpmath | SciPy | GSL | NumPy |
| --- | ---: | --- | --- | --- | --- | --- |
| `get_window` | 2006 | 10.577 / 0.98053 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `windows.barthann` | 2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `windows.bartlett` | 2008 | 3.875 / 0.10895 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | 7.75 / 0.47198 |
| `windows.blackman` | 2008 | 49.978 / 1.3382 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | 45.112 / 1.6828 |
| `windows.blackmanharris` | 2008 | 261.57 / 2.7067 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `windows.bohman` | 2008 | 349.14 / 2.8582 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `windows.boxcar` | 2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `windows.chebwin` | 2008 | 176.17 / 5.6427 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `windows.cosine` | 2008 | 39.487 / 0.99122 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `windows.dpss` | 2006 | 4.1693e+17 / 7.1674e+14 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `windows.exponential` | 2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `windows.flattop` | 2008 | 283.86 / 3.2804 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `windows.gaussian` | 2008 | 32.026 / 3.6175 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `windows.general_cosine` | 2008 | 66.002 / 1.5512 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `windows.general_gaussian` | 2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `windows.general_hamming` | 2008 | 4.1065 / 0.30727 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `windows.get_window` | 2006 | 10.577 / 0.96974 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `windows.hamming` | 2008 | 9.9559 / 1.1893 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | 8.3494 / 1.237 |
| `windows.hann` | 2008 | 48.852 / 1.1563 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | 34.414 / 1.382 |
| `windows.kaiser` | 2011 | 8.1057 / 0.9771 | —; unavailable 2011/2011 | 0.0 / 0.0 | —; unavailable 2011/2011 | —; unavailable 2011/2011 |
| `windows.kaiser_bessel_derived` | 2004 | 2.7043 / 0.22147 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `windows.lanczos` | 2008 | 9.3543 / 0.32059 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `windows.nuttall` | 2008 | 110.28 / 1.8934 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `windows.parzen` | 2008 | 5.2204 / 0.23339 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `windows.taylor` | 2008 | 0.98472 / 0.11782 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `windows.triang` | 2008 | 4.125 / 0.39496 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `windows.tukey` | 2012 | 48.629 / 0.75994 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | —; unavailable 2012/2012 |

## Evaluation failures

- `signal.fftconvolve` (NumPy): 1 failed; shape mismatch: expected 2 scalar outputs, received 3. Example: `signal.fftconvolve/recipe-0231`.
- `signal.oaconvolve` (NumPy): 1 failed; shape mismatch: expected 3 scalar outputs, received 10. Example: `signal.oaconvolve/recipe-1258`.


## Conventions and limitations

The pinned signal fixture adapters cover waveform, filtering, design, spectral, window, peak, spline, convolution and system functions through the public Sci facade. Output tuples flatten in the documented fixture order. Parameters, modes, shapes and normalization are retained in each function CSV; original real convolution/correlation also call NumPy. The limited GSL statistical-filter competitors compose public GSL linear-regression, mean/variance, sorting and median primitives for the exact saved 1-D/2-D detrend, Wiener and median-filter profiles. These are derived DSP implementations, not GSL exports of Sci's named signal-processing algorithms.

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
