# signal: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-signal)

Measured binary64 (`f64`) results using **Sci-signal as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 358507 cases and 179 public API names;
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

### `sci::signal::advanced_design` (10 functions, 20151 cases)

| Function | Cases | Sci vs reference | SciPy | NumPy |
| --- | ---: | --- | --- | --- |
| `band_stop_obj` | 2008 | 318.36 / 83.42 | 318.36 / 83.42 | — |
| `bessel` | 2024 | 5.5871 / 0.043311 | 5.5871 / 0.043311 | — |
| `besselap` | 2039 | 137.88 / 3.1479 | 137.88 / 3.1479 | — |
| `ellip` | 2024 | 4.954e+16 / 8.6359e+12 | 4.954e+16 / 8.6359e+12 | — |
| `ellipap` | 2012 | 29.688 / 2.4795 | 29.688 / 2.4795 | — |
| `ellipord` | 2005 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `findfreqs` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `gammatone` | 2008 | 43.952 / 2.3769 | 43.952 / 2.3769 | — |
| `iircomb` | 2008 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `iirdesign` | 2020 | 10.876 / 3.6349 | 10.876 / 3.6349 | — |

### `sci::signal::convolution` (8 functions, 16168 cases)

| Function | Cases | Sci vs reference | SciPy | NumPy |
| --- | ---: | --- | --- | --- |
| `convolve` | 2057 | 3.3777e+15 / 2.3581e+11 | 11.253 / 1.3636e-3 | 11.253 / 1.3636e-3 |
| `correlate` | 2057 | 580.78 / 0.5823 | 0.0 / 0.0 | 0.0 / 0.0 |
| `correlation_lags` | 2012 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `decimate` | 2012 | 89.245 / 1.3162 | 89.245 / 1.3162 | — |
| `deconvolve` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `fftconvolve` | 2012 | 4.5036e+15 / 1.6821e+12 | 4.5036e+15 / 1.6821e+12 | — |
| `resample_poly` | 2010 | 34.154 / 4.0757e-3 | 34.154 / 4.0757e-3 | — |
| `upfirdn` | 2005 | 0.0 / 0.0 | 0.0 / 0.0 | — |

### `sci::signal::design` (28 functions, 56165 cases)

| Function | Cases | Sci vs reference | SciPy | NumPy |
| --- | ---: | --- | --- | --- |
| `bilinear` | 2002 | 0.65691 / 0.26587 | 0.65691 / 0.26587 | — |
| `bilinear_zpk` | 2002 | 1.5627 / 0.12371 | 1.5627 / 0.12371 | — |
| `bode` | 2002 | 28.571 / 1.9145 | 28.571 / 1.9145 | — |
| `buttap` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `butter` | 2032 | 5.4464 / 0.033037 | 5.4464 / 0.033037 | — |
| `buttord` | 2004 | 0.68007 / 3.3769e-4 | 0.68007 / 3.3769e-4 | — |
| `cheb1ap` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `cheb1ord` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `cheb2ap` | 2004 | 1.0416 / 4.2034e-3 | 1.0416 / 4.2034e-3 | — |
| `cheb2ord` | 2004 | 1.0251 / 2.5575e-4 | 1.0251 / 2.5575e-4 | — |
| `cheby1` | 2032 | 15.347 / 0.040186 | 15.347 / 0.040186 | — |
| `cheby2` | 2032 | 1.2981e+16 / 3.0529e+12 | 1.2981e+16 / 3.0529e+12 | — |
| `dbode` | 2003 | 14.983 / 0.69538 | 14.983 / 0.69538 | — |
| `freqs_zpk` | 2002 | 1.0165 / 0.17286 | 1.0165 / 0.17286 | — |
| `freqz_zpk` | 2002 | 3.2497 / 0.48942 | 3.2497 / 0.48942 | — |
| `iirfilter` | 2006 | 1.2981e+16 / 8.0798e+11 | 1.2981e+16 / 8.0798e+11 | — |
| `lp2bp` | 2002 | 7.0 / 2.3208e-3 | 7.0 / 2.3208e-3 | — |
| `lp2bp_zpk` | 2002 | 0.53333 / 0.056134 | 0.53333 / 0.056134 | — |
| `lp2bs` | 2002 | 14.0 / 0.010471 | 14.0 / 0.010471 | — |
| `lp2bs_zpk` | 2002 | 0.90212 / 5.8197e-5 | 0.90212 / 5.8197e-5 | — |
| `lp2hp` | 2002 | 14.0 / 7.3015e-3 | 14.0 / 7.3015e-3 | — |
| `lp2hp_zpk` | 2002 | 0.8 / 6.1477e-5 | 0.8 / 6.1477e-5 | — |
| `lp2lp` | 2002 | 14.0 / 3.8276e-3 | 14.0 / 3.8276e-3 | — |
| `lp2lp_zpk` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `normalize` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `sos2zpk` | 2002 | 1.0 / 0.30555 | 1.0 / 0.30555 | — |
| `tf2zpk` | 2003 | 11.0 / 4.1576e-3 | 11.0 / 4.1576e-3 | — |
| `zpk2tf` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — |

### `sci::signal::filtering` (24 functions, 48113 cases)

| Function | Cases | Sci vs reference | SciPy | NumPy |
| --- | ---: | --- | --- | --- |
| `filtfilt` | 2008 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `firwin` | 2010 | 7.4555 / 1.7444 | 7.4555 / 1.7444 | — |
| `freqs` | 2002 | 1.9994 / 0.43274 | 1.9994 / 0.43274 | — |
| `freqz` | 2002 | 1.7141 / 0.37285 | 1.7141 / 0.37285 | — |
| `freqz_sos` | 2002 | 6.6952 / 0.42815 | 6.6952 / 0.42815 | — |
| `group_delay` | 2002 | 33.903 / 4.5991 | 33.903 / 4.5991 | — |
| `iirnotch` | 2002 | 2.0423 / 1.6387e-3 | 2.0423 / 1.6387e-3 | — |
| `iirpeak` | 2002 | 96.547 / 0.049031 | 96.547 / 0.049031 | — |
| `kaiser_atten` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `kaiser_beta` | 2006 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `kaiserord` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `lfilter` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `lfilter_with_state` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `lfilter_zi` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `lfiltic` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `medfilt` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `order_filter` | 2009 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `savgol_coeffs` | 2016 | 1.6413e+16 / 8.0553e+12 | 1.6413e+16 / 8.0553e+12 | — |
| `savgol_filter` | 2012 | 1.8014e+16 / 7.4612e+11 | 1.8014e+16 / 7.4612e+11 | — |
| `sos2tf` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `sosfilt` | 2002 | 8.3429 / 0.26274 | 8.3429 / 0.26274 | — |
| `sosfilt_zi` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `sosfiltfilt` | 2004 | 3.3461 / 0.33608 | 3.3461 / 0.33608 | — |
| `wiener` | 2005 | 0.0 / 0.0 | 0.0 / 0.0 | — |

### `sci::signal::fir` (5 functions, 10038 cases)

| Function | Cases | Sci vs reference | SciPy | NumPy |
| --- | ---: | --- | --- | --- |
| `firls` | 2009 | 4313.1 / 258.56 | 4313.1 / 258.56 | — |
| `firwin2` | 2009 | 1.3324e+15 / 2.7324e+13 | 1.3324e+15 / 2.7324e+13 | — |
| `firwin_2d` | 2002 | 4.8098 / 1.1811 | 4.8098 / 1.1811 | — |
| `minimum_phase` | 2012 | 229.3 / 67.865 | 229.3 / 67.865 | — |
| `remez` | 2006 | 3.8424e+8 / 1.5651e+7 | 3.8424e+8 / 1.5651e+7 | — |

### `sci::signal::multidimensional` (5 functions, 10045 cases)

| Function | Cases | Sci vs reference | SciPy | NumPy |
| --- | ---: | --- | --- | --- |
| `convolve2d` | 2018 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `correlate2d` | 2018 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `hilbert2` | 2004 | 2.3685e+16 / 2.7049e+12 | 2.3685e+16 / 2.7049e+12 | — |
| `medfilt2d` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `sepfir2d` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — |

### `sci::signal::peaks` (6 functions, 12082 cases)

| Function | Cases | Sci vs reference | SciPy | NumPy |
| --- | ---: | --- | --- | --- |
| `argrelextrema` | 2016 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `argrelmax` | 2016 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `argrelmin` | 2016 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `find_peaks` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `peak_prominences` | 2006 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `peak_widths` | 2024 | 0.0 / 0.0 | 0.0 / 0.0 | — |

### `sci::signal::representations` (7 functions, 14041 cases)

| Function | Cases | Sci vs reference | SciPy | NumPy |
| --- | ---: | --- | --- | --- |
| `invres` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `invresz` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `residue` | 2008 | 19.0 / 8.2256e-3 | 19.0 / 8.2256e-3 | — |
| `residuez` | 2008 | 19.0 / 8.8041e-3 | 19.0 / 8.8041e-3 | — |
| `tf2sos` | 2004 | 94.688 / 0.31332 | 94.688 / 0.31332 | — |
| `unique_roots` | 2009 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `zpk2sos` | 2004 | 110.68 / 1.0752 | 110.68 / 1.0752 | — |

### `sci::signal::spectral` (15 functions, 30122 cases)

| Function | Cases | Sci vs reference | SciPy | NumPy |
| --- | ---: | --- | --- | --- |
| `check_cola` | 2009 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `check_nola` | 2009 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `coherence` | 2003 | 34.906 / 2.6874 | 34.906 / 2.6874 | — |
| `csd` | 2012 | 50.358 / 1.7422 | 50.358 / 1.7422 | — |
| `czt` | 2004 | 22.931 / 1.9899 | 22.931 / 1.9899 | — |
| `czt_points` | 2006 | 2.4271 / 4.3636e-3 | 2.4271 / 4.3636e-3 | — |
| `hilbert` | 2003 | 4.5036e+15 / 1.1203e+12 | 4.5036e+15 / 1.1203e+12 | — |
| `istft` | 2008 | 1.0058e+17 / 2.1743e+15 | 1.0058e+17 / 2.1743e+15 | — |
| `lombscargle` | 2004 | 64.822 / 8.0263 | 64.822 / 8.0263 | — |
| `periodogram` | 2012 | 6.7554e+15 / 2.5913e+14 | 6.7554e+15 / 2.5913e+14 | — |
| `resample` | 2018 | 4.5036e+15 / 2.1767e+12 | 4.5036e+15 / 2.1767e+12 | — |
| `spectrogram` | 2012 | 5.5014 / 0.020721 | 5.5014 / 0.020721 | — |
| `stft` | 2008 | 4.5036e+15 / 1.2082e+11 | 4.5036e+15 / 1.2082e+11 | — |
| `welch` | 2012 | 6.0447 / 0.39974 | 6.0447 / 0.39974 | — |
| `zoom_fft` | 2002 | 220.68 / 9.9508 | 220.68 / 9.9508 | — |

### `sci::signal::splines` (10 functions, 20053 cases)

| Function | Cases | Sci vs reference | SciPy | NumPy |
| --- | ---: | --- | --- | --- |
| `cspline1d` | 2008 | 5.0063 / 0.039147 | 5.0063 / 0.039147 | — |
| `cspline1d_eval` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `cspline2d` | 2004 | 809.5 / 0.061668 | 809.5 / 0.061668 | — |
| `qspline1d` | 2004 | 1.8378 / 5.9731e-3 | 1.8378 / 5.9731e-3 | — |
| `qspline1d_eval` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `qspline2d` | 2002 | 12.471 / 0.03535 | 12.471 / 0.03535 | — |
| `spline_filter` | 2004 | 194.76 / 0.049308 | 194.76 / 0.049308 | — |
| `symiirorder1` | 2003 | 1.9875 / 3.6143e-4 | 1.9875 / 3.6143e-4 | — |
| `symiirorder2` | 2006 | 961.31 / 0.85295 | 961.31 / 0.85295 | — |
| `whittaker_henderson` | 2018 | 51.0 / 0.029234 | 51.0 / 0.029234 | — |

### `sci::signal::systems` (13 functions, 26058 cases)

| Function | Cases | Sci vs reference | SciPy | NumPy |
| --- | ---: | --- | --- | --- |
| `abcd_normalize` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `cont2discrete` | 2012 | 2.0687 / 0.22102 | 2.0687 / 0.22102 | — |
| `dimpulse` | 2004 | 0.625 / 8.3112e-5 | 0.625 / 8.3112e-5 | — |
| `dlsim` | 2003 | 0.625 / 0.1041 | 0.625 / 0.1041 | — |
| `dstep` | 2004 | 0.625 / 8.3112e-5 | 0.625 / 8.3112e-5 | — |
| `impulse` | 2002 | 3.3738 / 0.27516 | 3.3738 / 0.27516 | — |
| `lsim` | 2008 | 19.703 / 0.24737 | 19.703 / 0.24737 | — |
| `place_poles` | 2002 | 6.5062 / 1.0856 | 6.5062 / 1.0856 | — |
| `ss2tf` | 2004 | 4.5036e+15 / 1.1216e+12 | 4.5036e+15 / 1.1216e+12 | — |
| `ss2zpk` | 2007 | 4.5036e+15 / 7.4452e+11 | 4.5036e+15 / 7.4452e+11 | — |
| `step` | 2004 | 2.8941 / 0.22663 | 2.8941 / 0.22663 | — |
| `tf2ss` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `zpk2ss` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — |

### `sci::signal::utilities` (5 functions, 9120 cases)

| Function | Cases | Sci vs reference | SciPy | NumPy |
| --- | ---: | --- | --- | --- |
| `envelope` | 2072 | 6.0 / 0.48479 | 6.0 / 0.48479 | — |
| `find_peaks_cwt` | 2008 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `max_len_seq` | 2013 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `oaconvolve` | 2012 | 9.0072e+15 / 7.6696e+13 | 9.0072e+15 / 7.6696e+13 | — |
| `utilities.closest_stft_dual_window` | 1015 | 4.5036e+15 / 2.6721e+11 | 4.5036e+15 / 2.6721e+11 | — |

### `sci::signal::waveforms` (9 functions, 18092 cases)

| Function | Cases | Sci vs reference | SciPy | NumPy |
| --- | ---: | --- | --- | --- |
| `chirp` | 2016 | 1.9247 / 9.5472e-4 | 1.9247 / 9.5472e-4 | — |
| `detrend` | 2006 | 4.5036e+15 / 4.3606e+15 | 4.5036e+15 / 4.3606e+15 | — |
| `gauss_spline` | 2009 | 0.60225 / 0.1313 | 0.60225 / 0.1313 | — |
| `gausspulse` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `sawtooth` | 2024 | 2.2254 / 0.13805 | 2.2254 / 0.13805 | — |
| `square` | 2024 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `sweep_poly` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `unit_impulse` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `vectorstrength` | 2003 | 9.5247e+15 / 2.3776e+12 | 9.5247e+15 / 2.3776e+12 | — |

### `sci::signal::windows` (27 functions, 54213 cases)

| Function | Cases | Sci vs reference | SciPy | NumPy |
| --- | ---: | --- | --- | --- |
| `get_window` | 2006 | 10.577 / 0.98053 | 10.577 / 0.98053 | — |
| `windows.barthann` | 2008 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `windows.bartlett` | 2008 | 3.875 / 0.10895 | 3.875 / 0.10895 | — |
| `windows.blackman` | 2008 | 49.978 / 1.3382 | 49.978 / 1.3382 | — |
| `windows.blackmanharris` | 2008 | 261.57 / 2.7067 | 261.57 / 2.7067 | — |
| `windows.bohman` | 2008 | 349.14 / 2.8582 | 349.14 / 2.8582 | — |
| `windows.boxcar` | 2008 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `windows.chebwin` | 2008 | 176.17 / 5.6427 | 176.17 / 5.6427 | — |
| `windows.cosine` | 2008 | 39.487 / 0.99122 | 39.487 / 0.99122 | — |
| `windows.dpss` | 2006 | 8.3416e+17 / 1.6317e+15 | 8.3416e+17 / 1.6317e+15 | — |
| `windows.exponential` | 2008 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `windows.flattop` | 2008 | 283.86 / 3.2804 | 283.86 / 3.2804 | — |
| `windows.gaussian` | 2008 | 32.026 / 3.6175 | 32.026 / 3.6175 | — |
| `windows.general_cosine` | 2008 | 66.002 / 1.5512 | 66.002 / 1.5512 | — |
| `windows.general_gaussian` | 2008 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `windows.general_hamming` | 2008 | 4.1065 / 0.30727 | 4.1065 / 0.30727 | — |
| `windows.get_window` | 2006 | 10.577 / 0.96974 | 10.577 / 0.96974 | — |
| `windows.hamming` | 2008 | 9.9559 / 1.1893 | 9.9559 / 1.1893 | — |
| `windows.hann` | 2008 | 48.852 / 1.1563 | 48.852 / 1.1563 | — |
| `windows.kaiser` | 2011 | 8.1057 / 0.9771 | 8.1057 / 0.9771 | — |
| `windows.kaiser_bessel_derived` | 2004 | 2.7043 / 0.22147 | 2.7043 / 0.22147 | — |
| `windows.lanczos` | 2008 | 9.3543 / 0.32059 | 9.3543 / 0.32059 | — |
| `windows.nuttall` | 2008 | 110.28 / 1.8934 | 110.28 / 1.8934 | — |
| `windows.parzen` | 2008 | 5.2204 / 0.23339 | 5.2204 / 0.23339 | — |
| `windows.taylor` | 2008 | 0.98472 / 0.11782 | 0.98472 / 0.11782 | — |
| `windows.triang` | 2008 | 4.125 / 0.39496 | 4.125 / 0.39496 | — |
| `windows.tukey` | 2012 | 48.629 / 0.75994 | 48.629 / 0.75994 | — |

### `sci::signal` (7 functions, 14046 cases)

| Function | Cases | Sci vs reference | SciPy | NumPy |
| --- | ---: | --- | --- | --- |
| `check_COLA` | 2009 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `check_NOLA` | 2009 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `choose_conv_method` | 2012 | 0.0 / 0.0 | 0.0 / 0.0 | — |
| `closest_STFT_dual_window` | 2008 | 58.508 / 2.6203 | 58.508 / 2.6203 | — |
| `dfreqresp` | 2003 | 2.5684 / 0.20998 | 2.5684 / 0.20998 | — |
| `freqresp` | 2003 | 0.99656 / 0.14675 | 0.99656 / 0.14675 | — |
| `sosfreqz` | 2002 | 6.6952 / 0.42457 | 6.6952 / 0.42457 | — |


## Conventions and limitations

The pinned signal fixture adapters cover waveform, filtering, design, spectral, window, peak, spline, convolution and system functions through the public Sci facade. Output tuples flatten in the documented fixture order. Parameters, modes, shapes and normalization are retained in each function CSV; original real convolution/correlation also call NumPy.

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
