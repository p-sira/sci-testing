# signal: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-signal)

Measured binary64 (`f64`) results using **Sci-signal as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 357492 cases and 178 public API names;
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
| `abcd_normalize` | 2003 | 0.0 / 0.0 | — |
| `argrelextrema` | 2016 | 0.0 / 0.0 | — |
| `argrelmax` | 2016 | 0.0 / 0.0 | — |
| `argrelmin` | 2016 | 0.0 / 0.0 | — |
| `band_stop_obj` | 2008 | 318.36 / 83.42 | — |
| `bessel` | 2024 | 5.5871 / 0.043311 | — |
| `besselap` | 2039 | 137.88 / 3.1479 | — |
| `bilinear` | 2002 | 0.65691 / 0.26587 | — |
| `bilinear_zpk` | 2002 | 1.5627 / 0.12371 | — |
| `bode` | 2002 | 28.571 / 1.9145 | — |
| `buttap` | 2004 | 0.0 / 0.0 | — |
| `butter` | 2032 | 5.4464 / 0.033037 | — |
| `buttord` | 2004 | 0.68007 / 3.3769e-4 | — |
| `cheb1ap` | 2004 | 0.0 / 0.0 | — |
| `cheb1ord` | 2004 | 0.0 / 0.0 | — |
| `cheb2ap` | 2004 | 1.0416 / 4.2034e-3 | — |
| `cheb2ord` | 2004 | 1.0251 / 2.5575e-4 | — |
| `cheby1` | 2032 | 15.347 / 0.040186 | — |
| `cheby2` | 2032 | 1.2981e+16 / 3.0529e+12 | — |
| `check_COLA` | 2009 | 0.0 / 0.0 | — |
| `check_NOLA` | 2009 | 0.0 / 0.0 | — |
| `check_cola` | 2009 | 0.0 / 0.0 | — |
| `check_nola` | 2009 | 0.0 / 0.0 | — |
| `chirp` | 2016 | 1.9247 / 9.5472e-4 | — |
| `choose_conv_method` | 2012 | 0.0 / 0.0 | — |
| `closest_STFT_dual_window` | 2008 | 58.508 / 2.6203 | — |
| `coherence` | 2003 | 34.906 / 2.6874 | — |
| `cont2discrete` | 2012 | 2.0687 / 0.22102 | — |
| `convolve` | 2057 | 11.253 / 1.3636e-3 | 11.253 / 1.3636e-3 |
| `convolve2d` | 2018 | 0.0 / 0.0 | — |
| `correlate` | 2057 | 0.0 / 0.0 | 0.0 / 0.0 |
| `correlate2d` | 2018 | 0.0 / 0.0 | — |
| `correlation_lags` | 2012 | 0.0 / 0.0 | — |
| `csd` | 2012 | 50.358 / 1.7422 | — |
| `cspline1d` | 2008 | 5.0063 / 0.039147 | — |
| `cspline1d_eval` | 2002 | 0.0 / 0.0 | — |
| `cspline2d` | 2004 | 809.5 / 0.061668 | — |
| `czt` | 2004 | 22.931 / 1.9899 | — |
| `czt_points` | 2006 | 2.4271 / 4.3636e-3 | — |
| `dbode` | 2003 | 14.983 / 0.69538 | — |
| `decimate` | 2012 | 89.245 / 1.3162 | — |
| `deconvolve` | 2003 | 0.0 / 0.0 | — |
| `detrend` | 2006 | 4.5036e+15 / 4.3606e+15 | — |
| `dfreqresp` | 2003 | 2.5684 / 0.20998 | — |
| `dimpulse` | 2004 | 0.625 / 8.3112e-5 | — |
| `dlsim` | 2003 | 0.625 / 0.1041 | — |
| `dstep` | 2004 | 0.625 / 8.3112e-5 | — |
| `ellip` | 2024 | 4.954e+16 / 8.6359e+12 | — |
| `ellipap` | 2012 | 29.688 / 2.4795 | — |
| `ellipord` | 2005 | 0.0 / 0.0 | — |
| `envelope` | 2072 | 6.0 / 0.48479 | — |
| `fftconvolve` | 2012 | 4.5036e+15 / 1.6821e+12 | — |
| `filtfilt` | 2008 | 0.0 / 0.0 | — |
| `find_peaks` | 2004 | 0.0 / 0.0 | — |
| `find_peaks_cwt` | 2008 | 0.0 / 0.0 | — |
| `findfreqs` | 2003 | 0.0 / 0.0 | — |
| `firls` | 2009 | 4313.1 / 258.56 | — |
| `firwin` | 2010 | 7.4555 / 1.7444 | — |
| `firwin2` | 2009 | 1.3324e+15 / 2.7324e+13 | — |
| `firwin_2d` | 2002 | 4.8098 / 1.1811 | — |
| `freqresp` | 2003 | 0.99656 / 0.14675 | — |
| `freqs` | 2002 | 1.9994 / 0.43274 | — |
| `freqs_zpk` | 2002 | 1.0165 / 0.17286 | — |
| `freqz` | 2002 | 1.7141 / 0.37285 | — |
| `freqz_sos` | 2002 | 6.6952 / 0.42815 | — |
| `freqz_zpk` | 2002 | 3.2497 / 0.48942 | — |
| `gammatone` | 2008 | 43.952 / 2.3769 | — |
| `gauss_spline` | 2009 | 0.60225 / 0.1313 | — |
| `gausspulse` | 2003 | 0.0 / 0.0 | — |
| `get_window` | 2006 | 10.577 / 0.98053 | — |
| `group_delay` | 2002 | 33.903 / 4.5991 | — |
| `hilbert` | 2003 | 4.5036e+15 / 1.1203e+12 | — |
| `hilbert2` | 2004 | 2.3685e+16 / 2.7049e+12 | — |
| `iircomb` | 2008 | 0.0 / 0.0 | — |
| `iirdesign` | 2020 | 10.876 / 3.6349 | — |
| `iirfilter` | 2006 | 1.2981e+16 / 8.0798e+11 | — |
| `iirnotch` | 2002 | 2.0423 / 1.6387e-3 | — |
| `iirpeak` | 2002 | 96.547 / 0.049031 | — |
| `impulse` | 2002 | 3.3738 / 0.27516 | — |
| `invres` | 2004 | 0.0 / 0.0 | — |
| `invresz` | 2004 | 0.0 / 0.0 | — |
| `istft` | 2008 | 1.0058e+17 / 2.1743e+15 | — |
| `kaiser_atten` | 2003 | 0.0 / 0.0 | — |
| `kaiser_beta` | 2006 | 0.0 / 0.0 | — |
| `kaiserord` | 2003 | 0.0 / 0.0 | — |
| `lfilter` | 2004 | 0.0 / 0.0 | — |
| `lfilter_with_state` | 2004 | 0.0 / 0.0 | — |
| `lfilter_zi` | 2004 | 0.0 / 0.0 | — |
| `lfiltic` | 2004 | 0.0 / 0.0 | — |
| `lombscargle` | 2004 | 64.822 / 8.0263 | — |
| `lp2bp` | 2002 | 7.0 / 2.3208e-3 | — |
| `lp2bp_zpk` | 2002 | 0.53333 / 0.056134 | — |
| `lp2bs` | 2002 | 14.0 / 0.010471 | — |
| `lp2bs_zpk` | 2002 | 0.90212 / 5.8197e-5 | — |
| `lp2hp` | 2002 | 14.0 / 7.3015e-3 | — |
| `lp2hp_zpk` | 2002 | 0.8 / 6.1477e-5 | — |
| `lp2lp` | 2002 | 14.0 / 3.8276e-3 | — |
| `lp2lp_zpk` | 2002 | 0.0 / 0.0 | — |
| `lsim` | 2008 | 19.703 / 0.24737 | — |
| `max_len_seq` | 2013 | 0.0 / 0.0 | — |
| `medfilt` | 2003 | 0.0 / 0.0 | — |
| `medfilt2d` | 2003 | 0.0 / 0.0 | — |
| `minimum_phase` | 2012 | 229.3 / 67.865 | — |
| `normalize` | 2003 | 0.0 / 0.0 | — |
| `oaconvolve` | 2012 | 9.0072e+15 / 7.6696e+13 | — |
| `order_filter` | 2009 | 0.0 / 0.0 | — |
| `peak_prominences` | 2006 | 0.0 / 0.0 | — |
| `peak_widths` | 2024 | 0.0 / 0.0 | — |
| `periodogram` | 2012 | 6.7554e+15 / 2.5913e+14 | — |
| `place_poles` | 2002 | 6.5062 / 1.0856 | — |
| `qspline1d` | 2004 | 1.8378 / 5.9731e-3 | — |
| `qspline1d_eval` | 2002 | 0.0 / 0.0 | — |
| `qspline2d` | 2002 | 12.471 / 0.03535 | — |
| `remez` | 2006 | 3.8424e+8 / 1.5651e+7 | — |
| `resample` | 2018 | 4.5036e+15 / 2.1767e+12 | — |
| `resample_poly` | 2010 | 34.154 / 4.0757e-3 | — |
| `residue` | 2008 | 19.0 / 8.2256e-3 | — |
| `residuez` | 2008 | 19.0 / 8.8041e-3 | — |
| `savgol_coeffs` | 2016 | 1.6413e+16 / 8.0553e+12 | — |
| `savgol_filter` | 2012 | 1.8014e+16 / 7.4612e+11 | — |
| `sawtooth` | 2024 | 2.2254 / 0.13805 | — |
| `sepfir2d` | 2002 | 0.0 / 0.0 | — |
| `sos2tf` | 2002 | 0.0 / 0.0 | — |
| `sos2zpk` | 2002 | 1.0 / 0.30555 | — |
| `sosfilt` | 2002 | 8.3429 / 0.26274 | — |
| `sosfilt_zi` | 2002 | 0.0 / 0.0 | — |
| `sosfiltfilt` | 2004 | 3.3461 / 0.33608 | — |
| `sosfreqz` | 2002 | 6.6952 / 0.42457 | — |
| `spectrogram` | 2012 | 5.5014 / 0.020721 | — |
| `spline_filter` | 2004 | 194.76 / 0.049308 | — |
| `square` | 2024 | 0.0 / 0.0 | — |
| `ss2tf` | 2004 | 4.5036e+15 / 1.1216e+12 | — |
| `ss2zpk` | 2007 | 4.5036e+15 / 7.4452e+11 | — |
| `step` | 2004 | 2.8941 / 0.22663 | — |
| `stft` | 2008 | 4.5036e+15 / 1.2082e+11 | — |
| `sweep_poly` | 2003 | 0.0 / 0.0 | — |
| `symiirorder1` | 2003 | 1.9875 / 3.6143e-4 | — |
| `symiirorder2` | 2006 | 961.31 / 0.85295 | — |
| `tf2sos` | 2004 | 94.688 / 0.31332 | — |
| `tf2ss` | 2003 | 0.0 / 0.0 | — |
| `tf2zpk` | 2003 | 11.0 / 4.1576e-3 | — |
| `unique_roots` | 2009 | 0.0 / 0.0 | — |
| `unit_impulse` | 2004 | 0.0 / 0.0 | — |
| `upfirdn` | 2005 | 0.0 / 0.0 | — |
| `vectorstrength` | 2003 | 9.5247e+15 / 2.3776e+12 | — |
| `welch` | 2012 | 6.0447 / 0.39974 | — |
| `whittaker_henderson` | 2018 | 51.0 / 0.029234 | — |
| `wiener` | 2005 | 0.0 / 0.0 | — |
| `windows.barthann` | 2008 | 0.0 / 0.0 | — |
| `windows.bartlett` | 2008 | 3.875 / 0.10895 | — |
| `windows.blackman` | 2008 | 49.978 / 1.3382 | — |
| `windows.blackmanharris` | 2008 | 261.57 / 2.7067 | — |
| `windows.bohman` | 2008 | 349.14 / 2.8582 | — |
| `windows.boxcar` | 2008 | 0.0 / 0.0 | — |
| `windows.chebwin` | 2008 | 176.17 / 5.6427 | — |
| `windows.cosine` | 2008 | 39.487 / 0.99122 | — |
| `windows.dpss` | 2006 | 8.3416e+17 / 1.6317e+15 | — |
| `windows.exponential` | 2008 | 0.0 / 0.0 | — |
| `windows.flattop` | 2008 | 283.86 / 3.2804 | — |
| `windows.gaussian` | 2008 | 32.026 / 3.6175 | — |
| `windows.general_cosine` | 2008 | 66.002 / 1.5512 | — |
| `windows.general_gaussian` | 2008 | 0.0 / 0.0 | — |
| `windows.general_hamming` | 2008 | 4.1065 / 0.30727 | — |
| `windows.get_window` | 2006 | 10.577 / 0.96974 | — |
| `windows.hamming` | 2008 | 9.9559 / 1.1893 | — |
| `windows.hann` | 2008 | 48.852 / 1.1563 | — |
| `windows.kaiser` | 2011 | 8.1057 / 0.9771 | — |
| `windows.kaiser_bessel_derived` | 2004 | 2.7043 / 0.22147 | — |
| `windows.lanczos` | 2008 | 9.3543 / 0.32059 | — |
| `windows.nuttall` | 2008 | 110.28 / 1.8934 | — |
| `windows.parzen` | 2008 | 5.2204 / 0.23339 | — |
| `windows.taylor` | 2008 | 0.98472 / 0.11782 | — |
| `windows.triang` | 2008 | 4.125 / 0.39496 | — |
| `windows.tukey` | 2012 | 48.629 / 0.75994 | — |
| `zoom_fft` | 2002 | 220.68 / 9.9508 | — |
| `zpk2sos` | 2004 | 110.68 / 1.0752 | — |
| `zpk2ss` | 2002 | 0.0 / 0.0 | — |
| `zpk2tf` | 2002 | 0.0 / 0.0 | — |


## Conventions and limitations

The pinned signal fixture adapters cover waveform, filtering, design, spectral, window, peak, spline, convolution and system functions through the public Sci facade. Output tuples flatten in the documented fixture order. Parameters, modes, shapes and normalization are retained in each function CSV; original real convolution/correlation also call NumPy.

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
