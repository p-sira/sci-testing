# stats: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-stats)

Measured binary64 (`f64`) results compared against the **saved per-case reference**.
Each implementation, including Sci-stats, is evaluated against the same
reference values stored in the raw CSVs. This page
covers 336398 cases and 166 public API names;
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

### `sci::stats::association`

9 functions · 18131 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | GSL | NumPy | statrs |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `chatterjeexi` | 2013 | 26.351 / 0.010368 | —; unavailable 2013/2013 | 0.0 / 0.0 | —; unavailable 2013/2013 | —; unavailable 2013/2013 | —; unavailable 2013/2013 | —; unavailable 2013/2013 |
| `kendalltau` | 2024 | 0.0 / 0.0 | —; unavailable 2024/2024 | 0.0 / 0.0 | —; unavailable 2024/2024 | —; unavailable 2024/2024 | —; unavailable 2024/2024 | —; unavailable 2024/2024 |
| `multiscale_graphcorr` | 2015 | 1.1259e+15 / 1.8279e+12 | —; unavailable 2015/2015 | 0.0 / 0.0 | —; unavailable 2015/2015 | —; unavailable 2015/2015 | —; unavailable 2015/2015 | —; unavailable 2015/2015 |
| `pearsonr` | 2012 | 18.9 / 3.6841 | —; unavailable 2012/2012 | 0.0 / 0.0 | 22.679 / 4.2061 | 4.5036e+15 / 2.2384e+12 | —; unavailable 2012/2012 | —; unavailable 2012/2012 |
| `pointbiserialr` | 2002 | 9.4721 / 2.0126 | —; unavailable 2002/2002 | 0.0 / 0.0 | 14.208 / 2.5395 | 10.064 / 2.4319 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `somersd` | 2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | —; unavailable 2012/2012 | —; unavailable 2012/2012 | —; unavailable 2012/2012 |
| `spearmanr` | 2009 | 16.313 / 4.0574 | —; unavailable 2009/2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | 4.5036e+15 / 2.2417e+12 | —; unavailable 2009/2009 | —; unavailable 2009/2009 |
| `spearmanrho` | 2012 | 3.6053 / 0.50006 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | 4.5036e+15 / 2.2384e+12 | —; unavailable 2012/2012 | —; unavailable 2012/2012 |
| `weightedtau` | 2032 | 6.7659 / 0.01214 | —; unavailable 2032/2032 | 0.0 / 0.0 | —; unavailable 2032/2032 | —; unavailable 2032/2032 | —; unavailable 2032/2032 | —; unavailable 2032/2032 |

### `sci::stats::binning`

5 functions · 10027 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | GSL | NumPy | statrs |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `binned_statistic` | 2014 | 0.71971 / 2.7511e-5 | —; unavailable 2014/2014 | 0.0 / 0.0 | —; unavailable 2014/2014 | —; unavailable 2014/2014 | —; unavailable 2014/2014 | —; unavailable 2014/2014 |
| `binned_statistic_2d` | 2007 | 0.0 / 0.0 | —; unavailable 2007/2007 | 0.0 / 0.0 | —; unavailable 2007/2007 | —; unavailable 2007/2007 | —; unavailable 2007/2007 | —; unavailable 2007/2007 |
| `binned_statistic_dd` | 2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `cumfreq` | 2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `relfreq` | 2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |

### `sci::stats::contingency`

10 functions · 20115 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | GSL | NumPy | statrs |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `barnard_exact` | 2054 | 5.2267e+5 / 9008.9 | —; unavailable 2054/2054 | 0.0 / 0.0 | —; unavailable 2054/2054 | —; unavailable 2054/2054 | —; unavailable 2054/2054 | —; unavailable 2054/2054 |
| `boschloo_exact` | 2027 | 5.2267e+5 / 25373.0 | —; unavailable 2027/2027 | 0.0 / 0.0 | —; unavailable 2027/2027 | —; unavailable 2027/2027 | —; unavailable 2027/2027 | —; unavailable 2027/2027 |
| `chi2_contingency` | 2006 | 99.866 / 1.6943 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `contingency.association` | 2006 | 1.0572 / 0.48521 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `contingency.crosstab` | 2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `contingency.expected_freq` | 2003 | 0.92857 / 0.13877 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `contingency.margins` | 2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `contingency.odds_ratio` | 2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `contingency.relative_risk` | 2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `fisher_exact` | 2009 | 1192.0 / 170.67 | —; unavailable 2009/2009 | 0.0 / 0.0 | 1.0072 / 0.324 | —; unavailable 2009/2009 | —; unavailable 2009/2009 | —; unavailable 2009/2009 |

### `sci::stats::descriptive`

42 functions · 84446 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | GSL | NumPy | statrs |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `circmean` | 2003 | 2.6782 / 0.56516 | 14.319 / 0.31836 | 0.0 / 0.0 | 2.1458 / 0.4252 | 5.5074 / 0.50608 | 1.0731 / 0.23548 | —; unavailable 2003/2003 |
| `circstd` | 2006 | 97.198 / 21.938 | 69.622 / 17.273 | 0.0 / 0.0 | 97.189 / 11.734 | 65.571 / 15.339 | 97.14 / 23.615 | —; unavailable 2006/2006 |
| `circvar` | 2003 | 193.53 / 43.992 | 122.14 / 34.491 | 0.0 / 0.0 | 193.54 / 23.636 | 129.08 / 30.42 | 193.54 / 45.856 | —; unavailable 2003/2003 |
| `describe` | 2004 | 28.797 / 1.3885 | 22.627 / 1.1969 | 0.0 / 0.0 | 25.369 / 1.2738 | 29.483 / 1.118 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `differential_entropy` | 2008 | 1.6289 / 0.33512 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `directional_stats` | 2004 | 1.202 / 0.1994 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `entropy` | 2006 | 3.7856 / 5.1559e-3 | 0.94639 / 1.1537e-3 | 0.0 / 0.0 | —; unavailable 2006/2006 | 0.0 / 0.0 | 0.0 / 0.0 | —; unavailable 2006/2006 |
| `expectile` | 2015 | 0.92308 / 4.581e-4 | —; unavailable 2015/2015 | 0.0 / 0.0 | —; unavailable 2015/2015 | 0.78788 / 6.68e-4 | —; unavailable 2015/2015 | —; unavailable 2015/2015 |
| `gmean` | 2003 | 1.4156 / 0.083665 | 0.70781 / 0.27006 | 0.0 / 0.0 | 0.0 / 0.0 | 1.4156 / 0.083665 | 0.0 / 0.0 | —; unavailable 2003/2003 |
| `gstd` | 2006 | 0.92144 / 0.053096 | 0.92144 / 0.77986 | 0.0 / 0.0 | 0.92144 / 0.29108 | 0.0 / 0.0 | 0.0 / 0.0 | —; unavailable 2006/2006 |
| `gzscore` | 2004 | 3.9131 / 0.81393 | 2.2361 / 0.57927 | 0.0 / 0.0 | —; unavailable 2004/2004 | 2.7951 / 0.21289 | 0.0 / 0.0 | —; unavailable 2004/2004 |
| `hmean` | 2003 | 0.93844 / 0.26447 | 0.93844 / 0.26447 | 0.0 / 0.0 | 0.0 / 0.0 | 0.93843 / 0.1638 | 0.0 / 0.0 | —; unavailable 2003/2003 |
| `iqr` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | 0.0 / 0.0 | —; unavailable 2004/2004 |
| `kstat` | 2008 | 8.9074 / 7.0074e-3 | 1.0 / 0.18398 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `kstatvar` | 2004 | 3.9016 / 0.67103 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `kurtosis` | 2010 | 283.25 / 0.79221 | 50.896 / 0.66526 | 0.0 / 0.0 | 283.25 / 0.93819 | 283.25 / 1.2134 | —; unavailable 2010/2010 | —; unavailable 2010/2010 |
| `lmoment` | 2020 | 26.091 / 0.40889 | —; unavailable 2020/2020 | 0.0 / 0.0 | —; unavailable 2020/2020 | —; unavailable 2020/2020 | —; unavailable 2020/2020 | —; unavailable 2020/2020 |
| `mean` | 2004 | 0.83416 / 0.13556 | 0.83417 / 0.15761 | 0.0 / 0.0 | 0.83417 / 0.168 | 0.83417 / 0.15761 | 0.0 / 0.0 | —; unavailable 2004/2004 |
| `median` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | —; unavailable 2004/2004 |
| `median_abs_deviation` | 2008 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | —; unavailable 2008/2008 |
| `mode` | 2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | —; unavailable 2004/2004 |
| `moment` | 2024 | 1.9207 / 3.348e-3 | 4.5036e+15 / 2.2384e+12; unavailable 12/2024 | 0.0 / 0.0 | —; unavailable 2024/2024 | 0.93238 / 9.1311e-4 | 0.0 / 0.0 | —; unavailable 2024/2024 |
| `pmean` | 2018 | 1.0667 / 0.16185 | 0.70711 / 0.14656 | 0.0 / 0.0 | —; unavailable 2018/2018 | 0.65924 / 0.18858 | —; unavailable 2018/2018 | —; unavailable 2018/2018 |
| `quantile` | 2020 | 1.6667 / 1.3751e-3 | —; unavailable 2020/2020 | 0.0 / 0.0 | —; unavailable 2020/2020 | 1.6667 / 1.3751e-3 | 1.6667 / 1.3751e-3 | —; unavailable 2020/2020 |
| `scoreatpercentile` | 2020 | 0.0 / 0.0 | —; unavailable 2020/2020 | 0.0 / 0.0 | —; unavailable 2020/2020 | 0.0 / 0.0 | 0.0 / 0.0 | —; unavailable 2020/2020 |
| `sem` | 2006 | 1.4042 / 0.26019 | 1.4044 / 0.2903 | 0.0 / 0.0 | 1.4044 / 0.22347 | 1.4044 / 0.25494 | 0.0 / 0.0 | —; unavailable 2006/2006 |
| `sigmaclip` | 2006 | 2.4997 / 0.11594 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `skew` | 2006 | 28.797 / 8.445 | 22.627 / 7.4002 | 0.0 / 0.0 | 25.369 / 7.6733 | 29.483 / 6.1775 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `std` | 2007 | 1.2558 / 0.26731 | 0.62808 / 0.18196 | 0.0 / 0.0 | 1.256 / 0.22728 | 1.2562 / 0.26731 | 0.0 / 0.0 | —; unavailable 2007/2007 |
| `tmax` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `tmean` | 2016 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | —; unavailable 2016/2016 | —; unavailable 2016/2016 |
| `tmin` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `trim1` | 2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | —; unavailable 2012/2012 | —; unavailable 2012/2012 | —; unavailable 2012/2012 |
| `trim_mean` | 2016 | 0.83416 / 0.14467 | 0.83417 / 0.16122 | 0.0 / 0.0 | —; unavailable 2016/2016 | 0.83416 / 0.15456 | 0.0 / 0.0 | —; unavailable 2016/2016 |
| `trimboth` | 2016 | 0.0 / 0.0 | —; unavailable 2016/2016 | 0.0 / 0.0 | —; unavailable 2016/2016 | —; unavailable 2016/2016 | 0.0 / 0.0 | —; unavailable 2016/2016 |
| `tsem` | 2015 | 1.0 / 0.15831 | 1.0 / 0.26968 | 0.0 / 0.0 | 1.7329 / 0.27398 | 1.7335 / 0.31784 | —; unavailable 2015/2015 | —; unavailable 2015/2015 |
| `tstd` | 2015 | 0.7071 / 0.20201 | 0.7071 / 0.10205 | 0.0 / 0.0 | 1.0009 / 0.050647 | 0.5005 / 0.063804 | —; unavailable 2015/2015 | —; unavailable 2015/2015 |
| `tvar` | 2015 | 0.6 / 0.07645 | 0.6 / 0.07645 | 0.0 / 0.0 | 2.0038 / 0.10109 | 1.002 / 0.1222 | —; unavailable 2015/2015 | —; unavailable 2015/2015 |
| `variance` | 2079 | 2.3139 / 0.36574 | 0.48267 / 0.18462 | 2.4628e+8 / 3.0635e+5 | 6.1915e+11 / 5.1229e+8 | 6.1571e+7 / 62244.0 | 2.4628e+8 / 3.0635e+5 | —; unavailable 2079/2079 |
| `variation` | 2007 | 1.5059 / 0.3672 | 1.5059 / 0.32171 | 0.0 / 0.0 | 2.2588 / 0.30755 | 2.2588 / 0.37433 | 0.58482 / 2.9139e-4 | —; unavailable 2007/2007 |
| `zmap` | 2005 | 1.0 / 4.2184e-4 | 0.95257 / 0.060441 | 0.0 / 0.0 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | —; unavailable 2005/2005 |
| `zscore` | 2004 | 7.9687 / 0.85237 | 4.9804 / 0.66786 | 0.0 / 0.0 | 5.9765 / 0.53594 | 5.9765 / 0.58614 | 0.0 / 0.0 | —; unavailable 2004/2004 |

### `sci::stats::diagnostics`

7 functions · 14072 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | GSL | NumPy | statrs |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `boxcox_normmax` | 2016 | 5.0407e+23 / 2.0359e+21 | —; unavailable 2016/2016 | 0.0 / 0.0 | —; unavailable 2016/2016 | —; unavailable 2016/2016 | —; unavailable 2016/2016 | —; unavailable 2016/2016 |
| `boxcox_normplot` | 2006 | 2.2429 / 0.24871 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `ppcc_max` | 2015 | 1.786e+9 / 4.311e+8 | —; unavailable 2015/2015 | 0.0 / 0.0 | —; unavailable 2015/2015 | —; unavailable 2015/2015 | —; unavailable 2015/2015 | —; unavailable 2015/2015 |
| `ppcc_plot` | 2015 | 13.208 / 0.29083 | —; unavailable 2015/2015 | 0.0 / 0.0 | —; unavailable 2015/2015 | —; unavailable 2015/2015 | —; unavailable 2015/2015 | —; unavailable 2015/2015 |
| `probplot` | 2003 | 1.2205 / 6.0934e-5 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `yeojohnson_normmax` | 2008 | 2.2946e+8 / 4.7734e+7 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `yeojohnson_normplot` | 2009 | 2.6227 / 0.22653 | —; unavailable 2009/2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | —; unavailable 2009/2009 | —; unavailable 2009/2009 | —; unavailable 2009/2009 |

### `sci::stats::distances`

3 functions · 6016 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | GSL | NumPy | statrs |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `energy_distance` | 2005 | 1.1558 / 0.3196 | 1.1558 / 0.3263 | 0.0 / 0.0 | —; unavailable 2005/2005 | 1.1558 / 0.3196 | 2.8891 / 0.92232 | —; unavailable 2005/2005 |
| `wasserstein_distance` | 2005 | 0.0 / 0.0 | 0.99999 / 0.11744 | 0.0 / 0.0 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | —; unavailable 2005/2005 |
| `wasserstein_distance_nd` | 2006 | 1.4142 / 0.67917 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |

### `sci::stats::distributions`

12 functions · 25995 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | GSL | NumPy | statrs |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `cdf` | 2379 | 4.2666e+5 / 8782.7 | 2031.9 / 1.2188 | 0.0 / 0.0 | 1499.9 / 1.2083; unavailable 290/2379 | 1392.5 / 2.9578; unavailable 363/2379 | —; unavailable 2379/2379 | 2.2078e+5 / 3346.8; unavailable 587/2379 |
| `fit` | 2036 | 2.8955 / 2.0096e-3 | —; unavailable 2036/2036 | 0.0 / 0.0 | —; unavailable 2036/2036 | —; unavailable 2036/2036 | —; unavailable 2036/2036 | —; unavailable 2036/2036 |
| `logpdf` | 2163 | 832.64 / 10.889 | 76.599 / 1.1185; unavailable 7/2163 | 0.0 / 0.0 | 73.615 / 1.1669; unavailable 336/2163 | 82.601 / 3.1892; unavailable 387/2163 | —; unavailable 2163/2163 | 832.64 / 11.352; unavailable 374/2163 |
| `make_distribution` | 2015 | 5.6207 / 0.35452 | —; unavailable 2015/2015 | 0.0 / 0.0 | —; unavailable 2015/2015 | —; unavailable 2015/2015 | —; unavailable 2015/2015 | —; unavailable 2015/2015 |
| `normal_cdf` | 2111 | 940.89 / 110.9 | 0.49075 / 0.15686 | 940.89 / 110.9 | 844.83 / 109.78 | 3.1714 / 0.39997 | —; unavailable 2111/2111 | 4.6178e+5 / 1.3474e+5 |
| `normal_pdf` | 2000 | 1.5723 / 0.78551 | 8.4529 / 1.8876 | 0.0 / 0.0 | 0.0 / 0.0 | 0.99061 / 0.15999 | 0.0 / 0.0 | 1.5723 / 0.78551 |
| `normal_ppf` | 2000 | 0.0 / 0.0 | 2.785 / 0.41242 | 0.0 / 0.0 | 2.785 / 0.39778 | 3.7134 / 0.67276 | —; unavailable 2000/2000 | 3.7134 / 0.52153 |
| `order_statistic` | 2160 | 2107.7 / 3.1201 | —; unavailable 2160/2160 | 0.0 / 0.0 | —; unavailable 2160/2160 | —; unavailable 2160/2160 | —; unavailable 2160/2160 | —; unavailable 2160/2160 |
| `pdf` | 2377 | 490.26 / 9.7963 | 17.636 / 0.63913 | 0.0 / 0.0 | 18.616 / 0.62024; unavailable 291/2377 | 84.153 / 3.9631; unavailable 360/2377 | —; unavailable 2377/2377 | 490.26 / 9.274; unavailable 530/2377 |
| `ppf` | 2300 | 4.6057e+6 / 6617.6 | 2.9826e+7 / 13051.0; unavailable 14/2300 | 0.0 / 0.0 | 2464.5 / 3.9201; unavailable 448/2300 | 1.5289e+9 / 1.7337e+6; failed 1/2300; unavailable 1065/2300 | —; unavailable 2300/2300 | 8.1955e+14 / 2.6974e+12; unavailable 477/2300 |
| `sf` | 2379 | 2.9859e+5 / 4114.0 | 33.012 / 0.38789 | 0.0 / 0.0 | 20.88 / 0.35776; unavailable 289/2379 | 69.361 / 1.7378; unavailable 357/2379 | —; unavailable 2379/2379 | 3.639e+5 / 5135.7; unavailable 598/2379 |
| `truncate` | 2075 | 4.5036e+15 / 2.6045e+12 | —; unavailable 2075/2075 | 0.0 / 0.0 | —; unavailable 2075/2075 | —; unavailable 2075/2075 | —; unavailable 2075/2075 | —; unavailable 2075/2075 |

### `sci::stats::estimation`

2 functions · 4030 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | GSL | NumPy | statrs |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `bayes_mvs` | 2012 | 83.096 / 3.6539 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | —; unavailable 2012/2012 | —; unavailable 2012/2012 | —; unavailable 2012/2012 |
| `mvsdist` | 2018 | 516.22 / 1.083 | —; unavailable 2018/2018 | 0.0 / 0.0 | —; unavailable 2018/2018 | —; unavailable 2018/2018 | —; unavailable 2018/2018 | —; unavailable 2018/2018 |

### `sci::stats::goodness`

10 functions · 20137 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | GSL | NumPy | statrs |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `anderson` | 2011 | 99.687 / 0.4487 | —; unavailable 2011/2011 | 0.0 / 0.0 | —; unavailable 2011/2011 | —; unavailable 2011/2011 | —; unavailable 2011/2011 | —; unavailable 2011/2011 |
| `anderson_ksamp` | 2014 | 6.1166 / 2.0784 | —; unavailable 2014/2014 | 0.0 / 0.0 | —; unavailable 2014/2014 | —; unavailable 2014/2014 | —; unavailable 2014/2014 | —; unavailable 2014/2014 |
| `chisquare` | 2006 | 5.2329 / 1.0128 | —; unavailable 2006/2006 | 0.0 / 0.0 | 3.4945 / 0.38534 | 4.9907 / 0.15494 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `cramervonmises` | 2014 | 2.6969e+5 / 67.074 | —; unavailable 2014/2014 | 0.0 / 0.0 | —; unavailable 2014/2014 | —; unavailable 2014/2014 | —; unavailable 2014/2014 | —; unavailable 2014/2014 |
| `cramervonmises_2samp` | 2013 | 1.2769e+5 / 31.818 | —; unavailable 2013/2013 | 0.0 / 0.0 | —; unavailable 2013/2013 | —; unavailable 2013/2013 | —; unavailable 2013/2013 | —; unavailable 2013/2013 |
| `goodness_of_fit` | 2012 | 4.7452 / 0.14129 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | —; unavailable 2012/2012 | —; unavailable 2012/2012 | —; unavailable 2012/2012 |
| `ks_1samp` | 2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | 0.0 / 0.0 | —; unavailable 2009/2009 | —; unavailable 2009/2009 | —; unavailable 2009/2009 | —; unavailable 2009/2009 |
| `ks_2samp` | 2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | —; unavailable 2012/2012 | —; unavailable 2012/2012 | —; unavailable 2012/2012 |
| `kstest` | 2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | —; unavailable 2012/2012 | —; unavailable 2012/2012 | —; unavailable 2012/2012 |
| `power_divergence` | 2034 | 34.056 / 0.81664 | —; unavailable 2034/2034 | 0.0 / 0.0 | 10.586 / 0.56997 | —; unavailable 2034/2034 | —; unavailable 2034/2034 | —; unavailable 2034/2034 |

### `sci::stats::inference`

27 functions · 54625 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | GSL | NumPy | statrs |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `ansari` | 2015 | 0.0 / 0.0 | —; unavailable 2015/2015 | 0.0 / 0.0 | —; unavailable 2015/2015 | —; unavailable 2015/2015 | —; unavailable 2015/2015 | —; unavailable 2015/2015 |
| `bartlett` | 2003 | 14.111 / 3.3287 | —; unavailable 2003/2003 | 0.0 / 0.0 | 19.007 / 2.692 | 19.007 / 2.5758 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `binomtest` | 2048 | 138.29 / 3.8558 | —; unavailable 2048/2048 | 0.0 / 0.0 | 2.4031 / 0.080525 | —; unavailable 2048/2048 | —; unavailable 2048/2048 | —; unavailable 2048/2048 |
| `brunnermunzel` | 2030 | 87.763 / 8.0938 | —; unavailable 2030/2030 | 0.0 / 0.0 | —; unavailable 2030/2030 | —; unavailable 2030/2030 | —; unavailable 2030/2030 | —; unavailable 2030/2030 |
| `combine_pvalues` | 2011 | 1.9401e+5 / 100.03 | —; unavailable 2011/2011 | 0.0 / 0.0 | —; unavailable 2011/2011 | —; unavailable 2011/2011 | —; unavailable 2011/2011 | —; unavailable 2011/2011 |
| `dunnett` | 2015 | 1.3013e+10 / 7.1396e+6 | —; unavailable 2015/2015 | 0.0 / 0.0 | —; unavailable 2015/2015 | —; unavailable 2015/2015 | —; unavailable 2015/2015 | —; unavailable 2015/2015 |
| `epps_singleton_2samp` | 2020 | 5.5464e+6 / 36092.0 | —; unavailable 2020/2020 | 0.0 / 0.0 | —; unavailable 2020/2020 | —; unavailable 2020/2020 | —; unavailable 2020/2020 | —; unavailable 2020/2020 |
| `f_oneway` | 2003 | 36.246 / 11.131 | —; unavailable 2003/2003 | 0.0 / 0.0 | 8.9747 / 1.6891 | 13.774 / 4.0099 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `false_discovery_control` | 2006 | 0.93844 / 0.062734 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `fligner` | 2008 | 48.635 / 20.564 | —; unavailable 2008/2008 | 0.0 / 0.0 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 | —; unavailable 2008/2008 |
| `friedmanchisquare` | 2002 | 2.1827 / 1.0911 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `kruskal` | 2003 | 31.544 / 0.98758 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `levene` | 2002 | 7.7383 / 2.3693 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `mannwhitneyu` | 2050 | 0.73669 / 1.7968e-4 | —; unavailable 2050/2050 | 0.0 / 0.0 | —; unavailable 2050/2050 | —; unavailable 2050/2050 | —; unavailable 2050/2050 | —; unavailable 2050/2050 |
| `median_test` | 2022 | 17.422 / 5.0827e-3 | —; unavailable 2022/2022 | 0.0 / 0.0 | —; unavailable 2022/2022 | —; unavailable 2022/2022 | —; unavailable 2022/2022 | —; unavailable 2022/2022 |
| `mood` | 2021 | 0.0 / 0.0 | —; unavailable 2021/2021 | 0.0 / 0.0 | —; unavailable 2021/2021 | —; unavailable 2021/2021 | —; unavailable 2021/2021 | —; unavailable 2021/2021 |
| `page_trend_test` | 2010 | 0.0 / 0.0 | —; unavailable 2010/2010 | 0.0 / 0.0 | —; unavailable 2010/2010 | —; unavailable 2010/2010 | —; unavailable 2010/2010 | —; unavailable 2010/2010 |
| `poisson_means_test` | 2030 | 3.9349e+5 / 6892.3 | —; unavailable 2030/2030 | 0.0 / 0.0 | —; unavailable 2030/2030 | —; unavailable 2030/2030 | —; unavailable 2030/2030 | —; unavailable 2030/2030 |
| `power` | 2020 | 14.694 / 7.1974 | —; unavailable 2020/2020 | 0.0 / 0.0 | —; unavailable 2020/2020 | —; unavailable 2020/2020 | —; unavailable 2020/2020 | —; unavailable 2020/2020 |
| `quantile_test` | 2108 | 42.896 / 0.23919 | —; unavailable 2108/2108 | 0.0 / 0.0 | —; unavailable 2108/2108 | —; unavailable 2108/2108 | —; unavailable 2108/2108 | —; unavailable 2108/2108 |
| `ranksums` | 2015 | 4.1473e+5 / 1.0352e+5 | —; unavailable 2015/2015 | 0.0 / 0.0 | —; unavailable 2015/2015 | —; unavailable 2015/2015 | —; unavailable 2015/2015 | —; unavailable 2015/2015 |
| `ttest_1samp` | 2012 | 17.687 / 3.3458 | —; unavailable 2012/2012 | 0.0 / 0.0 | 7.5744 / 0.80021 | 14.648 / 2.3606 | —; unavailable 2012/2012 | —; unavailable 2012/2012 |
| `ttest_ind` | 2030 | 63.924 / 5.5843 | —; unavailable 2030/2030 | 0.0 / 0.0 | 6.0414 / 0.50132 | 20.152 / 2.0879 | —; unavailable 2030/2030 | —; unavailable 2030/2030 |
| `ttest_ind_from_stats` | 2002 | 554.04 / 50.071 | —; unavailable 2002/2002 | 0.0 / 0.0 | 9.0391 / 0.68986 | 35.094 / 2.8201 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `ttest_rel` | 2015 | 21.796 / 2.3601 | —; unavailable 2015/2015 | 0.0 / 0.0 | 10.929 / 1.057 | 17.854 / 1.6086 | —; unavailable 2015/2015 | —; unavailable 2015/2015 |
| `tukey_hsd` | 2010 | 3.1942e+15 / 3.8324e+11 | —; unavailable 2010/2010 | 0.0 / 0.0 | —; unavailable 2010/2010 | —; unavailable 2010/2010 | —; unavailable 2010/2010 | —; unavailable 2010/2010 |
| `wilcoxon` | 2114 | 0.52036 / 1.2307e-4 | —; unavailable 2114/2114 | 0.0 / 0.0 | —; unavailable 2114/2114 | —; unavailable 2114/2114 | —; unavailable 2114/2114 | —; unavailable 2114/2114 |

### `sci::stats::multivariate`

10 functions · 20083 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | GSL | NumPy | statrs |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `categorical_cdf` | 2024 | 0.0 / 0.0 | —; unavailable 2024/2024 | 0.0 / 0.0 | —; unavailable 2024/2024 | —; unavailable 2024/2024 | —; unavailable 2024/2024 | —; unavailable 2024/2024 |
| `categorical_pmf` | 2024 | 0.0 / 0.0 | —; unavailable 2024/2024 | 0.0 / 0.0 | —; unavailable 2024/2024 | —; unavailable 2024/2024 | —; unavailable 2024/2024 | —; unavailable 2024/2024 |
| `dirichlet_logpdf` | 2003 | 4.5036e+15 / 2.2484e+12 | 2.6432 / 0.70186 | 0.0 / 0.0 | —; unavailable 2003/2003 | 10.573 / 10.562 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `dirichlet_pdf` | 2003 | 17.983 / 4.0248 | 1.7008 / 0.42734 | 0.0 / 0.0 | —; unavailable 2003/2003 | 6.2365 / 5.9934 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |
| `multinomial_logpmf` | 2006 | 2.9398e+5 / 147.45 | 1.5731e+5 / 79.514 | 0.0 / 0.0 | —; unavailable 2006/2006 | 2056.5 / 1.6472 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `multinomial_pmf` | 2007 | 2.0972e+6 / 1048.6 | 1.1222e+6 / 564.27 | 0.0 / 0.0 | —; unavailable 2007/2007 | 24576.0 / 12.849 | —; unavailable 2007/2007 | —; unavailable 2007/2007 |
| `multivariate_normal_logpdf` | 2002 | 1.4343 / 0.64566 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | 1.4341 / 0.27413 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `multivariate_normal_pdf` | 2002 | 3.5594 / 1.001 | —; unavailable 2002/2002 | 0.0 / 0.0 | —; unavailable 2002/2002 | 8.0098 / 0.73923 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `multivariate_t_logpdf` | 2006 | 27.701 / 2.4638 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `multivariate_t_pdf` | 2006 | 200.81 / 8.5053 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |

### `sci::stats::normality`

7 functions · 14077 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | GSL | NumPy | statrs |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `alexandergovern` | 2005 | 88.121 / 4.8069 | —; unavailable 2005/2005 | 0.0 / 0.0 | —; unavailable 2005/2005 | —; unavailable 2005/2005 | —; unavailable 2005/2005 | —; unavailable 2005/2005 |
| `jarque_bera` | 2002 | 5.7494 / 1.0431 | —; unavailable 2002/2002 | 0.0 / 0.0 | 7.392 / 1.1445 | 7.392 / 0.96753 | —; unavailable 2002/2002 | —; unavailable 2002/2002 |
| `kurtosistest` | 2021 | 576.73 / 2.4259 | —; unavailable 2021/2021 | 0.0 / 0.0 | —; unavailable 2021/2021 | —; unavailable 2021/2021 | —; unavailable 2021/2021 | —; unavailable 2021/2021 |
| `normaltest` | 2006 | 6.78e+15 / 1.2996e+15 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `obrientransform` | 2004 | 21.538 / 0.72203 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `shapiro` | 2021 | 6.0048e+15 / 1.4856e+12 | —; unavailable 2021/2021 | 0.0 / 0.0 | —; unavailable 2021/2021 | —; unavailable 2021/2021 | —; unavailable 2021/2021 | —; unavailable 2021/2021 |
| `skewtest` | 2018 | 2.7302e+31 / 3.5875e+30 | —; unavailable 2018/2018 | 0.0 / 0.0 | —; unavailable 2018/2018 | —; unavailable 2018/2018 | —; unavailable 2018/2018 | —; unavailable 2018/2018 |

### `sci::stats::qmc`

4 functions · 8039 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | GSL | NumPy | statrs |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `qmc.discrepancy` | 2024 | 37.775 / 2.7787 | —; unavailable 2024/2024 | 0.0 / 0.0 | —; unavailable 2024/2024 | —; unavailable 2024/2024 | —; unavailable 2024/2024 | —; unavailable 2024/2024 |
| `qmc.geometric_discrepancy` | 2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `qmc.scale` | 2006 | 2.6666 / 0.537 | —; unavailable 2006/2006 | 0.0 / 0.0 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 | —; unavailable 2006/2006 |
| `qmc.update_discrepancy` | 2003 | 105.64 / 29.38 | —; unavailable 2003/2003 | 0.0 / 0.0 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 | —; unavailable 2003/2003 |

### `sci::stats::ranking`

5 functions · 10332 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | GSL | NumPy | statrs |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `ecdf` | 2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |
| `estimated_cdf` | 2224 | 0.0 / 0.0 | —; unavailable 2224/2224 | 0.0 / 0.0 | —; unavailable 2224/2224 | —; unavailable 2224/2224 | —; unavailable 2224/2224 | —; unavailable 2224/2224 |
| `percentileofscore` | 2080 | 0.0 / 0.0 | —; unavailable 2080/2080 | 0.0 / 0.0 | —; unavailable 2080/2080 | —; unavailable 2080/2080 | 0.0 / 0.0 | —; unavailable 2080/2080 |
| `rankdata` | 2020 | 0.0 / 0.0 | 0.0 / 0.0 | 0.0 / 0.0 | —; unavailable 2020/2020 | 0.0 / 0.0 | 0.0 / 0.0 | —; unavailable 2020/2020 |
| `tiecorrect` | 2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | 0.0 / 0.0 | —; unavailable 2004/2004 | —; unavailable 2004/2004 |

### `sci::stats::regression`

3 functions · 6071 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | GSL | NumPy | statrs |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `linregress` | 2009 | 4.5036e+15 / 7.4724e+11 | —; unavailable 2009/2009 | 0.0 / 0.0 | 4.5036e+15 / 7.4724e+11 | 4.5036e+15 / 7.4724e+11 | —; unavailable 2009/2009 | —; unavailable 2009/2009 |
| `siegelslopes` | 2014 | 0.5 / 0.015268 | —; unavailable 2014/2014 | 0.0 / 0.0 | —; unavailable 2014/2014 | —; unavailable 2014/2014 | —; unavailable 2014/2014 | —; unavailable 2014/2014 |
| `theilslopes` | 2048 | 0.5 / 0.015808 | —; unavailable 2048/2048 | 0.0 / 0.0 | —; unavailable 2048/2048 | —; unavailable 2048/2048 | —; unavailable 2048/2048 | —; unavailable 2048/2048 |

### `sci::stats::resampling`

4 functions · 8111 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | GSL | NumPy | statrs |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `bootstrap` | 2036 | 2.0385 / 0.096895 | —; unavailable 2036/2036 | 0.0 / 0.0 | —; unavailable 2036/2036 | —; unavailable 2036/2036 | —; unavailable 2036/2036 | —; unavailable 2036/2036 |
| `bws_test` | 2012 | 25.412 / 7.2321e-3 | —; unavailable 2012/2012 | 0.0 / 0.0 | —; unavailable 2012/2012 | —; unavailable 2012/2012 | —; unavailable 2012/2012 | —; unavailable 2012/2012 |
| `monte_carlo_test` | 2027 | 0.0 / 0.0 | —; unavailable 2027/2027 | 0.0 / 0.0 | —; unavailable 2027/2027 | —; unavailable 2027/2027 | —; unavailable 2027/2027 | —; unavailable 2027/2027 |
| `permutation_test` | 2036 | 0.0 / 0.0 | —; unavailable 2036/2036 | 0.0 / 0.0 | —; unavailable 2036/2036 | —; unavailable 2036/2036 | —; unavailable 2036/2036 | —; unavailable 2036/2036 |

### `sci::stats::sensitivity`

1 function · 2007 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | GSL | NumPy | statrs |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `sobol_indices` | 2007 | 9.0076e+15 / 1.683e+12 | —; unavailable 2007/2007 | 0.0 / 0.0 | —; unavailable 2007/2007 | —; unavailable 2007/2007 | —; unavailable 2007/2007 | —; unavailable 2007/2007 |

### `sci::stats::survival`

1 function · 2024 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | GSL | NumPy | statrs |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `logrank` | 2024 | 3.9598 / 4.3075e-3 | —; unavailable 2024/2024 | 0.0 / 0.0 | —; unavailable 2024/2024 | —; unavailable 2024/2024 | —; unavailable 2024/2024 | —; unavailable 2024/2024 |

### `sci::stats::transformations`

4 functions · 8060 cases

| Function | Cases | Sci | mpmath | SciPy | Boost.Math | GSL | NumPy | statrs |
| --- | ---: | --- | --- | --- | --- | --- | --- | --- |
| `boxcox` | 2015 | 0.0 / 0.0 | —; unavailable 2015/2015 | 0.0 / 0.0 | —; unavailable 2015/2015 | —; unavailable 2015/2015 | —; unavailable 2015/2015 | —; unavailable 2015/2015 |
| `boxcox_llf` | 2010 | 1.3356 / 0.33411 | —; unavailable 2010/2010 | 0.0 / 0.0 | —; unavailable 2010/2010 | 1.3519 / 0.3068 | —; unavailable 2010/2010 | —; unavailable 2010/2010 |
| `yeojohnson` | 2021 | 0.0 / 0.0 | —; unavailable 2021/2021 | 0.0 / 0.0 | —; unavailable 2021/2021 | —; unavailable 2021/2021 | —; unavailable 2021/2021 | —; unavailable 2021/2021 |
| `yeojohnson_llf` | 2014 | 1.9043 / 0.11034 | —; unavailable 2014/2014 | 0.0 / 0.0 | —; unavailable 2014/2014 | —; unavailable 2014/2014 | —; unavailable 2014/2014 | —; unavailable 2014/2014 |

## Evaluation failures

- `stats.ppf` (GSL): 1 failed; distribution returned NaN. Example: `stats.ppf/recipe-5112`.


## Conventions and limitations

Distribution parameters, tail conventions, alternatives, ddof, explicit resampling draws and output ordering follow the pinned statistics fixtures. SciPy supplies binary64 references for added APIs, including documented adaptations and independently seeded probabilistic references. Weighted tau compares its statistic; its always-NaN p-value is left to the ordinary API tests. Normality, inverse and joint-tail approximations have reference and algorithm accuracy limits. Public mpmath calculations and the independent Rust statrs distribution routines provide further numerical comparisons on supported fixtures. Finite results do not imply identical tail precision or work limits; the full results and errors remain in the saved measurements. Additional GSL-labelled entropy, ranks, tie corrections, mode, Wasserstein/energy distances, Box-Cox likelihood, z-map and expectiles are derived reductions using public GSL sorting, statistics or special functions. These are not identically named exports from GSL. Additional Boost.Math-derived statistics use its public sample-moment primitives and chi-square, binomial and hypergeometric distributions; alternatives, ties, degrees of freedom and paired outputs follow the specific saved Sci contracts. GSL-derived statistics combine public GSL sample moments with t, F and chi-square distribution tails and explicit circular/logarithmic transformations. These are mathematical compositions with independent SciPy validation, not claims that the libraries export identical named Sci statistical functions.

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
