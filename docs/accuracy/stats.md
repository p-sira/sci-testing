# stats: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-stats)

Measured binary64 (`f64`) results using **Sci-stats as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 342414 cases and 169 public API names;
coverage is limited to the operations and inputs below.

## Agreement with Sci

`Sci vs reference` is **maximum / mean difference between Sci and the saved
independent mathematical reference in ε**. Other library cells are **maximum /
mean difference from Sci in ε** over cases
where both evaluations succeed. `Cases` counts input cases; `Sci failed N/T`
identifies missing baseline evaluations. `failed N/T` counts that library's own
failures; `Sci unavailable N/T` counts its successful evaluations excluded
because Sci failed. A dash means no saved comparison or no paired finite results.
Means include all compared scalar components, including exact zeros.
Full outputs, references, warnings and errors remain in the per-function CSVs.
Tables are grouped by Sci 0.3 public submodule using the compiler-resolved
canonical API inventory. Historical flat recipe aliases appear under their
current public namespace, with their measured recipe names preserved.

### `sci::stats::association` (9 functions, 18131 cases)

| Function | Cases | Sci vs reference | SciPy | Boost.Math | GSL | NumPy | stdlib |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `chatterjeexi` | 2013 | 26.351 / 0.010368 | 26.351 / 0.010368 | — | — | — | — |
| `kendalltau` | 2024 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `multiscale_graphcorr` | 2015 | 1.1259e+15 / 1.8279e+12 | 1.1259e+15 / 1.8279e+12 | — | — | — | — |
| `pearsonr` | 2012 | 18.9 / 3.6841 | 18.9 / 3.6841 | — | — | — | — |
| `pointbiserialr` | 2002 | 9.4721 / 2.0126 | 9.4721 / 2.0126 | — | — | — | — |
| `somersd` | 2012 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `spearmanr` | 2009 | 16.313 / 4.0574 | 16.313 / 4.0574 | — | — | — | — |
| `spearmanrho` | 2012 | 3.6053 / 0.50006 | 3.6053 / 0.50006 | — | — | — | — |
| `weightedtau` | 2032 | 6.7659 / 0.01214 | 6.7659 / 0.01214 | — | — | — | — |

### `sci::stats::binning` (5 functions, 10027 cases)

| Function | Cases | Sci vs reference | SciPy | Boost.Math | GSL | NumPy | stdlib |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `binned_statistic` | 2014 | 0.71971 / 2.7511e-5 | 0.71971 / 2.7511e-5 | — | — | — | — |
| `binned_statistic_2d` | 2007 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `binned_statistic_dd` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `cumfreq` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `relfreq` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |

### `sci::stats::contingency` (10 functions, 20115 cases)

| Function | Cases | Sci vs reference | SciPy | Boost.Math | GSL | NumPy | stdlib |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `barnard_exact` | 2054 | 5.2267e+5 / 9008.9 | 5.2267e+5 / 9008.9 | — | — | — | — |
| `boschloo_exact` | 2027 | 5.2267e+5 / 25373.0 | 5.2267e+5 / 25373.0 | — | — | — | — |
| `chi2_contingency` | 2006 | 99.866 / 1.6943 | 99.866 / 1.6943 | — | — | — | — |
| `contingency.association` | 2006 | 1.0572 / 0.48521 | 1.0572 / 0.48521 | — | — | — | — |
| `contingency.crosstab` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `contingency.expected_freq` | 2003 | 0.92857 / 0.13877 | 0.92857 / 0.13877 | — | — | — | — |
| `contingency.margins` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `contingency.odds_ratio` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `contingency.relative_risk` | 2002 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `fisher_exact` | 2009 | 1192.0 / 170.67 | 1192.0 / 170.67 | — | — | — | — |

### `sci::stats::descriptive` (42 functions, 84446 cases)

| Function | Cases | Sci vs reference | SciPy | Boost.Math | GSL | NumPy | stdlib |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `circmean` | 2003 | 2.6782 / 0.56516 | 2.6782 / 0.56516 | — | — | — | — |
| `circstd` | 2006 | 97.198 / 21.938 | 97.198 / 21.938 | — | — | — | — |
| `circvar` | 2003 | 193.53 / 43.992 | 193.53 / 43.992 | — | — | — | — |
| `describe` | 2004 | 28.797 / 1.3885 | 28.797 / 1.3885 | — | — | — | — |
| `differential_entropy` | 2008 | 1.6289 / 0.33512 | 1.6289 / 0.33512 | — | — | — | — |
| `directional_stats` | 2004 | 1.202 / 0.1994 | 1.202 / 0.1994 | — | — | — | — |
| `entropy` | 2006 | 3.7856 / 5.1559e-3 | 3.7856 / 5.1559e-3 | — | — | — | — |
| `expectile` | 2015 | 0.92308 / 4.581e-4 | 0.92308 / 4.581e-4 | — | — | — | — |
| `gmean` | 2003 | 1.4156 / 0.083665 | 1.4156 / 0.083665 | — | — | — | — |
| `gstd` | 2006 | 0.92144 / 0.053096 | 0.92144 / 0.053096 | — | — | — | — |
| `gzscore` | 2004 | 3.9131 / 0.81393 | 3.9131 / 0.81393 | — | — | — | — |
| `hmean` | 2003 | 0.93844 / 0.26447 | 0.93844 / 0.26447 | — | — | — | — |
| `iqr` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `kstat` | 2008 | 8.9074 / 7.0074e-3 | 8.9074 / 7.0074e-3 | — | — | — | — |
| `kstatvar` | 2004 | 3.9016 / 0.67103 | 3.9016 / 0.67103 | — | — | — | — |
| `kurtosis` | 2010 | 283.25 / 0.79221 | 283.25 / 0.79221 | — | — | — | — |
| `lmoment` | 2020 | 26.091 / 0.40889 | 26.091 / 0.40889 | — | — | — | — |
| `mean` | 2004 | 0.83416 / 0.13556 | 0.83416 / 0.13556 | — | — | — | — |
| `median` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `median_abs_deviation` | 2008 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `mode` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `moment` | 2024 | 1.9207 / 3.348e-3 | 1.9207 / 3.348e-3 | — | — | — | — |
| `pmean` | 2018 | 1.0667 / 0.16185 | 1.0667 / 0.16185 | — | — | — | — |
| `quantile` | 2020 | 1.6667 / 1.3751e-3 | 1.6667 / 1.3751e-3 | — | — | — | — |
| `scoreatpercentile` | 2020 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `sem` | 2006 | 1.4042 / 0.26019 | 1.4042 / 0.26019 | — | — | — | — |
| `sigmaclip` | 2006 | 2.4997 / 0.11594 | 2.4997 / 0.11594 | — | — | — | — |
| `skew` | 2006 | 28.797 / 8.445 | 28.797 / 8.445 | — | — | — | — |
| `std` | 2007 | 1.2558 / 0.26731 | 1.2558 / 0.26731 | — | — | — | — |
| `tmax` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `tmean` | 2016 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `tmin` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `trim1` | 2012 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `trim_mean` | 2016 | 0.83416 / 0.14467 | 0.83416 / 0.14467 | — | — | — | — |
| `trimboth` | 2016 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `tsem` | 2015 | 1.0 / 0.15831 | 1.0 / 0.15831 | — | — | — | — |
| `tstd` | 2015 | 0.7071 / 0.20201 | 0.7071 / 0.20201 | — | — | — | — |
| `tvar` | 2015 | 0.6 / 0.07645 | 0.6 / 0.07645 | — | — | — | — |
| `variance` | 2079 | 2.3139 / 0.36574 | 2.4628e+8 / 3.0635e+5 | — | — | 2.4628e+8 / 3.0635e+5 | 2.5116 / 0.29218 |
| `variation` | 2007 | 1.5059 / 0.3672 | 1.5059 / 0.3672 | — | — | — | — |
| `zmap` | 2005 | 1.0 / 4.2184e-4 | 1.0 / 4.2184e-4 | — | — | — | — |
| `zscore` | 2004 | 7.9687 / 0.85237 | 7.9687 / 0.85237 | — | — | — | — |

### `sci::stats::diagnostics` (7 functions, 14072 cases)

| Function | Cases | Sci vs reference | SciPy | Boost.Math | GSL | NumPy | stdlib |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `boxcox_normmax` | 2016 | 5.0407e+23 / 2.0359e+21 | 5.0407e+23 / 2.0359e+21 | — | — | — | — |
| `boxcox_normplot` | 2006 | 2.2429 / 0.24871 | 2.2429 / 0.24871 | — | — | — | — |
| `ppcc_max` | 2015 | 1.786e+9 / 4.311e+8 | 1.786e+9 / 4.311e+8 | — | — | — | — |
| `ppcc_plot` | 2015 | 13.208 / 0.29083 | 13.208 / 0.29083 | — | — | — | — |
| `probplot` | 2003 | 1.2205 / 6.0934e-5 | 1.2205 / 6.0934e-5 | — | — | — | — |
| `yeojohnson_normmax` | 2008 | 2.2946e+8 / 4.7734e+7 | 2.2946e+8 / 4.7734e+7 | — | — | — | — |
| `yeojohnson_normplot` | 2009 | 2.6227 / 0.22653 | 2.6227 / 0.22653 | — | — | — | — |

### `sci::stats::distances` (3 functions, 6016 cases)

| Function | Cases | Sci vs reference | SciPy | Boost.Math | GSL | NumPy | stdlib |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `energy_distance` | 2005 | 1.1558 / 0.3196 | 1.1558 / 0.3196 | — | — | — | — |
| `wasserstein_distance` | 2005 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `wasserstein_distance_nd` | 2006 | 1.4142 / 0.67917 | 1.4142 / 0.67917 | — | — | — | — |

### `sci::stats::distributions` (12 functions, 25995 cases)

| Function | Cases | Sci vs reference | SciPy | Boost.Math | GSL | NumPy | stdlib |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `cdf` | 2379 | 4.2666e+5 / 8782.7 | 4.2666e+5 / 8782.7 | — | — | — | — |
| `fit` | 2036 | 2.8955 / 2.0096e-3 | 2.8955 / 2.0096e-3 | — | — | — | — |
| `logpdf` | 2163 | 832.64 / 10.889 | 832.64 / 10.889 | — | — | — | — |
| `make_distribution` | 2015 | 5.6207 / 0.35452 | 5.6207 / 0.35452 | — | — | — | — |
| `normal_cdf` | 2111 | 940.89 / 110.9 | 0.0 / 0.0 | 1488.8 / 216.15 | 940.82 / 110.99 | — | 1488.8 / 216.15 |
| `normal_pdf` | 2000 | 1.5723 / 0.78551 | 1.5723 / 0.78551 | 1.5723 / 0.78551 | 1.1589 / 0.62552 | — | — |
| `normal_ppf` | 2000 | 0.0 / 0.0 | 0.0 / 0.0 | 2.785 / 0.39778 | 3.7134 / 0.67276 | — | — |
| `order_statistic` | 2160 | 2107.7 / 3.1201 | 2107.7 / 3.1201 | — | — | — | — |
| `pdf` | 2377 | 490.26 / 9.7963 | 490.26 / 9.7963 | — | — | — | — |
| `ppf` | 2300 | 4.6057e+6 / 6617.6 | 4.6057e+6 / 6617.6 | — | — | — | — |
| `sf` | 2379 | 2.9859e+5 / 4114.0 | 2.9859e+5 / 4114.0 | — | — | — | — |
| `truncate` | 2075 | 4.5036e+15 / 2.6045e+12 | 4.5036e+15 / 2.6045e+12 | — | — | — | — |

### `sci::stats::estimation` (2 functions, 4030 cases)

| Function | Cases | Sci vs reference | SciPy | Boost.Math | GSL | NumPy | stdlib |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `bayes_mvs` | 2012 | 83.096 / 3.6539 | 83.096 / 3.6539 | — | — | — | — |
| `mvsdist` | 2018 | 516.22 / 1.083 | 516.22 / 1.083 | — | — | — | — |

### `sci::stats::goodness` (9 functions, 18125 cases)

| Function | Cases | Sci vs reference | SciPy | Boost.Math | GSL | NumPy | stdlib |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `anderson` | 2011 | 99.687 / 0.4487 | 99.687 / 0.4487 | — | — | — | — |
| `anderson_ksamp` | 2014 | 6.1166 / 2.0784 | 6.1166 / 2.0784 | — | — | — | — |
| `chisquare` | 2006 | 5.2329 / 1.0128 | 5.2329 / 1.0128 | — | — | — | — |
| `cramervonmises` | 2014 | 2.6969e+5 / 67.074 | 2.6969e+5 / 67.074 | — | — | — | — |
| `cramervonmises_2samp` | 2013 | 15961.0 / 3.9766 | 15961.0 / 3.9766 | — | — | — | — |
| `goodness_of_fit` | 2012 | 4.7452 / 0.14129 | 4.7452 / 0.14129 | — | — | — | — |
| `ks_1samp` | 2009 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `ks_2samp` | 2012 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `power_divergence` | 2034 | 34.056 / 0.81664 | 34.056 / 0.81664 | — | — | — | — |

### `sci::stats::inference` (27 functions, 54625 cases)

| Function | Cases | Sci vs reference | SciPy | Boost.Math | GSL | NumPy | stdlib |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `ansari` | 2015 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `bartlett` | 2003 | 14.111 / 3.3287 | 14.111 / 3.3287 | — | — | — | — |
| `binomtest` | 2048 | 138.29 / 3.8558 | 138.29 / 3.8558 | — | — | — | — |
| `brunnermunzel` | 2030 | 87.763 / 8.0938 | 87.763 / 8.0938 | — | — | — | — |
| `combine_pvalues` | 2011 | 1.9401e+5 / 100.03 | 1.9401e+5 / 100.03 | — | — | — | — |
| `dunnett` | 2015 | 1.3013e+10 / 7.1396e+6 | 1.3013e+10 / 7.1396e+6 | — | — | — | — |
| `epps_singleton_2samp` | 2020 | 5.5464e+6 / 36092.0 | 5.5464e+6 / 36092.0 | — | — | — | — |
| `f_oneway` | 2003 | 36.246 / 11.131 | 36.246 / 11.131 | — | — | — | — |
| `false_discovery_control` | 2006 | 0.93844 / 0.062734 | 0.93844 / 0.062734 | — | — | — | — |
| `fligner` | 2008 | 48.635 / 20.564 | 48.635 / 20.564 | — | — | — | — |
| `friedmanchisquare` | 2002 | 2.1827 / 1.0911 | 2.1827 / 1.0911 | — | — | — | — |
| `kruskal` | 2003 | 31.544 / 0.98758 | 31.544 / 0.98758 | — | — | — | — |
| `levene` | 2002 | 7.7383 / 2.3693 | 7.7383 / 2.3693 | — | — | — | — |
| `mannwhitneyu` | 2050 | 0.73669 / 1.7968e-4 | 0.73669 / 1.7968e-4 | — | — | — | — |
| `median_test` | 2022 | 17.422 / 5.0827e-3 | 17.422 / 5.0827e-3 | — | — | — | — |
| `mood` | 2021 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `page_trend_test` | 2010 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `poisson_means_test` | 2030 | 3.9349e+5 / 6892.3 | 3.9349e+5 / 6892.3 | — | — | — | — |
| `power` | 2020 | 14.694 / 7.1974 | 14.694 / 7.1974 | — | — | — | — |
| `quantile_test` | 2108 | 42.896 / 0.23919 | 42.896 / 0.23919 | — | — | — | — |
| `ranksums` | 2015 | 4.1473e+5 / 1.0352e+5 | 4.1473e+5 / 1.0352e+5 | — | — | — | — |
| `ttest_1samp` | 2012 | 17.687 / 3.3458 | 17.687 / 3.3458 | — | — | — | — |
| `ttest_ind` | 2030 | 63.924 / 5.5843 | 63.924 / 5.5843 | — | — | — | — |
| `ttest_ind_from_stats` | 2002 | 554.04 / 50.071 | 554.04 / 50.071 | — | — | — | — |
| `ttest_rel` | 2015 | 21.796 / 2.3601 | 21.796 / 2.3601 | — | — | — | — |
| `tukey_hsd` | 2010 | 3.1942e+15 / 3.8324e+11 | 3.1942e+15 / 3.8324e+11 | — | — | — | — |
| `wilcoxon` | 2114 | 0.52036 / 1.2307e-4 | 0.52036 / 1.2307e-4 | — | — | — | — |

### `sci::stats::multivariate` (10 functions, 20083 cases)

| Function | Cases | Sci vs reference | SciPy | Boost.Math | GSL | NumPy | stdlib |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `categorical_cdf` | 2024 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `categorical_pmf` | 2024 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `dirichlet_logpdf` | 2003 | 4.5036e+15 / 2.2484e+12 | 4.5036e+15 / 2.2484e+12 | — | — | — | — |
| `dirichlet_pdf` | 2003 | 17.983 / 4.0248 | 17.983 / 4.0248 | — | — | — | — |
| `multinomial_logpmf` | 2006 | 2.9398e+5 / 147.45 | 2.9398e+5 / 147.45 | — | — | — | — |
| `multinomial_pmf` | 2007 | 2.0972e+6 / 1048.6 | 2.0972e+6 / 1048.6 | — | — | — | — |
| `multivariate_normal_logpdf` | 2002 | 1.4343 / 0.64566 | 1.4343 / 0.64566 | — | — | — | — |
| `multivariate_normal_pdf` | 2002 | 3.5594 / 1.001 | 3.5594 / 1.001 | — | — | — | — |
| `multivariate_t_logpdf` | 2006 | 27.701 / 2.4638 | 27.701 / 2.4638 | — | — | — | — |
| `multivariate_t_pdf` | 2006 | 200.81 / 8.5053 | 200.81 / 8.5053 | — | — | — | — |

### `sci::stats::normality` (7 functions, 14077 cases)

| Function | Cases | Sci vs reference | SciPy | Boost.Math | GSL | NumPy | stdlib |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `alexandergovern` | 2005 | 88.121 / 4.8069 | 88.121 / 4.8069 | — | — | — | — |
| `jarque_bera` | 2002 | 5.7494 / 1.0431 | 5.7494 / 1.0431 | — | — | — | — |
| `kurtosistest` | 2021 | 576.73 / 2.4259 | 576.73 / 2.4259 | — | — | — | — |
| `normaltest` | 2006 | 6.78e+15 / 1.2996e+15 | 6.78e+15 / 1.2996e+15 | — | — | — | — |
| `obrientransform` | 2004 | 21.538 / 0.72203 | 21.538 / 0.72203 | — | — | — | — |
| `shapiro` | 2021 | 6.0048e+15 / 1.4856e+12 | 6.0048e+15 / 1.4856e+12 | — | — | — | — |
| `skewtest` | 2018 | 2.7302e+31 / 3.5875e+30 | 2.7302e+31 / 3.5875e+30 | — | — | — | — |

### `sci::stats::qmc` (4 functions, 8039 cases)

| Function | Cases | Sci vs reference | SciPy | Boost.Math | GSL | NumPy | stdlib |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `qmc.discrepancy` | 2024 | 37.775 / 2.7787 | 37.775 / 2.7787 | — | — | — | — |
| `qmc.geometric_discrepancy` | 2006 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `qmc.scale` | 2006 | 2.6666 / 0.537 | 2.6666 / 0.537 | — | — | — | — |
| `qmc.update_discrepancy` | 2003 | 105.64 / 29.38 | 105.64 / 29.38 | — | — | — | — |

### `sci::stats::ranking` (5 functions, 10332 cases)

| Function | Cases | Sci vs reference | SciPy | Boost.Math | GSL | NumPy | stdlib |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `ecdf` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `estimated_cdf` | 2224 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `percentileofscore` | 2080 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `rankdata` | 2020 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `tiecorrect` | 2004 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |

### `sci::stats::regression` (3 functions, 6071 cases)

| Function | Cases | Sci vs reference | SciPy | Boost.Math | GSL | NumPy | stdlib |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `linregress` | 2009 | 4.5036e+15 / 7.4724e+11 | 4.5036e+15 / 7.4724e+11 | — | — | — | — |
| `siegelslopes` | 2014 | 0.5 / 0.015268 | 0.5 / 0.015268 | — | — | — | — |
| `theilslopes` | 2048 | 0.5 / 0.015808 | 0.5 / 0.015808 | — | — | — | — |

### `sci::stats::resampling` (4 functions, 8111 cases)

| Function | Cases | Sci vs reference | SciPy | Boost.Math | GSL | NumPy | stdlib |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `bootstrap` | 2036 | 2.0385 / 0.096895 | 2.0385 / 0.096895 | — | — | — | — |
| `bws_test` | 2012 | 25.412 / 7.2321e-3 | 25.412 / 7.2321e-3 | — | — | — | — |
| `monte_carlo_test` | 2027 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `permutation_test` | 2036 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |

### `sci::stats::sensitivity` (1 functions, 2007 cases)

| Function | Cases | Sci vs reference | SciPy | Boost.Math | GSL | NumPy | stdlib |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `sobol_indices` | 2007 | 9.0076e+15 / 1.683e+12 | 9.0076e+15 / 1.683e+12 | — | — | — | — |

### `sci::stats::survival` (1 functions, 2024 cases)

| Function | Cases | Sci vs reference | SciPy | Boost.Math | GSL | NumPy | stdlib |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `logrank` | 2024 | 3.9598 / 4.3075e-3 | 3.9598 / 4.3075e-3 | — | — | — | — |

### `sci::stats::transformations` (4 functions, 8060 cases)

| Function | Cases | Sci vs reference | SciPy | Boost.Math | GSL | NumPy | stdlib |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `boxcox` | 2015 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `boxcox_llf` | 2010 | 1.3356 / 0.33411 | 1.3356 / 0.33411 | — | — | — | — |
| `yeojohnson` | 2021 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `yeojohnson_llf` | 2014 | 1.9043 / 0.11034 | 1.9043 / 0.11034 | — | — | — | — |

### `sci::stats` (4 functions, 8028 cases)

| Function | Cases | Sci vs reference | SciPy | Boost.Math | GSL | NumPy | stdlib |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| `abs` | 2006 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `exp` | 2007 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `kstest` | 2012 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |
| `log` | 2003 | 0.0 / 0.0 | 0.0 / 0.0 | — | — | — | — |


## Conventions and limitations

Distribution parameters, tail conventions, alternatives, ddof, explicit resampling draws and output ordering follow the pinned statistics fixtures. SciPy supplies binary64 references for added APIs, including documented adaptations and independently seeded probabilistic references. Weighted tau compares its statistic; its always-NaN p-value is left to the ordinary API tests. Normality, inverse and joint-tail approximations have reference and algorithm accuracy limits.

Inputs are generated deterministically into `sci-testing/cases/generated/<module>/<function>.csv` (optionally `.gz`); raw outputs, high-precision
references and the environment are stored in that repository's `results/`.
NaN, infinity, signed-zero contracts, invalid domains and work-limit enforcement
are outside this finite-input accuracy sample. Existing Sci API tests cover
additional behavior. The [coverage audit](../../results/coverage.csv) checks every public function
against the compiler-resolved inventory. These samples do not establish
worst-case bounds or full domain coverage. SciPy aliases, compatibility namespaces and Cython
variants are not counted as independent implementations.

## Reproduction

- Sci source: `c781e8259b175ebeb3c7d27a522ef3f3a7e24998` (numerical sources clean: true).
- Testing source SHA-256: `3d1f538d12ef58c6738188ac9a526fb0459c4c5cb1023a3c7edaa3c23d5a3452`.
- NumPy 2.5.4; SciPy 1.18.1; FilterPy 1.4.5; mpmath 1.3.0.
- Native comparison versions: {'Boost.Math': '1_92', 'GSL': '2.8'}.
- rustc 1.99.1 (58cf8f9c3 2026-10-13); Python 3.13.3; dev; opt-level=0; debug=line-tables-only; no custom RUSTFLAGS.
- Linux-7.2.11-arch1-1-x86_64-with-glibc2.44; AMD Ryzen 5 4600H with Radeon Graphics.

See the [shared metric, provenance and reproduction instructions](README.md).

Sci-only refresh: independent references and other library outputs are retained from the saved snapshot. Their original environment is recorded in `results/environment.json` under `reference_environment`.

Sci results were refreshed by namespace. Per-function Sci revisions are recorded in `results/environment.json` under `sci_measurement_commits`.
