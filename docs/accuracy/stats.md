# stats: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-stats)

Measured binary64 (`f64`) results using **Sci-stats as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 342414 cases and 169 public API names;
coverage is limited to the operations and inputs below.

## Agreement with Sci

Each library cell is **maximum / mean difference from Sci in ε** over cases
where both evaluations succeed. `Cases` counts input cases; `Sci failed N/T`
identifies missing baseline evaluations. `failed N/T` counts that library's own
failures; `Sci unavailable N/T` counts its successful evaluations excluded
because Sci failed. A dash means unavailable or no paired finite results.
Means include all compared scalar components, including exact zeros.
Full outputs, references, warnings and errors remain in the per-function CSVs.

| Function | Cases | SciPy | Boost.Math | GSL | NumPy | stdlib |
| --- | ---: | --- | --- | --- | --- | --- |
| `abs` | 2006 | 0.0 / 0.0 | — | — | — | — |
| `alexandergovern` | 2005 | 88.121 / 4.8069 | — | — | — | — |
| `anderson` | 2011 | 99.687 / 0.4487 | — | — | — | — |
| `anderson_ksamp` | 2014 | 6.1166 / 2.0784 | — | — | — | — |
| `ansari` | 2015 | 0.0 / 0.0 | — | — | — | — |
| `barnard_exact` | 2054 | 5.2267e+5 / 9008.9 | — | — | — | — |
| `bartlett` | 2003 | 14.111 / 3.3287 | — | — | — | — |
| `bayes_mvs` | 2012 | 83.096 / 3.6539 | — | — | — | — |
| `binned_statistic` | 2014 | 0.71971 / 2.7511e-5 | — | — | — | — |
| `binned_statistic_2d` | 2007 | 0.0 / 0.0 | — | — | — | — |
| `binned_statistic_dd` | 2002 | 0.0 / 0.0 | — | — | — | — |
| `binomtest` | 2048 | 138.29 / 3.8558 | — | — | — | — |
| `bootstrap` | 2036 | 2.0385 / 0.096895 | — | — | — | — |
| `boschloo_exact` | 2027 | 5.2267e+5 / 25373.0 | — | — | — | — |
| `boxcox` | 2015 | 0.0 / 0.0 | — | — | — | — |
| `boxcox_llf` | 2010 | 1.3356 / 0.33411 | — | — | — | — |
| `boxcox_normmax` | 2016 | 5.0407e+23 / 2.0359e+21 | — | — | — | — |
| `boxcox_normplot` | 2006 | 2.2429 / 0.24871 | — | — | — | — |
| `brunnermunzel` | 2030 | 87.763 / 8.0938 | — | — | — | — |
| `bws_test` | 2012 | 25.412 / 7.2321e-3 | — | — | — | — |
| `categorical_cdf` | 2024 | 0.0 / 0.0 | — | — | — | — |
| `categorical_pmf` | 2024 | 0.0 / 0.0 | — | — | — | — |
| `cdf` | 2379 | 4.2666e+5 / 8782.7 | — | — | — | — |
| `chatterjeexi` | 2013 | 26.351 / 0.010368 | — | — | — | — |
| `chi2_contingency` | 2006 | 99.866 / 1.6943 | — | — | — | — |
| `chisquare` | 2006 | 5.2329 / 1.0128 | — | — | — | — |
| `circmean` | 2003 | 2.6782 / 0.56516 | — | — | — | — |
| `circstd` | 2006 | 97.198 / 21.938 | — | — | — | — |
| `circvar` | 2003 | 193.53 / 43.992 | — | — | — | — |
| `combine_pvalues` | 2011 | 1.9401e+5 / 100.03 | — | — | — | — |
| `contingency.association` | 2006 | 1.0572 / 0.48521 | — | — | — | — |
| `contingency.crosstab` | 2002 | 0.0 / 0.0 | — | — | — | — |
| `contingency.expected_freq` | 2003 | 0.92857 / 0.13877 | — | — | — | — |
| `contingency.margins` | 2003 | 0.0 / 0.0 | — | — | — | — |
| `contingency.odds_ratio` | 2003 | 0.0 / 0.0 | — | — | — | — |
| `contingency.relative_risk` | 2002 | 0.0 / 0.0 | — | — | — | — |
| `cramervonmises` | 2014 | 2.6969e+5 / 67.074 | — | — | — | — |
| `cramervonmises_2samp` | 2013 | 15961.0 / 3.9766 | — | — | — | — |
| `cumfreq` | 2002 | 0.0 / 0.0 | — | — | — | — |
| `describe` | 2004 | 28.797 / 1.3885 | — | — | — | — |
| `differential_entropy` | 2008 | 1.6289 / 0.33512 | — | — | — | — |
| `directional_stats` | 2004 | 1.202 / 0.1994 | — | — | — | — |
| `dirichlet_logpdf` | 2003 | 4.5036e+15 / 2.2484e+12 | — | — | — | — |
| `dirichlet_pdf` | 2003 | 17.983 / 4.0248 | — | — | — | — |
| `dunnett` | 2015 | 1.3013e+10 / 7.1396e+6 | — | — | — | — |
| `ecdf` | 2004 | 0.0 / 0.0 | — | — | — | — |
| `energy_distance` | 2005 | 1.1558 / 0.3196 | — | — | — | — |
| `entropy` | 2006 | 3.7856 / 5.1559e-3 | — | — | — | — |
| `epps_singleton_2samp` | 2020 | 5.5464e+6 / 36092.0 | — | — | — | — |
| `estimated_cdf` | 2224 | 0.0 / 0.0 | — | — | — | — |
| `exp` | 2007 | 0.0 / 0.0 | — | — | — | — |
| `expectile` | 2015 | 0.92308 / 4.581e-4 | — | — | — | — |
| `f_oneway` | 2003 | 36.246 / 11.131 | — | — | — | — |
| `false_discovery_control` | 2006 | 0.93844 / 0.062734 | — | — | — | — |
| `fisher_exact` | 2009 | 1192.0 / 170.67 | — | — | — | — |
| `fit` | 2036 | 2.8955 / 2.0096e-3 | — | — | — | — |
| `fligner` | 2008 | 48.635 / 20.564 | — | — | — | — |
| `friedmanchisquare` | 2002 | 2.1827 / 1.0911 | — | — | — | — |
| `gmean` | 2003 | 1.4156 / 0.083665 | — | — | — | — |
| `goodness_of_fit` | 2012 | 4.7452 / 0.14129 | — | — | — | — |
| `gstd` | 2006 | 0.92144 / 0.053096 | — | — | — | — |
| `gzscore` | 2004 | 3.9131 / 0.81393 | — | — | — | — |
| `hmean` | 2003 | 0.93844 / 0.26447 | — | — | — | — |
| `iqr` | 2004 | 0.0 / 0.0 | — | — | — | — |
| `jarque_bera` | 2002 | 5.7494 / 1.0431 | — | — | — | — |
| `kendalltau` | 2024 | 0.0 / 0.0 | — | — | — | — |
| `kruskal` | 2003 | 31.544 / 0.98758 | — | — | — | — |
| `ks_1samp` | 2009 | 0.0 / 0.0 | — | — | — | — |
| `ks_2samp` | 2012 | 0.0 / 0.0 | — | — | — | — |
| `kstat` | 2008 | 8.9074 / 7.0074e-3 | — | — | — | — |
| `kstatvar` | 2004 | 3.9016 / 0.67103 | — | — | — | — |
| `kstest` | 2012 | 0.0 / 0.0 | — | — | — | — |
| `kurtosis` | 2010 | 283.25 / 0.79221 | — | — | — | — |
| `kurtosistest` | 2021 | 576.73 / 2.4259 | — | — | — | — |
| `levene` | 2002 | 7.7383 / 2.3693 | — | — | — | — |
| `linregress` | 2009 | 4.5036e+15 / 7.4724e+11 | — | — | — | — |
| `lmoment` | 2020 | 26.091 / 0.40889 | — | — | — | — |
| `log` | 2003 | 0.0 / 0.0 | — | — | — | — |
| `logpdf` | 2163 | 832.64 / 10.889 | — | — | — | — |
| `logrank` | 2024 | 3.9598 / 4.3075e-3 | — | — | — | — |
| `make_distribution` | 2015 | 5.6207 / 0.35452 | — | — | — | — |
| `mannwhitneyu` | 2050 | 0.73669 / 1.7968e-4 | — | — | — | — |
| `mean` | 2004 | 0.83416 / 0.13556 | — | — | — | — |
| `median` | 2004 | 0.0 / 0.0 | — | — | — | — |
| `median_abs_deviation` | 2008 | 0.0 / 0.0 | — | — | — | — |
| `median_test` | 2022 | 17.422 / 5.0827e-3 | — | — | — | — |
| `mode` | 2004 | 0.0 / 0.0 | — | — | — | — |
| `moment` | 2024 | 1.9207 / 3.348e-3 | — | — | — | — |
| `monte_carlo_test` | 2027 | 0.0 / 0.0 | — | — | — | — |
| `mood` | 2021 | 0.0 / 0.0 | — | — | — | — |
| `multinomial_logpmf` | 2006 | 2.9398e+5 / 147.45 | — | — | — | — |
| `multinomial_pmf` | 2007 | 2.0972e+6 / 1048.6 | — | — | — | — |
| `multiscale_graphcorr` | 2015 | 1.1259e+15 / 1.8279e+12 | — | — | — | — |
| `multivariate_normal_logpdf` | 2002 | 1.4343 / 0.64566 | — | — | — | — |
| `multivariate_normal_pdf` | 2002 | 3.5594 / 1.001 | — | — | — | — |
| `multivariate_t_logpdf` | 2006 | 27.701 / 2.4638 | — | — | — | — |
| `multivariate_t_pdf` | 2006 | 200.81 / 8.5053 | — | — | — | — |
| `mvsdist` | 2018 | 516.22 / 1.083 | — | — | — | — |
| `normal_cdf` | 2111 | 0.0 / 0.0 | 1488.8 / 216.15 | 940.82 / 110.99 | — | 1488.8 / 216.15 |
| `normal_pdf` | 2000 | 1.5723 / 0.78551 | 1.5723 / 0.78551 | 1.1589 / 0.62552 | — | — |
| `normal_ppf` | 2000 | 0.0 / 0.0 | 2.785 / 0.39778 | 3.7134 / 0.67276 | — | — |
| `normaltest` | 2006 | 6.78e+15 / 1.2996e+15 | — | — | — | — |
| `obrientransform` | 2004 | 21.538 / 0.72203 | — | — | — | — |
| `order_statistic` | 2160 | 2107.7 / 3.1201 | — | — | — | — |
| `page_trend_test` | 2010 | 0.0 / 0.0 | — | — | — | — |
| `pdf` | 2377 | 490.26 / 9.7963 | — | — | — | — |
| `pearsonr` | 2012 | 18.9 / 3.6841 | — | — | — | — |
| `percentileofscore` | 2080 | 0.0 / 0.0 | — | — | — | — |
| `permutation_test` | 2036 | 0.0 / 0.0 | — | — | — | — |
| `pmean` | 2018 | 1.0667 / 0.16185 | — | — | — | — |
| `pointbiserialr` | 2002 | 9.4721 / 2.0126 | — | — | — | — |
| `poisson_means_test` | 2030 | 3.9349e+5 / 6892.3 | — | — | — | — |
| `power` | 2020 | 14.694 / 7.1974 | — | — | — | — |
| `power_divergence` | 2034 | 34.056 / 0.81664 | — | — | — | — |
| `ppcc_max` | 2015 | 1.786e+9 / 4.311e+8 | — | — | — | — |
| `ppcc_plot` | 2015 | 13.208 / 0.29083 | — | — | — | — |
| `ppf` | 2300 | 4.6057e+6 / 6617.6 | — | — | — | — |
| `probplot` | 2003 | 1.2205 / 6.0934e-5 | — | — | — | — |
| `qmc.discrepancy` | 2024 | 37.775 / 2.7787 | — | — | — | — |
| `qmc.geometric_discrepancy` | 2006 | 0.0 / 0.0 | — | — | — | — |
| `qmc.scale` | 2006 | 2.6666 / 0.537 | — | — | — | — |
| `qmc.update_discrepancy` | 2003 | 105.64 / 29.38 | — | — | — | — |
| `quantile` | 2020 | 1.6667 / 1.3751e-3 | — | — | — | — |
| `quantile_test` | 2108 | 42.896 / 0.23919 | — | — | — | — |
| `rankdata` | 2020 | 0.0 / 0.0 | — | — | — | — |
| `ranksums` | 2015 | 4.1473e+5 / 1.0352e+5 | — | — | — | — |
| `relfreq` | 2002 | 0.0 / 0.0 | — | — | — | — |
| `scoreatpercentile` | 2020 | 0.0 / 0.0 | — | — | — | — |
| `sem` | 2006 | 1.4042 / 0.26019 | — | — | — | — |
| `sf` | 2379 | 2.9859e+5 / 4114.0 | — | — | — | — |
| `shapiro` | 2021 | 6.0048e+15 / 1.4856e+12 | — | — | — | — |
| `siegelslopes` | 2014 | 0.5 / 0.015268 | — | — | — | — |
| `sigmaclip` | 2006 | 2.4997 / 0.11594 | — | — | — | — |
| `skew` | 2006 | 28.797 / 8.445 | — | — | — | — |
| `skewtest` | 2018 | 2.7302e+31 / 3.5875e+30 | — | — | — | — |
| `sobol_indices` | 2007 | 9.0076e+15 / 1.683e+12 | — | — | — | — |
| `somersd` | 2012 | 0.0 / 0.0 | — | — | — | — |
| `spearmanr` | 2009 | 16.313 / 4.0574 | — | — | — | — |
| `spearmanrho` | 2012 | 3.6053 / 0.50006 | — | — | — | — |
| `std` | 2007 | 1.2558 / 0.26731 | — | — | — | — |
| `theilslopes` | 2048 | 0.5 / 0.015808 | — | — | — | — |
| `tiecorrect` | 2004 | 0.0 / 0.0 | — | — | — | — |
| `tmax` | 2004 | 0.0 / 0.0 | — | — | — | — |
| `tmean` | 2016 | 0.0 / 0.0 | — | — | — | — |
| `tmin` | 2004 | 0.0 / 0.0 | — | — | — | — |
| `trim1` | 2012 | 0.0 / 0.0 | — | — | — | — |
| `trim_mean` | 2016 | 0.83416 / 0.14467 | — | — | — | — |
| `trimboth` | 2016 | 0.0 / 0.0 | — | — | — | — |
| `truncate` | 2075 | 4.5036e+15 / 2.6045e+12 | — | — | — | — |
| `tsem` | 2015 | 1.0 / 0.15831 | — | — | — | — |
| `tstd` | 2015 | 0.7071 / 0.20201 | — | — | — | — |
| `ttest_1samp` | 2012 | 17.687 / 3.3458 | — | — | — | — |
| `ttest_ind` | 2030 | 63.924 / 5.5843 | — | — | — | — |
| `ttest_ind_from_stats` | 2002 | 554.04 / 50.071 | — | — | — | — |
| `ttest_rel` | 2015 | 21.796 / 2.3601 | — | — | — | — |
| `tukey_hsd` | 2010 | 3.1942e+15 / 3.8324e+11 | — | — | — | — |
| `tvar` | 2015 | 0.6 / 0.07645 | — | — | — | — |
| `variance` | 2079 | 2.4628e+8 / 3.0635e+5 | — | — | 2.4628e+8 / 3.0635e+5 | 2.5116 / 0.29218 |
| `variation` | 2007 | 1.5059 / 0.3672 | — | — | — | — |
| `wasserstein_distance` | 2005 | 0.0 / 0.0 | — | — | — | — |
| `wasserstein_distance_nd` | 2006 | 1.4142 / 0.67917 | — | — | — | — |
| `weightedtau` | 2032 | 6.7659 / 0.01214 | — | — | — | — |
| `wilcoxon` | 2114 | 0.52036 / 1.2307e-4 | — | — | — | — |
| `yeojohnson` | 2021 | 0.0 / 0.0 | — | — | — | — |
| `yeojohnson_llf` | 2014 | 1.9043 / 0.11034 | — | — | — | — |
| `yeojohnson_normmax` | 2008 | 2.2946e+8 / 4.7734e+7 | — | — | — | — |
| `yeojohnson_normplot` | 2009 | 2.6227 / 0.22653 | — | — | — | — |
| `zmap` | 2005 | 1.0 / 4.2184e-4 | — | — | — | — |
| `zscore` | 2004 | 7.9687 / 0.85237 | — | — | — | — |


## Conventions and limitations

Distribution parameters, tail conventions, alternatives, ddof, explicit resampling draws and output ordering follow the pinned statistics fixtures. SciPy supplies binary64 references for added APIs, including documented adaptations and independently seeded probabilistic references. Weighted tau compares its statistic; its always-NaN p-value is left to the ordinary API tests. Normality, inverse and joint-tail approximations have reference and algorithm accuracy limits.

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
