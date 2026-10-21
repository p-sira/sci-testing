# special: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-special)

Measured binary64 (`f64`) results using **Sci-special as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 728841 cases and 364 public API names;
coverage is limited to the operations and inputs below.

## Agreement with Sci

Each library cell is **maximum / mean difference from Sci in ε** over cases
where both evaluations succeed. `Cases` counts input cases; `Sci failed N/T`
identifies missing baseline evaluations. `failed N/T` counts that library's own
failures; `Sci unavailable N/T` counts its successful evaluations excluded
because Sci failed. A dash means unavailable or no paired finite results.
Means include all compared scalar components, including exact zeros.
Full outputs, references, warnings and errors remain in the per-function CSVs.

| Function | Cases | SciPy | Boost.Math | GSL | stdlib |
| --- | ---: | --- | --- | --- | --- |
| `agm` | 2001 | 0.0 / 0.0 | — | — | — |
| `ai_zeros` | 2001 | 43772.0 / 1367.2 | — | — | — |
| `airy` | 2001 | 0.0 / 0.0 | — | — | — |
| `airye` | 2001 | 0.0 / 0.0 | — | — | — |
| `assoc_laguerre` | 2001 | 6.8158 / 0.39437 | — | — | — |
| `assoc_legendre_p` | 2001 | 0.0 / 0.0 | — | — | — |
| `assoc_legendre_p_all` | 2001 | 0.0 / 0.0 | — | — | — |
| `bdtr` | 2001 | 0.0 / 0.0 | — | — | — |
| `bdtrc` | 2001 | 0.0 / 0.0 | — | — | — |
| `bdtri` | 2001 | 0.0 / 0.0 | — | — | — |
| `bdtrik` | 2001 | 7.1304 / 1.6731 | — | — | — |
| `bdtrin` | 2001 | 3.1597 / 0.71833 | — | — | — |
| `bei` | 2001 | 0.0 / 0.0 | — | — | — |
| `bei_zeros` | 2001 | 0.0 / 0.0 | — | — | — |
| `beip` | 2001 | 0.0 / 0.0 | — | — | — |
| `beip_zeros` | 2001 | 0.0 / 0.0 | — | — | — |
| `ber` | 2001 | 0.0 / 0.0 | — | — | — |
| `ber_zeros` | 2001 | 0.0 / 0.0 | — | — | — |
| `bernoulli` | 2001 | 7760.6 / 691.42 | — | — | — |
| `berp` | 2001 | 0.0 / 0.0 | — | — | — |
| `berp_zeros` | 2001 | 0.0 / 0.0 | — | — | — |
| `besselpoly` | 2001 | 0.0 / 0.0 | — | — | — |
| `beta` | 2001 | 5.9302 / 1.5594 | 8.4711 / 2.2767 | 22.06 / 7.0893 | — |
| `betainc` | 2001 | 7.6968 / 2.2647 | 6.1639 / 1.334 | 16.414 / 6.0764 | — |
| `betaincc` | 2001 | 3.2445 / 1.4697 | — | — | — |
| `betainccinv` | 2001 | 3.8609 / 0.77062 | — | — | — |
| `betaincinv` | 2001 | 5.6373 / 1.4373 | — | — | — |
| `betaln` | 2001 | 1.8183 / 0.45892 | 2.4244 / 0.6917 | 7.2662 / 1.6828 | — |
| `bi_zeros` | 2001 | 71439.0 / 3266.4 | — | — | — |
| `binom` | 2001 | 0.0 / 0.0 | — | — | — |
| `boxcox` | 2001 | 0.0 / 0.0 | — | — | — |
| `boxcox1p` | 2001 | 0.0 / 0.0 | — | — | — |
| `btdtria` | 2001 | 7.5267 / 1.4761 | — | — | — |
| `btdtrib` | 2001 | 5.7996 / 1.3467 | — | — | — |
| `cbrt` | 2001 | 0.63017 / 0.064224 | — | — | — |
| `cel` | 2001 | 2.7746 / 0.61479 | — | — | — |
| `cel1` | 2001 | 2.1483 / 0.41855 | — | — | — |
| `cel2` | 2001 | 2.4585 / 0.51361 | — | — | — |
| `cel3` | 2001 | 2.8543 / 0.60068 | — | — | — |
| `chdtr` | 2001 | 0.0 / 0.0 | — | — | — |
| `chdtrc` | 2001 | 0.0 / 0.0 | — | — | — |
| `chdtri` | 2001 | 0.0 / 0.0 | — | — | — |
| `chdtriv` | 2001 | 2.738 / 0.60018 | — | — | — |
| `chebyc` | 2003 | 4.5036e+15 / 1.6478e+15 | — | — | — |
| `chebys` | 2003 | 4.5036e+15 / 2.0928e+15 | — | — | — |
| `chebyt` | 2005 | 4.5036e+15 / 1.6591e+15 | — | — | — |
| `chebyu` | 2005 | 4.5036e+15 / 1.6831e+15 | — | — | — |
| `chndtr` | 2001 | 2.7537 / 0.56996 | — | — | — |
| `chndtridf` | 2001 | 4.9693 / 1.0191 | — | — | — |
| `chndtrinc` | 2001 | 7.9331 / 1.8815 | — | — | — |
| `chndtrix` | 2001 | 4.7771 / 1.2404 | — | — | — |
| `comb` | 2001 | 0.0 / 0.0 | — | — | — |
| `cosdg` | 2001 | 1.155 / 0.25795 | — | — | — |
| `cosm1` | 2001 | 1.8215 / 0.53599 | — | — | — |
| `cotdg` | 2001 | 2.8884 / 0.69191 | — | — | — |
| `dawsn` | 2001 | 0.0 / 0.0 | — | — | — |
| `digamma` | 2005 | 6.409 / 1.4938 | 6.409 / 1.4866 | 6.409 / 1.4836 | — |
| `diric` | 2001 | 0.0 / 0.0 | — | — | — |
| `el1` | 2001 | 2.5548 / 0.53071 | — | — | — |
| `el2` | 2001 | 2.1592 / 0.49421 | — | — | — |
| `el3` | 2001 | 1.6937 / 0.5146 | — | — | — |
| `ellip_harm` | 2037 | 1.8991 / 0.029334 | — | — | — |
| `ellip_harm_2` | 2037 | 56527.0 / 729.09 | — | — | — |
| `ellip_normal` | 2037 | 2669.1 / 5.9574 | — | — | — |
| `ellipd` | 2001 | 2.9801 / 0.61298 | — | — | — |
| `ellipdinc` | 2001 | 3.4883 / 1.0207 | — | — | — |
| `ellipe` | 2001 | 1.4808 / 0.34596 | 0.74053 / 0.069931 | 3.7018 / 0.83622 | — |
| `ellipeinc` | 2001 | 3.3583 / 0.8555 | — | — | — |
| `ellipf` | 2001 | 3.3085 / 0.74569 | — | — | — |
| `ellipj` | 2001 | 0.0 / 0.0 | — | — | — |
| `ellipk` | 2001 | 1.0789 / 0.34259 | 0.53947 / 0.083827 | 2.697 / 0.61294 | — |
| `ellipke` | 2001 | 1.4808 / 0.33868 | — | — | — |
| `ellipkinc` | 2001 | 4.1356 / 0.76222 | — | — | — |
| `ellipkm1` | 2001 | 1.618 / 0.3097 | — | — | — |
| `ellippi` | 2001 | 1.9116 / 0.42646 | — | — | — |
| `ellippiinc` | 2001 | 2.4666 / 0.63483 | — | — | — |
| `ellippiinc_bulirsch` | 2001 | 2.4666 / 0.65742 | — | — | — |
| `elliprc` | 2001 | 2.0939 / 0.49619 | — | — | — |
| `elliprd` | 2001 | 2.277 / 0.46848 | — | — | — |
| `elliprf` | 2001 | 2.9002 / 0.57425 | — | — | — |
| `elliprg` | 2001 | 2.4904 / 0.45294 | — | — | — |
| `elliprj` | 2001 | 3.586 / 0.90923 | — | — | — |
| `entr` | 2001 | 0.0 / 0.0 | — | — | — |
| `erf` | 2143 | 1.4674 / 0.049532 | 0.96062 / 0.017472 | 2.2205 / 0.15874 | 0.91085 / 9.6726e-3 |
| `erf_zeros` | 2001 | 27.467 / 8.5619 | — | — | — |
| `erfc` | 2001 | 1.5644 / 0.36296 | 1.0436 / 0.3034 | 1.0435 / 0.40725 | — |
| `erfcinv` | 2001 | 0.0 / 0.0 | — | — | — |
| `erfcx` | 2001 | 0.0 / 0.0 | — | — | — |
| `erfi` | 2001 | 0.0 / 0.0 | — | — | — |
| `erfinv` | 2001 | 1.5722 / 0.37406 | — | — | — |
| `euler` | 2001 | 5.1888 / 0.64984 | — | — | — |
| `eval_chebyc` | 2001 | 1.1465 / 0.16057 | — | — | — |
| `eval_chebys` | 2001 | 0.0 / 0.0 | — | — | — |
| `eval_chebyt` | 2001 | 0.88106 / 0.13637 | — | — | — |
| `eval_chebyu` | 2001 | 0.0 / 0.0 | — | — | — |
| `eval_gegenbauer` | 2001 | 3.574 / 0.983 | — | — | — |
| `eval_genlaguerre` | 2001 | 2.173 / 0.15088 | — | — | — |
| `eval_hermite` | 2001 | 0.0 / 0.0 | — | — | — |
| `eval_hermitenorm` | 2001 | 0.0 / 0.0 | — | — | — |
| `eval_jacobi` | 2001 | 4.5826 / 0.94428 | — | — | — |
| `eval_laguerre` | 2001 | 0.54563 / 0.09915 | — | — | — |
| `eval_legendre` | 2001 | 3.574 / 0.9598 | — | — | — |
| `eval_sh_chebyt` | 2001 | 0.53003 / 0.13579 | — | — | — |
| `eval_sh_chebyu` | 2001 | 0.0 / 0.0 | — | — | — |
| `eval_sh_jacobi` | 2001 | 10.352 / 6.3012 | — | — | — |
| `eval_sh_legendre` | 2001 | 3.408 / 0.69851 | — | — | — |
| `exp1` | 2001 | 0.0 / 0.0 | 1.7884 / 0.47006 | 1.7884 / 0.437 | — |
| `exp10` | 2001 | 0.63318 / 0.11726 | — | — | — |
| `exp2` | 2001 | 0.70735 / 0.18092 | — | — | — |
| `expi` | 2001 | 0.0 / 0.0 | 5.5144 / 0.99096 | 4.9629 / 0.96261 | — |
| `expit` | 2001 | 0.0 / 0.0 | — | — | — |
| `expm1` | 2001 | 1.9093 / 0.55906 | — | — | — |
| `expn` | 2001 | 0.0 / 0.0 | — | — | — |
| `exprel` | 2001 | 0.0 / 0.0 | — | — | — |
| `factorial` | 2001 | 4.2785 / 0.93799 | — | — | — |
| `factorial2` | 2001 | 1.0667 / 0.23223 | — | — | — |
| `factorialk` | 2001 | 5.4857 / 1.1432 | — | — | — |
| `fdtr` | 2001 | 4.2371 / 0.78464 | — | — | — |
| `fdtrc` | 2001 | 4.3164 / 1.1371 | — | — | — |
| `fdtri` | 2001 | 5.5325 / 1.3301 | — | — | — |
| `fdtridfd` | 2001 | 1.0855e+5 / 954.14 | — | — | — |
| `fresnel` | 2001 | 0.0 / 0.0 | — | — | — |
| `fresnel_zeros` | 2001 | 1.977 / 0.21988 | — | — | — |
| `fresnelc_zeros` | 2001 | 1.977 / 0.2654 | — | — | — |
| `fresnels_zeros` | 2001 | 1.8671 / 0.17796 | — | — | — |
| `gamma` | 2087 | 4.1771 / 0.84064 | 3.6888 / 0.65275 | 17.917 / 1.0401; failed 1/2087 | 2.9747 / 0.54139 |
| `gammainc` | 2001 | 0.0 / 0.0 | 4.1696 / 0.9599 | 3.3464 / 2.1961 | — |
| `gammaincc` | 2001 | 0.0 / 0.0 | 17.937 / 4.34 | 23.662 / 5.3599 | — |
| `gammainccinv` | 2001 | 0.0 / 0.0 | — | — | — |
| `gammaincinv` | 2001 | 0.0 / 0.0 | — | — | — |
| `gammaln` | 2001 | 1.7588 / 0.29268 | 0.88356 / 0.085589 | 48.921 / 20.063 | — |
| `gammasgn` | 2001 | 0.0 / 0.0 | — | — | — |
| `gdtr` | 2001 | 0.0 / 0.0 | — | — | — |
| `gdtrc` | 2001 | 0.0 / 0.0 | — | — | — |
| `gdtria` | 2001 | 0.0 / 0.0 | — | — | — |
| `gdtrib` | 2001 | 0.0 / 0.0 | — | — | — |
| `gdtrix` | 2001 | 0.0 / 0.0 | — | — | — |
| `gegenbauer` | 2003 | 7.3411 / 0.91612 | — | — | — |
| `gen_harmonic` | 2010 | 0.0 / 0.0 | — | — | — |
| `genlaguerre` | 2003 | 3.6596 / 0.41779 | — | — | — |
| `geterr` | 2000 | — | — | — | — |
| `h1vp` | 2001 | 7.5141 / 2.4821 | — | — | — |
| `h2vp` | 2001 | 7.5141 / 2.5024 | — | — | — |
| `hankel1` | 2001 | 2.1375 / 0.27458 | — | — | — |
| `hankel1e` | 2001 | 4.5036e+15 / 1.4404e+14 | — | — | — |
| `hankel2` | 2001 | 2.1354 / 0.27584 | — | — | — |
| `hankel2e` | 2001 | 4.5036e+15 / 1.3729e+14 | — | — | — |
| `harmonic` | 2006 | 0.81633 / 0.36478 | — | — | — |
| `hermite` | 2005 | 4.5036e+15 / 1.5805e+15 | — | — | — |
| `hermitenorm` | 2005 | 4.5036e+15 / 1.5584e+15 | — | — | — |
| `heuman_lambda` | 2001 | 4.9099 / 1.3411 | — | — | — |
| `huber` | 2001 | 0.0 / 0.0 | — | — | — |
| `hyp0f1` | 2001 | 0.0 / 0.0 | — | — | — |
| `hyp1f1` | 2001 | 0.73483 / 0.010646 | — | — | — |
| `hyp2f1` | 2001 | 0.0 / 0.0 | — | — | — |
| `hyperu` | 2001 | 0.0 / 0.0 | — | — | — |
| `i0` | 2001 | 0.0 / 0.0 | — | — | — |
| `i0e` | 2001 | 0.0 / 0.0 | — | — | — |
| `i1` | 2001 | 0.0 / 0.0 | — | — | — |
| `i1e` | 2001 | 0.0 / 0.0 | — | — | — |
| `inv_boxcox` | 2001 | 0.0 / 0.0 | — | — | — |
| `inv_boxcox1p` | 2001 | 0.0 / 0.0 | — | — | — |
| `it2i0k0` | 2001 | 0.0 / 0.0 | — | — | — |
| `it2j0y0` | 2001 | 0.0 / 0.0 | — | — | — |
| `it2struve0` | 2001 | 0.0 / 0.0 | — | — | — |
| `itairy` | 2001 | 0.0 / 0.0 | — | — | — |
| `iti0k0` | 2001 | 0.0 / 0.0 | — | — | — |
| `itj0y0` | 2001 | 0.0 / 0.0 | — | — | — |
| `itmodstruve0` | 2001 | 0.0 / 0.0 | — | — | — |
| `itstruve0` | 2001 | 0.0 / 0.0 | — | — | — |
| `iv` | 2001 | 0.0 / 0.0 | 8.8057 / 1.4094 | 13.698 / 4.9097 | — |
| `ive` | 2001 | 0.0 / 0.0 | — | — | — |
| `ivp` | 2001 | 0.0 / 0.0 | — | — | — |
| `j0` | 2001 | 0.0 / 0.0 | — | — | — |
| `j1` | 2001 | 0.0 / 0.0 | — | — | — |
| `jacobi` | 2003 | 38.39 / 5.3536 | — | — | — |
| `jacobi_zeta` | 2001 | 10.482 / 2.5691 | — | — | — |
| `jn` | 2001 | 0.0 / 0.0 | — | — | — |
| `jn_zeros` | 2004 | 2.7694 / 0.82563 | — | — | — |
| `jnjnp_zeros` | 2001 | 1.9042 / 0.20715 | — | — | — |
| `jnp_zeros` | 2004 | 2.1725 / 0.81284 | — | — | — |
| `jnyn_zeros` | 2001 | 2.1725 / 0.75254 | — | — | — |
| `jv` | 2001 | 0.0 / 0.0 | 7.8081 / 5.6399 | 12.687 / 10.404 | — |
| `jve` | 2001 | 0.0 / 0.0 | — | — | — |
| `jvp` | 2001 | 1.6738 / 0.0989 | — | — | — |
| `k0` | 2001 | 0.0 / 0.0 | — | — | — |
| `k0e` | 2001 | 0.0 / 0.0 | — | — | — |
| `k1` | 2001 | 0.0 / 0.0 | — | — | — |
| `k1e` | 2001 | 0.0 / 0.0 | — | — | — |
| `kei` | 2001 | 0.0 / 0.0 | — | — | — |
| `kei_zeros` | 2001 | 0.0 / 0.0 | — | — | — |
| `keip` | 2001 | 0.0 / 0.0 | — | — | — |
| `keip_zeros` | 2005 | 2.8127 / 1.0218 | — | — | — |
| `kelvin` | 2001 | 0.0 / 0.0 | — | — | — |
| `kelvin_zeros` | 2002 | 2.8127 / 0.14347 | — | — | — |
| `ker` | 2001 | 0.0 / 0.0 | — | — | — |
| `ker_zeros` | 2001 | 0.0 / 0.0 | — | — | — |
| `kerp` | 2001 | 0.0 / 0.0 | — | — | — |
| `kerp_zeros` | 2001 | 0.0 / 0.0 | — | — | — |
| `kl_div` | 2001 | 8.705 / 2.5692 | — | — | — |
| `kn` | 2001 | 2.6488 / 0.31225 | — | — | — |
| `kolmogi` | 2001 | 0.0 / 0.0 | — | — | — |
| `kolmogorov` | 2001 | 0.0 / 0.0 | — | — | — |
| `kv` | 2001 | 1.5672 / 0.24594 | 2.0891 / 0.46986 | 30.182 / 4.157 | — |
| `kve` | 2001 | 1.1289 / 0.24278 | — | — | — |
| `kvp` | 2001 | 2.4777 / 0.3032 | — | — | — |
| `laguerre` | 2005 | 8.0 / 1.2267 | — | — | — |
| `lambertw` | 2001 | 2.8438 / 0.61451 | — | — | — |
| `lambertwm1` | 2001 | 2.7966 / 0.73794 | — | — | — |
| `legendre` | 2005 | 4.5036e+15 / 1.7391e+15 | — | — | — |
| `legendre_p` | 2001 | 0.0 / 0.0 | — | — | — |
| `legendre_p_all` | 2001 | 0.0 / 0.0 | — | — | — |
| `lmbda` | 2001 | 0.0 / 0.0 | — | — | — |
| `ln_binomial` | 2007 | 47081.0 / 23.464 | — | — | — |
| `ln_factorial` | 2007 | 1.4189 / 0.34424 | — | — | — |
| `ln_multinomial` | 2006 | 699.67 / 111.84 | — | — | — |
| `log1p` | 2001 | 0.0 / 0.0 | — | — | — |
| `log_expit` | 2001 | 0.0 / 0.0 | — | — | — |
| `log_ndtr` | 2001 | 0.0 / 0.0 | — | — | — |
| `log_softmax` | 2001 | 0.0 / 0.0 | — | — | — |
| `log_wright_bessel` | 2001 | 0.0 / 0.0 | — | — | — |
| `loggamma` | 2001 | 0.0 / 0.0 | — | — | — |
| `logit` | 2001 | 1.7691 / 0.40078 | — | — | — |
| `logsumexp` | 2001 | 0.51156 / 2.5565e-4 | — | — | — |
| `lpmv` | 2001 | 0.0 / 0.0 | — | — | — |
| `lqmn` | 2004 | 0.0 / 0.0 | — | — | — |
| `lqn` | 2001 | 0.0 / 0.0 | — | — | — |
| `mathieu_a` | 2001 | 0.0 / 0.0 | — | — | — |
| `mathieu_b` | 2001 | 0.0 / 0.0 | — | — | — |
| `mathieu_cem` | 2001 | 0.0 / 0.0 | — | — | — |
| `mathieu_even_coef` | 2001 | 0.0 / 0.0 | — | — | — |
| `mathieu_modcem1` | 2001 | 0.0 / 0.0 | — | — | — |
| `mathieu_modcem2` | 2001 | 0.0 / 0.0 | — | — | — |
| `mathieu_modsem1` | 2001 | 0.0 / 0.0 | — | — | — |
| `mathieu_modsem2` | 2001 | 0.0 / 0.0 | — | — | — |
| `mathieu_odd_coef` | 2001 | 0.0 / 0.0 | — | — | — |
| `mathieu_sem` | 2001 | 0.0 / 0.0 | — | — | — |
| `modfresnelm` | 2001 | 0.0 / 0.0 | — | — | — |
| `modfresnelp` | 2001 | 0.0 / 0.0 | — | — | — |
| `modstruve` | 2001 | 0.0 / 0.0 | — | — | — |
| `multigammaln` | 2001 | 0.0 / 0.0 | — | — | — |
| `multinomial` | 2006 | 8776.0 / 1371.5 | — | — | — |
| `nbdtr` | 2001 | 0.0 / 0.0 | — | — | — |
| `nbdtrc` | 2001 | 0.0 / 0.0 | — | — | — |
| `nbdtri` | 2001 | 0.0 / 0.0 | — | — | — |
| `nbdtrik` | 2001 | 12.598 / 3.0836 | — | — | — |
| `nbdtrin` | 2001 | 4.6169 / 1.052 | — | — | — |
| `ncfdtr` | 2001 | 3.8155 / 0.79819 | — | — | — |
| `ncfdtri` | 2001 | 5.9921 / 1.2441 | — | — | — |
| `ncfdtridfd` | 2001 | 1.9925e+11 / 1.9815e+11 | — | — | — |
| `ncfdtridfn` | 2001 | 1.524e+10 / 1.5088e+10 | — | — | — |
| `ncfdtrinc` | 2001 | 4.6295e+10 / 4.562e+10 | — | — | — |
| `nctdtr` | 2037 | 4.2427e+10 / 2.1083e+7 | — | — | — |
| `nctdtridf` | 2001 | 2.7026e+7 / 2.6696e+7 | — | — | — |
| `nctdtrinc` | 2001 | 1.2499e+6 / 1.2391e+6 | — | — | — |
| `nctdtrit` | 2001 | 43.516 / 10.793 | — | — | — |
| `ndtr` | 2001 | 0.0 / 0.0 | — | — | — |
| `ndtri` | 2001 | 0.0 / 0.0 | — | — | — |
| `ndtri_exp` | 2001 | 0.57617 / 1.7271e-3 | — | — | — |
| `nrdtrimn` | 2001 | 0.0 / 0.0 | — | — | — |
| `nrdtrisd` | 2001 | 0.0 / 0.0 | — | — | — |
| `obl_ang1` | 2001 | 0.0 / 0.0 | — | — | — |
| `obl_ang1_cv` | 2001 | 0.0 / 0.0 | — | — | — |
| `obl_cv` | 2001 | 0.0 / 0.0 | — | — | — |
| `obl_cv_seq` | 2001 | 37.334 / 5.1976 | — | — | — |
| `obl_rad1` | 2001 | 0.0 / 0.0 | — | — | — |
| `obl_rad1_cv` | 2001 | 0.0 / 0.0 | — | — | — |
| `obl_rad2` | 2001 | 0.0 / 0.0 | — | — | — |
| `obl_rad2_cv` | 2001 | 0.0 / 0.0 | — | — | — |
| `owens_t` | 2001 | 0.0 / 0.0 | — | — | — |
| `pbdn_seq` | 2001 | 0.7873 / 0.09299 | — | — | — |
| `pbdv` | 2001 | 0.0 / 0.0 | — | — | — |
| `pbdv_seq` | 2001 | 0.0 / 0.0 | — | — | — |
| `pbvv` | 2001 | 0.0 / 0.0 | — | — | — |
| `pbvv_seq` | 2001 | 4.7193 / 0.76015 | — | — | — |
| `pbwa` | 2001 | 0.0 / 0.0 | — | — | — |
| `pdtr` | 2001 | 0.0 / 0.0 | — | — | — |
| `pdtrc` | 2001 | 0.0 / 0.0 | — | — | — |
| `pdtri` | 2001 | 0.0 / 0.0 | — | — | — |
| `pdtrik` | 2001 | 6.2386 / 1.3086 | — | — | — |
| `perm` | 2001 | 0.0 / 0.0 | — | — | — |
| `poch` | 2001 | 0.0 / 0.0 | — | — | — |
| `polygamma` | 2001 | 0.0 / 0.0 | — | — | — |
| `powm1` | 2001 | 0.0 / 0.0 | — | — | — |
| `pro_ang1` | 2001 | 0.0 / 0.0 | — | — | — |
| `pro_ang1_cv` | 2001 | 0.0 / 0.0 | — | — | — |
| `pro_cv` | 2001 | 0.0 / 0.0 | — | — | — |
| `pro_cv_seq` | 2001 | 35.371 / 4.9756 | — | — | — |
| `pro_rad1` | 2001 | 0.0 / 0.0 | — | — | — |
| `pro_rad1_cv` | 2001 | 0.0 / 0.0 | — | — | — |
| `pro_rad2` | 2001 | 0.0 / 0.0 | — | — | — |
| `pro_rad2_cv` | 2001 | 0.0 / 0.0 | — | — | — |
| `pseudo_huber` | 2001 | 1.6202 / 0.40591 | — | — | — |
| `psi` | 2001 | 6.406 / 1.5338 | — | — | — |
| `radian` | 2001 | 0.94435 / 0.18972 | — | — | — |
| `rel_entr` | 2001 | 0.0 / 0.0 | — | — | — |
| `rgamma` | 2001 | 0.0 / 0.0 | — | — | — |
| `riccati_jn` | 2001 | 15.919 / 2.8225 | — | — | — |
| `riccati_yn` | 2001 | 3.3448 / 0.54361 | — | — | — |
| `roots_chebyc` | 2001 | 4.5036e+15 / 1.1458e+14 | — | — | — |
| `roots_chebys` | 2001 | 1.6474e+16 / 2.6089e+14 | — | — | — |
| `roots_chebyt` | 2001 | 4.5036e+15 / 1.1502e+14 | — | — | — |
| `roots_chebyu` | 2001 | 1.6474e+16 / 2.6038e+14 | — | — | — |
| `roots_gegenbauer` | 2001 | 4.5036e+15 / 4.9477e+14 | — | — | — |
| `roots_genlaguerre` | 2001 | 21.72 / 2.1716 | — | — | — |
| `roots_hermite` | 2001 | 4.5036e+15 / 1.027e+14 | — | — | — |
| `roots_hermitenorm` | 2001 | 4.5036e+15 / 1.1496e+14 | — | — | — |
| `roots_jacobi` | 2001 | 49.914 / 3.9343 | — | — | — |
| `roots_laguerre` | 2001 | 4379.2 / 51.388 | — | — | — |
| `roots_legendre` | 2001 | 4.5036e+15 / 1.1386e+14 | — | — | — |
| `roots_sh_chebyt` | 2001 | 830.69 / 20.158 | — | — | — |
| `roots_sh_chebyu` | 2001 | 176.68 / 11.645 | — | — | — |
| `roots_sh_jacobi` | 2001 | 17.07 / 2.1388 | — | — | — |
| `roots_sh_legendre` | 2001 | 399.94 / 17.093 | — | — | — |
| `round` | 2001 | 0.0 / 0.0 | — | — | — |
| `seterr` | 2001 | — | — | — | — |
| `sh_chebyt` | 2003 | 58.0 / 2.1017 | — | — | — |
| `sh_chebyu` | 2003 | 9.6 / 0.64961 | — | — | — |
| `sh_jacobi` | 2003 | 15.985 / 6.281 | — | — | — |
| `sh_legendre` | 2003 | 510.0 / 8.8597 | — | — | — |
| `shichi` | 2001 | 0.0 / 0.0 | — | — | — |
| `sici` | 2001 | 0.0 / 0.0 | — | — | — |
| `sinc` | 2001 | 0.0 / 0.0 | — | — | — |
| `sindg` | 2001 | 0.99993 / 0.031487 | — | — | — |
| `smirnov` | 2001 | 0.0 / 0.0 | — | — | — |
| `smirnovi` | 2001 | 0.0 / 0.0 | — | — | — |
| `softmax` | 2001 | 0.0 / 0.0 | — | — | — |
| `softplus` | 2001 | 0.0 / 0.0 | — | — | — |
| `spence` | 2001 | 0.0 / 0.0 | — | — | — |
| `sph_harm_y` | 2001 | 0.0 / 0.0 | — | — | — |
| `sph_harm_y_all` | 2001 | 0.0 / 0.0 | — | — | — |
| `sph_legendre_p` | 2001 | 0.0 / 0.0 | — | — | — |
| `sph_legendre_p_all` | 2001 | 0.0 / 0.0 | — | — | — |
| `spherical_in` | 2001 | 0.0 / 0.0 | — | — | — |
| `spherical_jn` | 2001 | 0.0 / 0.0 | — | — | — |
| `spherical_kn` | 2001 | 3.5448 / 0.45645 | — | — | — |
| `spherical_yn` | 2001 | 0.0 / 0.0 | — | — | — |
| `stdtr` | 2001 | 4.4493 / 1.4397 | — | — | — |
| `stdtridf` | 2001 | 179.73 / 55.412 | — | — | — |
| `stdtrit` | 2001 | 9.4201 / 2.2158 | — | — | — |
| `stirling2` | 2001 | 0.0 / 0.0 | — | — | — |
| `struve` | 2001 | 0.0 / 0.0 | — | — | — |
| `tandg` | 2001 | 1.7314 / 0.31981 | — | — | — |
| `tklmbda` | 2001 | 0.0 / 0.0 | — | — | — |
| `trigamma` | 2005 | 3.5692 / 1.3033 | 2.0418 / 0.95601 | 3.0612 / 1.1846 | — |
| `voigt_profile` | 2001 | 0.0 / 0.0 | — | — | — |
| `wofz` | 2001 | 0.0 / 0.0 | — | — | — |
| `wright_bessel` | 2001 | 0.0 / 0.0 | — | — | — |
| `wrightomega` | 2001 | 0.0 / 0.0 | — | — | — |
| `xlog1py` | 2001 | 0.0 / 0.0 | — | — | — |
| `xlogy` | 2001 | 0.0 / 0.0 | — | — | — |
| `y0` | 2001 | 0.0 / 0.0 | — | — | — |
| `y0_zeros` | 2001 | 1.7537e+7 / 1.2955e+6 | — | — | — |
| `y1` | 2001 | 0.0 / 0.0 | — | — | — |
| `y1_zeros` | 2001 | 1.0344e+7 / 7.5412e+5 | — | — | — |
| `y1p_zeros` | 2001 | 3.2631e+7 / 2.2931e+6 | — | — | — |
| `yn` | 2001 | 5.1374 / 0.93603 | — | — | — |
| `yn_zeros` | 2004 | 2.2382 / 0.69188 | — | — | — |
| `ynp_zeros` | 2004 | 1.9465 / 0.87439 | — | — | — |
| `yv` | 2001 | 2.1341 / 0.29492 | 3.7236 / 1.1416 | 10.148 / 2.1179 | — |
| `yve` | 2001 | 2.1372 / 0.30725 | — | — | — |
| `yvp` | 2001 | 1.6774 / 0.22764 | — | — | — |
| `zeta` | 2001 | 0.0 / 0.0 | — | — | — |
| `zetac` | 2001 | 0.0 / 0.0 | 4.3971 / 1.0808 | 23.403 / 5.9072 | — |

## Evaluation failures

- `special.gamma` (GSL): 1 failed; overflow. Example: `special.gamma/expanded-0076`.


## Conventions and limitations

All public special-function implementations have typed registrations or an explicit policy recipe. Original erf/gamma recipes use high-precision mpmath; added recipes use independent public SciPy calls, identities or quadrature. Polynomial constructors return descending coefficients; sequence and tuple ordering follows the pinned fixtures. Boost.Math/GSL are called where the bridge has matching scalar APIs, converting elliptic parameter m to modulus sqrt(m).

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
