// C ABI for independent Boost.Math and GSL comparisons. No exceptions cross it.
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/airy.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/jacobi_elliptic.hpp>
#include <boost/math/special_functions/jacobi_zeta.hpp>
#include <boost/math/special_functions/heuman_lambda.hpp>
#include <boost/math/special_functions/ellint_d.hpp>
#include <boost/math/special_functions/ellint_rc.hpp>
#include <boost/math/special_functions/ellint_rd.hpp>
#include <boost/math/special_functions/ellint_rf.hpp>
#include <boost/math/special_functions/ellint_rg.hpp>
#include <boost/math/special_functions/ellint_rj.hpp>
#include <boost/math/special_functions/gegenbauer.hpp>
#include <boost/math/special_functions/chebyshev.hpp>
#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/special_functions/laguerre.hpp>
#include <boost/math/special_functions/jacobi.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/special_functions/hypergeometric_0F1.hpp>
#include <boost/math/special_functions/hypergeometric_1F1.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>
#include <boost/math/special_functions/lambert_w.hpp>
#include <boost/math/special_functions/owens_t.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/distributions/logistic.hpp>
#include <boost/math/special_functions/powm1.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/zeta.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/version.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_elljac.h>
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_sf_hermite.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_sf_mathieu.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_version.h>
#include <cmath>
#include <array>
#include <cstring>
#include <stdexcept>
#include <string>

namespace {
thread_local std::string error;
// One concise registration gives both implementations, their arity and Sci's
// parameter convention. Elliptic arguments convert parameter m to modulus k.
#define SCALAR_FUNCTIONS(X) \
 X(erf, 1, boost::math::erf(a[0]), gsl_sf_erf_e(a[0], &r)) \
 X(erfc, 1, boost::math::erfc(a[0]), gsl_sf_erfc_e(a[0], &r)) \
 X(gamma, 1, boost::math::tgamma(a[0]), gsl_sf_gamma_e(a[0], &r)) \
 X(gammaln, 1, boost::math::lgamma(a[0]), gsl_sf_lngamma_e(a[0], &r)) \
 X(loggamma, 1, boost::math::lgamma(a[0]), gsl_sf_lngamma_e(a[0], &r)) \
 X(gammasgn, 1, ([&] {int s=0; boost::math::lgamma(a[0], &s); return static_cast<double>(s);}()), ([&] {double s=0.; int status=gsl_sf_lngamma_sgn_e(a[0], &r, &s); r.val=s; return status;}())) \
 X(factorial, 1, boost::math::tgamma(a[0]+1.), gsl_sf_gamma_e(a[0]+1., &r)) \
 X(factorial2, 1, boost::math::double_factorial<double>(static_cast<unsigned>(a[0])), gsl_sf_doublefact_e(static_cast<unsigned>(a[0]), &r)) \
 X(ln_factorial, 1, boost::math::lgamma(a[0]+1.), gsl_sf_lnfact_e(static_cast<unsigned>(a[0]), &r)) \
 X(ln_binomial, 2, boost::math::lgamma(a[0]+1.)-boost::math::lgamma(a[1]+1.)-boost::math::lgamma(a[0]-a[1]+1.), gsl_sf_lnchoose_e(static_cast<unsigned>(a[0]), static_cast<unsigned>(a[1]), &r)) \
 X(binom, 2, boost::math::tgamma(a[0]+1.)/(boost::math::tgamma(a[1]+1.)*boost::math::tgamma(a[0]-a[1]+1.)), (r.val = std::exp(gsl_sf_lngamma(a[0]+1.)-gsl_sf_lngamma(a[1]+1.)-gsl_sf_lngamma(a[0]-a[1]+1.)), GSL_SUCCESS)) \
 X(comb, 2, boost::math::tgamma(a[0]+1.)/(boost::math::tgamma(a[1]+1.)*boost::math::tgamma(a[0]-a[1]+1.)), (r.val = std::exp(gsl_sf_lngamma(a[0]+1.)-gsl_sf_lngamma(a[1]+1.)-gsl_sf_lngamma(a[0]-a[1]+1.)), GSL_SUCCESS)) \
 X(poch, 2, boost::math::tgamma_ratio(a[0]+a[1], a[0]), gsl_sf_poch_e(a[0], a[1], &r)) \
 X(gamma_ratio, 2, boost::math::tgamma_ratio(a[0], a[1]), gsl_sf_poch_e(a[1], a[0]-a[1], &r)) \
 X(digamma, 1, boost::math::digamma(a[0]), gsl_sf_psi_e(a[0], &r)) \
 X(trigamma, 1, boost::math::trigamma(a[0]), gsl_sf_psi_1_e(a[0], &r)) \
 X(beta, 2, boost::math::beta(a[0], a[1]), gsl_sf_beta_e(a[0], a[1], &r)) \
 X(betaln, 2, std::log(boost::math::beta(a[0], a[1])), gsl_sf_lnbeta_e(a[0], a[1], &r)) \
 X(gammainc, 2, boost::math::gamma_p(a[0], a[1]), gsl_sf_gamma_inc_P_e(a[0], a[1], &r)) \
 X(gammaincc, 2, boost::math::gamma_q(a[0], a[1]), gsl_sf_gamma_inc_Q_e(a[0], a[1], &r)) \
 X(betainc, 3, boost::math::ibeta(a[0], a[1], a[2]), gsl_sf_beta_inc_e(a[0], a[1], a[2], &r)) \
 X(jv, 2, boost::math::cyl_bessel_j(a[0], a[1]), gsl_sf_bessel_Jnu_e(a[0], a[1], &r)) \
 X(yv, 2, boost::math::cyl_neumann(a[0], a[1]), gsl_sf_bessel_Ynu_e(a[0], a[1], &r)) \
 X(iv, 2, boost::math::cyl_bessel_i(a[0], a[1]), gsl_sf_bessel_Inu_e(a[0], a[1], &r)) \
 X(kv, 2, boost::math::cyl_bessel_k(a[0], a[1]), gsl_sf_bessel_Knu_e(a[0], a[1], &r)) \
 X(ellipk, 1, boost::math::ellint_1(std::sqrt(a[0])), gsl_sf_ellint_Kcomp_e(std::sqrt(a[0]), GSL_PREC_DOUBLE, &r)) \
 X(ellipe, 1, boost::math::ellint_2(std::sqrt(a[0])), gsl_sf_ellint_Ecomp_e(std::sqrt(a[0]), GSL_PREC_DOUBLE, &r)) \
 X(exp1, 1, boost::math::expint(1, a[0]), gsl_sf_expint_E1_e(a[0], &r)) \
 X(expi, 1, boost::math::expint(a[0]), gsl_sf_expint_Ei_e(a[0], &r)) \
 X(zetac, 1, boost::math::zeta(a[0])-1, gsl_sf_zetam1_e(a[0], &r)) \
 X(i0, 1, boost::math::cyl_bessel_i(0, a[0]), gsl_sf_bessel_I0_e(a[0], &r)) \
 X(i1, 1, boost::math::cyl_bessel_i(1, a[0]), gsl_sf_bessel_I1_e(a[0], &r)) \
 X(j0, 1, boost::math::cyl_bessel_j(0, a[0]), gsl_sf_bessel_J0_e(a[0], &r)) \
 X(j1, 1, boost::math::cyl_bessel_j(1, a[0]), gsl_sf_bessel_J1_e(a[0], &r)) \
 X(k0, 1, boost::math::cyl_bessel_k(0, a[0]), gsl_sf_bessel_K0_e(a[0], &r)) \
 X(k1, 1, boost::math::cyl_bessel_k(1, a[0]), gsl_sf_bessel_K1_e(a[0], &r)) \
 X(y0, 1, boost::math::cyl_neumann(0, a[0]), gsl_sf_bessel_Y0_e(a[0], &r)) \
 X(y1, 1, boost::math::cyl_neumann(1, a[0]), gsl_sf_bessel_Y1_e(a[0], &r)) \
 X(kn, 2, boost::math::cyl_bessel_k(static_cast<int>(a[0]), a[1]), gsl_sf_bessel_Kn_e(static_cast<int>(a[0]), a[1], &r)) \
 X(jn, 2, boost::math::cyl_bessel_j(static_cast<int>(a[0]), a[1]), gsl_sf_bessel_Jn_e(static_cast<int>(a[0]), a[1], &r)) \
 X(yn, 2, boost::math::cyl_neumann(static_cast<int>(a[0]), a[1]), gsl_sf_bessel_Yn_e(static_cast<int>(a[0]), a[1], &r)) \
 X(i0e, 1, static_cast<double>(boost::math::cyl_bessel_i(0, static_cast<long double>(a[0]))*std::exp(-std::abs(static_cast<long double>(a[0])))), gsl_sf_bessel_I0_scaled_e(a[0], &r)) \
 X(i1e, 1, static_cast<double>(boost::math::cyl_bessel_i(1, static_cast<long double>(a[0]))*std::exp(-std::abs(static_cast<long double>(a[0])))), gsl_sf_bessel_I1_scaled_e(a[0], &r)) \
 X(k0e, 1, static_cast<double>(boost::math::cyl_bessel_k(0, static_cast<long double>(a[0]))*std::exp(static_cast<long double>(a[0]))), gsl_sf_bessel_K0_scaled_e(a[0], &r)) \
 X(k1e, 1, static_cast<double>(boost::math::cyl_bessel_k(1, static_cast<long double>(a[0]))*std::exp(static_cast<long double>(a[0]))), gsl_sf_bessel_K1_scaled_e(a[0], &r)) \
 X(ive, 2, static_cast<double>(boost::math::cyl_bessel_i(static_cast<long double>(a[0]), static_cast<long double>(a[1]))*std::exp(-std::abs(static_cast<long double>(a[1])))), gsl_sf_bessel_Inu_scaled_e(a[0], a[1], &r)) \
 X(kve, 2, static_cast<double>(boost::math::cyl_bessel_k(static_cast<long double>(a[0]), static_cast<long double>(a[1]))*std::exp(static_cast<long double>(a[1]))), gsl_sf_bessel_Knu_scaled_e(a[0], a[1], &r)) \
 X(jve, 2, boost::math::cyl_bessel_j(a[0], a[1]), gsl_sf_bessel_Jnu_e(a[0], a[1], &r)) \
 X(yve, 2, boost::math::cyl_neumann(a[0], a[1]), gsl_sf_bessel_Ynu_e(a[0], a[1], &r)) \
 X(spherical_jn, 3, ([&] {int n=static_cast<int>(a[0]); double x=a[1]; if (!a[2]) return boost::math::sph_bessel(n,x); if(n==0) return -boost::math::sph_bessel(1,x); if(x==0.) return n==1 ? 1./3. : 0.; return boost::math::sph_bessel(n-1,x)-(n+1.)/x*boost::math::sph_bessel(n,x);}()), ([&] {int n=static_cast<int>(a[0]); int code=gsl_sf_bessel_jl_e(n,a[1],&r); if (code || !a[2]) return code; if(n==0) {gsl_sf_result next{}; code=gsl_sf_bessel_jl_e(1,a[1],&next); r.val=-next.val; return code;} if(a[1]==0.) {r.val=n==1 ? 1./3. : 0.; return static_cast<int>(GSL_SUCCESS);} gsl_sf_result prev{}; code=gsl_sf_bessel_jl_e(n-1,a[1],&prev); r.val=prev.val-(n+1.)/a[1]*r.val; return code;}())) \
 X(spherical_yn, 3, ([&] {int n=static_cast<int>(a[0]); double x=a[1]; if (!a[2]) return boost::math::sph_neumann(n,x); if(n==0) return -boost::math::sph_neumann(1,x); return boost::math::sph_neumann(n-1,x)-(n+1.)/x*boost::math::sph_neumann(n,x);}()), ([&] {int n=static_cast<int>(a[0]); int code=gsl_sf_bessel_yl_e(n,a[1],&r); if (code || !a[2]) return code; if(n==0) {gsl_sf_result next{}; code=gsl_sf_bessel_yl_e(1,a[1],&next); r.val=-next.val; return code;} gsl_sf_result prev{}; code=gsl_sf_bessel_yl_e(n-1,a[1],&prev); r.val=prev.val-(n+1.)/a[1]*r.val; return code;}())) \
 X(elliprc, 2, boost::math::ellint_rc(a[0], a[1]), gsl_sf_ellint_RC_e(a[0], a[1], GSL_PREC_DOUBLE, &r)) \
 X(elliprd, 3, boost::math::ellint_rd(a[0], a[1], a[2]), gsl_sf_ellint_RD_e(a[0], a[1], a[2], GSL_PREC_DOUBLE, &r)) \
 X(elliprf, 3, boost::math::ellint_rf(a[0], a[1], a[2]), gsl_sf_ellint_RF_e(a[0], a[1], a[2], GSL_PREC_DOUBLE, &r)) \
 X(elliprj, 4, boost::math::ellint_rj(a[0], a[1], a[2], a[3]), gsl_sf_ellint_RJ_e(a[0], a[1], a[2], a[3], GSL_PREC_DOUBLE, &r)) \
 X(ellipkinc, 2, boost::math::ellint_1(std::sqrt(a[1]), a[0]), gsl_sf_ellint_F_e(a[0], std::sqrt(a[1]), GSL_PREC_DOUBLE, &r)) \
 X(ellipf, 2, boost::math::ellint_1(std::sqrt(a[1]), a[0]), gsl_sf_ellint_F_e(a[0], std::sqrt(a[1]), GSL_PREC_DOUBLE, &r)) \
 X(ellipeinc, 2, boost::math::ellint_2(std::sqrt(a[1]), a[0]), gsl_sf_ellint_E_e(a[0], std::sqrt(a[1]), GSL_PREC_DOUBLE, &r)) \
 X(ellippi, 2, boost::math::ellint_3(std::sqrt(a[1]), a[0]), gsl_sf_ellint_Pcomp_e(std::sqrt(a[1]), -a[0], GSL_PREC_DOUBLE, &r)) \
 X(ellippiinc, 3, boost::math::ellint_3(std::sqrt(a[2]), a[1], a[0]), gsl_sf_ellint_P_e(a[0], std::sqrt(a[2]), -a[1], GSL_PREC_DOUBLE, &r)) \
 X(ellippiinc_bulirsch, 3, boost::math::ellint_3(std::sqrt(a[2]), a[1], a[0]), gsl_sf_ellint_P_e(a[0], std::sqrt(a[2]), -a[1], GSL_PREC_DOUBLE, &r)) \
 X(ellipd, 1, boost::math::ellint_d(std::sqrt(a[0])), gsl_sf_ellint_Dcomp_e(std::sqrt(a[0]), GSL_PREC_DOUBLE, &r)) \
 X(ellipdinc, 2, boost::math::ellint_d(std::sqrt(a[1]), a[0]), gsl_sf_ellint_D_e(a[0], std::sqrt(a[1]), GSL_PREC_DOUBLE, &r)) \
 X(betaincc, 3, boost::math::ibetac(a[0], a[1], a[2]), (r.val = gsl_cdf_beta_Q(a[2], a[0], a[1]), GSL_SUCCESS)) \
 X(betaincinv, 3, boost::math::ibeta_inv(a[0], a[1], a[2]), (r.val = gsl_cdf_beta_Pinv(a[2], a[0], a[1]), GSL_SUCCESS)) \
 X(betainccinv, 3, boost::math::ibetac_inv(a[0], a[1], a[2]), (r.val = gsl_cdf_beta_Qinv(a[2], a[0], a[1]), GSL_SUCCESS)) \
 X(gammaincinv, 2, boost::math::gamma_p_inv(a[0], a[1]), (r.val = gsl_cdf_gamma_Pinv(a[1], a[0], 1.), GSL_SUCCESS)) \
 X(gammainccinv, 2, boost::math::gamma_q_inv(a[0], a[1]), (r.val = gsl_cdf_gamma_Qinv(a[1], a[0], 1.), GSL_SUCCESS)) \
 X(hyp0f1, 2, boost::math::hypergeometric_0F1(a[0], a[1]), gsl_sf_hyperg_0F1_e(a[0], a[1], &r)) \
 X(hyp1f1, 3, boost::math::hypergeometric_1F1(a[0], a[1], a[2]), gsl_sf_hyperg_1F1_e(a[0], a[1], a[2], &r)) \
 X(hyp2f1, 4, boost::math::hypergeometric_pFq({a[0], a[1]}, {a[2]}, a[3]), gsl_sf_hyperg_2F1_e(a[0], a[1], a[2], a[3], &r)) \
 X(eval_legendre, 2, boost::math::legendre_p(static_cast<int>(a[0]), a[1]), gsl_sf_legendre_Pl_e(static_cast<int>(a[0]), a[1], &r)) \
 X(eval_sh_legendre, 2, boost::math::legendre_p(static_cast<int>(a[0]), 2.*a[1]-1.), gsl_sf_legendre_Pl_e(static_cast<int>(a[0]), 2.*a[1]-1., &r)) \
 X(eval_hermite, 2, boost::math::hermite(static_cast<unsigned>(a[0]), a[1]), gsl_sf_hermite_phys_e(static_cast<int>(a[0]), a[1], &r)) \
 X(eval_hermitenorm, 2, std::pow(2., -a[0]/2.)*boost::math::hermite(static_cast<unsigned>(a[0]), a[1]/std::sqrt(2.)), gsl_sf_hermite_prob_e(static_cast<int>(a[0]), a[1], &r)) \
 X(eval_gegenbauer, 3, boost::math::gegenbauer(static_cast<unsigned>(a[0]), a[1], a[2]), gsl_sf_gegenpoly_n_e(static_cast<int>(a[0]), a[1], a[2], &r)) \
 X(eval_laguerre, 2, boost::math::laguerre(static_cast<unsigned>(a[0]), a[1]), gsl_sf_laguerre_n_e(static_cast<int>(a[0]), 0., a[1], &r)) \
 X(eval_genlaguerre, 3, boost::math::hypergeometric_1F1(-a[0], a[1]+1., a[2])*boost::math::tgamma_ratio(a[1]+a[0]+1., a[1]+1.)/boost::math::tgamma(a[0]+1.), gsl_sf_laguerre_n_e(static_cast<int>(a[0]), a[1], a[2], &r)) \
 X(polygamma, 2, boost::math::polygamma(static_cast<int>(a[0]), a[1]), gsl_sf_psi_n_e(static_cast<int>(a[0]), a[1], &r)) \
 X(rgamma, 1, 1./boost::math::tgamma(a[0]), gsl_sf_gammainv_e(a[0], &r)) \
 X(expn, 2, boost::math::expint(static_cast<unsigned>(a[0]), a[1]), gsl_sf_expint_En_e(static_cast<int>(a[0]), a[1], &r)) \
 X(lambertw, 1, boost::math::lambert_w0(a[0]), gsl_sf_lambert_W0_e(a[0], &r)) \
 X(lambertwm1, 1, boost::math::lambert_wm1(a[0]), gsl_sf_lambert_Wm1_e(a[0], &r)) \
 X(ellipkm1, 1, boost::math::ellint_rf(0., a[0], 1.), gsl_sf_ellint_RF_e(0., a[0], 1., GSL_PREC_DOUBLE, &r)) \
 X(log1pmx, 1, boost::math::log1pmx(a[0]), gsl_sf_log_1plusx_mx_e(a[0], &r)) \
 X(expit, 1, boost::math::cdf(boost::math::logistic_distribution<double>(), a[0]), (r.val = gsl_cdf_logistic_P(a[0], 1.), GSL_SUCCESS)) \
 X(logit, 1, boost::math::quantile(boost::math::logistic_distribution<double>(), a[0]), (r.val = gsl_cdf_logistic_Pinv(a[0], 1.), GSL_SUCCESS)) \
 X(chdtr, 2, boost::math::cdf(boost::math::chi_squared_distribution<double>(a[0]), a[1]), (r.val = gsl_cdf_chisq_P(a[1], a[0]), GSL_SUCCESS)) \
 X(chdtrc, 2, boost::math::cdf(boost::math::complement(boost::math::chi_squared_distribution<double>(a[0]), a[1])), (r.val = gsl_cdf_chisq_Q(a[1], a[0]), GSL_SUCCESS)) \
 X(chdtri, 2, boost::math::quantile(boost::math::complement(boost::math::chi_squared_distribution<double>(a[0]), a[1])), (r.val = gsl_cdf_chisq_Qinv(a[1], a[0]), GSL_SUCCESS)) \
 X(stdtr, 2, boost::math::cdf(boost::math::students_t_distribution<double>(a[0]), a[1]), (r.val = gsl_cdf_tdist_P(a[1], a[0]), GSL_SUCCESS)) \
 X(stdtrit, 2, boost::math::quantile(boost::math::students_t_distribution<double>(a[0]), a[1]), (r.val = gsl_cdf_tdist_Pinv(a[1], a[0]), GSL_SUCCESS)) \
 X(fdtr, 3, boost::math::cdf(boost::math::fisher_f_distribution<double>(a[0], a[1]), a[2]), (r.val = gsl_cdf_fdist_P(a[2], a[0], a[1]), GSL_SUCCESS)) \
 X(fdtrc, 3, boost::math::cdf(boost::math::complement(boost::math::fisher_f_distribution<double>(a[0], a[1]), a[2])), (r.val = gsl_cdf_fdist_Q(a[2], a[0], a[1]), GSL_SUCCESS)) \
 X(fdtri, 3, boost::math::quantile(boost::math::fisher_f_distribution<double>(a[0], a[1]), a[2]), (r.val = gsl_cdf_fdist_Pinv(a[2], a[0], a[1]), GSL_SUCCESS)) \
 X(gdtr, 3, boost::math::cdf(boost::math::gamma_distribution<double>(a[1], 1./a[0]), a[2]), (r.val = gsl_cdf_gamma_P(a[2], a[1], 1./a[0]), GSL_SUCCESS)) \
 X(gdtrc, 3, boost::math::cdf(boost::math::complement(boost::math::gamma_distribution<double>(a[1], 1./a[0]), a[2])), (r.val = gsl_cdf_gamma_Q(a[2], a[1], 1./a[0]), GSL_SUCCESS)) \
 X(gdtrix, 3, boost::math::quantile(boost::math::gamma_distribution<double>(a[1], 1./a[0]), a[2]), (r.val = gsl_cdf_gamma_Pinv(a[2], a[1], 1./a[0]), GSL_SUCCESS)) \
 X(pdtr, 2, boost::math::cdf(boost::math::poisson_distribution<double>(a[1]), a[0]), (r.val = gsl_cdf_poisson_P(static_cast<unsigned>(a[0]), a[1]), GSL_SUCCESS)) \
 X(pdtrc, 2, boost::math::cdf(boost::math::complement(boost::math::poisson_distribution<double>(a[1]), a[0])), (r.val = gsl_cdf_poisson_Q(static_cast<unsigned>(a[0]), a[1]), GSL_SUCCESS)) \
 X(pdtri, 2, boost::math::gamma_q_inv(a[0]+1., a[1]), (r.val = gsl_cdf_gamma_Qinv(a[1], a[0]+1., 1.), GSL_SUCCESS)) \
 X(bdtr, 3, boost::math::cdf(boost::math::binomial_distribution<double>(a[1], a[2]), std::floor(a[0])), (r.val = gsl_cdf_binomial_P(static_cast<unsigned>(a[0]), a[2], static_cast<unsigned>(a[1])), GSL_SUCCESS)) \
 X(bdtrc, 3, boost::math::cdf(boost::math::complement(boost::math::binomial_distribution<double>(a[1], a[2]), std::floor(a[0]))), (r.val = gsl_cdf_binomial_Q(static_cast<unsigned>(a[0]), a[2], static_cast<unsigned>(a[1])), GSL_SUCCESS)) \
 X(bdtri, 3, boost::math::ibetac_inv(a[0]+1., a[1]-a[0], a[2]), (r.val = gsl_cdf_beta_Qinv(a[2], a[0]+1., a[1]-a[0]), GSL_SUCCESS)) \
 X(ndtr, 1, boost::math::cdf(boost::math::normal_distribution<double>(), a[0]), (r.val = gsl_cdf_ugaussian_P(a[0]), GSL_SUCCESS)) \
 X(ndtri, 1, boost::math::quantile(boost::math::normal_distribution<double>(), a[0]), (r.val = gsl_cdf_ugaussian_Pinv(a[0]), GSL_SUCCESS))

// Separate registries prevent presenting a missing implementation as a failure.
// Every entry is an independently callable public library API.
#define BOOST_ONLY_FUNCTIONS(X) \
 X(erfinv, 1, boost::math::erf_inv(a[0])) \
 X(erfcinv, 1, boost::math::erfc_inv(a[0])) \
 X(elliprg, 3, boost::math::ellint_rg(a[0], a[1], a[2])) \
 X(eval_chebyt, 2, boost::math::chebyshev_t(static_cast<unsigned>(a[0]), a[1])) \
 X(eval_chebyu, 2, boost::math::chebyshev_u(static_cast<unsigned>(a[0]), a[1])) \
 X(eval_jacobi, 4, boost::math::jacobi(static_cast<unsigned>(a[0]), a[1], a[2], a[3])) \
 X(eval_sh_chebyt, 2, boost::math::chebyshev_t(static_cast<unsigned>(a[0]), 2.*a[1]-1.)) \
 X(eval_sh_chebyu, 2, boost::math::chebyshev_u(static_cast<unsigned>(a[0]), 2.*a[1]-1.)) \
 X(jacobi_zeta, 2, boost::math::jacobi_zeta(std::sqrt(a[1]), a[0])) \
 X(heuman_lambda, 2, boost::math::heuman_lambda(std::sqrt(a[1]), a[0])) \
 X(owens_t, 2, boost::math::owens_t(a[0], a[1])) \
 X(powm1, 2, boost::math::powm1(a[0], a[1]))

#define GSL_ONLY_FUNCTIONS(X) \
 X(hyperu, 3, gsl_sf_hyperg_U_e(a[0], a[1], a[2], &r)) \
 X(zeta, 2, gsl_sf_hzeta_e(a[0], a[1], &r)) \
 X(mathieu_a, 2, gsl_sf_mathieu_a_e(static_cast<int>(a[0]), a[1], &r)) \
 X(mathieu_b, 2, gsl_sf_mathieu_b_e(static_cast<int>(a[0]), a[1], &r)) \
 X(scaled_exp1, 1, a[0] == 0. ? (r.val = 0., GSL_SUCCESS) : (gsl_sf_expint_E1_scaled_e(a[0], &r) == GSL_SUCCESS ? (r.val *= a[0], GSL_SUCCESS) : GSL_EFAILED)) \
 X(dawsn, 1, gsl_sf_dawson_e(a[0], &r)) \
 X(spence, 1, gsl_sf_dilog_e(1.-a[0], &r))
}
extern "C" {
const char* comparison_version(int library) { return library == 0 ? BOOST_LIB_VERSION : gsl_version; }
const char* comparison_error() { return error.c_str(); }
int comparison_supported(const char* name) {
    if (std::strcmp(name, "airy") == 0 || std::strcmp(name, "airye") == 0 ||
        std::strcmp(name, "sici") == 0 || std::strcmp(name, "shichi") == 0) return 1;
    if (std::strcmp(name, "ellipke") == 0 || std::strcmp(name, "ellipj") == 0) return std::strcmp(name, "ellipj") == 0 ? 2 : 1;
#define SUPPORTED(name_, arity, boost_, gsl_) if (std::strcmp(name, #name_) == 0) return arity;
    SCALAR_FUNCTIONS(SUPPORTED)
#undef SUPPORTED
#define SUPPORTED_BOOST(name_, arity, boost_) if (std::strcmp(name, #name_) == 0) return arity;
    BOOST_ONLY_FUNCTIONS(SUPPORTED_BOOST)
#undef SUPPORTED_BOOST
#define SUPPORTED_GSL(name_, arity, gsl_) if (std::strcmp(name, #name_) == 0) return arity;
    GSL_ONLY_FUNCTIONS(SUPPORTED_GSL)
#undef SUPPORTED_GSL
    if (std::strcmp(name, "normal_cdf") == 0 || std::strcmp(name, "normal_pdf") == 0 || std::strcmp(name, "normal_ppf") == 0) return 1;
    return 0;
}
int comparison_library_supported(int library, const char* name) {
    if (library != 0 && library != 1) return 0;
    if (std::strcmp(name, "airy") == 0 || std::strcmp(name, "airye") == 0) return 1;
    if (std::strcmp(name, "ellipke") == 0 || std::strcmp(name, "ellipj") == 0) return 1;
    if (library == 1 && (std::strcmp(name, "sici") == 0 || std::strcmp(name, "shichi") == 0)) return 1;
    if (std::strcmp(name, "normal_cdf") == 0 || std::strcmp(name, "normal_pdf") == 0 || std::strcmp(name, "normal_ppf") == 0) return 1;
#define HAS_BOTH(name_, arity, boost_, gsl_) if (std::strcmp(name, #name_) == 0) return 1;
    SCALAR_FUNCTIONS(HAS_BOTH)
#undef HAS_BOTH
#define HAS_BOOST(name_, arity, boost_) if (library == 0 && std::strcmp(name, #name_) == 0) return 1;
    BOOST_ONLY_FUNCTIONS(HAS_BOOST)
#undef HAS_BOOST
#define HAS_GSL(name_, arity, gsl_) if (library == 1 && std::strcmp(name, #name_) == 0) return 1;
    GSL_ONLY_FUNCTIONS(HAS_GSL)
#undef HAS_GSL
    return 0;
}
int comparison_output_count(const char* name) {
    if (std::strcmp(name, "airy") == 0 || std::strcmp(name, "airye") == 0) return 4;
    if (std::strcmp(name, "ellipke") == 0) return 2;
    if (std::strcmp(name, "ellipj") == 0) return 4;
    if (std::strcmp(name, "sici") == 0 || std::strcmp(name, "shichi") == 0) return 2;
    return comparison_supported(name) ? 1 : 0;
}
// Public GSL complex special functions. The real scalar ABI above is kept
// separate because complex inputs and outputs both contain two doubles.
int comparison_complex_supported(const char* name) {
    if (!name) return 0;
    return (std::strcmp(name, "digamma") == 0 ||
            std::strcmp(name, "spence") == 0) ? 2 : 0;
}
int comparison_complex_evaluate(int library, const char* name, const double* a, int count,
                                double* result, int capacity) {
    error.clear();
    gsl_set_error_handler_off();
    try {
        if (library != 1) throw std::invalid_argument("complex function unavailable for this library");
        if (!name || !a || !result) throw std::invalid_argument("null complex argument");
        if (!comparison_complex_supported(name)) throw std::invalid_argument("unsupported complex function");
        if (count != 2 || capacity < 2) throw std::invalid_argument("incorrect complex arity or output capacity");
        if (!std::isfinite(a[0]) || !std::isfinite(a[1]))
            throw std::invalid_argument("nonfinite complex input");
        gsl_sf_result re{}, im{};
        const int status = std::strcmp(name, "digamma") == 0
            ? gsl_sf_complex_psi_e(a[0], a[1], &re, &im)
            : gsl_sf_complex_spence_xy_e(a[0], a[1], &re, &im);
        if (status != GSL_SUCCESS) { error = gsl_strerror(status); return status; }
        if (!std::isfinite(re.val) || !std::isfinite(im.val))
            throw std::runtime_error("nonfinite complex output");
        result[0] = re.val;
        result[1] = im.val;
        return 0;
    } catch (const std::exception& e) { error = e.what(); return -1; }
      catch (...) { error = "unknown native exception"; return -1; }
}
// Multi-component real Airy comparisons. The scalar ABI remains unchanged.
int comparison_evaluate_many(int library, const char* name, const double* a, int count,
                             double* result, int capacity) {
    error.clear();
    gsl_set_error_handler_off();
    try {
        if (library != 0 && library != 1) throw std::invalid_argument("unknown library");
        if (std::strcmp(name, "ellipke") == 0) {
            if (count != 1 || capacity < 2) throw std::invalid_argument("incorrect elliptic K/E capacity");
            const double k = std::sqrt(a[0]);
            if (library == 0) {
                result[0] = boost::math::ellint_1(k);
                result[1] = boost::math::ellint_2(k);
            } else {
                gsl_sf_result ek{}, ee{};
                int code = gsl_sf_ellint_Kcomp_e(k, GSL_PREC_DOUBLE, &ek);
                if (code) { error = gsl_strerror(code); return code; }
                code = gsl_sf_ellint_Ecomp_e(k, GSL_PREC_DOUBLE, &ee);
                if (code) { error = gsl_strerror(code); return code; }
                result[0] = ek.val;
                result[1] = ee.val;
            }
            if (!std::isfinite(result[0]) || !std::isfinite(result[1]))
                throw std::runtime_error("nonfinite output");
            return 0;
        }
        if (std::strcmp(name, "ellipj") == 0) {
            if (count != 2 || capacity < 4) throw std::invalid_argument("incorrect Jacobi output capacity");
            const double u = a[0], m = a[1];
            if (!(m >= 0. && m < 1.)) throw std::invalid_argument("real Jacobi parameter requires 0 <= m < 1");
            const double k = std::sqrt(m);
            const double quarter = library == 0 ? boost::math::ellint_1(k)
                : gsl_sf_ellint_Kcomp(k, GSL_PREC_DOUBLE);
            const double cycles = std::floor((u + 2.*quarter)/(4.*quarter));
            const double reduced_u = u - 4.*quarter*cycles;
            double sn{}, cn{}, dn{};
            if (library == 0) {
                sn = boost::math::jacobi_elliptic(k, reduced_u, &cn, &dn);
                // Boost can return a noticeably inconsistent dn precisely at
                // quarter-period arguments. Preserve the Jacobi identity.
                dn = std::sqrt(std::max(0., 1. - m * sn * sn));
            } else {
                const int code = gsl_sf_elljac_e(reduced_u, m, &sn, &cn, &dn);
                if (code) { error = gsl_strerror(code); return code; }
            }
            const double amplitude = std::atan2(sn, cn) + 2.*std::acos(-1.)*cycles;
            const std::array<double, 4> values{sn, cn, dn, amplitude};
            for (size_t i = 0; i < values.size(); ++i) {
                if (!std::isfinite(values[i])) throw std::runtime_error("nonfinite Jacobi output");
                result[i] = values[i];
            }
            return 0;
        }
        if (std::strcmp(name, "sici") == 0 || std::strcmp(name, "shichi") == 0) {
            if (library != 1 || count != 1 || capacity < 2)
                throw std::invalid_argument("incorrect GSL integral arity or output capacity");
            gsl_sf_result first{}, second{};
            const bool hyperbolic = std::strcmp(name, "shichi") == 0;
            const int first_status = hyperbolic ? gsl_sf_Shi_e(a[0], &first) : gsl_sf_Si_e(a[0], &first);
            if (first_status != GSL_SUCCESS) { error = gsl_strerror(first_status); return first_status; }
            const int second_status = hyperbolic ? gsl_sf_Chi_e(a[0], &second) : gsl_sf_Ci_e(a[0], &second);
            if (second_status != GSL_SUCCESS) { error = gsl_strerror(second_status); return second_status; }
            if (!std::isfinite(first.val) || !std::isfinite(second.val))
                throw std::runtime_error("nonfinite output");
            result[0] = first.val;
            result[1] = second.val;
            return 0;
        }
        if (count != 1 || comparison_output_count(name) != 4 || capacity < 4)
            throw std::invalid_argument("incorrect arity or output capacity");
        const double x = a[0];
        if (!std::isfinite(x)) throw std::invalid_argument("nonfinite Airy input");
        const bool scaled = std::strcmp(name, "airye") == 0;
        // Sci's real airye has no real scaled Ai or Ai' for x < 0.
        // GSL instead returns their unscaled real values on this domain.
        if (scaled && x < 0) throw std::invalid_argument("real airye requires x >= 0");

        std::array<double, 4> values{}; // Ai, Ai', Bi, Bi'
        if (library == 0) {
            if (scaled) {
                // Boost provides only unscaled functions. Use extended-range
                // intermediates so exp(zeta) does not overflow before the
                // exponentially small Ai values can cancel its magnitude.
                const long double y = x;
                const long double zeta = (2.0L / 3.0L) * y * std::sqrt(y);
                const long double ai_scale = std::exp(zeta);
                const long double bi_scale = std::exp(-zeta);
                values = {static_cast<double>(boost::math::airy_ai(y) * ai_scale),
                          static_cast<double>(boost::math::airy_ai_prime(y) * ai_scale),
                          static_cast<double>(boost::math::airy_bi(y) * bi_scale),
                          static_cast<double>(boost::math::airy_bi_prime(y) * bi_scale)};
            } else {
                values = {boost::math::airy_ai(x), boost::math::airy_ai_prime(x),
                          boost::math::airy_bi(x), boost::math::airy_bi_prime(x)};
            }
        } else {
            using AiryFunction = int (*)(double, gsl_mode_t, gsl_sf_result*);
            const std::array<AiryFunction, 4> functions = scaled
                ? std::array<AiryFunction, 4>{gsl_sf_airy_Ai_scaled_e,
                                              gsl_sf_airy_Ai_deriv_scaled_e,
                                              gsl_sf_airy_Bi_scaled_e,
                                              gsl_sf_airy_Bi_deriv_scaled_e}
                : std::array<AiryFunction, 4>{gsl_sf_airy_Ai_e,
                                              gsl_sf_airy_Ai_deriv_e,
                                              gsl_sf_airy_Bi_e,
                                              gsl_sf_airy_Bi_deriv_e};
            for (size_t i = 0; i < values.size(); ++i) {
                gsl_sf_result component{};
                const int status = functions[i](x, GSL_PREC_DOUBLE, &component);
                if (status != GSL_SUCCESS) { error = gsl_strerror(status); return status; }
                values[i] = component.val;
            }
        }
        for (double value : values)
            if (!std::isfinite(value)) throw std::runtime_error("nonfinite output");
        for (size_t i = 0; i < values.size(); ++i) result[i] = values[i];
        return 0;
    } catch (const std::exception& e) { error = e.what(); return -1; }
      catch (...) { error = "unknown native exception"; return -1; }
}
int comparison_evaluate(int library, const char* name, const double* a, int count, double* result) {
    error.clear();
    gsl_set_error_handler_off();
    try {
        if (library != 0 && library != 1) throw std::invalid_argument("unknown library");
        if (count == 0 || count != comparison_supported(name)) throw std::invalid_argument("incorrect arity");
        if (!comparison_library_supported(library, name)) throw std::invalid_argument("unsupported library function");
        gsl_sf_result r{};
        int status = GSL_SUCCESS;
#define EVALUATE(name_, arity, boost_, gsl_) \
        if (std::strcmp(name, #name_) == 0) { \
            if (library == 0) *result = (boost_); \
            else { status = (gsl_); *result = r.val; } \
        } else
        SCALAR_FUNCTIONS(EVALUATE)
#define EVALUATE_BOOST(name_, arity, boost_) \
        if (library == 0 && std::strcmp(name, #name_) == 0) { *result = (boost_); } else
        BOOST_ONLY_FUNCTIONS(EVALUATE_BOOST)
#undef EVALUATE_BOOST
#define EVALUATE_GSL(name_, arity, gsl_) \
        if (library == 1 && std::strcmp(name, #name_) == 0) { status = (gsl_); *result = r.val; } else
        GSL_ONLY_FUNCTIONS(EVALUATE_GSL)
#undef EVALUATE_GSL
        if (std::strcmp(name, "normal_cdf") == 0) {
            *result = library == 0 ? boost::math::cdf(boost::math::normal_distribution<double>(), a[0]) : gsl_cdf_ugaussian_P(a[0]);
        } else if (std::strcmp(name, "normal_pdf") == 0) {
            *result = library == 0 ? boost::math::pdf(boost::math::normal_distribution<double>(), a[0]) : gsl_ran_ugaussian_pdf(a[0]);
        } else if (std::strcmp(name, "normal_ppf") == 0) {
            *result = library == 0 ? boost::math::quantile(boost::math::normal_distribution<double>(), a[0]) : gsl_cdf_ugaussian_Pinv(a[0]);
        } else throw std::invalid_argument("unsupported function");
#undef EVALUATE
        if (status != GSL_SUCCESS) { error = gsl_strerror(status); return status; }
        if (!std::isfinite(*result)) throw std::runtime_error("nonfinite output");
        return 0;
    } catch (const std::exception& e) { error = e.what(); return -1; }
      catch (...) { error = "unknown native exception"; return -1; }
}
}
