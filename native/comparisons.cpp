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
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/special_functions/zeta.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/version.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_ellint.h>
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
 X(zetac, 1, boost::math::zeta(a[0])-1, gsl_sf_zetam1_e(a[0], &r))
}
extern "C" {
const char* comparison_version(int library) { return library == 0 ? BOOST_LIB_VERSION : gsl_version; }
const char* comparison_error() { return error.c_str(); }
int comparison_supported(const char* name) {
    if (std::strcmp(name, "airy") == 0 || std::strcmp(name, "airye") == 0) return 1;
#define SUPPORTED(name_, arity, boost_, gsl_) if (std::strcmp(name, #name_) == 0) return arity;
    SCALAR_FUNCTIONS(SUPPORTED)
#undef SUPPORTED
    if (std::strcmp(name, "normal_cdf") == 0 || std::strcmp(name, "normal_pdf") == 0 || std::strcmp(name, "normal_ppf") == 0) return 1;
    return 0;
}
int comparison_output_count(const char* name) {
    if (std::strcmp(name, "airy") == 0 || std::strcmp(name, "airye") == 0) return 4;
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
        gsl_sf_result r{};
        int status = GSL_SUCCESS;
#define EVALUATE(name_, arity, boost_, gsl_) \
        if (std::strcmp(name, #name_) == 0) { \
            if (library == 0) *result = (boost_); \
            else { status = (gsl_); *result = r.val; } \
        } else
        SCALAR_FUNCTIONS(EVALUATE)
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
