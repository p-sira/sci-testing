// Public Boost.Math and GSL special-function oracles missing from the base ABI.
// The argument order follows the saved SciPy-style special.* fixture inputs.
#include <boost/math/special_functions/cbrt.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/hypergeometric_1F1.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_psi.h>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <string>

namespace {
thread_local std::string last_error;

bool named(const char* actual, const char* expected) {
    return actual && std::strcmp(actual, expected) == 0;
}

void require(bool valid, const char* message) {
    if (!valid) throw std::domain_error(message);
}

int arity(int provider, const char* name) {
    if (provider != 0 && provider != 1) return 0;
    if (named(name, "psi") || named(name, "expm1") || named(name, "log1p"))
        return 1;
    if (named(name, "assoc_laguerre")) return 3;
    if (provider == 0 && named(name, "cbrt")) return 1;
    return 0;
}

double boost_value(const char* name, const double* a) {
    if (named(name, "psi")) return boost::math::digamma(a[0]);
    if (named(name, "expm1")) return boost::math::expm1(a[0]);
    if (named(name, "log1p")) {
        require(a[0] > -1., "log1p requires x > -1 for finite results");
        return boost::math::log1p(a[0]);
    }
    if (named(name, "cbrt")) return boost::math::cbrt(a[0]);
    if (named(name, "assoc_laguerre")) {
        // SciPy assoc_laguerre(x,n,k) = L_n^(k)(x). Boost.Math's
        // three-argument laguerre accepts an integer k, whereas the saved
        // fixtures use noninteger k, so use Boost's public 1F1 and gamma.
        const int n = static_cast<int>(a[1]);
        const double k = a[2];
        const double normalization = boost::math::tgamma_ratio(k + n + 1., k + 1.)
                                   / boost::math::tgamma(n + 1.);
        return normalization * boost::math::hypergeometric_1F1(-n, k + 1., a[0]);
    }
    throw std::invalid_argument("unknown Boost scalar special function");
}

double gsl_value(const char* name, const double* a) {
    gsl_sf_result result{};
    int status = GSL_EINVAL;
    if (named(name, "psi")) status = gsl_sf_psi_e(a[0], &result);
    else if (named(name, "expm1")) status = gsl_sf_expm1_e(a[0], &result);
    else if (named(name, "log1p")) {
        require(a[0] > -1., "log1p requires x > -1 for finite results");
        status = gsl_sf_log_1plusx_e(a[0], &result);
    } else if (named(name, "assoc_laguerre"))
        status = gsl_sf_laguerre_n_e(static_cast<int>(a[1]), a[2], a[0], &result);
    require(status == GSL_SUCCESS, gsl_strerror(status));
    return result.val;
}
} // namespace

extern "C" {
const char* extra_special_scalar_error() { return last_error.c_str(); }

int extra_special_scalar_supported(int provider, const char* name) {
    return arity(provider, name);
}

int extra_special_scalar_evaluate(int provider, const char* name,
                                  const double* arguments, int count,
                                  double* output, int capacity) {
    last_error.clear();
    try {
        require(arity(provider, name) == count && count > 0,
                "unsupported scalar special function or arity");
        require(arguments && output && capacity >= 1,
                "invalid scalar special input or output buffer");
        for (int i = 0; i < count; ++i)
            require(std::isfinite(arguments[i]), "nonfinite special input");
        if (named(name, "assoc_laguerre")) {
            require(arguments[1] >= 0. && arguments[1] <= 128. &&
                    std::trunc(arguments[1]) == arguments[1],
                    "Laguerre degree must be an integer in [0,128]");
            require(arguments[2] > -1., "Laguerre parameter must be > -1");
        }
        // Handle GSL domain errors via explicit return codes and preserve
        // the same exception boundary used by other independent adapters.
        if (provider == 1) gsl_set_error_handler_off();
        const double value = provider == 0
            ? boost_value(name, arguments) : gsl_value(name, arguments);
        require(std::isfinite(value), "nonfinite scalar special result");
        output[0] = value;
        return 0;
    } catch (const std::exception& e) { last_error = e.what(); }
      catch (...) { last_error = "unknown scalar special failure"; }
    return -1;
}
} // extern C
