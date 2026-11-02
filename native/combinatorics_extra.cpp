// Independent public Boost.Math/GSL comparisons for Sci combinatorial series.
// Bernoulli B_0..B_n is vector-valued; other functions return one real scalar.
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/bernoulli.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_zeta.h>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
thread_local std::string error;
constexpr int max_harmonic_n = 1000000;
// B_258 is finite in IEEE binary64; B_260 exceeds DBL_MAX.
constexpr int max_bernoulli_n = 258;

bool is_named(const char* name, const char* expected) {
    return name && std::strcmp(name, expected) == 0;
}

int arity(int library, const char* name) {
    if (library < 0 || library > 1) return 0;
    if (is_named(name, "bernoulli")) return library == 0 ? 1 : 0;
    if (is_named(name, "harmonic")) return 1;
    // gsl_sf_hzeta(s,q) is public only for s > 1 and q > 0. The identity
    // H_n^(s) = zeta(s,1) - zeta(s,n+1) covers that domain, not all real s.
    if (is_named(name, "gen_harmonic")) return library == 1 ? 2 : 0;
    return 0;
}

int checked_n(double value, int limit) {
    if (!std::isfinite(value) || value < 0 || value > limit ||
        std::floor(value) != value)
        throw std::invalid_argument("n must be a nonnegative integer in the finite supported range");
    return static_cast<int>(value);
}

void check_gsl(int status) {
    if (status != GSL_SUCCESS) throw std::runtime_error(gsl_strerror(status));
}

double gsl_hurwitz(double exponent, double q) {
    gsl_sf_result result{};
    const int status = gsl_sf_hzeta_e(exponent, q, &result);
    // A positive tail smaller than the smallest subnormal can be represented
    // as zero in binary64; GSL reports that expected underflow explicitly.
    if (status == GSL_EUNDRFLW && result.val == 0.) return 0.;
    check_gsl(status);
    return result.val;
}

double gsl_harmonic(int n) {
    if (n == 0) return 0.;
    gsl_sf_result value{};
    check_gsl(gsl_sf_psi_int_e(n + 1, &value));
    return value.val + boost::math::constants::euler<double>();
}
}  // namespace

extern "C" {
const char* combinatorics_extra_error() { return error.c_str(); }

int combinatorics_extra_supported(int library, const char* name) {
    return arity(library, name);
}

int combinatorics_extra_output_count(const char* name, const double* args, int count) {
    if (!name || !args || count < 1) return 0;
    if (is_named(name, "bernoulli")) {
        if (count != 1 || !std::isfinite(args[0]) || args[0] < 0 ||
            args[0] > max_bernoulli_n || std::floor(args[0]) != args[0]) return 0;
        return static_cast<int>(args[0]) + 1;
    }
    return 1;
}

int combinatorics_extra_evaluate(int library, const char* name,
                                 const double* args, int count,
                                 double* output, int capacity) {
    error.clear();
    gsl_set_error_handler_off();
    try {
        if (!args || !name || !output || !arity(library, name) ||
            count != arity(library, name))
            throw std::invalid_argument("unsupported library, function, or argument count");
        const int expected = combinatorics_extra_output_count(name, args, count);
        if (!expected || capacity < expected)
            throw std::invalid_argument("invalid output shape or capacity");
        if (is_named(name, "bernoulli")) {
            const int n = checked_n(args[0], max_bernoulli_n);
            // Boost returns even Bernoulli numbers B_(2k); odd B_(2k+1)
            // vanish except for B_1 = -1/2 under SciPy's convention.
            for (int i = 0; i <= n; ++i) {
                double value = i == 0 ? 1. : i == 1 ? -0.5 :
                    (i & 1) ? 0. : boost::math::bernoulli_b2n<double>(i / 2);
                if (!std::isfinite(value)) throw std::runtime_error("nonfinite Bernoulli coefficient");
                output[i] = value;
            }
        } else if (is_named(name, "harmonic")) {
            const int n = checked_n(args[0], max_harmonic_n);
            if (n == 0) {
                output[0] = 0.;
            } else if (library == 0) {
                output[0] = static_cast<double>(
                    boost::math::digamma(static_cast<long double>(n + 1)) +
                    boost::math::constants::euler<long double>());
            } else {
                output[0] = gsl_harmonic(n);
            }
        } else if (is_named(name, "gen_harmonic")) {
            const int n = checked_n(args[0], max_harmonic_n);
            const double s = args[1];
            if (!std::isfinite(s)) throw std::invalid_argument("exponent must be finite");
            if (n == 0) output[0] = 0.;
            else if (s == 1.) output[0] = gsl_harmonic(n);
            else {
                if (s < 1.) throw std::invalid_argument(
                    "GSL Hurwitz zeta representation requires exponent s > 1");
                const double total = gsl_hurwitz(s, 1.);
                const double tail = gsl_hurwitz(s, static_cast<double>(n) + 1.);
                output[0] = total - tail;
                // Two enormous zeta values can almost cancel near s = 1.
                // Reject when binary64 subtraction cannot retain the accuracy
                // needed by the comparison report.
                const double cancellation_bound = 8. * std::numeric_limits<double>::epsilon() *
                                                  (std::abs(total) + std::abs(tail));
                if (cancellation_bound > 1e-10 * std::abs(output[0]) + 1e-13)
                    throw std::runtime_error("loss of precision near Hurwitz zeta pole");
            }
        } else {
            throw std::invalid_argument("unsupported function");
        }
        for (int i = 0; i < expected; ++i)
            if (!std::isfinite(output[i])) throw std::runtime_error("nonfinite output");
        return 0;
    } catch (const std::exception& exception) {
        error = exception.what(); return -1;
    } catch (...) { error = "unknown native exception"; return -1; }
}
}
