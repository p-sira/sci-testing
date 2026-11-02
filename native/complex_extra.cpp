// Standalone GSL complex Gamma / analytic log-Gamma / reciprocal Gamma bridge.
//
// Link this translation unit with -lgsl -lgslcblas.  This C ABI is separate
// from comparison_complex_evaluate, so integration does not require changes
// to another worker's code.
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>

#include <cmath>
#include <complex>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
thread_local std::string last_error;

using complex_ld = std::complex<long double>;
constexpr long double pi = 3.141592653589793238462643383279502884L;
constexpr long double two_pi = 2.0L * pi;
constexpr long double half_log_two_pi = 0.9189385332046727417803297364056176399L;

bool is_function(const char* name) {
    return name && (!std::strcmp(name, "gamma") ||
                    !std::strcmp(name, "loggamma") ||
                    !std::strcmp(name, "rgamma"));
}

bool pole(double x, double y) {
    return y == 0.0 && x <= 0.0 && std::floor(x) == x;
}

struct gamma_data {
    complex_ld loggamma;
    complex_ld shifted_gamma;
    complex_ld product;
    bool product_finite;
};

// The Stirling imaginary part is only used to identify which multiple of
// 2*pi GSL discarded. It is never substituted for GSL's log magnitude.
// Shifting the input to Re >= 16 bounds the omitted expansion at <1e-13
// and makes branch selection robust even when Im(log(Gamma)) is very large.
long double stirling_imag(const complex_ld& z) {
    const complex_ld reciprocal = 1.0L / z;
    const complex_ld r2 = reciprocal * reciprocal;
    const complex_ld correction =
        reciprocal * (1.0L / 12.0L +
            r2 * (-1.0L / 360.0L +
            r2 * (1.0L / 1260.0L +
            r2 * (-1.0L / 1680.0L +
            r2 * (1.0L / 1188.0L)))));
    const complex_ld asymptotic =
        (z - 0.5L) * std::log(z) - z + half_log_two_pi + correction;
    return asymptotic.imag();
}

gamma_data analytic_loggamma(double x, double y) {
    // Work within a stated, bounded common domain; outside it, do not claim
    // a comparison from this GSL-based adapter.
    if (x < -2048.0 || x > 1.0e6 || std::abs(y) > 1.0e6)
        throw std::domain_error("complex Gamma outside validated GSL comparison range");
    if (pole(x, y))
        throw std::domain_error("Gamma pole at a nonpositive integer");

    const int shifts = x < 16.0 ? static_cast<int>(std::ceil(16.0 - x)) : 0;
    const complex_ld input(static_cast<long double>(x), static_cast<long double>(y));
    const complex_ld shifted = input + static_cast<long double>(shifts);
    gsl_sf_result log_magnitude{}, wrapped_phase{};
    const int code = gsl_sf_lngamma_complex_e(
        static_cast<double>(shifted.real()), static_cast<double>(shifted.imag()),
        &log_magnitude, &wrapped_phase);
    if (code != GSL_SUCCESS)
        throw std::runtime_error(std::string("gsl_sf_lngamma_complex_e: ") +
                                 gsl_strerror(code));
    if (!std::isfinite(log_magnitude.val) || !std::isfinite(wrapped_phase.val))
        throw std::runtime_error("GSL returned nonfinite log-Gamma components");

    const long double phase = static_cast<long double>(wrapped_phase.val);
    const long double winding = std::round((stirling_imag(shifted) - phase) / two_pi);
    complex_ld logarithm(
        static_cast<long double>(log_magnitude.val), phase + two_pi * winding);

    // Gamma(z+n) = Gamma(z) * product(z+j), j=0..n-1.
    // Principal logarithms of each factor preserve continuous analytic
    // log-Gamma and opposite sides of the negative-real branch cut.
    // Crucially, a tiny nonzero imaginary part near a negative integer is
    // retained in log(z+j), rather than being lost inside GSL's reflection.
    complex_ld product(1.0L, 0.0L);
    for (int j = 0; j < shifts; ++j) {
        const complex_ld factor = input + static_cast<long double>(j);
        logarithm -= std::log(factor);
        product *= factor;
    }
    // Signed-zero is not consistently preserved by std::complex arithmetic
    // on negative-real factors. The two boundary values of analytic
    // log-Gamma differ by the cut's winding number. Ordinary Gamma is
    // meromorphic and has no corresponding branch distinction.
    if (y == 0.0 && x < 0.0)
        logarithm.imag(std::copysign(pi * std::ceil(-static_cast<long double>(x)),
                                    std::signbit(y) ? 1.0L : -1.0L));
    else if (y == 0.0 && x > 0.0)
        logarithm.imag(0.0L);
    if (!std::isfinite(logarithm.real()) || !std::isfinite(logarithm.imag()))
        throw std::runtime_error("nonfinite reconstructed log-Gamma");
    const complex_ld shifted_gamma = std::exp(complex_ld(
        static_cast<long double>(log_magnitude.val), phase));
    const bool safe_product = std::isfinite(product.real()) &&
                              std::isfinite(product.imag());
    return {logarithm, shifted_gamma, product, safe_product};
}
}

extern "C" {

// Mirrors the existing complex ABI's input/output scalar width, not the
// count of complex numbers: 2 doubles in and 2 doubles out.
int complex_extra_supported(const char* name) {
    return is_function(name) ? 2 : 0;
}

const char* complex_extra_error() {
    return last_error.c_str();
}

// library = 1 denotes GSL, consistent with native/comparisons.cpp.
// Returns zero on success, negative one on invalid/unsupported/nonnumeric
// input, and the GSL error in the diagnostic string when relevant.
// The output buffer remains untouched on failure.
int complex_extra_evaluate(int library, const char* name,
                           const double* input, int count,
                           double* output, int capacity) {
    last_error.clear();
    try {
        if (library != 1 || !is_function(name))
            throw std::invalid_argument("complex extra API requires supported GSL function");
        if (!input || !output || count != 2 || capacity < 2)
            throw std::invalid_argument("incorrect complex input or output capacity");
        const double x = input[0], y = input[1];
        if (!std::isfinite(x) || !std::isfinite(y))
            throw std::invalid_argument("nonfinite complex input");
        if (pole(x, y)) {
            if (!std::strcmp(name, "rgamma")) {
                output[0] = 0.0;
                output[1] = 0.0;
                return 0;
            }
            throw std::domain_error("Gamma pole at a nonpositive integer");
        }

        // GSL's default error handler aborts the host process. Disable it
        // before any GSL special-function evaluation, as in the existing ABI.
        gsl_set_error_handler_off();
        const gamma_data computed = analytic_loggamma(x, y);
        complex_ld value;
        if (!std::strcmp(name, "loggamma")) {
            value = computed.loggamma;
        } else if (computed.product_finite &&
                   std::isfinite(computed.shifted_gamma.real()) &&
                   std::isfinite(computed.shifted_gamma.imag())) {
            // Direct recurrence preserves tiny components close to a pole;
            // exp(logGamma) would multiply a huge magnitude by the tiny
            // phase error of an integer multiple of pi.
            value = !std::strcmp(name, "rgamma")
                ? computed.product / computed.shifted_gamma
                : computed.shifted_gamma / computed.product;
        } else {
            // A large intermediate factorial may exceed long double even
            // when the final value is representable. The log form bounds it.
            value = std::exp(!std::strcmp(name, "rgamma")
                ? -computed.loggamma : computed.loggamma);
        }
        if (y == 0.0 && std::strcmp(name, "loggamma") != 0)
            value.imag(0.0L);
        const double real = static_cast<double>(value.real());
        const double imag = static_cast<double>(value.imag());
        if (!std::isfinite(real) || !std::isfinite(imag))
            throw std::overflow_error("complex Gamma result outside finite double range");
        output[0] = real;
        output[1] = imag;
        return 0;
    } catch (const std::exception& error) {
        last_error = error.what();
    } catch (...) {
        last_error = "unknown complex Gamma error";
    }
    return -1;
}
}
