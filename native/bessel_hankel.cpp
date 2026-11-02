// Complex-valued cylindrical Hankel functions for real positive arguments,
// assembled solely from the public real-argument Boost.Math / GSL Bessel J,Y.
//
// H^(1)_v(x) = J_v(x) + i Y_v(x); H^(2)_v(x) = J_v(x) - i Y_v(x).
// The e-suffixed functions multiply by exp(-i*x) or exp(+i*x), respectively.
// The n-th x derivative is 2^-n sum_{k=0}^n (-1)^k binom(n,k)
// H_{v-n+2*k}(x). No complex-argument support is implied by this bridge.
//
// The C ABI uses two or three real input scalars and returns two doubles
// [real, imaginary]. Library 0 is Boost.Math; library 1 is GSL.

#include <boost/math/special_functions/bessel.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>

#include <cmath>
#include <complex>
#include <cstring>
#include <stdexcept>
#include <string>

namespace {
thread_local std::string last_error;
constexpr double pi = 3.141592653589793238462643383279502884;

int arity(int library, const char* name) {
    if ((library != 0 && library != 1) || !name) return 0;
    if (!std::strcmp(name, "hankel1") || !std::strcmp(name, "hankel1e") ||
        !std::strcmp(name, "hankel2") || !std::strcmp(name, "hankel2e")) return 2;
    if (!std::strcmp(name, "h1vp") || !std::strcmp(name, "h2vp")) return 3;
    return 0;
}

std::complex<double> hankel(int library, bool first, double order, double x) {
    if (library == 0) {
        return {boost::math::cyl_bessel_j(order, x),
                (first ? 1.0 : -1.0) * boost::math::cyl_neumann(order, x)};
    }

    // Public GSL Jnu/Ynu require nonnegative orders. For a negative order
    // use H^(1)_{-a} = exp(i*pi*a) H^(1)_a and its H^(2) conjugate.
    const double positive = std::abs(order);
    gsl_sf_result j{}, y{};
    const int j_status = gsl_sf_bessel_Jnu_e(positive, x, &j);
    if (j_status != GSL_SUCCESS)
        throw std::runtime_error(std::string("GSL Jnu: ") + gsl_strerror(j_status));
    const int y_status = gsl_sf_bessel_Ynu_e(positive, x, &y);
    if (y_status != GSL_SUCCESS)
        throw std::runtime_error(std::string("GSL Ynu: ") + gsl_strerror(y_status));

    std::complex<double> value(j.val, (first ? 1.0 : -1.0) * y.val);
    if (order < 0.0) {
        // Remainder keeps trig phases accurate for ordinary finite orders.
        const double phase = pi * std::remainder(positive, 2.0);
        const double signed_phase = first ? phase : -phase;
        value *= std::complex<double>(std::cos(signed_phase), std::sin(signed_phase));
    }
    return value;
}

std::complex<double> derivative(int library, bool first, double order,
                                double x, int n) {
    std::complex<double> sum(0.0, 0.0);
    double coefficient = 1.0;
    for (int k = 0; k <= n; ++k) {
        const auto term = hankel(library, first, order - n + 2.0 * k, x);
        sum += (k & 1 ? -coefficient : coefficient) * term;
        if (k < n) coefficient *= double(n - k) / double(k + 1);
    }
    return std::ldexp(1.0, -n) * sum;
}
}  // namespace

extern "C" {

int bessel_hankel_supported(int library, const char* name) {
    return arity(library, name);
}

const char* bessel_hankel_error() {
    return last_error.c_str();
}

int bessel_hankel_evaluate(int library, const char* name,
                           const double* args, int count,
                           double* output, int capacity) {
    last_error.clear();
    try {
        const int expected = arity(library, name);
        if (!expected || !args || count != expected || !output || capacity < 2)
            throw std::invalid_argument("unsupported Hankel function or incorrect scalar input/output shape");
        const double order = args[0];
        const double x = args[1];
        // The comparison is defined on the shared real Bessel J/Y domain.
        // Bound order and x to the validated finite region; no analytic
        // continuation through the complex x branch cut is attempted.
        if (!std::isfinite(order) || !std::isfinite(x) ||
            std::abs(order) > 256.0 || x <= 0.0 || x > 1.0e6)
            throw std::domain_error("real Hankel comparison requires -256 <= order <= 256 and 0 < x <= 1e6");
        int n = 0;
        if (expected == 3) {
            if (!std::isfinite(args[2]) || args[2] != std::trunc(args[2]) ||
                args[2] < 0.0 || args[2] > 12.0)
                throw std::domain_error("Hankel derivative order must be an integer in [0,12]");
            n = static_cast<int>(args[2]);
        }

        gsl_set_error_handler_off();
        const bool first = !std::strcmp(name, "hankel1") ||
                           !std::strcmp(name, "hankel1e") ||
                           !std::strcmp(name, "h1vp");
        std::complex<double> value = derivative(library, first, order, x, n);
        if (!std::strcmp(name, "hankel1e") || !std::strcmp(name, "hankel2e")) {
            const double phase = std::remainder(x, 2.0 * pi);
            const double angle = first ? -phase : phase;
            value *= std::complex<double>(std::cos(angle), std::sin(angle));
        }
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
            throw std::overflow_error("nonfinite real-argument Hankel result");
        output[0] = value.real();
        output[1] = value.imag();
        return 0;
    } catch (const std::exception& error) {
        last_error = error.what();
    } catch (...) {
        last_error = "unknown real-argument Hankel evaluation error";
    }
    return -1;
}
}
