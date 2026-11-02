// Public Boost.Math/GSL spherical Bessel functions adapted to the Riccati
// functions returned by scipy.special.riccati_jn and riccati_yn.
//
// ψ_k(x) = x*j_k(x), χ_k(x) = x*y_k(x) [SciPy uses the PLUS sign].
// For k>=1: ψ'_k = ψ_{k-1} - k*ψ_k/x, similarly for χ'_k.
// For order zero: ψ_0=sin(x), ψ'_0=cos(x),
//                 χ_0=-cos(x), χ'_0=sin(x).
//
// Returned vector: [value(0..n), derivative(0..n)] -- length 2*(n+1).
// Boost/GSL spherical functions are defined for nonnegative/positive x.
// For ψ and x<0, integer-order parity gives exact real continuation.
// For χ and x<=0, SciPy uses non-mathematical 1e300 sentinels; reject
// rather than claiming these come from Boost or GSL.
#include <boost/math/special_functions/bessel.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_result.h>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
thread_local std::string last_error;

bool supported(int library, const char* name) {
    if (library != 0 && library != 1) return false;
    if (!name) return false;
    return std::strcmp(name, "riccati_jn") == 0 ||
           std::strcmp(name, "riccati_yn") == 0;
}

int order(const double* parameters, int length) {
    if (!parameters || length != 2)
        throw std::invalid_argument("expected scalar (maximum_order, argument)");
    const double n = parameters[0], x = parameters[1];
    if (!std::isfinite(n) || n < 0 || n > 512 || std::trunc(n) != n)
        throw std::invalid_argument("maximum order must be an integer in [0, 512]");
    if (!std::isfinite(x))
        throw std::invalid_argument("argument must be a finite real scalar");
    return static_cast<int>(n);
}

double spherical(int lib, bool first_kind, int n, double x) {
    if (lib == 0)
        return first_kind ? boost::math::sph_bessel(static_cast<unsigned>(n), x)
                          : boost::math::sph_neumann(static_cast<unsigned>(n), x);

    gsl_sf_result result{};
    const int status = first_kind
        ? gsl_sf_bessel_jl_e(n, x, &result)
        : gsl_sf_bessel_yl_e(n, x, &result);
    if (status != GSL_SUCCESS)
        throw std::runtime_error(std::string("GSL spherical Bessel: ") +
                                 gsl_strerror(status));
    return result.val;
}
}

extern "C" {

const char* bessel_riccati_error() {
    return last_error.c_str();
}

int bessel_riccati_supported(int library, const char* name) {
    return supported(library, name) ? 2 : 0;
}

int bessel_riccati_output_count(const char* name,
                                const double* parameters, int length) {
    if (!supported(0, name)) return 0;
    try {
        return 2 * (order(parameters, length) + 1);
    } catch (...) {
        return 0;
    }
}

int bessel_riccati_evaluate(int library, const char* name,
                            const double* parameters, int length,
                            double* output, int capacity) {
    last_error.clear();
    gsl_set_error_handler_off();
    try {
        if (!supported(library, name))
            throw std::invalid_argument("unsupported library or Riccati function");
        const int n = order(parameters, length);
        const int count = n + 1;
        if (!output || capacity < 2 * count)
            throw std::invalid_argument("insufficient output capacity");
        const double x = parameters[1];
        const bool first = std::strcmp(name, "riccati_jn") == 0;

        if (!first && x <= 0)
            throw std::invalid_argument(
                "Riccati y at x<=0 has no finite real Boost/GSL equivalent; "
                "SciPy's +/-1e300 outputs are sentinel values");

        if (x == 0) {
            // ψ_0(x) = sin(x), and ψ_k(0)=0 for all k.
            // ψ'_0(0)=1, ψ'_k(0)=0 for k>0.
            for (int k = 0; k < count; ++k) {
                output[k] = 0.;
                output[count+k] = k == 0 ? 1. : 0.;
            }
            return 0;
        }

        const double argument = std::abs(x);
        std::vector<double> values(count);
        values[0] = first ? std::sin(argument) : -std::cos(argument);
        output[count] = first ? std::cos(argument) : std::sin(argument);
        for (int k = 1; k < count; ++k) {
            const double sph = spherical(library, first, k, argument);
            values[k] = argument * sph;
            // Recurrence follows directly from the defining spherical
            // Bessel derivative. This avoids an additional differentiation.
            output[count+k] = values[k-1] - (k/argument)*values[k];
        }

        if (x < 0) {
            // ψ_k(-a)=(-1)^(k+1) ψ_k(a), hence
            // ψ'_k(-a)=(-1)^k ψ'_k(a).
            values[0] = -values[0];
            for (int k = 1; k < count; ++k) {
                if ((k & 1) == 0) values[k] = -values[k];
                if (k & 1) output[count+k] = -output[count+k];
            }
        }
        for (int k = 0; k < count; ++k) {
            output[k] = values[k];
            if (!std::isfinite(output[k]) || !std::isfinite(output[count+k]))
                throw std::runtime_error("nonfinite spherical Bessel value or derivative");
        }
        return 0;
    } catch (const std::exception& error) {
        last_error = error.what();
        return -1;
    } catch (...) {
        last_error = "unknown native Riccati evaluation failure";
        return -1;
    }
}
}
