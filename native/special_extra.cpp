// Additional *public* Boost.Math and GSL special-function comparisons.
// Keep this ABI independent of native/comparisons.cpp while the main bridge
// is being extended in parallel.
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/chebyshev.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/jacobi.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/sinc.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_trig.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
thread_local std::string last_error;

bool named(const char* name, const char* expected) {
    return name && std::strcmp(name, expected) == 0;
}

int arity(const char* name, int library) {
    if (library != 0 && library != 1) return 0;
    if (named(name, "jvp") || named(name, "yvp") ||
        named(name, "ivp") || named(name, "kvp") ||
        named(name, "assoc_legendre_p") || named(name, "spherical_in") ||
        named(name, "spherical_kn")) return 3;
    if (named(name, "eval_sh_jacobi")) return library == 0 ? 4 : 0;
    if (named(name, "eval_chebyc") || named(name, "eval_chebys") ||
        named(name, "perm")) return library == 0 ? 2 : 0;
    if (named(name, "jn_zeros") || named(name, "legendre_p")) return 2;
    if (named(name, "yn_zeros")) return library == 0 ? 2 : 0;
    if (named(name, "sinc")) return 1;
    if (named(name, "exprel")) return library == 1 ? 1 : 0;
    // No direct public Boost.Math/GSL softplus API; do not ascribe a generic
    // std::log1p/std::exp implementation to either comparison library.
    return 0;
}

int exact_integer(double value, int max_value, const char* label) {
    if (!std::isfinite(value) || value != std::floor(value) ||
        value < 0 || value > max_value)
        throw std::invalid_argument(std::string("invalid ") + label);
    return static_cast<int>(value);
}

int signed_integer(double value, int max_absolute, const char* label) {
    if (!std::isfinite(value) || value != std::floor(value) ||
        value < -max_absolute || value > max_absolute)
        throw std::invalid_argument(std::string("invalid ") + label);
    return static_cast<int>(value);
}

void gsl_check(int status) {
    if (status != GSL_SUCCESS) throw std::runtime_error(gsl_strerror(status));
}

// GSL cylindrical functions have a nonnegative-order contract. The Bessel
// reflection identities extend the comparison to the negative orders needed
// by SciPy's order-n derivative finite-sum recurrence.
double gsl_cylinder(char kind, double order, double x) {
    gsl_sf_result result{};
    if (order >= 0) {
        switch (kind) {
            case 'j': gsl_check(gsl_sf_bessel_Jnu_e(order, x, &result)); break;
            case 'y': gsl_check(gsl_sf_bessel_Ynu_e(order, x, &result)); break;
            case 'i': gsl_check(gsl_sf_bessel_Inu_e(order, x, &result)); break;
            case 'k': gsl_check(gsl_sf_bessel_Knu_e(order, x, &result)); break;
        }
        return result.val;
    }
    const double positive = -order;
    if (kind == 'k') return gsl_cylinder('k', positive, x);
    if (kind == 'i') {
        return gsl_cylinder('i', positive, x) +
               (2. / std::acos(-1.)) * std::sin(std::acos(-1.) * positive) *
               gsl_cylinder('k', positive, x);
    }
    const double j = gsl_cylinder('j', positive, x);
    const double y = gsl_cylinder('y', positive, x);
    const double angle = std::acos(-1.) * positive;
    return kind == 'j' ? std::cos(angle)*j - std::sin(angle)*y
                       : std::sin(angle)*j + std::cos(angle)*y;
}

double boost_cylinder(char kind, double order, double x) {
    switch (kind) {
        case 'j': return boost::math::cyl_bessel_j(order, x);
        case 'y': return boost::math::cyl_neumann(order, x);
        case 'i': return boost::math::cyl_bessel_i(order, x);
        default: return boost::math::cyl_bessel_k(order, x);
    }
}

double cylinder_derivative(int library, char kind, double order, double x, int derivative) {
    double result = 0.;
    double choose = 1.;
    for (int k = 0; k <= derivative; ++k) {
        const double shifted = order - derivative + 2*k;
        const double term = library == 0 ? boost_cylinder(kind, shifted, x)
                                         : gsl_cylinder(kind, shifted, x);
        const double sign = ((kind == 'j' || kind == 'y') && (k & 1)) ? -1. : 1.;
        result += choose * sign * term;
        if (k < derivative) choose *= double(derivative-k) / double(k+1);
    }
    result = std::ldexp(result, -derivative);
    if (kind == 'k' && (derivative & 1)) result = -result;
    return result;
}

double spherical(int library, char kind, int n, double x) {
    if (library == 0) {
        const double factor = std::sqrt(std::acos(-1.)/(2.*x));
        return factor * (kind == 'i' ? boost::math::cyl_bessel_i(n+.5,x)
                                     : boost::math::cyl_bessel_k(n+.5,x));
    }
    gsl_sf_result result{};
    const int status = kind == 'i' ? gsl_sf_bessel_il_scaled_e(n, x, &result)
                                   : gsl_sf_bessel_kl_scaled_e(n, x, &result);
    gsl_check(status);
    // GSL gives exp(-|x|) i_n and exp(x) k_n, respectively.
    return result.val * std::exp(kind == 'i' ? std::abs(x) : -x);
}

double spherical_value(int library, char kind, int n, double x, bool derivative) {
    if (kind == 'i' && x == 0.) {
        if (!derivative) return n == 0 ? 1. : 0.;
        return n == 1 ? 1./3. : 0.;
    }
    if (kind == 'i' && x < 0.) {
        const double value = spherical_value(library, kind, n, -x, derivative);
        return ((n + int(derivative)) & 1) ? -value : value;
    }
    if (x <= 0.) throw std::invalid_argument("spherical K requires positive x");
    const double value = spherical(library, kind, n, x);
    if (!derivative) return value;
    if (n == 0) return (kind == 'i' ? 1. : -1.) * spherical(library, kind, 1, x);
    const double prev = spherical(library, kind, n-1, x);
    return (kind == 'i' ? prev : -prev) - (n+1.)*value/x;
}

double scalar(int library, const char* name, const double* a) {
    if (named(name, "jvp") || named(name, "yvp") ||
        named(name, "ivp") || named(name, "kvp")) {
        const int derivative = exact_integer(a[2], 12, "derivative order");
        if (a[1] <= 0 || !std::isfinite(a[0]))
            throw std::invalid_argument("Bessel derivative requires x > 0 and finite order");
        return cylinder_derivative(library, name[0], a[0], a[1], derivative);
    }
    if (named(name, "spherical_in") || named(name, "spherical_kn")) {
        const int n = exact_integer(a[0], 100, "spherical order");
        if (a[2] != 0 && a[2] != 1) throw std::invalid_argument("invalid derivative flag");
        return spherical_value(library, name[10], n, a[1], a[2] != 0);
    }
    if (named(name, "legendre_p")) {
        const int n = exact_integer(a[0], 10000, "Legendre degree");
        if (library == 0) return boost::math::legendre_p(n, a[1]);
        gsl_sf_result r{};
        gsl_check(gsl_sf_legendre_Pl_e(n, a[1], &r));
        return r.val;
    }
    if (named(name, "assoc_legendre_p")) {
        const int n = exact_integer(a[0], 10000, "Legendre degree");
        const int m = signed_integer(a[1], 10000, "Legendre order");
        if (std::abs(m) > n) return 0.;
        if (library == 0) return boost::math::legendre_p(n, m, a[2]);
        gsl_sf_result r{};
        const int positive_order = std::abs(m);
        gsl_check(gsl_sf_legendre_Plm_e(n, positive_order, a[2], &r));
        if (m >= 0) return r.val;
        // P_n^{-m} = (-1)^m (n-m)!/(n+m)! P_n^m.
        gsl_sf_result lower{}, upper{};
        gsl_check(gsl_sf_lnfact_e(static_cast<unsigned>(n-positive_order), &lower));
        gsl_check(gsl_sf_lnfact_e(static_cast<unsigned>(n+positive_order), &upper));
        return (positive_order & 1 ? -1. : 1.) *
               std::exp(lower.val-upper.val) * r.val;
    }
    if (named(name, "eval_chebyc") || named(name, "eval_chebys")) {
        const auto degree = static_cast<unsigned>(exact_integer(a[0], 10000, "polynomial degree"));
        return named(name, "eval_chebyc") ? 2.*boost::math::chebyshev_t(degree, a[1]/2.)
                                           : boost::math::chebyshev_u(degree, a[1]/2.);
    }
    if (named(name, "eval_sh_jacobi")) {
        const auto degree = static_cast<unsigned>(exact_integer(a[0], 1000, "polynomial degree"));
        // SciPy G_n^(p,q)(x) = P_n^(p-q,q-1)(2x-1) / binom(2n+p-1,n).
        const double normalization = boost::math::factorial<double>(degree) /
           boost::math::falling_factorial(2.*degree+a[1]-1., degree);
        return boost::math::jacobi(degree, a[1]-a[2], a[2]-1., 2.*a[3]-1.)*normalization;
    }
    if (named(name, "perm")) {
        const unsigned k = static_cast<unsigned>(exact_integer(a[1], 1000, "permutation length"));
        if (!std::isfinite(a[0])) throw std::invalid_argument("permutation requires finite N");
        if (a[0] < 0 || a[0] < k) return 0.;
        if (k == 0) return 1.;
        return boost::math::falling_factorial(a[0], k);
    }
    if (named(name, "sinc")) {
        if (library == 0) return boost::math::sinc_pi(a[0]*std::acos(-1.));
        gsl_sf_result r{};
        gsl_check(gsl_sf_sinc_e(a[0], &r));
        return r.val;
    }
    if (named(name, "exprel")) {
        gsl_sf_result r{};
        gsl_check(gsl_sf_exprel_e(a[0], &r));
        return r.val;
    }
    throw std::invalid_argument("unsupported function");
}
} // namespace

extern "C" {
const char* special_extra_error() { return last_error.c_str(); }
int special_extra_supported(int library, const char* name) { return arity(name, library); }
int special_extra_output_count(const char* name, const double* args, int count) {
    if (!args || !name) return 0;
    if (named(name, "jn_zeros") || named(name, "yn_zeros")) {
        if (count != 2 || !std::isfinite(args[1]) ||
            args[1] != std::floor(args[1]) || args[1] < 1 || args[1] > 1000)
            return 0;
        return static_cast<int>(args[1]);
    }
    return 1;
}
int special_extra_evaluate(int library, const char* name, const double* args, int count,
                           double* output, int capacity) {
    last_error.clear();
    gsl_set_error_handler_off();
    try {
        if (!name || !args || !output || !arity(name, library) ||
            count != arity(name, library)) throw std::invalid_argument("unsupported function or arity");
        const int expected = special_extra_output_count(name, args, count);
        if (!expected || capacity < expected) throw std::invalid_argument("invalid output capacity");
        if (named(name, "jn_zeros") || named(name, "yn_zeros")) {
            const int order = exact_integer(args[0], 10000, "Bessel order");
            for (int k = 1; k <= expected; ++k) {
                if (library == 0) {
                    output[k-1] = named(name, "jn_zeros") ? boost::math::cyl_bessel_j_zero(double(order), k)
                                                           : boost::math::cyl_neumann_zero(double(order), k);
                } else {
                    gsl_sf_result result{};
                    gsl_check(gsl_sf_bessel_zero_Jnu_e(order, static_cast<unsigned>(k), &result));
                    // GSL's public Jnu-zero routine occasionally leaves a
                    // few e-9 relative error for larger order/rank. Refine
                    // its root with GSL's own Jnu function and recurrence.
                    double root = result.val;
                    for (int iteration = 0; iteration < 3; ++iteration) {
                        const double j = gsl_cylinder('j', order, root);
                        const double derivative = .5 * (
                            gsl_cylinder('j', order - 1., root) -
                            gsl_cylinder('j', order + 1., root));
                        if (derivative == 0.) throw std::runtime_error("zero Bessel derivative");
                        const double step = j / derivative;
                        root -= step;
                        if (std::abs(step) < 1e-14 * std::max(1., std::abs(root))) break;
                    }
                    output[k-1] = root;
                }
                if (!std::isfinite(output[k-1])) throw std::runtime_error("nonfinite output");
            }
        } else {
            output[0] = scalar(library, name, args);
            if (!std::isfinite(output[0])) throw std::runtime_error("nonfinite output");
        }
        return 0;
    } catch (const std::exception& e) { last_error = e.what(); return -1; }
      catch (...) { last_error = "unknown native exception"; return -1; }
}
}
