// Public Boost.Math / GSL spherical harmonics and normalized Legendre forms.
// Independent C ABI: link alongside native/comparisons.cpp without changing it.
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>
#include <cmath>
#include <complex>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
thread_local std::string last_error;

struct Specification {
    int n;
    int m;
    double theta;
    double phi;
};

int arity(const char* name) {
    if (!name) return 0;
    if (std::strcmp(name, "sph_harm_y") == 0) return 4;
    if (std::strcmp(name, "sph_legendre_p") == 0) return 3;
    return 0;
}

Specification validate(const double* arguments, int count) {
    if (!arguments || count < 3 || count > 4)
        throw std::invalid_argument("invalid spherical argument count");
    for (int i = 0; i < count; ++i)
        if (!std::isfinite(arguments[i]))
            throw std::invalid_argument("nonfinite spherical argument");

    // A practical degree/order limit bounds native work. The public Boost/GSL
    // routines may support still larger values; such inputs fail explicitly
    // rather than entering an unbounded special-function calculation.
    const double degree_limit = 1024;
    if (arguments[0] < 0 || arguments[0] > degree_limit ||
        std::trunc(arguments[0]) != arguments[0])
        throw std::invalid_argument("spherical degree must be a supported nonnegative integer");
    if (arguments[1] < -degree_limit || arguments[1] > degree_limit ||
        std::trunc(arguments[1]) != arguments[1])
        throw std::invalid_argument("spherical order must be an integer in native range");
    if (arguments[2] < 0 || arguments[2] > std::acos(-1.0))
        throw std::invalid_argument("polar angle theta outside [0, pi]");
    return {static_cast<int>(arguments[0]), static_cast<int>(arguments[1]),
            arguments[2], count == 4 ? arguments[3] : 0.0};
}

double gsl_amplitude(const Specification& args) {
    // GSL sphPlm already includes the same spherical normalization and the
    // Condon-Shortley phase as SciPy (verified for odd/even positive orders).
    // Y_n^(-m) = (-1)^m conjugate(Y_n^m); at phi=0 its normalized Legendre
    // amplitude therefore differs by (-1)^m for negative orders.
    if (std::abs(args.m) > args.n) return 0.0;
    gsl_sf_result result{};
    const int positive_order = std::abs(args.m);
    const double cosine = std::cos(args.theta);
    // GSL accepts x=cos(theta), which cannot resolve a small nonzero polar
    // angle once x rounds to +/-1. With m>0 this can silently replace an
    // actual nonzero Y_n^m by zero. Reject precision-sensitive inputs rather
    // than reporting a false GSL equivalence; exact poles are well-defined.
    if (positive_order != 0 && args.theta != 0.0 &&
        args.theta != std::acos(-1.0) &&
        1.0 - std::abs(cosine) < 1e-10)
        throw std::domain_error("GSL spherical Legendre loses polar-angle precision near a pole");
    const int status = gsl_sf_legendre_sphPlm_e(
        args.n, positive_order, cosine, &result);
    if (status != GSL_SUCCESS)
        throw std::runtime_error(std::string("gsl_sf_legendre_sphPlm_e: ") +
                                 gsl_strerror(status));
    double amplitude = result.val;
    if (args.m < 0 && positive_order % 2 != 0) amplitude = -amplitude;
    return amplitude;
}

double boost_amplitude(const Specification& args) {
    if (std::abs(args.m) > args.n) return 0.0;
    // Public real spherical-harmonic API, with azimuth=0, equals SciPy's
    // sph_legendre_p(n,m,theta) for either sign of m.
    return boost::math::spherical_harmonic_r(
        static_cast<unsigned>(args.n), args.m, args.theta, 0.0);
}

std::complex<double> boost_harmonic(const Specification& args) {
    if (std::abs(args.m) > args.n) return {0.0, 0.0};
    return boost::math::spherical_harmonic(
        static_cast<unsigned>(args.n), args.m, args.theta, args.phi);
}

std::complex<double> gsl_harmonic(const Specification& args) {
    const double amplitude = gsl_amplitude(args);
    // Multiplication by e^(i m phi) sets the signed complex azimuthal phase.
    const double angle = static_cast<double>(args.m) * args.phi;
    return {amplitude * std::cos(angle), amplitude * std::sin(angle)};
}
}  // namespace

extern "C" {
const char* spherical_extra_error() { return last_error.c_str(); }

int spherical_extra_supported(int library, const char* name) {
    return (library == 0 || library == 1) ? arity(name) : 0;
}

int spherical_extra_output_count(const char* name) {
    const int count = arity(name);
    return count == 4 ? 2 : count == 3 ? 1 : 0;
}

int spherical_extra_evaluate(int library, const char* name, const double* args,
                             int count, double* out, int capacity) {
    last_error.clear();
    try {
        if (spherical_extra_supported(library, name) != count || count == 0)
            throw std::invalid_argument("unsupported spherical API or library");
        const int output_count = spherical_extra_output_count(name);
        if (!out || capacity < output_count)
            throw std::invalid_argument("incorrect spherical output capacity");
        const Specification spec = validate(args, count);
        if (library == 1) gsl_set_error_handler_off();

        if (count == 3) {
            out[0] = library == 0 ? boost_amplitude(spec) : gsl_amplitude(spec);
            if (!std::isfinite(out[0]))
                throw std::runtime_error("nonfinite normalized Legendre result");
        } else {
            const std::complex<double> y = library == 0 ?
                boost_harmonic(spec) : gsl_harmonic(spec);
            if (!std::isfinite(y.real()) || !std::isfinite(y.imag()))
                throw std::runtime_error("nonfinite complex spherical harmonic");
            out[0] = y.real();
            out[1] = y.imag();
        }
        return 0;
    } catch (const std::exception& e) {
        last_error = e.what();
    } catch (...) {
        last_error = "unknown spherical function failure";
    }
    return -1;
}
}
