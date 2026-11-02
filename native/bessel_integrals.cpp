// Bessel values supplied by public Boost.Math/GSL APIs and integrated using
// Boost's adaptive Gauss-Kronrod quadrature. No Sci/Cephes code is linked.
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <string>

namespace {
thread_local std::string last_error;
bool named(const char* name, const char* expected) {
    return name && std::strcmp(name, expected) == 0;
}
int arity(int library, const char* name) {
    if (library != 0 && library != 1) return 0;
    if (named(name, "besselpoly")) return 3;
    if (named(name, "it2i0k0") || named(name, "it2j0y0") ||
        named(name, "iti0k0") || named(name, "itj0y0")) return 1;
    return 0;
}
double bessel(int library, char kind, double nu, double x) {
    if (library == 0) {
        switch (kind) {
        case 'j': return boost::math::cyl_bessel_j(nu, x);
        case 'y': return boost::math::cyl_neumann(nu, x);
        case 'i': return boost::math::cyl_bessel_i(nu, x);
        default: return boost::math::cyl_bessel_k(nu, x);
        }
    }
    gsl_sf_result result{};
    int status = GSL_EUNIMPL;
    if (nu == 0.) {
        switch (kind) {
        case 'j': status = gsl_sf_bessel_J0_e(x, &result); break;
        case 'y': status = gsl_sf_bessel_Y0_e(x, &result); break;
        case 'i': status = gsl_sf_bessel_I0_e(x, &result); break;
        default: status = gsl_sf_bessel_K0_e(x, &result); break;
        }
    } else if (kind == 'j' && nu >= 0.) {
        status = gsl_sf_bessel_Jnu_e(nu, x, &result);
    }
    if (status != GSL_SUCCESS)
        throw std::runtime_error(std::string("GSL Bessel: ") + gsl_strerror(status));
    return result.val;
}
template <typename F> double integrate(const F& f, double a, double b) {
    if (a == b) return 0.;
    double error = 0.;
    double result = boost::math::quadrature::gauss_kronrod<double, 31>::integrate(
        f, a, b, 13, 2e-13, &error);
    if (!std::isfinite(result) || !std::isfinite(error) ||
        error > 3e-10 * std::max(1e-8, std::abs(result)))
        throw std::runtime_error("Bessel quadrature did not converge");
    return result;
}
double one_minus_j0(double t) {
    const double z = t*t/4.;
    return z * (1 - z*(1./4. - z*(1./36. - z/576.)));
}
double i0_minus_one(double t) {
    const double z = t*t/4.;
    return z * (1 + z*(1./4. + z*(1./36. + z/576.)));
}
// Analytic integral from 0 to a of Y0 or K0. Both have a logarithmic
// singularity; isolating its integral before quadrature avoids a misleading
// Gauss-Kronrod error estimate at the endpoint.
double small_log_bessel_integral(char kind, double a) {
    if (a == 0.) return 0.;
    constexpr long double pi = 3.141592653589793238462643383279502884L;
    constexpr long double gamma = 0.577215664901532860606512090082402431L;
    const long double x = a, L = std::log(x/2) + gamma;
    const long double z = x*x/4;
    long double term = 1., harmonic = 0.;
    long double sum = kind == 'y' ? L-1 : 1-L;
    for (int k=1; k<80; ++k) {
        term *= z / (static_cast<long double>(k)*k);
        harmonic += 1./static_cast<long double>(k);
        const long double m = 2*k+1;
        const long double factor = (L-harmonic)/m - 1/(m*m);
        const long double coefficient = kind == 'y'
            ? (k % 2 ? -term : term)*factor : -term*factor;
        sum += coefficient;
        if (std::abs(coefficient) < 1e-25L) break;
    }
    return static_cast<double>(x*sum*(kind=='y' ? 2/pi : 1));
}
// The fixed endpoint 10 avoids numerical integration to infinity for Y0/t.
// Its tail comes from the Mellin transform
// -2^(s-1) Gamma(s/2)^2 cos(pi*s/2)/pi, expanded at s=0, and the
// integrated Taylor series of Y0. The library's Y0 integrates [x,10].
double y0_tail_at_ten() {
    constexpr long double pi = 3.141592653589793238462643383279502884L;
    constexpr long double gamma = 0.577215664901532860606512090082402431L;
    const long double L = std::log(5.L) + gamma, z = 25.L;
    long double term = 1., harmonic = 0., sum = 0.;
    for (int k = 1; k <= 100; ++k) {
        term *= -z / (static_cast<long double>(k)*k);
        harmonic += 1./static_cast<long double>(k);
        long double add = term * ((L-harmonic)/(2*k)-1./(4.L*k*k));
        sum += add;
        if (k > 30 && std::abs(add) < 1e-23L) break;
    }
    return static_cast<double>(pi/6 - L*L/pi - 2*sum/pi);
}
void evaluate(int library, const char* name, const double* args, double* out) {
    const double x = args[0];
    if (!std::isfinite(x)) throw std::domain_error("nonfinite input");
    if (named(name, "besselpoly")) {
        const double lambda = args[1], nu = args[2];
        if (!std::isfinite(lambda) || !std::isfinite(nu) || x < 0. ||
            nu < 0. || lambda + nu <= -1. || x > 20.)
            throw std::domain_error("besselpoly outside validated domain");
        if (x == 0.) { out[0] = nu == 0. ? 1./(lambda+1.) : 0.; return; }
        out[0] = integrate([&](double t) {
            return std::pow(t,lambda)*bessel(library,'j',nu,2*x*t);
        }, 0., 1.);
        return;
    }
    const bool modified = named(name, "iti0k0") || named(name, "it2i0k0");
    const bool divided = named(name, "it2i0k0") || named(name, "it2j0y0");
    if (x == 0. && !divided) {
        out[0] = out[1] = 0.;
        return;
    }
    if (x <= 0. || x > 10.)
        throw std::domain_error("validated Bessel integral range: 0 < x <= 10");
    const char first = modified ? 'i' : 'j';
    const char second = modified ? 'k' : 'y';
    if (divided) {
        // t=x*u^2: removable endpoint; small-t series prevents cancellation.
        out[0] = integrate([&](double u) {
            if (u == 0.) return 0.;
            const double t = x*u*u;
            const double numerator = t < .02 ?
                (modified ? i0_minus_one(t) : one_minus_j0(t)) :
                (modified ? bessel(library,first,0.,t)-1.
                          : 1.-bessel(library,first,0.,t));
            return 2*numerator/u;
        }, 0., 1.);
        if (modified) {
            // Exponential K0 decay: tail past x+52 < 1e-24 absolute.
            out[1] = integrate([&](double t) {
                return bessel(library,second,0.,t)/t;
            },x,x+52.);
        } else {
            static const double tail = y0_tail_at_ten();
            out[1] = tail + integrate([&](double t) {
                return bessel(library,second,0.,t)/t;
            },x,10.);
        }
    } else {
        out[0] = integrate([&](double t) {
            return bessel(library,first,0.,t);
        }, 0., x);
        const double cutoff = std::min(x,0.1);
        out[1] = small_log_bessel_integral(second,cutoff) +
            integrate([&](double t) {
                return bessel(library,second,0.,t);
            },cutoff,x);
    }
}
}
extern "C" {
int bessel_integrals_supported(int library, const char* name) {
    return arity(library,name);
}
int bessel_integrals_output_count(const char* name) {
    if (named(name,"besselpoly")) return 1;
    return named(name,"it2i0k0") || named(name,"it2j0y0") ||
           named(name,"iti0k0") || named(name,"itj0y0") ? 2 : 0;
}
const char* bessel_integrals_error() { return last_error.c_str(); }
int bessel_integrals_evaluate(int library, const char* name,
                              const double* args, int count,
                              double* output, int capacity) {
    last_error.clear();
    gsl_set_error_handler_off();
    try {
        int ar = arity(library,name), n = bessel_integrals_output_count(name);
        if (!ar || !args || !output || count != ar || capacity < n)
            throw std::invalid_argument("invalid Bessel integral ABI");
        double values[2]{};
        evaluate(library,name,args,values);
        for (int i=0; i<n; ++i)
            if (!std::isfinite(values[i]))
                throw std::runtime_error("nonfinite integral result");
        std::copy_n(values,n,output);
        return 0;
    } catch (const std::exception& error) {
        last_error = error.what();
    } catch (...) { last_error = "unexpected Bessel error"; }
    return -1;
}
}
