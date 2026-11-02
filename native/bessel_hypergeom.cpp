// Independent real-domain Bessel-family comparisons through PUBLIC Boost.Math
// gamma and generalized hypergeometric pFq functions (Boost >= 1.83).
// This bridge does not call Sci/SciPy/Cephes or private Boost detail APIs.
// Library 0 = Boost.Math; library 1 = GSL (not supported: no public 0F2/0F3).
// Validated domains are intentionally narrower than the full Sci APIs.
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>
#include <cmath>
#include <cstring>
#include <initializer_list>
#include <stdexcept>
#include <string>

namespace {
thread_local std::string last_error;

bool named(const char* actual, const char* expected) {
    return actual && std::strcmp(actual, expected) == 0;
}

int arity(int library, const char* name) {
    if (library != 0) return 0;
    if (named(name, "struve") || named(name, "modstruve")) return 2;
    if (named(name, "ber") || named(name, "bei") ||
        named(name, "berp") || named(name, "beip")) return 1;
    if (named(name, "wright_bessel") || named(name, "log_wright_bessel")) return 3;
    return 0;
}

double pfq(std::initializer_list<double> numerators,
           std::initializer_list<double> denominators, double argument) {
    return boost::math::hypergeometric_pFq(numerators, denominators, argument);
}

double struve_family(const char* name, double nu, double x) {
    // The principal real branch of the 1F2 representation, away from
    // negative-order poles and large-x cancellation / overflow regimes.
    if (!(0. <= nu && nu <= 4. && 0. <= x && x <= 8.))
        throw std::domain_error("Struve comparison validated only for 0<=nu<=4, 0<=x<=8");
    if (x == 0.) return 0.;
    const double b = nu + 1.5;
    const double prefactor = std::pow(x / 2., nu + 1.) /
        (boost::math::tgamma(1.5) * boost::math::tgamma(b));
    const double argument = x*x/4. * (named(name, "struve") ? -1. : 1.);
    return prefactor * pfq({1.}, {1.5, b}, argument);
}

double kelvin_family(const char* name, double x) {
    // 0F3 series of Re/Im J0(x exp(3 pi i/4)), and analytic derivatives.
    // Entire real solutions; |x|<=8 limits large-argument cancellation.
    if (std::abs(x) > 8.)
        throw std::domain_error("Kelvin comparison validated only for |x|<=8");
    const double x2 = x*x;
    const double z = -x2*x2/256.;
    if (named(name, "ber"))
        return pfq({}, {0.5, 0.5, 1.}, z);
    if (named(name, "bei"))
        return x2/4. * pfq({}, {1., 1.5, 1.5}, z);
    if (named(name, "berp"))
        return -(x*x2)/16. * pfq({}, {1.5, 1.5, 2.}, z);
    // d/dx (x^2/4 0F3(;1,3/2,3/2;-x^4/256))
    return x/2. * pfq({}, {1., 1.5, 1.5}, z) -
           (x*x2*x2)/576. * pfq({}, {2., 2.5, 2.5}, z);
}

double wright_half(const char* name, double a, double b, double x) {
    // ONLY the exact a=1/2 branch. Splitting even and odd terms of
    // Phi(a,b,x) = sum_{k>=0} x^k/[k! Gamma(a*k+b)] gives two public 0F2s:
    // 0F2(;1/2,b;x²/4)/Gamma(b)
    //   + x*0F2(;3/2,b+1/2;x²/4)/Gamma(b+1/2).
    // Other a are not representable by this single reduction.
    if (a != 0.5 || !(0.5 <= b && b <= 4.) || !(0. <= x && x <= 8.))
        throw std::domain_error(
            "Wright comparison only covers exact a=1/2, 0.5<=b<=4, 0<=x<=8");
    const double z = x*x/4.;
    const double value = pfq({}, {0.5, b}, z) / boost::math::tgamma(b) +
                         x * pfq({}, {1.5, b+0.5}, z) / boost::math::tgamma(b+0.5);
    if (!(value > 0.) || !std::isfinite(value))
        throw std::runtime_error("nonpositive or nonfinite Wright series");
    return named(name, "log_wright_bessel") ? std::log(value) : value;
}
}  // namespace

extern "C" {

int bessel_hypergeom_supported(int library, const char* name) {
    return arity(library, name);
}

const char* bessel_hypergeom_error() {
    return last_error.c_str();
}

int bessel_hypergeom_evaluate(int library, const char* name,
                             const double* inputs, int count,
                             double* outputs, int capacity) {
    last_error.clear();
    try {
        const int expected = arity(library, name);
        if (!expected || !inputs || count != expected || !outputs || capacity < 1)
            throw std::invalid_argument("unsupported library/API or incorrect ABI shape");
        for (int i=0; i<count; ++i)
            if (!std::isfinite(inputs[i]))
                throw std::invalid_argument("Bessel hypergeometric inputs must be finite real numbers");

        double value;
        if (named(name, "struve") || named(name, "modstruve"))
            value = struve_family(name, inputs[0], inputs[1]);
        else if (expected == 1)
            value = kelvin_family(name, inputs[0]);
        else
            value = wright_half(name, inputs[0], inputs[1], inputs[2]);

        if (!std::isfinite(value))
            throw std::runtime_error("nonfinite Boost hypergeometric result");
        outputs[0] = value;
        return 0;
    } catch (const std::exception& exc) {
        last_error = exc.what();
    } catch (...) {
        last_error = "unknown native hypergeometric exception";
    }
    return -1;
}

}  // extern "C"
