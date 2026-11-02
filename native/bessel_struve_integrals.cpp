// Real Struve integral comparisons using Boost.Math's PUBLIC generalized
// hypergeometric_pFq API (present in Boost 1.83).
//
//   itstruve0(x) = integral_0^x H_0(t) dt
//                = x^2/pi * 2F3(1,1;2,3/2,3/2;-x^2/4)
//   itmodstruve0(x) = integral_0^x L_0(t) dt
//                   = x^2/pi * 2F3(1,1;2,3/2,3/2;+x^2/4)
//   it2struve0(x) = integral_x^infinity H_0(t)/t dt
//                 = pi/2 - 2*x/pi * 2F3(1,1/2;3/2,3/2,3/2;-x^2/4)
//
// The first two are even in real x; the third retains the sign of x.
// All formulas follow by termwise integration of the real Struve series.
// The bounded |x|<=8 domain avoids poorly conditioned subtraction in
// it2struve0 and eventual loss of precision in alternating pFq series.
// No GSL provider is advertised: GSL has no corresponding public 2F3 API.

#include <boost/math/special_functions/hypergeometric_pFq.hpp>

#include <cmath>
#include <cstring>
#include <stdexcept>
#include <string>

namespace {
thread_local std::string last_error;
constexpr double pi = 3.141592653589793238462643383279502884;

bool named(const char* name, const char* expected) {
    return name && std::strcmp(name, expected) == 0;
}

int arity(int library, const char* name) {
    if (library != 0) return 0;
    return named(name, "itstruve0") || named(name, "itmodstruve0") ||
           named(name, "it2struve0") ? 1 : 0;
}

double integral(const char* name, double x) {
    if (!std::isfinite(x) || std::abs(x) > 8.0)
        throw std::domain_error("Struve integral comparison requires finite real |x|<=8");
    if (x == 0.0) return named(name, "it2struve0") ? pi/2.0 : 0.0;

    const double z = (x/2.0)*(x/2.0);
    if (named(name, "it2struve0")) {
        const double f = boost::math::hypergeometric_pFq(
            {1.0, 0.5}, {1.5, 1.5, 1.5}, -z);
        return pi/2.0 - (2.0*x/pi)*f;
    }

    const double f = boost::math::hypergeometric_pFq(
        {1.0, 1.0}, {2.0, 1.5, 1.5}, named(name, "itstruve0") ? -z : z);
    return (x*x/pi)*f;
}
}  // namespace

extern "C" {
int bessel_struve_integrals_supported(int library, const char* name) {
    return arity(library, name);
}

const char* bessel_struve_integrals_error() {
    return last_error.c_str();
}

// library 0 = Boost.Math. Exactly one real scalar input and one real output.
// Returns 0 on success; -1 with a human-readable error on failure. Output is
// left untouched on failure (including unsupported library, invalid shape,
// nonfinite input, or an out-of-domain real argument).
int bessel_struve_integrals_evaluate(int library, const char* name,
                                    const double* inputs, int count,
                                    double* outputs, int capacity) {
    last_error.clear();
    try {
        if (!arity(library, name) || !inputs || count != 1 ||
            !outputs || capacity < 1)
            throw std::invalid_argument("unsupported Struve integral library/API or incorrect ABI shape");
        const double value = integral(name, inputs[0]);
        if (!std::isfinite(value))
            throw std::overflow_error("nonfinite Struve integral result");
        outputs[0] = value;
        return 0;
    } catch (const std::exception& exc) {
        last_error = exc.what();
    } catch (...) {
        last_error = "unknown Boost Struve integral evaluation error";
    }
    return -1;
}
}  // extern "C"
