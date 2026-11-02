// Real adjacent modified-Bessel ratios using Boost.Math's public I_v API.
//
// Boost.Math 1.83 exposes cyl_bessel_i(v,x) in bessel.hpp, but does not
// expose a public scaled-I ratio API. For arguments beyond the safe range
// of long double I_v, use an asymptotic Riccati expansion derived from the
// modified-Bessel differential equation, not a private Boost/Sci function.
// The inverse is a bracketed, safeguarded Newton solver for that ratio.
#include <boost/math/special_functions/bessel.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace {
thread_local std::string error_message;
constexpr double max_order = 128.0;
constexpr double max_argument = 1.0e12;
constexpr long double boost_limit = 4000.0L;

bool named(const char* name, const char* expected) {
    return name && std::strcmp(name, expected) == 0;
}
bool supported(const char* name) {
    return named(name, "iv_ratio") || named(name, "iv_ratio_c")
        || named(name, "iv_ratioinv");
}

void check_order(double v) {
    if (!std::isfinite(v) || v < 0.5 || v > max_order)
        throw std::domain_error("order outside validated Boost Bessel-ratio range [0.5,128]");
}

struct Ratio {
    long double value;
    long double complement;
};

Ratio ratio(double vd, double xd) {
    check_order(vd);
    if (!std::isfinite(xd) || xd < 0.0 || xd > max_argument)
        throw std::domain_error("argument outside validated Bessel-ratio range [0,1e12]");
    const long double v = vd, x = xd;
    // Sci preserves -0.0 for the direct ratio's zero limit.
    if (x == 0) return {std::signbit(xd) ? -0.0L : 0.0L, 1.0L};
    if (v == 0.5L) {
        // R_(1/2)=tanh(x); 1-R=2*exp(-2x)/(1+exp(-2x)).
        const long double decayed = std::exp(-2.0L * x);
        return {std::tanh(x), 2.0L * decayed / (1.0L + decayed)};
    }

    // I_v(x) ~ (x/2)^v/Gamma(v+1) for tiny x.  Both adjacent
    // Boost evaluations can underflow (e.g. v=128, x=1e-300),
    // though their ratio is perfectly representable.
    if (x < 1.0e-5L) {
        const long double leading = x / (2.0L * v);
        const long double x2 = x * x;
        const long double numerator = 1.0L + x2 / (4.0L * (v + 1.0L));
        const long double denominator = 1.0L + x2 / (4.0L * v);
        const long double value = leading * (numerator / denominator);
        return {value, 1.0L - value};
    }

    if (x <= boost_limit) {
        // Public Boost functions, in extended precision to avoid overflow
        // of I_v(x) at x ~ 1000 in binary64.
        const long double numerator = boost::math::cyl_bessel_i(v, x);
        const long double denominator = boost::math::cyl_bessel_i(v - 1.0L, x);
        if (numerator > 0 && denominator > 0 &&
            std::isfinite(numerator) && std::isfinite(denominator)) {
            const long double value = numerator / denominator;
            if (value >= 0.0L && value <= 1.0L)
                return {value, 1.0L - value};
        }
        throw std::runtime_error("Boost Bessel I values cannot form a finite reliable ratio");
    }

    // For large x use R'=1-R^2-(2v-1)R/x. If R=1+sum a_n/x^n,
    // a_1=-(2v-1)/2 and
    // 2a_n=(n-2v)a_{n-1}-sum_{k=1}^{n-1} a_k a_{n-k}.
    // Computing 1-R directly from the small series avoids cancellation
    // when x≫v, where an ordinary ratio rounds to one.
    const long double a = 2.0L*v - 1.0L;
    long double coefficients[25]{};
    coefficients[0] = 1.0L;
    coefficients[1] = -a / 2.0L;
    long double power = 1.0L/x;
    long double complement = -coefficients[1]*power;
    long double final_term = std::fabs(coefficients[1]*power);
    for (int n = 2; n <= 24; ++n) {
        long double convolution = 0.0L;
        for (int k = 1; k < n; ++k)
            convolution += coefficients[k]*coefficients[n-k];
        coefficients[n] = ((static_cast<long double>(n) - 1.0L - a)
                           * coefficients[n-1] - convolution) / 2.0L;
        power /= x;
        const long double term = coefficients[n]*power;
        complement -= term;
        final_term = std::fabs(term);
    }
    if (!(complement > 0.0L && complement < 1.0L) ||
        !std::isfinite(complement) ||
        final_term > 5e-13L * complement)
        throw std::runtime_error("Bessel-ratio asymptotic expansion not reliable");
    return {1.0L - complement, complement};
}

double inverse(double v, double r) {
    check_order(v);
    if (!std::isfinite(r) || r < 0.0 || r >= 1.0)
        throw std::domain_error("inverse target must be finite in [0,1)");
    if (r == 0.0) return 0.0;
    if (v == 0.5) return std::atanh(r);

    // Sharp upper bound from the first positive modified-Bessel ratio
    // inequalities. It is also the initial analytic estimate at small r.
    const long double target = r;
    const long double upper = 2.0L*v*target /
        ((1.0L-target)*(1.0L+target));
    if (upper > max_argument)
        throw std::domain_error("inverse target exceeds validated argument range");
    long double lo = 0.0L, hi = upper;
    long double x = upper;
    const long double epsilon = std::numeric_limits<double>::epsilon();
    for (int iter = 0; iter < 96; ++iter) {
        const Ratio value = ratio(v, static_cast<double>(x));
        // Near r=1, use the *complement* rather than subtracting a
        // near-unity, rounded I_v/I_(v-1) quotient.
        const long double difference = target > 0.5L ?
            (1.0L-target) - value.complement : value.value - target;
        if (difference < 0) lo = x;
        else hi = x;
        // The Riccati derivative is a public-identity derivative of I_v/I_(v-1).
        const long double derivative = 1.0L - value.value*value.value -
            (2.0L*v - 1.0L) * value.value / x;
        if (std::fabs(difference) <= 8.0L*epsilon*std::min(target, 1.0L-target)
            || hi-lo <= 8.0L*epsilon*x)
            return static_cast<double>(x);
        long double next = derivative > 0 ?
            x - difference/derivative : (lo + hi)/2.0L;
        if (!std::isfinite(next) || next <= lo || next >= hi)
            next = (lo + hi)/2.0L;
        if (next == x) return static_cast<double>(x);
        x = next;
    }
    throw std::runtime_error("Bessel-ratio inverse failed to converge");
}
}

extern "C" {

int bessel_ratio_supported(int library, const char* name) {
    return library == 0 && supported(name) ? 2 : 0;
}

const char* bessel_ratio_error() {
    return error_message.c_str();
}

// input[0]=v, input[1]=x or r; output[0] scalar.  Library 0 is Boost.Math.
// Error -1 keeps output unmodified; no misleading competitor result is emitted.
int bessel_ratio_evaluate(int library, const char* name,
                          const double* input, int count,
                          double* output, int capacity) {
    error_message.clear();
    try {
        if (bessel_ratio_supported(library, name) != 2)
            throw std::invalid_argument("unsupported Boost Bessel-ratio function");
        if (!input || !output || count != 2 || capacity < 1)
            throw std::invalid_argument("invalid Bessel-ratio input/output capacity");
        const double result = named(name, "iv_ratioinv") ?
            inverse(input[0], input[1]) :
            static_cast<double>(named(name, "iv_ratio_c") ?
                ratio(input[0], input[1]).complement :
                ratio(input[0], input[1]).value);
        if (!std::isfinite(result))
            throw std::runtime_error("nonfinite Bessel-ratio output");
        output[0] = result;
        return 0;
    } catch (const std::exception& error) {
        error_message = error.what();
    } catch (...) {
        error_message = "unknown Bessel-ratio error";
    }
    return -1;
}
}
