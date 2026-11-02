// Public Boost.Math Legendre functions adapted to SciPy's full array outputs.
//
// roots_legendre(n): [n ascending nodes, n Gauss-Legendre weights]. Boost's
// legendre_p_zeros returns only the nonnegative half; use parity to reconstruct
// the full quadrature rule, then the standard derivative weight identity.
// lqn(n,x): [Q_0..Q_n, Q'_0..Q'_n]. Boost exposes legendre_q only on (-1,1).
// Outside that interval, Q_0=sign(x)*atanh(1/abs(x)), Q_1=x*Q_0-1,
// and the public Legendre three-term and derivative recurrences give SciPy's
// real-valued exterior branch. For larger exterior orders, backward recurrence
// is used to prevent the growing first-kind solution from contaminating Q_n.
#include <boost/math/special_functions/legendre.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
thread_local std::string last_error;
constexpr double miller_threshold = 1.0001;

bool named(const char* name, const char* expected) {
    return name && std::strcmp(name, expected) == 0;
}

bool supported(int library, const char* name) {
    return library == 0 && (named(name, "roots_legendre") || named(name, "lqn"));
}

int validated_order(const char* name, const double* input, int length) {
    const bool roots = named(name, "roots_legendre");
    if (!input || length != (roots ? 1 : 2))
        throw std::invalid_argument("incorrect scalar input signature");
    const double order = input[0];
    if (!std::isfinite(order) || std::trunc(order) != order ||
        order < (roots ? 1 : 0) || order > 128)
        throw std::invalid_argument("Legendre order is outside integer [1,128] roots or [0,128] Q range");
    if (!roots) {
        const double x = input[1];
        if (!std::isfinite(x) || std::fabs(x) == 1.)
            throw std::invalid_argument("Legendre Q requires finite x != +/-1");
        if (std::fabs(x) > 1. && (std::fabs(x) > 1.0e12 ||
            (std::fabs(x) < miller_threshold && order > 8.)))
            throw std::invalid_argument("Legendre Q exterior argument outside validated stable recurrence range");
    }
    return static_cast<int>(order);
}

void roots_legendre(int n, double* output) {
    const auto half = boost::math::legendre_p_zeros<long double>(n);
    if (half.size() != static_cast<std::size_t>((n + 1) / 2))
        throw std::runtime_error("unexpected Boost nonnegative Legendre zero count");

    for (int k = 0; k < n; ++k) {
        // For odd n half[0]=0 and is kept only once at the center.
        const bool negative = k < n / 2;
        const std::size_t index = negative ? half.size() - 1 - k
                                           : static_cast<std::size_t>(k - n / 2);
        long double node = half[index];
        if (negative) node = -node;
        if (!std::isfinite(node) || !(node > -1. && node < 1.))
            throw std::runtime_error("Boost returned invalid Legendre node");
        const long double derivative = boost::math::legendre_p_prime(n, node);
        const long double weight = 2.0L / ((1.0L - node * node) * derivative * derivative);
        if (!(weight > 0.) || !std::isfinite(weight))
            throw std::runtime_error("invalid Boost-derived Gauss-Legendre weight");
        output[k] = static_cast<double>(node);
        output[n + k] = static_cast<double>(weight);
    }
}

void legendre_q_table(int n, double xd, double* output) {
    const long double x = static_cast<long double>(xd);
    const long double ax = std::fabs(x);
    long double previous;
    if (ax < 1.) {
        // The interior seed is the public Boost.Math Legendre Q function.
        previous = boost::math::legendre_q(0u, x);
    } else {
        // Stable for x large and x just outside the boundary; atanh(1/abs(x))
        // would round its argument to 1 near x=1.
        previous = std::copysign(0.5L * std::log1p(2.L / (ax - 1.L)), x);
    }
    std::vector<long double> q(n + 1);
    q[0] = previous;
    if (ax > 1. && ax >= static_cast<long double>(miller_threshold) && n >= 1) {
        // Q_n is recessive for |x|>1. Forward recurrence eventually loses
        // all of its significant digits to the growing P_n solution. Miller's
        // backward recurrence is stable; normalize at the known Q_0 value.
        // Extra iterations suppress admixture by exp(-2*extra*acosh(|x|)).
        const int extra = std::max(32, static_cast<int>(
            std::ceil(48.L / (2.L * std::acosh(ax)))));
        const int top = n + extra;
        std::vector<long double> table(top + 2, 0.L);
        table[top] = 1.L;
        for (int k = top; k > 0; --k) {
            table[k - 1] = ((2.L * k + 1.L) * x * table[k]
                             - (k + 1.L) * table[k + 1]) / k;
        }
        if (!std::isfinite(table[0]) || table[0] == 0.)
            throw std::runtime_error("Legendre Q backward recurrence overflow");
        const long double factor = previous / table[0];
        for (int k = 1; k <= n; ++k)
            q[k] = factor * table[k];
    } else if (n > 0) {
        q[1] = x * q[0] - 1.L;
        for (int k = 1; k < n; ++k)
            q[k + 1] = boost::math::legendre_next(
                static_cast<unsigned>(k), x, q[k], q[k - 1]);
    }
    for (int k = 0; k <= n; ++k) {
        const long double derivative = (k == 0) ? 1.L / (1.L - x*x) :
            k * (q[k - 1] - x * q[k]) / (1.L - x*x);
        if (!std::isfinite(q[k]) || !std::isfinite(derivative))
            throw std::runtime_error("nonfinite Legendre Q recurrence");
        output[k] = static_cast<double>(q[k]);
        output[n + 1 + k] = static_cast<double>(derivative);
    }
}
} // namespace

extern "C" {
const char* boost_special_arrays_error() {
    return last_error.c_str();
}

int boost_special_arrays_supported(int library, const char* name) {
    return supported(library, name) ? 2 : 0;
}

int boost_special_arrays_output_count(const char* name, const double* input, int length) {
    if (!supported(0, name)) return 0;
    try {
        const int n = validated_order(name, input, length);
        return 2 * (n + (named(name, "lqn") ? 1 : 0));
    } catch (...) {
        return 0;
    }
}

int boost_special_arrays_evaluate(int library, const char* name,
                                  const double* input, int length,
                                  double* output, int capacity) {
    last_error.clear();
    try {
        if (!supported(library, name))
            throw std::invalid_argument("unsupported Boost.Math Legendre array function");
        const int n = validated_order(name, input, length);
        const int required = 2 * (n + (named(name, "lqn") ? 1 : 0));
        if (!output || capacity < required)
            throw std::invalid_argument("insufficient Legendre output capacity");
        if (named(name, "roots_legendre"))
            roots_legendre(n, output);
        else
            legendre_q_table(n, input[1], output);
        for (int k = 0; k < required; ++k)
            if (!std::isfinite(output[k]))
                throw std::runtime_error("nonfinite Boost-derived Legendre array value");
        return 0;
    } catch (const std::exception& error) {
        last_error = error.what();
    } catch (...) {
        last_error = "unknown Boost Legendre array evaluation error";
    }
    return -1;
}
}
