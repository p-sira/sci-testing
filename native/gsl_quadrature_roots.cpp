// GSL public fixed Gaussian quadrature adapted to SciPy roots_* signatures.
//
// This uses the public GSL 2.7+ API gsl_integration_fixed_alloc(),
// gsl_integration_fixed_nodes(), and gsl_integration_fixed_weights().
// Each supported SciPy function returns two arrays, ordered as ascending
// abscissae followed by quadrature weights. A scalar degree n therefore
// produces an output buffer of exactly 2*n doubles; optional mu=True
// is not part of the saved Sci signatures and is deliberately not claimed.
//
// GSL's fixed Gaussian weights represent the full orthogonality measure,
// not probabilities. In particular:
//   Chebyshev C / S: x -> 2*x, w -> 2*w (relative to T / U).
//   Hermitenorm: physicists' Gaussian with b=1/2 (e^-x²/2), not b=1.
//   Shifted Jacobi(p,q): [0,1] Jacobi alpha=p-q, beta=q-1.
// All mappings are from public GSL rules; no SciPy/Sci internals are used.

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include <cmath>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
thread_local std::string last_error;
constexpr int MAX_DEGREE = 128;

bool is(const char* name, const char* expected) {
    return name && std::strcmp(name, expected) == 0;
}

int arity(const char* name) {
    if (is(name, "roots_jacobi") || is(name, "roots_sh_jacobi")) return 3;
    if (is(name, "roots_gegenbauer") || is(name, "roots_genlaguerre")) return 2;
    if (is(name, "roots_legendre") || is(name, "roots_laguerre") ||
        is(name, "roots_hermite") || is(name, "roots_hermitenorm") ||
        is(name, "roots_chebyt") || is(name, "roots_chebyu") ||
        is(name, "roots_chebyc") || is(name, "roots_chebys") ||
        is(name, "roots_sh_legendre") || is(name, "roots_sh_chebyt") ||
        is(name, "roots_sh_chebyu")) return 1;
    return 0;
}

struct Rule {
    const gsl_integration_fixed_type* type;
    double a;
    double b;
    double alpha;
    double beta;
    double node_scale = 1.;
    double weight_scale = 1.;
};

int check_inputs(const char* name, const double* args, int nargs) {
    const int expected = arity(name);
    if (!expected || !args || nargs != expected)
        throw std::invalid_argument("unsupported roots_* name or argument arity");
    for (int i = 0; i < expected; ++i)
        if (!std::isfinite(args[i]))
            throw std::invalid_argument("roots_* parameters must be finite real scalars");
    if (args[0] != std::trunc(args[0]) ||
        args[0] < 1. || args[0] > MAX_DEGREE)
        throw std::domain_error("quadrature degree must be integer in [1, 128]");
    if (is(name, "roots_jacobi") &&
        (!(args[1] > -1.) || !(args[2] > -1.)))
        throw std::domain_error("Jacobi alpha and beta must exceed -1");
    if (is(name, "roots_sh_jacobi") &&
        (!(args[2] > 0.) || !(args[1]-args[2] > -1.)))
        throw std::domain_error("shifted Jacobi requires q>0 and p-q>-1");
    if (is(name, "roots_gegenbauer") && !(args[1] > 0.))
        throw std::domain_error("Gegenbauer lambda must exceed 0");
    if (is(name, "roots_genlaguerre") && !(args[1] > -1.))
        throw std::domain_error("generalized Laguerre alpha must exceed -1");
    return static_cast<int>(args[0]);
}

Rule select_rule(const char* name, const double* args) {
    Rule rule{gsl_integration_fixed_legendre, -1., 1., 0., 0.};
    if (is(name, "roots_legendre")) return rule;
    if (is(name, "roots_sh_legendre")) {
        rule.a = 0.; rule.b = 1.;
        return rule;
    }
    if (is(name, "roots_jacobi") || is(name, "roots_sh_jacobi")) {
        rule.type = gsl_integration_fixed_jacobi;
        if (is(name, "roots_jacobi")) {
            rule.alpha = args[1]; rule.beta = args[2];
        } else {
            rule.a = 0.; rule.b = 1.;
            rule.alpha = args[1]-args[2]; // exponent at (1-x)
            rule.beta = args[2]-1.;        // exponent at x
        }
        return rule;
    }
    if (is(name, "roots_gegenbauer")) {
        rule.type = gsl_integration_fixed_gegenbauer;
        rule.alpha = args[1]-.5; // GSL's exponent is lambda - 1/2
        return rule;
    }
    if (is(name, "roots_laguerre") || is(name, "roots_genlaguerre")) {
        rule.type = gsl_integration_fixed_laguerre;
        rule.a = 0.; rule.b = 1.; // x^alpha exp(-x), x>=0
        rule.alpha = is(name, "roots_genlaguerre") ? args[1] : 0.;
        return rule;
    }
    if (is(name, "roots_hermite") || is(name, "roots_hermitenorm")) {
        rule.type = gsl_integration_fixed_hermite;
        rule.a = 0.;
        rule.b = is(name, "roots_hermitenorm") ? .5 : 1.;
        return rule;
    }
    if (is(name, "roots_chebyt") || is(name, "roots_chebyc") ||
        is(name, "roots_sh_chebyt")) {
        rule.type = gsl_integration_fixed_chebyshev;
    } else if (is(name, "roots_chebyu") || is(name, "roots_chebys") ||
               is(name, "roots_sh_chebyu")) {
        rule.type = gsl_integration_fixed_chebyshev2;
    } else {
        throw std::invalid_argument("unrecognized public GSL roots rule");
    }
    if (is(name, "roots_sh_chebyt") || is(name, "roots_sh_chebyu")) {
        rule.a = 0.; rule.b = 1.;
    }
    if (is(name, "roots_chebyc") || is(name, "roots_chebys")) {
        // Orthogonality measures for roots_chebyc and roots_chebys are
        // sqrt(1 - (x/2)^2) and its inverse, over [-2,2].
        rule.node_scale = 2.;
        rule.weight_scale = 2.;
    }
    return rule;
}
} // namespace

extern "C" {

const char* gsl_quadrature_roots_error() { return last_error.c_str(); }

int gsl_quadrature_roots_supported(int library, const char* name) {
    // 0 = Boost.Math, 1 = GSL. This translation unit is GSL ONLY.
    return library == 1 ? arity(name) : 0;
}

int gsl_quadrature_roots_output_count(const char* name,
                                       const double* args, int nargs) {
    try {
        return 2 * check_inputs(name, args, nargs);
    } catch (...) {
        return 0;
    }
}

int gsl_quadrature_roots_evaluate(int library, const char* name,
                                  const double* args, int nargs,
                                  double* out, int capacity) {
    last_error.clear();
    // GSL's default error handler aborts the process; disable it before
    // requesting any of the public fixed-Gaussian quadrature rules.
    static const bool handler_disabled = [] {
        gsl_set_error_handler_off();
        return true;
    }();
    (void)handler_disabled;
    try {
        if (gsl_quadrature_roots_supported(library, name) == 0)
            throw std::invalid_argument("unsupported GSL roots family/library");
        const int degree = check_inputs(name, args, nargs);
        if (!out || capacity < 2*degree)
            throw std::invalid_argument("insufficient roots output capacity");
        const Rule rule = select_rule(name, args);
        const auto deleter = &gsl_integration_fixed_free;
        std::unique_ptr<gsl_integration_fixed_workspace, decltype(deleter)> workspace(
            gsl_integration_fixed_alloc(rule.type, static_cast<size_t>(degree),
                                        rule.a, rule.b, rule.alpha, rule.beta),
            deleter);
        if (!workspace)
            throw std::runtime_error("public GSL fixed quadrature allocation failed");
        if (gsl_integration_fixed_n(workspace.get()) != static_cast<size_t>(degree))
            throw std::runtime_error("unexpected GSL quadrature node count");
        const double* nodes = gsl_integration_fixed_nodes(workspace.get());
        const double* weights = gsl_integration_fixed_weights(workspace.get());
        if (!nodes || !weights)
            throw std::runtime_error("GSL fixed quadrature returned null arrays");
        std::vector<double> checked(static_cast<size_t>(2*degree));
        for (int i=0; i<degree; ++i) {
            const double x = nodes[i] * rule.node_scale;
            const double w = weights[i] * rule.weight_scale;
            if (!std::isfinite(x) || !std::isfinite(w) || !(w > 0.))
                throw std::runtime_error("GSL returned invalid quadrature node/weight");
            if (i && x < checked[i-1])
                throw std::runtime_error("GSL returned unordered quadrature abscissae");
            checked[i] = x;
            checked[degree+i] = w;
        }
        // Never leave partially written output after an error.
        for (int i=0; i<2*degree; ++i) out[i] = checked[i];
        return 0;
    } catch (const std::exception& error) {
        last_error = error.what();
        return -1;
    } catch (...) {
        last_error = "unknown GSL fixed quadrature failure";
        return -1;
    }
}
}
