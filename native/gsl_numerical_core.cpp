#include <gsl/gsl_integration.h>

#include <algorithm>
#include <cmath>
#include <string>

namespace {
thread_local std::string last_error;
double exponential(double x, void *) { return std::exp(x); }
}

extern "C" {
const char *gsl_numerical_core_error() { return last_error.c_str(); }

int gsl_numerical_core_supported(int library, const char *name) {
    return library == 1 && name && std::string(name) == "fixed_quad" ? 2 : 0;
}

int gsl_numerical_core_evaluate(int library, const char *name,
                                const double *arguments, int arity,
                                double *output, int capacity) {
    last_error.clear();
    if (!gsl_numerical_core_supported(library, name) || !arguments || arity != 2 ||
        !output || capacity < 1) {
        last_error = "unsupported GSL numerical-core signature or output capacity";
        return -1;
    }
    const double a = arguments[0], b = arguments[1];
    if (!std::isfinite(a) || !std::isfinite(b) || std::max(a,b) > 700.) {
        last_error = "nonfinite or overflowing Gauss-Legendre interval";
        return -1;
    }
    if (a == b) { output[0] = 0.; return 0; }
    // The saved Sci operation is fixed_quad(exp, a, b, n=16), not adaptive
    // quadrature. glfixed is GSL's public fixed Gauss-Legendre algorithm.
    gsl_integration_glfixed_table *table = gsl_integration_glfixed_table_alloc(16);
    if (!table) { last_error = "GSL fixed integration table allocation failed"; return -1; }
    gsl_function function = {&exponential, nullptr};
    const double result = gsl_integration_glfixed(&function, a, b, table);
    gsl_integration_glfixed_table_free(table);
    if (!std::isfinite(result)) {
        last_error = "GSL fixed integration yielded nonfinite output";
        return -1;
    }
    output[0] = result;
    return 0;
}
}
