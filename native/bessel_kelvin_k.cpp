// GSL-DERIVED Kelvin comparisons, NOT native GSL Kelvin special functions.
//
// The public GSL integrator gsl_integration_qag evaluates DLMF 10.32.9 and
// 10.32.1 for z = x * exp(i*pi/4), with x restricted to 0.1 <= x <= 8:
//
//   K_0(z) = integral_0^infinity exp(-z*cosh(t)) dt, Re(z)>0.
//   I_0(z) = 1/pi * integral_0^pi exp(z*cos(t)) dt.
//
// Differentiation under these absolutely convergent integrals gives the
// Kelvin derivatives with respect to x. All integrands are evaluated using
// real exp/cos/sin, so this translation unit does NOT claim GSL provides a
// complex Bessel K API. The K upper cutoff has x*cosh(T)/sqrt(2)=70,
// leaving a tail far below double precision throughout the guarded domain.
//
// Library ID 1 = GSL; 0 = Boost.Math is deliberately unsupported.
// Scalar ker, kei, kerp, keip output one real double. Full kelvin output
// [ber,bei,ker,kei,berp,beip,kerp,keip] has eight real doubles, matching
// SciPy's four complex outputs flattened as pairs (real,imag).

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include <array>
#include <cmath>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>

namespace {
thread_local std::string last_error;
constexpr double ROOT_HALF = 0.707106781186547524400844362104849039;
constexpr double INV_PI = 0.318309886183790671537767526745028724;
constexpr double PI = 3.141592653589793238462643383279502884;
constexpr double MIN_X = 0.1;
constexpr double MAX_X = 8.;
constexpr size_t MAX_INTERVALS = 1024;

bool named(const char* name, const char* expected) {
    return name && std::strcmp(name, expected) == 0;
}

int output_count(const char* name) {
    if (named(name, "kelvin")) return 8;
    if (named(name, "ker") || named(name, "kei") ||
        named(name, "kerp") || named(name, "keip")) return 1;
    return 0;
}

struct IntegralContext {
    double x_over_root_two;
    int mode;  // 0,1 I0; 2,3 I0'; 4,5 K0; 6,7 K0'
};

double real_integrand(double t, void* context_pointer) {
    const auto& context = *static_cast<const IntegralContext*>(context_pointer);
    const bool first_kind = context.mode < 4;
    const double c = first_kind ? std::cos(t) : std::cosh(t);
    const double u = context.x_over_root_two * c;
    const double envelope = std::exp(first_kind ? u : -u);
    const double co = std::cos(u);
    const double si = std::sin(u);
    switch (context.mode) {
        case 0: return INV_PI * envelope * co;  // ber
        case 1: return INV_PI * envelope * si;  // bei
        case 2: return INV_PI * ROOT_HALF * c * envelope * (co-si); // berp
        case 3: return INV_PI * ROOT_HALF * c * envelope * (co+si); // beip
        case 4: return envelope * co;            // ker
        case 5: return -envelope * si;           // kei
        case 6: return -ROOT_HALF * c * envelope * (co+si); // kerp
        case 7: return -ROOT_HALF * c * envelope * (co-si); // keip
        default: return NAN;
    }
}

double integrate_one(int mode, IntegralContext& context,
                     gsl_integration_workspace* workspace, double upper_k) {
    context.mode = mode;
    gsl_function fn{&real_integrand, &context};
    double result = 0.;
    double absolute_error = 0.;
    const double upper = mode < 4 ? PI : upper_k;
    const int status = gsl_integration_qag(
        &fn, 0., upper, 1e-13, 1e-12, MAX_INTERVALS, 6,
        workspace, &result, &absolute_error);
    if (status != GSL_SUCCESS)
        throw std::runtime_error(
            std::string("GSL-DERIVED Kelvin quadrature: ") + gsl_strerror(status));
    if (!std::isfinite(result) || !std::isfinite(absolute_error))
        throw std::runtime_error("nonfinite GSL Kelvin quadrature result");
    return result;
}
} // namespace

extern "C" {

const char* bessel_kelvin_k_error() { return last_error.c_str(); }

int bessel_kelvin_k_supported(int library, const char* name) {
    return library == 1 ? output_count(name) : 0;
}

int bessel_kelvin_k_output_count(const char* name) {
    return output_count(name);
}

int bessel_kelvin_k_evaluate(int library, const char* name,
                             const double* inputs, int input_count,
                             double* output, int capacity) {
    last_error.clear();
    // GSL defaults to aborting on errors; explicitly request checked status.
    static const bool error_handler_disabled = [] {
        gsl_set_error_handler_off();
        return true;
    }();
    (void)error_handler_disabled;
    try {
        const int expected = bessel_kelvin_k_supported(library, name);
        if (!expected)
            throw std::invalid_argument("no matching GSL-DERIVED Kelvin function");
        if (!inputs || input_count != 1 || !output || capacity < expected)
            throw std::invalid_argument("invalid scalar input or output capacity");
        const double x = inputs[0];
        if (!std::isfinite(x) || x < MIN_X || x > MAX_X)
            throw std::domain_error(
                "GSL-DERIVED Kelvin comparison validated only on 0.1 <= x <= 8");
        auto workspace = std::unique_ptr<
            gsl_integration_workspace, decltype(&gsl_integration_workspace_free)>(
                gsl_integration_workspace_alloc(MAX_INTERVALS),
                &gsl_integration_workspace_free);
        if (!workspace)
            throw std::runtime_error("failed to allocate GSL integration workspace");
        IntegralContext context{x*ROOT_HALF, 0};
        const double upper_k = std::acosh(70. / context.x_over_root_two);
        if (expected == 1) {
            const int mode = named(name, "ker") ? 4 :
                             named(name, "kei") ? 5 :
                             named(name, "kerp") ? 6 : 7;
            const double value = integrate_one(mode, context, workspace.get(), upper_k);
            output[0] = value;
        } else {
            std::array<double, 8> values{};
            // Use a single GSL backend for all eight real components. This
            // avoids mixing Boost and GSL provenance for kelvin.
            for (int mode = 0; mode < 8; ++mode)
                values[mode] = integrate_one(mode, context, workspace.get(), upper_k);
            constexpr std::array<int,8> scipy_order{{0,1,4,5,2,3,6,7}};
            for (int i = 0; i < 8; ++i)
                output[i] = values[scipy_order[i]];
        }
        return 0;
    } catch (const std::exception& exception) {
        last_error = exception.what();
        return -1;
    } catch (...) {
        last_error = "unknown GSL-DERIVED Kelvin evaluation failure";
        return -1;
    }
}

} // extern "C"
