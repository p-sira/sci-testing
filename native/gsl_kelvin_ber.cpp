// GSL-DERIVED real Kelvin ber/bei and x-derivatives, NOT native GSL special
// functions. Reuse the existing GSL quadrature in bessel_kelvin_k.cpp through
// its public C ABI; do not duplicate any quadrature or complex-Bessel logic.
//
// The existing `kelvin` output is eight scalar lanes representing SciPy's
// four complex Kelvin values (Be, Ke, Bep, Kep), flattened as
// [ber, bei, ker, kei, berp, beip, kerp, keip].
// This adapter extracts only lanes 0, 1, 4, 5, respectively.

#include <cmath>
#include <cstring>
#include <stdexcept>
#include <string>

extern "C" {
int bessel_kelvin_k_supported(int library, const char* name);
int bessel_kelvin_k_evaluate(int library, const char* name,
                             const double* inputs, int input_count,
                             double* output, int capacity);
const char* bessel_kelvin_k_error();
}

namespace {
thread_local std::string last_error;
constexpr int gsl_library = 1;
constexpr int kelvin_lanes = 8;

int lane(const char* name) {
    if (!name) return -1;
    if (std::strcmp(name, "ber") == 0) return 0;
    if (std::strcmp(name, "bei") == 0) return 1;
    if (std::strcmp(name, "berp") == 0) return 4;
    if (std::strcmp(name, "beip") == 0) return 5;
    return -1;
}
}

extern "C" {

const char* gsl_kelvin_ber_error() { return last_error.c_str(); }

int gsl_kelvin_ber_supported(int library, const char* name) {
    return library == gsl_library && lane(name) >= 0 ? 1 : 0;
}

int gsl_kelvin_ber_output_count(const char* name) {
    return lane(name) >= 0 ? 1 : 0;
}

int gsl_kelvin_ber_evaluate(int library, const char* name,
                            const double* input, int count,
                            double* output, int capacity) {
    last_error.clear();
    try {
        if (!gsl_kelvin_ber_supported(library, name))
            throw std::invalid_argument("unsupported GSL-derived Kelvin function/library");
        if (!input || !output || count != 1 || capacity < 1)
            throw std::invalid_argument("invalid GSL Kelvin input/output capacity");

        // Keep precisely the domain validated by the shared GSL integrator.
        if (!std::isfinite(input[0]) || input[0] < 0.1 || input[0] > 8.0)
            throw std::domain_error("GSL Kelvin quadrature requires finite 0.1 <= x <= 8");
        if (bessel_kelvin_k_supported(gsl_library, "kelvin") != kelvin_lanes)
            throw std::runtime_error("GSL Kelvin tuple backend unavailable");

        double values[kelvin_lanes]{};
        if (bessel_kelvin_k_evaluate(gsl_library, "kelvin", input, 1,
                                    values, kelvin_lanes) != 0) {
            const char* reason = bessel_kelvin_k_error();
            throw std::runtime_error(reason && *reason ? reason :
                                     "GSL Kelvin tuple quadrature failed");
        }

        const double value = values[lane(name)];
        if (!std::isfinite(value))
            throw std::runtime_error("nonfinite GSL-derived Kelvin component");
        output[0] = value;
        return 0;
    } catch (const std::exception& exception) {
        last_error = exception.what();
    } catch (...) {
        last_error = "unknown GSL-derived Kelvin failure";
    }
    return -1;
}
}
