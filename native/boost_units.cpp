// Independent scientific unit conversions from PUBLIC Boost.Units (Boost >= 1.83).
// No Sci/SciPy, hardcoded temperature affine arithmetic, or bare c/x oracle:
//   - absolute Celsius -> absolute Fahrenheit via Boost.Units conversions
//   - wavelength [SI metres] <-> frequency [SI hertz] via Boost.Units
//     quantities and its codata::c physical_constant (exact SI c=299792458).
// CODATA key lookup has a deliberately tiny supported subset: the EXACT SI
// speed-of-light key, whose value and zero uncertainty are invariant between
// Boost.Units' earlier CODATA and SciPy CODATA 2022. Other dynamic text lookup
// keys are NOT advertised by this module.

#include <boost/units/absolute.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/systems/temperature/celsius.hpp>
#include <boost/units/systems/temperature/fahrenheit.hpp>
#include <boost/units/systems/si/codata/universal_constants.hpp>
#include <boost/units/systems/si/frequency.hpp>
#include <boost/units/systems/si/length.hpp>
#include <cmath>
#include <cstring>
#include <exception>
#include <string>

namespace {
thread_local std::string last_error;

bool eq(const char* actual, const char* expected) {
    return actual != nullptr && std::strcmp(actual, expected) == 0;
}

int input_count(const char* name) {
    if (eq(name, "convert_temperature") || eq(name, "lambda2nu") ||
        eq(name, "nu2lambda")) return 1;
    return 0;
}

bool supported_codataset_key(const char* name, const char* key) {
    return (eq(name, "value") || eq(name, "precision")) &&
           eq(key, "speed of light in vacuum");
}

int evaluate_units(const char* name, double value, double* output) {
    namespace u = boost::units;
    namespace si = boost::units::si;
    if (!std::isfinite(value)) {
        last_error = "nonfinite input outside tested Boost.Units numeric domain";
        return 1;
    }
    double result = 0.;
    if (eq(name, "convert_temperature")) {
        using celsius_abs = u::absolute<u::celsius::temperature>;
        using fahrenheit_abs = u::absolute<u::fahrenheit::temperature>;
        const u::quantity<celsius_abs, double> c(value * celsius_abs{});
        const u::quantity<fahrenheit_abs, double> f(c);
        result = f.value();
    } else if (eq(name, "lambda2nu")) {
        if (value == 0.) {
            last_error = "zero wavelength has nonfinite frequency";
            return 1;
        }
        const u::quantity<si::length, double> wavelength(value * si::meters);
        const u::quantity<si::frequency, double> freq(
            si::constants::codata::c.value() / wavelength);
        result = freq.value();
    } else if (eq(name, "nu2lambda")) {
        if (value == 0.) {
            last_error = "zero frequency has nonfinite wavelength";
            return 1;
        }
        const u::quantity<si::frequency, double> freq(value * si::hertz);
        const u::quantity<si::length, double> wavelength(
            si::constants::codata::c.value() / freq);
        result = wavelength.value();
    } else return -1;

    if (!std::isfinite(result)) {
        last_error = "converted result is not representable as finite binary64";
        return 1;
    }
    output[0] = result;
    return 0;
}
} // namespace

extern "C" {
int boost_units_supported(const char* name) { return input_count(name); }
int boost_units_key_supported(const char* name, const char* key) {
    return supported_codataset_key(name, key) ? 1 : 0;
}
const char* boost_units_error() { return last_error.c_str(); }

int boost_units_evaluate_key(const char* name, const char* key,
                             double* output, int capacity) {
    last_error.clear();
    if (!supported_codataset_key(name, key) || !output || capacity < 1) {
        last_error = "keyed CODATA entry is not independently verified";
        return -1;
    }
    try {
        // Values (and uncertainties) come directly from public Boost.Units
        // CODATA c quantity, not from a hard-coded Python/C++ table.
        const double value = eq(name, "value")
            ? boost::units::si::constants::codata::c.value().value()
            : boost::units::si::constants::codata::c.uncertainty().value();
        if (!std::isfinite(value)) {
            last_error = "nonfinite CODATA key result";
            return 1;
        }
        output[0] = value;
        return 0;
    } catch (const std::exception& e) {
        last_error = e.what();
    } catch (...) {
        last_error = "unknown Boost.Units CODATA key exception";
    }
    return -1;
}

int boost_units_evaluate(const char* name, const double* args, int count,
                         double* output, int capacity) {
    last_error.clear();
    if (!input_count(name) || !args || count != 1 || !output || capacity < 1) {
        last_error = "unsupported Boost.Units API or invalid input/output ABI";
        return -1;
    }
    try {
        double staged[1]{0.};
        const int status = evaluate_units(name, args[0], staged);
        if (status == 0) output[0] = staged[0];
        return status;
    } catch (const std::exception& e) {
        last_error = e.what();
    } catch (...) {
        last_error = "unknown Boost.Units exception";
    }
    return -1;
}
} // extern "C"
