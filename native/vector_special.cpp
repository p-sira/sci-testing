// Boost.Math public vector log-sum-exp and standard softmax transformations.
//
// Standalone compilation:
//   c++ -std=c++17 -O2 -shared -fPIC native/vector_special.cpp -o vector_special.so
//
// All exported symbols use a unique prefix, so this source may also be
// compiled into the existing native/comparisons.cpp shared object without
// changing the original scalar/complex ABI. Library index 0 means Boost.Math.
#include <boost/math/special_functions/logsumexp.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
thread_local std::string last_error;
constexpr int max_vector_length = 1'048'576;

enum class operation { unknown, logsumexp, softmax, log_softmax };

operation identify(const char* name) {
    if (!name) return operation::unknown;
    if (std::strcmp(name, "logsumexp") == 0) return operation::logsumexp;
    if (std::strcmp(name, "softmax") == 0) return operation::softmax;
    if (std::strcmp(name, "log_softmax") == 0) return operation::log_softmax;
    return operation::unknown;
}

bool valid_size(int n) {
    return n >= 0 && n <= max_vector_length;
}

std::vector<double> evaluate(operation which, const std::vector<double>& x) {
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const int size = static_cast<int>(x.size());

    // Boost's container implementation dereferences max_element on empty
    // input. SciPy's scalar logsumexp([]) is -inf. SciPy's vector softmax
    // variants reject empty vectors rather than returning a scalar.
    if (size == 0) {
        if (which == operation::logsumexp) return {-inf};
        throw std::invalid_argument("empty vector has no softmax reduction");
    }

    const bool any_nan = std::any_of(x.begin(), x.end(),
                                    [](double value) { return std::isnan(value); });
    if (any_nan) return std::vector<double>(
        which == operation::logsumexp ? 1 : size, nan);

    const double maximum = *std::max_element(x.begin(), x.end());

    // SciPy's logsumexp treats +inf as +inf, and an all -inf vector as -inf.
    // Boost's finite-container routine instead constructs inf-inf -> NaN.
    if (which == operation::logsumexp) {
        if (maximum == inf || maximum == -inf) return {maximum};
        // Ordinary arguments are computed by the actual public Boost API.
        return {boost::math::logsumexp(x)};
    }

    // The SciPy public contract propagates an undefined softmax reduction
    // when max is infinite; it is not a uniform distribution over +inf.
    if (maximum == -inf || (maximum == inf && which == operation::softmax))
        return std::vector<double>(size, nan);
    if (maximum == inf) {
        // SciPy log_softmax([+inf, finite, -inf]) == [NaN, -inf, -inf].
        std::vector<double> output;
        output.reserve(x.size());
        for (double value : x) output.push_back(value == inf ? nan : -inf);
        return output;
    }

    std::vector<double> shifted;
    shifted.reserve(x.size());
    for (double value : x) shifted.push_back(value - maximum);

    // Avoid subtracting logsumexp(x) from x directly: at x about 1e308
    // the logsumexp rounds to the maximum and loses the softmax offsets.
    // All shifted inputs are <= 0 and at least one is exactly 0.
    const double normalizer = boost::math::logsumexp(shifted);
    std::vector<double> output;
    output.reserve(x.size());
    for (double value : shifted)
        output.push_back(which == operation::softmax
            ? std::exp(value - normalizer)
            : value - normalizer);
    return output;
}
}

extern "C" {

// Return 1 for a Boost-supported vector name; 0 for other libraries/names.
int vector_special_supported(int library, const char* name) {
    return library == 0 && identify(name) != operation::unknown ? 1 : 0;
}

int vector_special_output_count(const char* name, int input_length) {
    const operation which = identify(name);
    if (which == operation::unknown || !valid_size(input_length)) return 0;
    if (which == operation::logsumexp) return 1;
    return input_length > 0 ? input_length : 0;
}

const char* vector_special_error() {
    return last_error.c_str();
}

// count and capacity represent numbers of doubles (not bytes).
// Returns 0 on success, -1 on any error. On failure, output is not modified.
// IEEE nonfinite results are preserved deliberately where SciPy defines them
// (e.g. logsumexp([inf]) = inf; softmax([-inf,-inf]) = [nan,nan]).
int vector_special_evaluate(int library, const char* name,
                            const double* input, int count,
                            double* output, int capacity) {
    last_error.clear();
    try {
        if (!vector_special_supported(library, name))
            throw std::invalid_argument("unsupported Boost vector-special function");
        const int needed = vector_special_output_count(name, count);
        if (needed == 0 || capacity < needed || !output)
            throw std::invalid_argument("invalid input size or output capacity");
        if (count > 0 && !input)
            throw std::invalid_argument("missing vector input");
        const std::vector<double> values = count > 0
            ? std::vector<double>(input, input + count)
            : std::vector<double>{};
        const auto computed = evaluate(identify(name), values);
        if (static_cast<int>(computed.size()) != needed)
            throw std::runtime_error("unexpected vector output length");
        std::copy(computed.begin(), computed.end(), output);
        return 0;
    } catch (const std::exception& exception) {
        last_error = exception.what();
    } catch (...) {
        last_error = "unknown vector-special exception";
    }
    return -1;
}
}
