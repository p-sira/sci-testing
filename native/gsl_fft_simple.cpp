// Public GSL 2.7+ one-dimensional complex DFTs for saved SciPy fft.fft/ifft.
// GSL's inverse normalizes by 1/n, as does scipy.fft.ifft (default norm).
// Interleaved [real, imaginary] input and output; provider GSL = 1.
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
thread_local std::string last_error;
constexpr int max_complex_length = 65536;

bool is_fft(const char* name) {
    return name && std::strcmp(name, "fft.fft") == 0;
}

bool is_ifft(const char* name) {
    return name && std::strcmp(name, "fft.ifft") == 0;
}

int checked_count(const char* name, int n) {
    return (is_fft(name) || is_ifft(name)) && n >= 1 &&
                   n <= max_complex_length
               ? 2 * n
               : 0;
}

void require(bool valid, const char* reason) {
    if (!valid) throw std::invalid_argument(reason);
}
}

extern "C" {
const char* gsl_fft_simple_error() { return last_error.c_str(); }

int gsl_fft_simple_supported(int library, const char* name) {
    return library == 1 && (is_fft(name) || is_ifft(name)) ? 1 : 0;
}

int gsl_fft_simple_output_count(const char* name, int n) {
    return checked_count(name, n);
}

int gsl_fft_simple_evaluate(int library, const char* name,
                            const double* data, int n,
                            double* output, int capacity) {
    last_error.clear();
    try {
        require(gsl_fft_simple_supported(library, name) == 1,
                "unsupported FFT function or GSL provider");
        const int count = checked_count(name, n);
        require(count > 0, "FFT length must be between 1 and 65536");
        require(data && output && capacity >= count,
                "invalid FFT data, output or capacity");
        for (int i = 0; i < count; ++i)
            require(std::isfinite(data[i]), "FFT requires finite complex input");

        std::vector<double> transformed(data, data + count);
        using table_ptr = std::unique_ptr<gsl_fft_complex_wavetable,
                                          decltype(&gsl_fft_complex_wavetable_free)>;
        using workspace_ptr = std::unique_ptr<gsl_fft_complex_workspace,
                                              decltype(&gsl_fft_complex_workspace_free)>;
        table_ptr table(gsl_fft_complex_wavetable_alloc(static_cast<size_t>(n)),
                        &gsl_fft_complex_wavetable_free);
        workspace_ptr workspace(gsl_fft_complex_workspace_alloc(static_cast<size_t>(n)),
                                &gsl_fft_complex_workspace_free);
        require(table && workspace, "GSL FFT workspace allocation failed");
        const int status = is_fft(name)
            ? gsl_fft_complex_forward(transformed.data(), 1, static_cast<size_t>(n),
                                      table.get(), workspace.get())
            : gsl_fft_complex_inverse(transformed.data(), 1, static_cast<size_t>(n),
                                      table.get(), workspace.get());
        require(status == GSL_SUCCESS, "GSL complex FFT failed");
        for (double value : transformed)
            require(std::isfinite(value), "GSL complex FFT returned nonfinite output");
        std::copy(transformed.begin(), transformed.end(), output);
        return 0;
    } catch (const std::exception& error) {
        last_error = error.what();
    } catch (...) {
        last_error = "unknown GSL FFT failure";
    }
    return -1;
}
}
