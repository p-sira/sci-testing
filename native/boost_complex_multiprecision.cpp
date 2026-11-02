// Independent complex arithmetic via the public Boost.Multiprecision
// cpp_complex_50 numerical backend, NOT std::complex, NumPy, Sci or SciPy.
// Header and functions are public in Boost 1.83 and later.
//
// 0: success with 1 or 2 scalar lanes; 1: unsupported numerical domain;
// -1: invalid input / exception. On failure the caller's output is untouched.
#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
using boost::multiprecision::cpp_complex_50;
using boost::multiprecision::cpp_bin_float_50;
thread_local std::string last_error;

bool named(const char* a, const char* b) {
    return a != nullptr && std::strcmp(a, b) == 0;
}

int output_lanes(const char* n) {
    if (named(n,"abs") || named(n,"arg") || named(n,"norm_sqr")
        || named(n,"real") || named(n,"imag")) return 1;
    if (named(n,"acos") || named(n,"acosh") || named(n,"add")
        || named(n,"asin") || named(n,"asinh") || named(n,"atan")
        || named(n,"atanh") || named(n,"cbrt") || named(n,"conjugate")
        || named(n,"cos") || named(n,"cosh") || named(n,"divide")
        || named(n,"exp") || named(n,"from_polar") || named(n,"ln")
        || named(n,"multiply") || named(n,"negate") || named(n,"powc")
        || named(n,"powf") || named(n,"powi") || named(n,"powu")
        || named(n,"reciprocal") || named(n,"sin") || named(n,"sinh")
        || named(n,"sqrt") || named(n,"subtract") || named(n,"tan")
        || named(n,"tanh") || named(n,"to_polar")) return 2;
    return 0;
}

int input_lanes(const char* n) {
    if (named(n,"add") || named(n,"divide") || named(n,"multiply")
        || named(n,"powc") || named(n,"subtract")) return 4;
    if (named(n,"powf") || named(n,"powi") || named(n,"powu")) return 3;
    return output_lanes(n) ? 2 : 0;
}

// Refuse signed-zero branches: cpp_bin_float's signed-zero preservation is
// not generally identical to Sci/IEEE binary64 on logarithmic branch cuts.
// This is a conservative, *input-only* gate independent of reference output.
bool sensitive_axis(const char* n) {
    return named(n,"acos") || named(n,"acosh") || named(n,"asin")
        || named(n,"asinh") || named(n,"atan") || named(n,"atanh")
        || named(n,"cbrt") || named(n,"ln") || named(n,"powc")
        || named(n,"powf") || named(n,"sqrt");
}

int compute(const char* n, const double* in, double* out) {
    const int arity = input_lanes(n);
    for (int i=0;i<arity;++i) {
        if (!std::isfinite(in[i])) {
            last_error = "complex input is nonfinite";
            return 1;
        }
    }
    // Both large input exponents and output under/overflow require validation
    // against binary64. We only advertise a numerically bounded real domain.
    for (int i=0;i<arity;++i) {
        if (std::abs(in[i]) > 1.0e12) {
            last_error = "input outside validated Boost complex magnitude range";
            return 1;
        }
    }
    if (sensitive_axis(n) && in[1] == 0.) {
        last_error = "principal-branch axis and signed-zero bank not validated";
        return 1;
    }
    // Sci's IEEE signed-zero branch direction on the imaginary axis does
    // not always coincide with Boost's high-precision complex adaptor.
    if ((named(n,"asin") || named(n,"asinh") ||
         named(n,"atan") || named(n,"atanh") ||
         named(n,"acos") || named(n,"acosh")) && in[0] == 0.) {
        last_error = "imaginary-axis branch and signed zero not validated";
        return 1;
    }
    if (named(n,"atan") && in[1] != 0. &&
        std::abs(in[1]) < std::numeric_limits<double>::min()) {
        last_error = "atan subnormal imaginary signed-zero underflow";
        return 1;
    }
    if (named(n,"powc") && in[3] == 0. && in[2] == 0.
        && in[0] < 0. && in[1] == 0.) {
        last_error = "complex power on negative real cut";
        return 1;
    }
    cpp_complex_50 z(in[0],in[1]), w;
    if (named(n,"abs")) {
        double value = abs(z).convert_to<double>();
        if (!std::isfinite(value)) return 1;
        out[0] = value; return 0;
    }
    if (named(n,"arg")) {
        // Boost arg uses the complex-adaptor backend atan2 implementation.
        auto value=arg(z).convert_to<double>();
        if (!std::isfinite(value)) return 1;
        out[0]=value; return 0;
    }
    if (named(n,"real") || named(n,"imag") || named(n,"norm_sqr")) {
        const cpp_bin_float_50 value =
            named(n,"real") ? z.real() :
            named(n,"imag") ? z.imag() :
            (z.real()*z.real() + z.imag()*z.imag()).convert_to<cpp_bin_float_50>();
        auto d=value.convert_to<double>();
        if (!std::isfinite(d)) return 1;
        out[0]=d; return 0;
    }
    if (named(n,"to_polar")) {
        const auto radius=abs(z).convert_to<double>();
        const auto phase=arg(z).convert_to<double>();
        if (!std::isfinite(radius) || !std::isfinite(phase)) return 1;
        out[0]=radius;out[1]=phase;return 0;
    }
    if (named(n,"from_polar")) {
        if (in[0] < 0.) return 1;
        // Radius × Boost-computed polar phase (not std::polar/std::complex).
        const cpp_bin_float_50 radius(in[0]), angle(in[1]);
        w=cpp_complex_50(radius*cos(angle),radius*sin(angle));
    } else if (named(n,"acos")) w=acos(z);
    else if (named(n,"acosh")) w=acosh(z);
    else if (named(n,"asin")) w=asin(z);
    else if (named(n,"asinh")) w=asinh(z);
    else if (named(n,"atan")) w=atan(z);
    else if (named(n,"atanh")) w=atanh(z);
    else if (named(n,"cbrt")) w=pow(z,cpp_complex_50(cpp_bin_float_50(1)/3));
    else if (named(n,"conjugate")) w=conj(z);
    else if (named(n,"cos")) w=cos(z);
    else if (named(n,"cosh")) w=cosh(z);
    else if (named(n,"exp")) w=exp(z);
    else if (named(n,"ln")) w=log(z);
    else if (named(n,"negate")) w=-z;
    else if (named(n,"reciprocal")) {
        if (in[0]==0. && in[1]==0.) return 1;
        w=cpp_complex_50(1)/z;
    } else if (named(n,"sin")) w=sin(z);
    else if (named(n,"sinh")) w=sinh(z);
    else if (named(n,"sqrt")) w=sqrt(z);
    else if (named(n,"tan")) w=tan(z);
    else if (named(n,"tanh")) w=tanh(z);
    else if (named(n,"add")) w=z+cpp_complex_50(in[2],in[3]);
    else if (named(n,"subtract")) w=z-cpp_complex_50(in[2],in[3]);
    else if (named(n,"multiply")) w=z*cpp_complex_50(in[2],in[3]);
    else if (named(n,"divide")) {
        if (in[2]==0. && in[3]==0.) return 1;
        w=z/cpp_complex_50(in[2],in[3]);
    } else if (named(n,"powc")) {
        if (in[0]==0. && in[1]==0.) return 1;
        w=pow(z,cpp_complex_50(in[2],in[3]));
    } else if (named(n,"powf")) {
        if (in[0]==0. && in[1]==0.) return 1;
        w=pow(z,cpp_complex_50(in[2]));
    } else if (named(n,"powi") || named(n,"powu")) {
        if (in[2] != std::trunc(in[2]) || std::abs(in[2]) > 1000.) return 1;
        if (named(n,"powu") && in[2] < 0.) return 1;
        if (in[0]==0. && in[1]==0. && in[2]<0.) return 1;
        w=pow(z,cpp_complex_50(in[2]));
    } else return -1;
    const auto re=w.real().convert_to<double>();
    const auto im=w.imag().convert_to<double>();
    if (!std::isfinite(re) || !std::isfinite(im)) {
        last_error="complex output does not fit finite binary64";
        return 1;
    }
    out[0]=re;out[1]=im;
    return 0;
}
} // namespace

extern "C" {
int boost_complex_mp_supported(const char* name) { return output_lanes(name); }
const char* boost_complex_mp_error() { return last_error.c_str(); }
int boost_complex_mp_evaluate(const char* name, const double* inputs,
                              int count, double* outputs, int capacity) {
    last_error.clear();
    try {
        const int expected=input_lanes(name), out_count=output_lanes(name);
        if (!expected || !inputs || count!=expected || !outputs || capacity<out_count) {
            last_error="unsupported operation or invalid C ABI shapes";
            return -1;
        }
        // Stage the results and commit only after success.
        double stage[2]={0.,0.};
        int status=compute(name,inputs,stage);
        if (status==0) for (int i=0;i<out_count;++i) outputs[i]=stage[i];
        return status;
    } catch (const std::exception& e) {
        last_error=e.what();
    } catch (...) {
        last_error="unknown Boost.Multiprecision exception";
    }
    return -1;
}
}
