// Independent Boost.Math-derived comparisons for real Sci special functions.
// Only public Boost.Math APIs are used. These are *compositions* around Boost
// gamma, log1p, expm1, normal quantile, inverse erfc, and complete elliptic K;
// they are NOT claimed to be directly-named Boost.Math implementations.
// Library 0 = Boost.Math; library 1 (GSL) is deliberately unsupported.
// A successful call returns one finite double. Failure never modifies output.
#include <boost/math/constants/constants.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
thread_local std::string last_error;

bool eq(const char* a, const char* b) { return a && std::strcmp(a,b) == 0; }

int arity(int library, const char* name) {
    if (library != 0 || !name) return 0;
    if (eq(name,"multigammaln") || eq(name,"boxcox") ||
        eq(name,"boxcox1p") || eq(name,"factorialk") ||
        eq(name,"agm")) return 2;
    if (eq(name,"nrdtrimn") || eq(name,"nrdtrisd")) return 3;
    if (eq(name,"softplus") || eq(name,"log_expit") ||
        eq(name,"ndtri_exp") || eq(name,"log_ndtr")) return 1;
    return 0;
}

int positive_integer(double x, int max, const char* label) {
    if (x < 1 || x > max || x != std::floor(x))
        throw std::domain_error(std::string(label) + " must be a positive integer within validated range");
    return static_cast<int>(x);
}

double multigammaln(double a, double d) {
    const int dimension = positive_integer(d,16,"multivariate gamma dimension");
    if (a <= (dimension-1.)/2. || a > 10000.)
        throw std::domain_error("multigammaln requires (d-1)/2<a<=10000");
    double sum = dimension*(dimension-1.)/4. *
        std::log(boost::math::constants::pi<double>());
    for (int j=1; j<=dimension; ++j)
        sum += boost::math::lgamma(a + (1.-j)/2.);
    return sum;
}

double boxcox(double x, double lambda, bool shifted) {
    const double z = shifted ? boost::math::log1p(x) : std::log(x);
    if (shifted ? x <= -1. : x <= 0.)
        throw std::domain_error("boxcox real finite comparison requires positive base");
    if (std::abs(lambda) > 128. || x > 1e100)
        throw std::domain_error("boxcox parameters outside validated finite domain");
    return lambda == 0. ? z : boost::math::expm1(lambda*z)/lambda;
}

double softplus(double x) {
    return std::max(x,0.) + boost::math::log1p(std::exp(-std::abs(x)));
}

double gaussian_quantile(double p) {
    if (!(p > 0. && p < 1.))
        throw std::domain_error("normal inverse requires 0<p<1");
    return boost::math::quantile(boost::math::normal_distribution<double>(),p);
}

double inverse_log_normal(double y) {
    // Quantile(Phi) in the log-probability domain. The erfc_inv route avoids
    // exp(y) rounding to 1 near zero, but an underflowing exp(y) cannot be
    // inverted reliably here: reject it rather than advertise full range.
    if (y > -1e-12 || y < -700.)
        throw std::domain_error("ndtri_exp comparison requires -700<=y<=-1e-12");
    const double inverse_erfc_arg = y < -boost::math::constants::ln_two<double>()
        ? 2.*std::exp(y)
        : -2.*boost::math::expm1(y);
    const double sign = y < -boost::math::constants::ln_two<double>() ? -1. : 1.;
    return sign * boost::math::constants::root_two<double>() *
        boost::math::erfc_inv(inverse_erfc_arg);
}

double log_normal_cdf(double x) {
    if (x < -35. || x > 35.)
        throw std::domain_error("log_ndtr comparison validated on [-35,35]");
    const double z = std::abs(x)/boost::math::constants::root_two<double>();
    const double tail = boost::math::erfc(z)/2.;
    return x <= 0. ? std::log(tail) : boost::math::log1p(-tail);
}

double factorialk(double n_input, double k_input) {
    const int k = positive_integer(k_input,20,"factorialk step");
    // The supplied integer Sci fixtures require only n>=0, k>=1. Noninteger
    // n and complex extensions have a different Sci contract, not offered.
    if (n_input < 0. || n_input > 100. || n_input != std::floor(n_input))
        throw std::domain_error("factorialk requires integer 0<=n<=100");
    const int n=static_cast<int>(n_input);
    if (n == 0) return 1.;
    const int quotient=n/k, rem=n%k;
    if (rem == 0)
        return std::pow(double(k),quotient)*boost::math::tgamma(double(quotient+1));
    const double r=double(rem)/k;
    return std::pow(double(k),quotient+1)*
        boost::math::tgamma_ratio(r+quotient+1,r);
}

double agm(double a, double b) {
    // Landen: AGM(a,b) = pi(a+b)/(4*K(|a-b|/(a+b))), where K uses
    // elliptic *modulus* (not the parameter m). Public Boost ellint_1(k).
    if (!(a > 0. && b > 0. && a <= 1e8 && b <= 1e8 && a >= 1e-8 && b >= 1e-8))
        throw std::domain_error("AGM comparison requires positive finite a,b in [1e-8,1e8]");
    // When a/b is very extreme, k rounds too close to one and the K(k)
    // transform loses precision; do not silently claim valid oracle output.
    if (std::max(a,b)/std::min(a,b) > 1e4)
        throw std::domain_error("AGM elliptic-K comparison requires max(a,b)/min(a,b)<=1e4");
    const double k=std::abs(a-b)/(a+b);
    return boost::math::constants::pi<double>()*(a+b)/
        (4.*boost::math::ellint_1(k));
}

double compute(const char* name, const double* x) {
    if (eq(name,"multigammaln")) return multigammaln(x[0],x[1]);
    if (eq(name,"boxcox")) return boxcox(x[0],x[1],false);
    if (eq(name,"boxcox1p")) return boxcox(x[0],x[1],true);
    if (eq(name,"softplus")) return softplus(x[0]);
    if (eq(name,"log_expit")) return -softplus(-x[0]);
    if (eq(name,"nrdtrimn")) {
        if (!(x[1]>0.) || x[1]>1e8)
            throw std::domain_error("nrdtrimn requires positive finite std<=1e8");
        return x[2]-x[1]*gaussian_quantile(x[0]);
    }
    if (eq(name,"nrdtrisd")) {
        const double q=gaussian_quantile(x[1]);
        if (std::abs(q)<1e-12 || (x[2]-x[0])/q<=0.)
            throw std::domain_error("nrdtrisd has no finite positive std for these inputs");
        return (x[2]-x[0])/q;
    }
    if (eq(name,"ndtri_exp")) return inverse_log_normal(x[0]);
    if (eq(name,"log_ndtr")) return log_normal_cdf(x[0]);
    if (eq(name,"factorialk")) return factorialk(x[0],x[1]);
    if (eq(name,"agm")) return agm(x[0],x[1]);
    throw std::invalid_argument("unknown Boost special transform name");
}
} // namespace

extern "C" {
int boost_special_transforms_supported(int library, const char* name) {
    return arity(library,name);
}

const char* boost_special_transforms_error() { return last_error.c_str(); }

int boost_special_transforms_evaluate(int library, const char* name,
                                      const double* args, int count,
                                      double* output, int capacity) {
    last_error.clear();
    try {
        const int expected=arity(library,name);
        if (!expected || !args || !output || count != expected || capacity < 1)
            throw std::invalid_argument("unsupported provider/name or ABI shape");
        for (int i=0;i<count;++i)
            if (!std::isfinite(args[i]))
                throw std::domain_error("Boost-derived transform requires finite real arguments");
        const double value=compute(name,args);
        if (!std::isfinite(value))
            throw std::runtime_error("nonfinite Boost-derived transform result");
        output[0]=value;
        return 0;
    } catch(const std::exception& err) { last_error=err.what(); }
      catch(...) { last_error="unknown native transform exception"; }
    return -1;
}
} // extern "C"
