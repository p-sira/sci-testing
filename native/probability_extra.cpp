// Additional Sci special probability functions exposed through public Boost.Math
// and GSL APIs. Separate symbols allow linking alongside native/comparisons.cpp.
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/distributions/non_central_f.hpp>
#include <boost/math/distributions/non_central_t.hpp>
#include <boost/math/distributions/kolmogorov_smirnov.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
thread_local std::string error;

struct Function {
    const char* name;
    int arity;
    int libraries;  // bit 0: Boost.Math, bit 1: GSL
};

constexpr Function functions[] = {
    {"nbdtr", 3, 3}, {"nbdtrc", 3, 3}, {"nbdtri", 3, 3},
    {"chndtr", 3, 1}, {"chndtrix", 3, 1},
    {"chndtridf", 3, 1}, {"chndtrinc", 3, 1},
    {"ncfdtr", 4, 1}, {"ncfdtri", 4, 1}, {"ncfdtrinc", 4, 1},
    {"ncfdtridfd", 4, 1}, {"ncfdtridfn", 4, 1},
    {"nctdtr", 3, 1}, {"nctdtrit", 3, 1}, {"nctdtrinc", 3, 1},
    {"nctdtridf", 3, 1},
    {"btdtria", 3, 1}, {"btdtrib", 3, 1}, {"bdtrin", 3, 1},
    {"nbdtrik", 3, 1}, {"nbdtrin", 3, 1},
    {"pdtrik", 2, 1}, {"gdtrib", 3, 1}, {"chdtriv", 2, 1},
    {"kolmogorov", 1, 1}, {"kolmogi", 1, 1},
};

const Function* find_function(const char* name) {
    if (!name) return nullptr;
    for (const Function& fn : functions)
        if (std::strcmp(name, fn.name) == 0) return &fn;
    return nullptr;
}

void require(bool condition, const char* reason) {
    if (!condition) throw std::domain_error(reason);
}

void positive(double x) { require(x > 0, "parameter must be positive"); }
void nonnegative(double x) { require(x >= 0, "parameter must be nonnegative"); }
void probability(double p) { require(p >= 0 && p <= 1, "invalid probability"); }
void interior(double p) { require(p > 0 && p < 1, "probability must be strictly between 0 and 1"); }

// Independently reject an inverse unless it satisfies the requested CDF,
// including in the smaller probability tail. Only public Boost.Math cdf()
// is needed, including on older Boost releases without df inverse methods.
template <typename F>
double verified_degree_of_freedom(double candidate, double target, F cdf) {
    positive(candidate);
    const double actual = cdf(candidate);
    const double tolerance = 1e-9 * std::max(std::min(target, 1-target), 1e-10);
    require(std::isfinite(actual) && std::abs(actual - target) <= tolerance,
            "degrees-of-freedom inverse failed CDF round trip");
    return candidate;
}

// The noncentral F CDF is *not* always monotone in either degree of freedom:
// for some inputs two distinct positive parameters have the same CDF. Scan a
// broad positive log-domain to locate every bracket before solving; do not
// silently report an arbitrary branch or claim SciPy compatibility where
// there are multiple roots. A fixed finite search range is preferable to an
// unbounded scan because the public Boost CDF can saturate in the far tails.
// The saved SciPy inverse fixtures have one well-conditioned positive root.
template <typename F>
double invert_degree_of_freedom(F cdf, double target) {
    interior(target);
    constexpr double lower_bound = 1e-10;
    constexpr double upper_bound = 1e10;
    constexpr double factor = 1.75;
    double low = lower_bound;
    double low_residual = cdf(low) - target;
    require(std::isfinite(low_residual), "invalid degrees-of-freedom CDF at search bound");
    int root_count = 0;
    double root_low = 0, root_high = 0;
    for (int i = 0; i < 128 && low < upper_bound; ++i) {
        const double high = std::min(upper_bound, low * factor);
        const double high_residual = cdf(high) - target;
        require(std::isfinite(high_residual), "nonfinite degrees-of-freedom CDF in search");
        // A bracket crossing is a root. A zero on a grid point is counted
        // exactly once (not once for each adjacent interval).
        if (low_residual == 0 ||
            (low_residual != 0 && high_residual != 0 &&
             (low_residual < 0) != (high_residual < 0))) {
            ++root_count;
            root_low = low;
            root_high = low_residual == 0 ? low : high;
            require(root_count == 1,
                    "ambiguous degrees-of-freedom inverse: multiple CDF roots");
        }
        low = high;
        low_residual = high_residual;
    }
    if (low_residual == 0) {
        ++root_count;
        root_low = root_high = upper_bound;
    }
    require(root_count == 1,
            root_count ? "ambiguous degrees-of-freedom inverse: multiple CDF roots" :
                         "degrees-of-freedom inverse has no bracketed root");

    // Bisection is reliable in either CDF direction and never assumes
    // monotonicity outside the verified local bracket. A tight numerical
    // round-trip check guards against tangencies/flat-tail false matches.
    double left = root_low, right = root_high;
    double left_residual = cdf(left) - target;
    for (int i = 0; i < 120 && left != right; ++i) {
        const double middle = left + (right - left) * .5;
        if (middle == left || middle == right) break;
        const double residual = cdf(middle) - target;
        require(std::isfinite(residual), "nonfinite CDF during df root solve");
        if (residual == 0) {
            left = right = middle;
            break;
        }
        if ((left_residual < 0) != (residual < 0)) right = middle;
        else {
            left = middle;
            left_residual = residual;
        }
    }
    return verified_degree_of_freedom(left + (right-left)*.5, target, cdf);
}

// Invert a continuous monotone *decreasing* public-library CDF with respect to
// a parameter. This is used only where the parameter has an ordered stochastic
// effect: chi-squared df/ncp, F ncp, and t ncp. A target outside the range is
// an error, rather than a falsely successful best-effort result.
template <typename F>
double invert_decreasing(F cdf, double target, double low, double high) {
    interior(target);
    double at_low = cdf(low), at_high = cdf(high);
    require(std::isfinite(at_low) && std::isfinite(at_high), "invalid CDF bracket");
    require(target <= at_low, "inverse has no solution in parameter domain");
    for (int i = 0; i < 100 && at_high > target; ++i) {
        high = high > 0 ? high * 2 : 1;
        require(std::isfinite(high) && high <= 1e12, "inverse solution outside supported range");
        at_high = cdf(high);
    }
    require(std::isfinite(at_high) && target >= at_high, "could not bracket inverse");
    for (int i = 0; i < 100; ++i) {
        double mid = low + (high - low) * .5;
        if (mid == low || mid == high) break;
        double p = cdf(mid);
        require(std::isfinite(p), "nonfinite CDF during inversion");
        if (p > target) low = mid;
        else high = mid;
    }
    double result = low + (high - low) * .5;
    // Relative accuracy in the smaller tail prevents a numerically plausible
    // but incorrect parameter from being accepted at extreme probabilities.
    require(std::abs(cdf(result) - target) <=
                2e-8 * std::max(std::min(target, 1-target), 1e-10),
            "inverse did not converge");
    return result;
}

double boost_value(const char* name, const double* a) {
    namespace bm = boost::math;
    if (!std::strcmp(name, "nbdtr") || !std::strcmp(name, "nbdtrc")) {
        nonnegative(a[0]); positive(a[1]); probability(a[2]);
        // SciPy's nbdtr family converts both count parameters to integers.
        positive(std::floor(a[1]));
        if (!std::strcmp(name, "nbdtr")) return bm::ibeta(std::floor(a[1]), std::floor(a[0])+1, a[2]);
        return bm::ibetac(std::floor(a[1]), std::floor(a[0])+1, a[2]);
    }
    if (!std::strcmp(name, "nbdtri")) {
        nonnegative(a[0]); positive(a[1]); probability(a[2]);
        positive(std::floor(a[1]));
        return bm::ibeta_inv(std::floor(a[1]), std::floor(a[0])+1, a[2]);
    }
    if (!std::strcmp(name, "nbdtrik")) {
        interior(a[0]); positive(a[1]); interior(a[2]);
        // The discrete quantile is zero below the first attainable CDF.
        if (a[0] <= std::pow(a[2], a[1])) return 0;
        return bm::ibeta_invb(a[1], a[2], a[0])-1;
    }
    if (!std::strcmp(name, "nbdtrin")) {
        nonnegative(a[0]); interior(a[1]); interior(a[2]);
        return bm::ibeta_inva(a[0]+1, a[2], a[1]);
    }
    if (!std::strcmp(name, "bdtrin")) {
        nonnegative(a[0]); interior(a[1]); interior(a[2]);
        // Binomial CDF(k,n,p) = I_(1-p)(n-k,k+1).
        return a[0] + bm::ibetac_invb(a[0]+1, a[2], a[1]);
    }
    if (!std::strcmp(name, "pdtrik")) {
        interior(a[0]); positive(a[1]);
        // Poisson CDF(k,m) = Q(k+1,m).
        return bm::gamma_q_inva(a[1], a[0])-1;
    }
    if (!std::strcmp(name, "btdtria")) {
        interior(a[0]); positive(a[1]); interior(a[2]);
        return bm::ibeta_inva(a[1], a[2], a[0]);
    }
    if (!std::strcmp(name, "btdtrib")) {
        positive(a[0]); interior(a[1]); interior(a[2]);
        return bm::ibeta_invb(a[0], a[2], a[1]);
    }
    if (!std::strcmp(name, "gdtrib")) {
        positive(a[0]); interior(a[1]); positive(a[2]);
        return bm::gamma_p_inva(a[0]*a[2], a[1]);
    }
    if (!std::strcmp(name, "chdtriv")) {
        interior(a[0]); positive(a[1]);
        return 2*bm::gamma_p_inva(a[1]*.5, a[0]);
    }
    if (!std::strcmp(name, "chndtr") || !std::strcmp(name, "chndtrix")) {
        positive(a[1]); nonnegative(a[2]);
        const auto dist = bm::non_central_chi_squared_distribution<double>(a[1],a[2]);
        if (!std::strcmp(name, "chndtr")) {
            nonnegative(a[0]);
            return bm::cdf(dist,a[0]);
        }
        probability(a[0]);
        return bm::quantile(dist,a[0]);
    }
    if (!std::strcmp(name, "chndtridf")) {
        positive(a[0]); interior(a[1]); nonnegative(a[2]);
        return invert_decreasing([&](double df) {
            return bm::cdf(bm::non_central_chi_squared_distribution<double>(df,a[2]),a[0]);
        },a[1],1e-7,2);
    }
    if (!std::strcmp(name, "chndtrinc")) {
        positive(a[0]); positive(a[1]); interior(a[2]);
        return invert_decreasing([&](double nc) {
            return bm::cdf(bm::non_central_chi_squared_distribution<double>(a[1],nc),a[0]);
        },a[2],0,2);
    }
    if (!std::strcmp(name, "ncfdtr") || !std::strcmp(name, "ncfdtri")) {
        positive(a[0]); positive(a[1]); nonnegative(a[2]);
        const auto dist = bm::non_central_f_distribution<double>(a[0],a[1],a[2]);
        if (!std::strcmp(name, "ncfdtr")) {
            nonnegative(a[3]);
            return bm::cdf(dist,a[3]);
        }
        probability(a[3]);
        return bm::quantile(dist,a[3]);
    }
    if (!std::strcmp(name, "ncfdtrinc")) {
        positive(a[0]); positive(a[1]); interior(a[2]); positive(a[3]);
        return invert_decreasing([&](double nc) {
            return bm::cdf(bm::non_central_f_distribution<double>(a[0],a[1],nc),a[3]);
        },a[2],0,2);
    }
    if (!std::strcmp(name, "ncfdtridfd")) {
        // SciPy: ncfdtridfd(dfn, p, nc, f), invert denominator df.
        positive(a[0]); interior(a[1]); nonnegative(a[2]); positive(a[3]);
        return invert_degree_of_freedom([&](double dfd) {
            return bm::cdf(bm::non_central_f_distribution<double>(a[0],dfd,a[2]),a[3]);
        }, a[1]);
    }
    if (!std::strcmp(name, "ncfdtridfn")) {
        // SciPy: ncfdtridfn(p, dfd, nc, f), invert numerator df.
        interior(a[0]); positive(a[1]); nonnegative(a[2]); positive(a[3]);
        return invert_degree_of_freedom([&](double dfn) {
            return bm::cdf(bm::non_central_f_distribution<double>(dfn,a[1],a[2]),a[3]);
        }, a[0]);
    }
    if (!std::strcmp(name, "nctdtr") || !std::strcmp(name, "nctdtrit")) {
        positive(a[0]);
        const auto dist = bm::non_central_t_distribution<double>(a[0],a[1]);
        if (!std::strcmp(name, "nctdtr")) return bm::cdf(dist,a[2]);
        probability(a[2]);
        return bm::quantile(dist,a[2]);
    }
    if (!std::strcmp(name, "nctdtrinc")) {
        positive(a[0]); interior(a[1]);
        // The noncentral-t CDF decreases continuously with noncentrality.
        // Bracket in both directions because nc can be negative.
        double low=-2, high=2;
        auto cdf=[&](double nc) {return bm::cdf(bm::non_central_t_distribution<double>(a[0],nc),a[2]);};
        while (cdf(low) < a[1] && low > -1e7) low*=2;
        require(cdf(low)>=a[1],"inverse solution outside supported range");
        return invert_decreasing(cdf,a[1],low,high);
    }
    if (!std::strcmp(name, "nctdtridf")) {
        // SciPy: nctdtridf(p, nc, t), invert noncentral t df.
        interior(a[0]);
        return invert_degree_of_freedom([&](double df) {
            return bm::cdf(bm::non_central_t_distribution<double>(df,a[1]),a[2]);
        }, a[0]);
    }
    if (!std::strcmp(name, "kolmogorov") || !std::strcmp(name, "kolmogi")) {
        const auto dist = bm::kolmogorov_smirnov_distribution<double>(1);
        if (!std::strcmp(name,"kolmogorov")) {
            if (a[0] <= 0) return 1;
            return bm::cdf(bm::complement(dist,a[0]));
        }
        probability(a[0]);
        if (a[0] == 1) return 0;
        return bm::quantile(bm::complement(dist,a[0]));
    }
    throw std::invalid_argument("unsupported Boost probability function");
}

double gsl_value(const char* name, const double* a) {
    if (!std::strcmp(name, "nbdtr") || !std::strcmp(name, "nbdtrc")) {
        nonnegative(a[0]); positive(a[1]); probability(a[2]);
        require(a[0] < static_cast<double>(std::numeric_limits<unsigned>::max()),"count too large");
        positive(std::floor(a[1]));
        unsigned k=static_cast<unsigned>(std::floor(a[0]));
        if (!std::strcmp(name,"nbdtr")) return gsl_cdf_negative_binomial_P(k,a[2],std::floor(a[1]));
        return gsl_cdf_negative_binomial_Q(k,a[2],std::floor(a[1]));
    }
    if (!std::strcmp(name,"nbdtri")) {
        nonnegative(a[0]); positive(a[1]); probability(a[2]);
        positive(std::floor(a[1]));
        return gsl_cdf_beta_Pinv(a[2],std::floor(a[1]),std::floor(a[0])+1);
    }
    throw std::invalid_argument("unsupported GSL probability function");
}
}  // namespace

extern "C" {
const char* probability_extra_error() { return error.c_str(); }
int probability_extra_supported(int library, const char* name) {
    const Function* fn=find_function(name);
    return fn && library>=0 && library<=1 && (fn->libraries & (1<<library)) ? fn->arity : 0;
}
int probability_extra_evaluate(int library, const char* name, const double* a,
                               int count, double* out) {
    error.clear();
    try {
        if (library == 1) gsl_set_error_handler_off();
        const int arity=probability_extra_supported(library,name);
        if (!arity || arity!=count || !a || !out)
            throw std::invalid_argument("unsupported function or incorrect arity");
        for (int i=0; i<count; ++i)
            if (!std::isfinite(a[i])) throw std::invalid_argument("nonfinite input");
        double result=library==0 ? boost_value(name,a) : gsl_value(name,a);
        if (!std::isfinite(result)) throw std::domain_error("nonfinite probability result");
        *out=result;
        return 0;
    } catch (const std::exception& e) {
        error=e.what();
    } catch (...) {
        error="unknown probability calculation failure";
    }
    return -1;
}
}
