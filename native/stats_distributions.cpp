// Independent public Boost.Math and GSL distribution comparisons.
// The ABI is separate from the special-function bridge; no C++ exception crosses it.
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/cauchy.hpp>
#include <boost/math/distributions/laplace.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/uniform.hpp>
#include <boost/math/distributions/weibull.hpp>
#include <boost/math/distributions/pareto.hpp>
#include <boost/math/distributions/inverse_gamma.hpp>
#include <boost/math/distributions/extreme_value.hpp>
#include <boost/math/distributions/triangular.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/geometric.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/policies/policy.hpp>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
thread_local std::string distribution_error;
using RoundDiscrete = boost::math::policies::policy<
    boost::math::policies::discrete_quantile<
        boost::math::policies::integer_round_up>>;
bool eq(const char* a, const char* b) { return a && std::strcmp(a, b) == 0; }
// Operations are 0 pdf, 1 logpdf, 2 cdf, 3 sf, 4 ppf.
int operation(const char* op) {
    if (eq(op, "pdf")) return 0;
    if (eq(op, "logpdf")) return 1;
    if (eq(op, "cdf")) return 2;
    if (eq(op, "sf")) return 3;
    if (eq(op, "ppf")) return 4;
    return -1;
}
int arity(const char* n) {
    if (eq(n, "Normal") || eq(n, "Beta") || eq(n, "Gamma") ||
        eq(n, "F") || eq(n, "Cauchy") || eq(n, "Laplace") ||
        eq(n, "LogNormal") || eq(n, "Uniform") || eq(n, "Weibull") ||
        eq(n, "Pareto") || eq(n, "InverseGamma") || eq(n, "Gumbel") ||
        eq(n, "Binomial") || eq(n, "NegativeBinomial") ||
        eq(n, "Erlang") || eq(n, "Levy") || eq(n, "DiscreteUniform"))
        return 2;
    if (eq(n, "Triangular") || eq(n, "Hypergeometric")) return 3;
    if (eq(n, "ChiSquared") || eq(n, "Chi") || eq(n, "StudentT") ||
        eq(n, "Exponential") || eq(n, "Bernoulli") || eq(n, "Poisson") ||
        eq(n, "Geometric")) return 1;
    return 0;
}
bool supported(int lib, const char* n, int op) {
    if (lib != 0 && lib != 1) return false;
    if (op < 0 || op > 4 || arity(n) == 0) return false;
    if (eq(n, "Levy")) return false; // No matching Boost.Math or GSL density distribution.
    if (lib == 0)
        return !eq(n, "Chi") && !eq(n, "DiscreteUniform");
    if (eq(n, "InverseGamma") || eq(n, "Triangular") ||
        eq(n, "DiscreteUniform"))
        return false;
    if (op == 4 && (eq(n, "Bernoulli") || eq(n, "Binomial") ||
        eq(n, "Poisson") || eq(n, "Geometric") ||
        eq(n, "NegativeBinomial") || eq(n, "Hypergeometric")))
        return false; // GSL has no public discrete inverse CDF.
    return true;
}
void validate(const char* n, const double* p, int count) {
    if (count != arity(n)) throw std::invalid_argument("incorrect parameter arity");
    for (int i = 0; i < count; ++i)
        if (!std::isfinite(p[i])) throw std::invalid_argument("nonfinite distribution parameter");
    auto positive = [](double x) { return x > 0; };
    auto probability = [](double x) { return x >= 0 && x <= 1; };
    auto integer = [](double x) {
        return x >= 0 && x <= 65536 && std::floor(x) == x;
    };
    if ((eq(n,"Normal") || eq(n,"Cauchy") || eq(n,"Laplace") ||
         eq(n,"Gumbel") || eq(n,"Levy") || eq(n,"LogNormal")) && !positive(p[1]))
        throw std::invalid_argument("invalid scale");
    if ((eq(n,"Gamma") || eq(n,"Erlang") || eq(n,"Weibull") ||
         eq(n,"Pareto") || eq(n,"InverseGamma")) &&
        (!positive(p[0]) || !positive(p[1])))
        throw std::invalid_argument("invalid shape or scale");
    if ((eq(n,"Beta") || eq(n,"F")) && (!positive(p[0]) || !positive(p[1])))
        throw std::invalid_argument("invalid shape or degrees of freedom");
    if ((eq(n,"ChiSquared") || eq(n,"Chi") || eq(n,"StudentT") ||
         eq(n,"Exponential")) && !positive(p[0]))
        throw std::invalid_argument("invalid degrees of freedom or scale");
    if ((eq(n,"Uniform") || eq(n,"DiscreteUniform")) && !(p[1] > p[0]))
        throw std::invalid_argument("invalid bounds");
    if (eq(n,"Triangular") && !(p[0] < p[1] && p[0] <= p[2] && p[2] <= p[1]))
        throw std::invalid_argument("invalid triangular mode");
    if (eq(n,"Erlang") && std::floor(p[0]) != p[0])
        throw std::invalid_argument("Erlang shape must be integral");
    if (eq(n,"Bernoulli") && !probability(p[0]))
        throw std::invalid_argument("invalid Bernoulli probability");
    if (eq(n,"Binomial") && (!integer(p[0]) || !probability(p[1])))
        throw std::invalid_argument("invalid binomial parameters");
    if (eq(n,"Poisson") && p[0] < 0)
        throw std::invalid_argument("invalid Poisson rate");
    if (eq(n,"Geometric") && !(p[0] > 0 && p[0] <= 1))
        throw std::invalid_argument("invalid geometric probability");
    if (eq(n,"NegativeBinomial") && !(p[0] > 0 && p[1] > 0 && p[1] <= 1))
        throw std::invalid_argument("invalid negative binomial parameters");
    if (eq(n,"Hypergeometric") && (!integer(p[0]) || !integer(p[1]) ||
         !integer(p[2]) || p[1] > p[0] || p[2] > p[0]))
        throw std::invalid_argument("invalid hypergeometric counts");
}
// Several public Boost/GSL functions reject observations outside mathematical
// support, where Sci's PDFs/CDFs have well-defined zero/one limits.
// Degenerate p=1 discrete laws and Poisson(rate=0) are likewise valid in Sci.
bool support_limit(const char* n, const double* p, int op, double x, double* result) {
    const bool cdf = op == 2, sf = op == 3, logpdf = op == 1;
    auto below = [&] { *result = cdf ? 0 : sf ? 1 : logpdf ? -INFINITY : 0; return true; };
    auto above = [&] { *result = cdf ? 1 : sf ? 0 : logpdf ? -INFINITY : 0; return true; };
    if (op == 4) {
        if (x == 0) {
            if (eq(n,"Bernoulli") || eq(n,"Binomial") ||
                eq(n,"Poisson") || eq(n,"NegativeBinomial")) {
                *result = -1; return true;
            }
            if (eq(n,"Geometric") || eq(n,"InverseGamma")) {
                *result = 0; return true;
            }
            if (eq(n,"Hypergeometric")) {
                *result = std::max(0., p[1]+p[2]-p[0])-1; return true;
            }
        }
        if (eq(n,"Poisson") && p[0] == 0) { *result = 0; return true; }
        if (eq(n,"NegativeBinomial") && p[1] == 1) { *result = 0; return true; }
        if (eq(n,"Geometric") && p[0] == 1) { *result = x == 0 ? 0 : 1; return true; }
        return false;
    }
    if (x == 0 && op <= 1 &&
        (eq(n,"Gamma") || eq(n,"Erlang") || eq(n,"ChiSquared"))) {
        const double shape = eq(n,"ChiSquared") ? p[0]/2 : p[0];
        const double scale = eq(n,"ChiSquared") ? 2 : p[1];
        if (shape < 1) *result = INFINITY;
        else if (shape > 1) *result = op == 1 ? -INFINITY : 0;
        else *result = op == 1 ? -std::log(scale) : 1./scale;
        return true;
    }
    if (eq(n,"Beta") && (x < 0 || x > 1)) return x < 0 ? below() : above();
    if ((eq(n,"Gamma") || eq(n,"Erlang") || eq(n,"ChiSquared") ||
         eq(n,"Chi") || eq(n,"F") || eq(n,"Exponential") ||
         eq(n,"Weibull") || eq(n,"InverseGamma") || eq(n,"LogNormal")) && x < 0)
        return below();
    if (eq(n,"Pareto") && x < p[0]) return below();
    if (eq(n,"Triangular") && (x < p[0] || x > p[1]))
        return x < p[0] ? below() : above();
    if ((eq(n,"Bernoulli") || eq(n,"Binomial")) &&
        x > (eq(n,"Bernoulli") ? 1 : p[0])) return above();
    if (eq(n,"Hypergeometric")) {
        const double low = std::max(0., p[1] + p[2] - p[0]);
        const double high = std::min(p[1], p[2]);
        if (x < low) return below();
        if (x > high) return above();
    }
    if (eq(n,"Geometric") && x < 1) return below();
    if (eq(n,"Poisson") && p[0] == 0) {
        if (x < 0) return below();
        if (cdf || sf) { *result = cdf ? 1 : 0; return true; }
        *result = x == 0 ? (logpdf ? 0 : 1) : (logpdf ? -INFINITY : 0);
        return true;
    }
    if (eq(n,"NegativeBinomial") && p[1] == 1) {
        if (x < 0) return below();
        if (cdf || sf) { *result = cdf ? 1 : 0; return true; }
        *result = x == 0 ? (logpdf ? 0 : 1) : (logpdf ? -INFINITY : 0);
        return true;
    }
    if (eq(n,"Geometric") && p[0] == 1) {
        if (cdf || sf) { *result = cdf ? 1 : 0; return true; }
        *result = x == 1 ? (logpdf ? 0 : 1) : (logpdf ? -INFINITY : 0);
        return true;
    }
    return false;
}
template <class Distribution>
double boost_value(const Distribution& distribution, int op, double x) {
    switch (op) {
        case 0: return boost::math::pdf(distribution, x);
        case 1: return std::log(boost::math::pdf(distribution, x));
        case 2: return boost::math::cdf(distribution, x);
        case 3: return boost::math::cdf(boost::math::complement(distribution, x));
        case 4: return boost::math::quantile(distribution, x);
        default: throw std::invalid_argument("invalid operation");
    }
}
double boost_eval(const char* n, const double* p, int op, double x) {
    // Boost geometric counts failures; Sci and SciPy count trials, starting at 1.
    if (eq(n,"Geometric")) {
        auto d = boost::math::geometric_distribution<double, RoundDiscrete>(p[0]);
        // Boost's geometric quantile is a continuous interpolation even with
        // the integer-rounding policy; round up to the first discrete trial.
        if (op == 4) return std::ceil(boost_value(d, op, x) + 1);
        if (op == 0 || op == 1) {
            if (x < 1 || std::floor(x) != x) return op == 0 ? 0 : -INFINITY;
        }
        if (op >= 2) x = std::floor(x);
        if (x < 1) return op == 2 ? 0 : op == 3 ? 1 : op == 1 ? -INFINITY : 0;
        return boost_value(d, op, x - 1);
    }
    // The discrete PMF vanishes off the integers, and CDF/SF use floor(x).
    const bool discrete = eq(n,"Bernoulli") || eq(n,"Binomial") ||
                          eq(n,"Poisson") || eq(n,"NegativeBinomial") ||
                          eq(n,"Hypergeometric");
    if (discrete && op != 4) {
        if (op < 2 && (std::floor(x) != x || x < 0))
            return op == 0 ? 0 : -INFINITY;
        if (op >= 2) x = std::floor(x);
        if (x < 0) return op == 2 ? 0 : op == 3 ? 1 : 0;
        // Hypergeometric has stricter support bounds than SciPy.
        if (eq(n,"Hypergeometric")) {
            const double low = std::max(0., p[1] + p[2] - p[0]);
            const double high = std::min(p[1], p[2]);
            if (x < low) return op == 2 ? 0 : op == 3 ? 1 : op == 1 ? -INFINITY : 0;
            if (x > high) return op == 2 ? 1 : op == 3 ? 0 : op == 1 ? -INFINITY : 0;
        }
    }
#define BOOST_DIST(NAME, TYPE, ...) \
    if (eq(n, NAME)) return boost_value(boost::math::TYPE##_distribution<double>(__VA_ARGS__), op, x)
    BOOST_DIST("Normal", normal, p[0], p[1]);
    BOOST_DIST("Beta", beta, p[0], p[1]);
    BOOST_DIST("Gamma", gamma, p[0], p[1]);
    BOOST_DIST("Erlang", gamma, p[0], p[1]);
    BOOST_DIST("ChiSquared", chi_squared, p[0]);
    BOOST_DIST("StudentT", students_t, p[0]);
    BOOST_DIST("F", fisher_f, p[0], p[1]);
    BOOST_DIST("Exponential", exponential, 1/p[0]);
    BOOST_DIST("Cauchy", cauchy, p[0], p[1]);
    BOOST_DIST("Laplace", laplace, p[0], p[1]);
    BOOST_DIST("LogNormal", lognormal, p[0], p[1]);
    BOOST_DIST("Uniform", uniform, p[0], p[1]);
    BOOST_DIST("Weibull", weibull, p[0], p[1]);
    BOOST_DIST("Pareto", pareto, p[0], p[1]);
    BOOST_DIST("InverseGamma", inverse_gamma, p[0], p[1]);
    BOOST_DIST("Gumbel", extreme_value, p[0], p[1]);
    BOOST_DIST("Triangular", triangular, p[0], p[2], p[1]);
#define BOOST_DISCRETE(NAME, TYPE, ...) \
    if (eq(n, NAME)) return boost_value(boost::math::TYPE##_distribution<double, RoundDiscrete>(__VA_ARGS__), op, x)
    BOOST_DISCRETE("Bernoulli", binomial, 1., p[0]);
    BOOST_DISCRETE("Binomial", binomial, p[0], p[1]);
    BOOST_DISCRETE("Poisson", poisson, p[0]);
    BOOST_DISCRETE("NegativeBinomial", negative_binomial, p[0], p[1]);
    BOOST_DISCRETE("Hypergeometric", hypergeometric,
                   static_cast<std::uint64_t>(p[1]), static_cast<std::uint64_t>(p[2]),
                   static_cast<std::uint64_t>(p[0]));
#undef BOOST_DISCRETE
#undef BOOST_DIST
    throw std::invalid_argument("unsupported Boost distribution");
}
double gsl_eval(const char* n, const double* p, int op, double x) {
    double y = x;
    if (eq(n,"Normal") || eq(n,"Cauchy") || eq(n,"Laplace")) y = x - p[0];
    // For all continuous distributions use public GSL PDF and P/Q/inverse P APIs.
#define GSL_CONT(NAME, PDF, CDF, SF, INV) \
    if (eq(n, NAME)) { \
        if (op == 4) return (INV); \
        if (op == 2) return (CDF); \
        if (op == 3) return (SF); \
        const double density = (PDF); \
        return op == 1 ? std::log(density) : density; \
    }
    GSL_CONT("Normal", gsl_ran_gaussian_pdf(y,p[1]),
             gsl_cdf_gaussian_P(y,p[1]), gsl_cdf_gaussian_Q(y,p[1]),
             p[0]+gsl_cdf_gaussian_Pinv(x,p[1]));
    GSL_CONT("Beta", gsl_ran_beta_pdf(x,p[0],p[1]),
             gsl_cdf_beta_P(x,p[0],p[1]), gsl_cdf_beta_Q(x,p[0],p[1]),
             gsl_cdf_beta_Pinv(x,p[0],p[1]));
    GSL_CONT("Gamma", gsl_ran_gamma_pdf(x,p[0],p[1]),
             gsl_cdf_gamma_P(x,p[0],p[1]), gsl_cdf_gamma_Q(x,p[0],p[1]),
             gsl_cdf_gamma_Pinv(x,p[0],p[1]));
    GSL_CONT("Erlang", gsl_ran_erlang_pdf(x,p[1],p[0]),
             gsl_cdf_gamma_P(x,p[0],p[1]), gsl_cdf_gamma_Q(x,p[0],p[1]),
             gsl_cdf_gamma_Pinv(x,p[0],p[1]));
    GSL_CONT("ChiSquared", gsl_ran_chisq_pdf(x,p[0]),
             gsl_cdf_chisq_P(x,p[0]), gsl_cdf_chisq_Q(x,p[0]),
             gsl_cdf_chisq_Pinv(x,p[0]));
    // Chi is the positive square root of a chi-squared random variable.
    GSL_CONT("Chi", x < 0 ? 0 : 2*x*gsl_ran_chisq_pdf(x*x,p[0]),
             x < 0 ? 0 : gsl_cdf_chisq_P(x*x,p[0]),
             x < 0 ? 1 : gsl_cdf_chisq_Q(x*x,p[0]),
             std::sqrt(gsl_cdf_chisq_Pinv(x,p[0])));
    GSL_CONT("StudentT", gsl_ran_tdist_pdf(x,p[0]),
             gsl_cdf_tdist_P(x,p[0]), gsl_cdf_tdist_Q(x,p[0]),
             gsl_cdf_tdist_Pinv(x,p[0]));
    GSL_CONT("F", gsl_ran_fdist_pdf(x,p[0],p[1]),
             gsl_cdf_fdist_P(x,p[0],p[1]), gsl_cdf_fdist_Q(x,p[0],p[1]),
             gsl_cdf_fdist_Pinv(x,p[0],p[1]));
    GSL_CONT("Exponential", gsl_ran_exponential_pdf(x,p[0]),
             gsl_cdf_exponential_P(x,p[0]), gsl_cdf_exponential_Q(x,p[0]),
             gsl_cdf_exponential_Pinv(x,p[0]));
    GSL_CONT("Cauchy", gsl_ran_cauchy_pdf(y,p[1]),
             gsl_cdf_cauchy_P(y,p[1]), gsl_cdf_cauchy_Q(y,p[1]),
             p[0]+gsl_cdf_cauchy_Pinv(x,p[1]));
    GSL_CONT("Laplace", gsl_ran_laplace_pdf(y,p[1]),
             gsl_cdf_laplace_P(y,p[1]), gsl_cdf_laplace_Q(y,p[1]),
             p[0]+gsl_cdf_laplace_Pinv(x,p[1]));
    GSL_CONT("LogNormal", gsl_ran_lognormal_pdf(x,p[0],p[1]),
             gsl_cdf_lognormal_P(x,p[0],p[1]), gsl_cdf_lognormal_Q(x,p[0],p[1]),
             gsl_cdf_lognormal_Pinv(x,p[0],p[1]));
    GSL_CONT("Uniform", gsl_ran_flat_pdf(x,p[0],p[1]),
             gsl_cdf_flat_P(x,p[0],p[1]), gsl_cdf_flat_Q(x,p[0],p[1]),
             gsl_cdf_flat_Pinv(x,p[0],p[1]));
    GSL_CONT("Weibull", gsl_ran_weibull_pdf(x,p[1],p[0]),
             gsl_cdf_weibull_P(x,p[1],p[0]), gsl_cdf_weibull_Q(x,p[1],p[0]),
             gsl_cdf_weibull_Pinv(x,p[1],p[0]));
    // GSL takes Pareto (shape, minimum); Sci exposes (minimum, shape).
    GSL_CONT("Pareto", gsl_ran_pareto_pdf(x,p[1],p[0]),
             gsl_cdf_pareto_P(x,p[1],p[0]), gsl_cdf_pareto_Q(x,p[1],p[0]),
             gsl_cdf_pareto_Pinv(x,p[1],p[0]));
    // gsl_gumbel1 uses (a=1/scale, b=exp(location/scale)), unlike
    // SciPy and Boost's (location, scale). The transformed public GSL
    // density/CDF/quantile is the same Gumbel maximum distribution.
    if (eq(n, "Gumbel")) {
        const double a = 1. / p[1];
        const double b = std::exp(p[0] / p[1]);
        if (!std::isfinite(a) || !std::isfinite(b) || b == 0)
            throw std::invalid_argument("GSL Gumbel parameter transform out of range");
        GSL_CONT("Gumbel", gsl_ran_gumbel1_pdf(x,a,b),
                 gsl_cdf_gumbel1_P(x,a,b), gsl_cdf_gumbel1_Q(x,a,b),
                 gsl_cdf_gumbel1_Pinv(x,a,b));
    }
#undef GSL_CONT
    if (x < 0 && op != 4 && (eq(n,"Bernoulli") || eq(n,"Binomial") ||
        eq(n,"Poisson") || eq(n,"Geometric") || eq(n,"NegativeBinomial") ||
        eq(n,"Hypergeometric")))
        return op == 2 ? 0 : op == 3 ? 1 : op == 1 ? -INFINITY : 0;
    const bool integer_x = x >= 0 && x <= static_cast<double>(UINT32_MAX) &&
                           std::floor(x) == x;
    if (op <= 1 && !integer_x)
        return op == 0 ? 0 : -INFINITY;
    // cdf and sf floor nonintegral values, matching the SciPy semantics.
    const unsigned int k = x >= 0 && x < static_cast<double>(UINT32_MAX)
        ? static_cast<unsigned int>(std::floor(x)) : UINT32_MAX;
#define GSL_DISC(NAME, PDF, CDF, SF) \
    if (eq(n, NAME)) { \
        if (op == 2) return (CDF); \
        if (op == 3) return (SF); \
        const double mass = (PDF); \
        return op == 1 ? std::log(mass) : mass; \
    }
    GSL_DISC("Bernoulli", gsl_ran_bernoulli_pdf(k,p[0]),
             gsl_cdf_binomial_P(k,p[0],1), gsl_cdf_binomial_Q(k,p[0],1));
    GSL_DISC("Binomial", gsl_ran_binomial_pdf(k,p[1],static_cast<unsigned int>(p[0])),
             gsl_cdf_binomial_P(k,p[1],static_cast<unsigned int>(p[0])),
             gsl_cdf_binomial_Q(k,p[1],static_cast<unsigned int>(p[0])));
    GSL_DISC("Poisson", gsl_ran_poisson_pdf(k,p[0]),
             gsl_cdf_poisson_P(k,p[0]), gsl_cdf_poisson_Q(k,p[0]));
    GSL_DISC("Geometric", gsl_ran_geometric_pdf(k,p[0]),
             gsl_cdf_geometric_P(k,p[0]), gsl_cdf_geometric_Q(k,p[0]));
    GSL_DISC("NegativeBinomial", gsl_ran_negative_binomial_pdf(k,p[1],p[0]),
             gsl_cdf_negative_binomial_P(k,p[1],p[0]),
             gsl_cdf_negative_binomial_Q(k,p[1],p[0]));
    GSL_DISC("Hypergeometric", gsl_ran_hypergeometric_pdf(k,
             static_cast<unsigned int>(p[1]), static_cast<unsigned int>(p[0]-p[1]),
             static_cast<unsigned int>(p[2])),
             gsl_cdf_hypergeometric_P(k,static_cast<unsigned int>(p[1]),
             static_cast<unsigned int>(p[0]-p[1]),static_cast<unsigned int>(p[2])),
             gsl_cdf_hypergeometric_Q(k,static_cast<unsigned int>(p[1]),
             static_cast<unsigned int>(p[0]-p[1]),static_cast<unsigned int>(p[2])));
#undef GSL_DISC
    throw std::invalid_argument("unsupported GSL distribution");
}
}
extern "C" {
const char* stats_distribution_error() { return distribution_error.c_str(); }
int stats_distribution_supported(int library, const char* name, const char* op) {
    return supported(library, name, operation(op)) ? arity(name) : 0;
}
int stats_distribution_evaluate(int library, const char* name, const char* op,
                                const double* p, int count, double x, double* output) {
    distribution_error.clear();
    gsl_set_error_handler_off();
    try {
        const int operation_id = operation(op);
        if (!name || !op || !p || !output) throw std::invalid_argument("null ABI argument");
        if (!supported(library, name, operation_id))
            throw std::invalid_argument("unsupported library distribution or operation");
        validate(name, p, count);
        if (std::isnan(x)) throw std::invalid_argument("NaN evaluation argument");
        if (operation_id == 4 && (x < 0 || x > 1))
            throw std::invalid_argument("probability outside [0,1]");
        // Mathematically exact nonfinite observation limits, without feeding
        // infinity into Boost APIs which usually reject nonfinite inputs.
        if (operation_id != 4 && std::isinf(x)) {
            if (operation_id <= 1) *output = operation_id == 0 ? 0 : -INFINITY;
            else if (operation_id == 2) *output = x > 0 ? 1 : 0;
            else *output = x > 0 ? 0 : 1;
            return 0;
        }
        if (support_limit(name, p, operation_id, x, output)) return 0;
        *output = library == 0 ? boost_eval(name,p,operation_id,x)
                               : gsl_eval(name,p,operation_id,x);
        if (std::isnan(*output)) throw std::runtime_error("distribution returned NaN");
        return 0;
    } catch (const std::exception& ex) {
        distribution_error = ex.what();
        return -1;
    } catch (...) {
        distribution_error = "unknown native exception";
        return -1;
    }
}
}
