// Independent public Boost.Math/GSL statistics comparator.
//
// ABI: stats_summary_supported(library, name) -> bool;
// stats_summary_output_count(name, sample_count) -> output count;
// stats_summary_evaluate(library, name, x, nx, y, ny, p, np, out, capacity)
// returns 0 on success, -1 on unsupported/invalid input. Library 0 is Boost,
// library 1 GSL. All arrays and outputs contain binary64 values.
//
// Sci fixture parameters: variance/std/sem/variation [ddof],
// skew [bias], kurtosis [fisher, bias], median_abs_deviation [scale],
// quantile [probability], scoreatpercentile [percentile],
// describe [ddof, bias], correlation [alternative: 0 two-sided, 1 less,
// 2 greater]. Correlation results are [statistic, pvalue].
#include <boost/math/statistics/univariate_statistics.hpp>
#include <boost/math/statistics/bivariate_statistics.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_cdf.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
thread_local std::string summary_error;
namespace bs = boost::math::statistics;

bool is(const std::string& name, const char* choice) { return name == choice; }
bool paired(const std::string& name) {
    return is(name, "pearsonr") || is(name, "spearmanr") ||
           is(name, "spearmanrho") || is(name, "pointbiserialr");
}
bool gsl_only(const std::string& name) {
    return is(name, "quantile") || is(name, "scoreatpercentile") ||
           is(name, "iqr") || is(name, "spearmanr") ||
           is(name, "spearmanrho");
}
bool known(const std::string& name) {
    return is(name, "mean") || is(name, "variance") || is(name, "std") ||
           is(name, "sem") || is(name, "variation") || is(name, "skew") ||
           is(name, "kurtosis") || is(name, "median") ||
           is(name, "median_abs_deviation") || is(name, "quantile") ||
           is(name, "scoreatpercentile") || is(name, "iqr") ||
           is(name, "describe") || is(name, "zscore") ||
           is(name, "pearsonr") || is(name, "pointbiserialr") ||
           is(name, "spearmanr") || is(name, "spearmanrho");
}
bool valid_count(const std::string& name, int np) {
    if (is(name, "variance") || is(name, "std") || is(name, "sem") ||
        is(name, "variation") || is(name, "skew") || is(name, "quantile") ||
        is(name, "scoreatpercentile") || is(name, "median_abs_deviation") ||
        is(name, "zscore"))
        return np == 1;
    if (is(name, "kurtosis") || is(name, "describe")) return np == 2;
    if (is(name, "pearsonr") || is(name, "spearmanr") ||
        is(name, "spearmanrho")) return np == 1;
    return np == 0;
}
int integral_parameter(double parameter, int low, int high) {
    if (!std::isfinite(parameter) || parameter < low || parameter > high ||
        std::floor(parameter) != parameter)
        throw std::invalid_argument("invalid integer parameter");
    return static_cast<int>(parameter);
}
int ddof(const double* p, size_t n) {
    return integral_parameter(p[0], 0, static_cast<int>(n - 1));
}
bool flag(double value) {
    return integral_parameter(value, 0, 1) != 0;
}
double raw_mean(int lib, const std::vector<double>& x) {
    return lib == 0 ? bs::mean(x) : gsl_stats_mean(x.data(), 1, x.size());
}
double raw_variance(int lib, const std::vector<double>& x, int correction) {
    if (x.size() == 1) return 0.;
    const double n = static_cast<double>(x.size());
    // Boost variance is population-normalized, GSL variance uses n-1.
    const double population = lib == 0 ? bs::variance(x)
        : gsl_stats_variance(x.data(), 1, x.size()) * ((n - 1.) / n);
    return population * (n / (n - correction));
}
double raw_skew(int lib, const std::vector<double>& x, bool bias) {
    const double n = static_cast<double>(x.size());
    if (x.size() < 2) return std::numeric_limits<double>::quiet_NaN();
    if (raw_variance(lib, x, 0) == 0.)
        return std::numeric_limits<double>::quiet_NaN();
    // Boost's public sequential skewness already matches m3/m2^(3/2).
    // GSL skew uses sample SD and must be converted to population moments.
    double skew = lib == 0 ? bs::skewness(x)
        : gsl_stats_skew(x.data(), 1, x.size());
    if (lib == 1) skew *= std::pow(n / (n - 1.), 1.5);
    if (!bias && n > 2.) skew *= std::sqrt(n * (n - 1.)) / (n - 2.);
    return skew;
}
double raw_kurtosis(int lib, const std::vector<double>& x, bool fisher, bool bias) {
    const double n = static_cast<double>(x.size());
    if (x.size() < 2 || raw_variance(lib, x, 0) == 0.)
        return std::numeric_limits<double>::quiet_NaN();
    // Boost kurtosis is population m4/m2^2, GSL uses sample SD and
    // returns excess, so both denominator and offset require correction.
    double k = lib == 0 ? bs::kurtosis(x)
        : (gsl_stats_kurtosis(x.data(), 1, x.size()) + 3.) *
          std::pow(n / (n - 1.), 2.);
    if (!bias && n > 3.)
        k = ((n * n - 1.) * k - 3. * (n - 1.) * (n - 1.)) /
            ((n - 2.) * (n - 3.)) + 3.;
    return fisher ? k - 3. : k;
}
double raw_median(int lib, const std::vector<double>& x) {
    auto copy = x; // Both public libraries may rearrange their input.
    return lib == 0 ? bs::median(copy) : gsl_stats_median(copy.data(), 1, copy.size());
}
double raw_quantile(const std::vector<double>& x, double probability) {
    if (!(probability >= 0. && probability <= 1.))
        throw std::invalid_argument("quantile probability must be in [0,1]");
    auto sorted = x;
    std::sort(sorted.begin(), sorted.end());
    return gsl_stats_quantile_from_sorted_data(sorted.data(), 1, sorted.size(), probability);
}
double raw_mad(int lib, const std::vector<double>& x) {
    if (lib == 0) {
        auto copy = x;
        return bs::median_absolute_deviation(copy);
    }
    std::vector<double> work(x.size());
    return gsl_stats_mad0(x.data(), 1, x.size(), work.data());
}
double raw_correlation(int lib, const std::vector<double>& x, const std::vector<double>& y,
                       bool ranks) {
    if (ranks) {
        if (lib == 0) throw std::invalid_argument("Boost.Math has no direct Spearman statistic");
        std::vector<double> work(2 * x.size());
        return gsl_stats_spearman(x.data(), 1, y.data(), 1, x.size(), work.data());
    }
    return lib == 0 ? bs::correlation_coefficient(x, y)
        : gsl_stats_correlation(x.data(), 1, y.data(), 1, x.size());
}
double pvalue(int lib, double r, size_t count, int alternative, bool approximate_spearman) {
    if (std::isnan(r)) return std::numeric_limits<double>::quiet_NaN();
    r = std::max(-1., std::min(1., r));
    if (count == 2) return approximate_spearman ?
        std::numeric_limits<double>::quiet_NaN() : 1.;
    if (approximate_spearman) {
        const double df = static_cast<double>(count - 2);
        const double t = r * std::sqrt(df / ((1. - r) * (1. + r)));
        if (lib == 0) {
            const boost::math::students_t_distribution<double> dist(df);
            if (alternative == 1) return boost::math::cdf(dist, t);
            if (alternative == 2) return boost::math::cdf(boost::math::complement(dist, t));
            return 2. * boost::math::cdf(boost::math::complement(dist, std::abs(t)));
        }
        if (alternative == 1) return gsl_cdf_tdist_P(t, df);
        if (alternative == 2) return gsl_cdf_tdist_Q(t, df);
        return 2. * gsl_cdf_tdist_Q(std::abs(t), df);
    }
    const double a = static_cast<double>(count) / 2. - 1.;
    // Under the zero-correlation null, (R+1)/2 ~ Beta(a,a).
    const double x = (alternative == 0) ? (1. - std::abs(r)) / 2. : (1. + r) / 2.;
    const double tail = lib == 0 ? (alternative == 2 ?
        boost::math::ibetac(a, a, x) : boost::math::ibeta(a, a, x))
        : (alternative == 2 ? gsl_cdf_beta_Q(x, a, a)
           : gsl_cdf_beta_P(x, a, a));
    return std::clamp(alternative == 0 ? 2. * tail : tail, 0., 1.);
}
}

extern "C" {
const char* stats_summary_error() { return summary_error.c_str(); }

int stats_summary_supported(int library, const char* name) {
    if ((library != 0 && library != 1) || !name) return 0;
    const std::string n(name);
    return known(n) && !(library == 0 && gsl_only(n));
}

int stats_summary_output_count(const char* name, int sample_count) {
    if (!name || sample_count <= 0 || !known(name)) return 0;
    const std::string n(name);
    if (is(n, "describe")) return 7;
    if (paired(n)) return 2;
    if (is(n, "zscore")) return sample_count;
    return 1;
}

int stats_summary_evaluate(int library, const char* name,
                           const double* x, int nx, const double* y, int ny,
                           const double* p, int np, double* output, int capacity) {
    summary_error.clear();
    try {
        if (!stats_summary_supported(library, name))
            throw std::invalid_argument("unsupported native statistics API");
        const std::string operation(name);
        if (!x || nx <= 0 || nx > 1048576 || !output || np < 0 ||
            (np && !p) || !valid_count(operation, np))
            throw std::invalid_argument("invalid input or parameter shape");
        if (capacity < stats_summary_output_count(name, nx))
            throw std::invalid_argument("insufficient output capacity");
        const std::vector<double> values(x, x + nx);
        if (std::any_of(values.begin(), values.end(), [](double v){ return !std::isfinite(v); }))
            throw std::invalid_argument("nonfinite observations");
        const bool is_pair = paired(operation);
        if (is_pair && (!y || ny != nx || nx < 2))
            throw std::invalid_argument("paired inputs must have equal lengths >= 2");
        if (!is_pair && ny != 0)
            throw std::invalid_argument("unexpected second sample");
        const std::vector<double> second = is_pair ? std::vector<double>(y, y + ny)
                                                   : std::vector<double>();
        if (std::any_of(second.begin(), second.end(), [](double v){ return !std::isfinite(v); }))
            throw std::invalid_argument("nonfinite paired observations");
        if (is(operation, "mean")) output[0] = raw_mean(library, values);
        else if (is(operation, "variance") || is(operation, "std") ||
                 is(operation, "sem") || is(operation, "variation") ||
                 is(operation, "zscore")) {
            const int correction = ddof(p, values.size());
            const double var = raw_variance(library, values, correction);
            const double standard = std::sqrt(var);
            if (is(operation, "variance")) output[0] = var;
            else if (is(operation, "std")) output[0] = standard;
            else if (is(operation, "sem")) {
                if (nx < 2) throw std::invalid_argument("sem requires at least two observations");
                output[0] = standard / std::sqrt(static_cast<double>(nx));
            } else if (is(operation, "variation")) {
                const double mean = raw_mean(library, values);
                output[0] = standard / mean;
            } else {
                const double mean = raw_mean(library, values);
                for (int i = 0; i < nx; ++i)
                    output[i] = (values[i] - mean) / standard;
            }
        } else if (is(operation, "skew")) output[0] = raw_skew(library, values, flag(p[0]));
        else if (is(operation, "kurtosis"))
            output[0] = raw_kurtosis(library, values, flag(p[0]), flag(p[1]));
        else if (is(operation, "median")) output[0] = raw_median(library, values);
        else if (is(operation, "median_abs_deviation")) {
            if (!(p[0] > 0.) || !std::isfinite(p[0]))
                throw std::invalid_argument("median absolute deviation scale must be positive");
            output[0] = raw_mad(library, values) / p[0];
        } else if (is(operation, "quantile")) output[0] = raw_quantile(values, p[0]);
        else if (is(operation, "scoreatpercentile")) output[0] = raw_quantile(values, p[0] / 100.);
        else if (is(operation, "iqr")) output[0] = raw_quantile(values, .75) - raw_quantile(values, .25);
        else if (is(operation, "describe")) {
            if (nx < 2) throw std::invalid_argument("describe requires at least two observations");
            auto bounds = std::minmax_element(values.begin(), values.end());
            output[0] = nx;
            output[1] = *bounds.first;
            output[2] = *bounds.second;
            output[3] = raw_mean(library, values);
            output[4] = raw_variance(library, values, ddof(p, values.size()));
            output[5] = raw_skew(library, values, flag(p[1]));
            output[6] = raw_kurtosis(library, values, true, flag(p[1]));
        } else if (is_pair) {
            const int alt = is(operation, "pointbiserialr") ? 0
                : integral_parameter(p[0], 0, 2);
            if (is(operation, "pointbiserialr")) {
                if (std::any_of(values.begin(), values.end(), [](double v){ return v != 0. && v != 1.; }))
                    throw std::invalid_argument("pointbiserialr expects binary first sample");
            }
            const bool ranks = is(operation, "spearmanr") || is(operation, "spearmanrho");
            const double r = raw_correlation(library, values, second, ranks);
            output[0] = r;
            output[1] = pvalue(library, r, values.size(), alt, is(operation, "spearmanr"));
        } else throw std::invalid_argument("unsupported operation");
        return 0;
    } catch (const std::exception& exception) {
        summary_error = exception.what();
    } catch (...) {
        summary_error = "unknown native statistics failure";
    }
    return -1;
}
}
