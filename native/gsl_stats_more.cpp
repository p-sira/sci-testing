// GSL-derived statistics for saved Sci 0.3 recipe signatures.
// Provider 1 is GSL. All quantities use public GSL 2.7/2.8 statistics
// (mean, sample variance, total sum of squares) and distribution CDFs.
// Composite tests and circular summaries are mathematical reductions;
// none is misrepresented as a standalone GSL function.
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_statistics_double.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
thread_local std::string last_error;
constexpr double pi = 3.141592653589793238462643383279502884;
constexpr double two_pi = 2.0 * pi;

bool same(const char* a, const char* b) { return a && std::strcmp(a, b) == 0; }

int width(int provider, const char* name) {
    if (provider != 1 || !name) return 0;
    if (same(name, "gmean") || same(name, "hmean") ||
        same(name, "circmean") || same(name, "circstd") ||
        same(name, "circvar")) return 1;
    if (same(name, "jarque_bera") || same(name, "chisquare") ||
        same(name, "f_oneway") || same(name, "bartlett")) return 2;
    if (same(name, "ttest_ind_from_stats")) return 3;
    return 0;
}

void check_finite(const double* data, int count, const char* label) {
    if (count < 0 || count > 100000 || (count && !data))
        throw std::invalid_argument(std::string("invalid ") + label + " buffer or length");
    for (int i = 0; i < count; ++i)
        if (!std::isfinite(data[i]))
            throw std::invalid_argument(std::string("nonfinite ") + label + " input");
}

double mean(const double* x, size_t n) { return gsl_stats_mean(x, 1, n); }
double variance(const double* x, size_t n) { return gsl_stats_variance(x, 1, n); }
double tss(const double* x, size_t n) { return gsl_stats_tss(x, 1, n); }

void require(bool condition, const char* message) {
    if (!condition) throw std::invalid_argument(message);
}

int whole(double x, int low, int high, const char* label) {
    if (!std::isfinite(x) || x != std::trunc(x) || x < low || x > high)
        throw std::invalid_argument(std::string("invalid ") + label);
    return static_cast<int>(x);
}

double p_t(double statistic, double df, int alt) {
    if (alt == 1) return gsl_cdf_tdist_P(statistic, df);
    if (alt == 2) return gsl_cdf_tdist_Q(statistic, df);
    return 2.0 * gsl_cdf_tdist_Q(std::abs(statistic), df);
}

void single_sample(const char* name, const double* x, int nx,
                   int ny, const double* p, int np,
                   double* out) {
    require(nx > 0 && !ny, "expected a nonempty sample without y");
    const bool geometric = same(name, "gmean");
    const bool harmonic = same(name, "hmean");
    const bool circular = same(name, "circmean") || same(name, "circstd") ||
                          same(name, "circvar");
    if (geometric || harmonic) {
        require(np == 0, "gmean/hmean expect no parameters");
        std::vector<double> transformed;
        transformed.reserve(nx);
        for (int i = 0; i < nx; ++i) {
            require(x[i] >= 0.0, "geometric/harmonic mean requires nonnegative observations");
            if (x[i] == 0.0) { out[0] = 0.0; return; }
            transformed.push_back(geometric ? std::log(x[i]) : 1.0 / x[i]);
        }
        // Both means are derived from public gsl_stats_mean.
        out[0] = geometric ? std::exp(mean(transformed.data(), transformed.size()))
                           : 1.0 / mean(transformed.data(), transformed.size());
        return;
    }
    if (circular) {
        const bool is_std = same(name, "circstd");
        require(np == (is_std ? 3 : 2), "circular statistic requires low/high and optional normalize flag");
        const double low = p[0], high = p[1];
        require(high > low && std::isfinite(high-low), "invalid circular interval");
        const int normal = is_std ? whole(p[2], 0, 1, "circular normalize flag") : 0;
        const double period = high-low;
        std::vector<double> s(nx), c(nx);
        for (int i = 0; i < nx; ++i) {
            const double angle = two_pi * ((x[i]-low) / period);
            s[i] = std::sin(angle);
            c[i] = std::cos(angle);
        }
        const double sbar = mean(s.data(), s.size());
        const double cbar = mean(c.data(), c.size());
        if (same(name, "circmean")) {
            double angle = std::atan2(sbar, cbar);
            if (angle < 0) angle += two_pi;
            out[0] = low + period * angle / two_pi;
        } else {
            const double r = std::min(1.0, std::hypot(sbar, cbar));
            if (same(name, "circvar")) out[0] = 1.0 - r;
            else {
                require(r > 0, "undefined circular standard deviation at zero mean resultant");
                const double sd = std::sqrt(-2.0 * std::log(r));
                out[0] = normal ? sd : sd*period/two_pi;
            }
        }
        return;
    }
    require(same(name, "jarque_bera") && np == 0 && nx >= 2,
            "Jarque-Bera requires at least two values and no parameters");
    const double center = mean(x, nx);
    const double central_second = tss(x, nx) / nx;
    require(central_second > 0 && std::isfinite(central_second),
            "Jarque-Bera undefined for constant sample");
    std::vector<double> cubes(nx), fourths(nx);
    for (int i = 0; i < nx; ++i) {
        const double delta = x[i]-center;
        cubes[i] = delta*delta*delta;
        fourths[i] = delta*delta*delta*delta;
    }
    const double skew = mean(cubes.data(), cubes.size()) /
                        std::pow(central_second, 1.5);
    const double kurt_excess = mean(fourths.data(), fourths.size()) /
                               (central_second*central_second) - 3.0;
    out[0] = nx/6.0 * (skew*skew + kurt_excess*kurt_excess/4.0);
    out[1] = gsl_cdf_chisq_Q(out[0], 2.0);
}

void chisquare(const double* x, int nx, const double* y, int ny,
               const double* p, int np, double* out) {
    require(nx >= 2 && np == 1 && (ny == 0 || ny == nx),
            "chisquare requires observations, optional equal-length expectations and ddof");
    const int ddof = whole(p[0], 0, nx-2, "chi-square ddof");
    const double sx = std::accumulate(x, x+nx, 0.0);
    require(sx > 0, "chi-square requires positive observation total");
    const double constant = ny ? 0.0 : mean(x, nx);
    const double sy = ny ? std::accumulate(y, y+ny, 0.0) : sx;
    require(std::abs(sx-sy) <= 1e-8 * std::max(sx, sy),
            "chi-square observed and expected totals disagree");
    double chi = 0.0;
    for (int i = 0; i < nx; ++i) {
        const double expected = ny ? y[i] : constant;
        require(x[i] >= 0 && expected > 0, "chi-square requires nonnegative observations and positive expectations");
        const double delta = x[i]-expected;
        chi += delta*delta/expected;
    }
    out[0] = chi;
    out[1] = gsl_cdf_chisq_Q(chi, nx-1-ddof);
}

void ttest_from_stats(const double* p, int np, double* out) {
    require(np == 8, "ttest_ind_from_stats needs mean1,std1,n1,mean2,std2,n2,equal_var,alternative");
    const int n1 = whole(p[2], 2, 100000000, "first t-test sample count");
    const int n2 = whole(p[5], 2, 100000000, "second t-test sample count");
    const int pooled = whole(p[6], 0, 1, "equal_var flag");
    const int alt = whole(p[7], 0, 2, "alternative flag");
    require(p[1] >= 0.0 && p[4] >= 0.0, "negative input standard deviation");
    const double s1 = p[1]*p[1], s2 = p[4]*p[4];
    double denom, df;
    if (pooled) {
        df = double(n1)+n2-2;
        denom = ((n1-1.)*s1 + (n2-1.)*s2)/df*(1./n1 + 1./n2);
    } else {
        const double a = s1/n1, b = s2/n2;
        denom = a+b;
        const double welch = a*a/(n1-1.)+b*b/(n2-1.);
        require(welch > 0, "undefined Welch degrees of freedom");
        df = denom*denom/welch;
    }
    require(denom > 0, "undefined t statistic due to zero pooled variance");
    out[0] = (p[0]-p[3])/std::sqrt(denom);
    out[1] = p_t(out[0],df,alt);
    out[2] = df;
}

void groups_test(const char* name, const double* g, int ng,
                 const int* lengths, int n_groups, double* out) {
    require(n_groups >= 2 && n_groups <= 256 && ng >= 4 && lengths,
            "group test expects two or more group length entries");
    int used = 0;
    std::vector<double> means, variances;
    means.reserve(n_groups);
    variances.reserve(n_groups);
    double within = 0.0, sum_reciprocal_df = 0.0;
    for (int j = 0; j < n_groups; ++j) {
        const int count = lengths[j];
        require(count >= 2 && count <= ng-used, "invalid group size");
        const double* v = g+used;
        const double m = mean(v, count);
        const double variance_j = variance(v, count);
        means.push_back(m);
        variances.push_back(variance_j);
        within += tss(v, count);
        sum_reciprocal_df += 1.0/(count-1);
        used += count;
    }
    require(used == ng && ng > n_groups, "group lengths do not match flattened input");
    if (same(name, "f_oneway")) {
        require(within > 0, "ANOVA requires positive within-group variation");
        const double grand = mean(g, ng);
        double between = 0.0;
        for (int i = 0; i < n_groups; ++i) {
            const double difference = means[i]-grand;
            between += lengths[i] * difference*difference;
        }
        out[0] = (between/(n_groups-1)) / (within/(ng-n_groups));
        out[1] = gsl_cdf_fdist_Q(out[0],n_groups-1,ng-n_groups);
    } else {
        require(within > 0, "Bartlett requires positive pooled sample variance");
        double weighted_logs = 0.0;
        for (int i = 0; i < n_groups; ++i) {
            require(variances[i] > 0, "Bartlett undefined for zero-variance group");
            weighted_logs += (lengths[i]-1) * std::log(variances[i]);
        }
        const double df = ng-n_groups;
        const double correction = 1.0 + (sum_reciprocal_df - 1.0/df) /
                                         (3.0*(n_groups-1));
        out[0] = (df*std::log(within/df)-weighted_logs)/correction;
        if (out[0] < 0 && out[0] > -1e-12) out[0] = 0.0;
        require(out[0] >= 0, "negative Bartlett statistic from numerical cancellation");
        out[1] = gsl_cdf_chisq_Q(out[0],n_groups-1);
    }
}
} // namespace

extern "C" {
int gsl_stats_more_supported(int library, const char* name) { return width(library,name); }
const char* gsl_stats_more_error() { return last_error.c_str(); }

int gsl_stats_more_evaluate(int library, const char* name,
                            const double* x, int nx, const double* y, int ny,
                            const double* p, int np, const double* g, int ng,
                            const int* group_sizes, int n_groups,
                            double* output, int capacity) {
    last_error.clear();
    try {
        const int result_size = width(library,name);
        require(result_size > 0 && output && capacity >= result_size,
                "unsupported library/API or insufficient output capacity");
        check_finite(x,nx,"x");
        check_finite(y,ny,"y");
        check_finite(p,np,"p");
        check_finite(g,ng,"g");
        require(n_groups >= 0 && n_groups <= 256,
                "invalid number of groups");
        gsl_set_error_handler_off();
        double computed[3] = {0.0, 0.0, 0.0};
        if (same(name,"ttest_ind_from_stats")) {
            require(nx == 0 && ny == 0 && ng == 0 && n_groups == 0,
                    "t-test summary accepts parameters only");
            ttest_from_stats(p,np,computed);
        } else if (same(name,"f_oneway") || same(name,"bartlett")) {
            require(nx == 0 && ny == 0 && np == 0,
                    "group test accepts group data only");
            groups_test(name,g,ng,group_sizes,n_groups,computed);
        } else {
            require(ng == 0 && n_groups == 0,
                    "scalar statistic does not accept group input");
            if (same(name,"chisquare")) chisquare(x,nx,y,ny,p,np,computed);
            else single_sample(name,x,nx,ny,p,np,computed);
        }
        for (int j=0; j<result_size; ++j)
            require(std::isfinite(computed[j]),"nonfinite GSL-derived statistic");
        for (int j=0; j<result_size; ++j) output[j]=computed[j];
        return 0;
    } catch (const std::exception& exc) {
        last_error = exc.what();
    } catch (...) {
        last_error = "unknown GSL-derived statistical failure";
    }
    return -1;
}
} // extern "C"
