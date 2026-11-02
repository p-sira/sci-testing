// Checked public-library counterparts to the remaining Sci statistics recipes.
// Libraries: 0 = Boost.Math, 1 = GSL.  Public library primitives perform
// estimates; transformations only account for Sci's saved signature.
#include <boost/math/statistics/univariate_statistics.hpp>
#include <boost/math/statistics/bivariate_statistics.hpp>
#include <boost/math/statistics/linear_regression.hpp>
#include <boost/math/statistics/t_test.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
thread_local std::string extra_error;
bool same(const char* a, const char* b) { return a && std::strcmp(a,b) == 0; }
bool trimmed(const char* name) {
    return same(name,"tmean") || same(name,"tvar") || same(name,"tstd") ||
           same(name,"tsem") || same(name,"tmin") || same(name,"tmax");
}
bool is_test(const char* name) {
    return same(name,"ttest_1samp") || same(name,"ttest_rel") ||
           same(name,"ttest_ind");
}
bool multivariate(const char* name) {
    return same(name,"dirichlet_pdf") || same(name,"dirichlet_logpdf") ||
           same(name,"multinomial_pmf") || same(name,"multinomial_logpmf");
}
bool gaussian(const char* name) {
    return same(name,"multivariate_normal_pdf") ||
           same(name,"multivariate_normal_logpdf");
}
bool available(int lib, const char* name) {
    if (lib != 0 && lib != 1) return false;
    if (same(name,"trim_mean")) return lib == 1; // GSL public trimmed mean
    if (same(name,"tmin") || same(name,"tmax")) return lib == 1;
    if (multivariate(name) || gaussian(name)) return lib == 1;
    return trimmed(name) || is_test(name) || same(name,"linregress");
}
int count(const char* name) { return same(name,"linregress") ? 6 : is_test(name) ? 3 : 1; }
int flag(double p) {
    if (p != 0. && p != 1.) throw std::invalid_argument("expected boolean 0 or 1");
    return static_cast<int>(p);
}
int alternative(double p) {
    if (p != 0. && p != 1. && p != 2.)
        throw std::invalid_argument("invalid alternative (0 two-sided, 1 less, 2 greater)");
    return static_cast<int>(p);
}
double mean(int lib, const std::vector<double>& x) {
    return lib == 0 ? boost::math::statistics::mean(x)
                    : gsl_stats_mean(x.data(),1,x.size());
}
double variance(int lib, const std::vector<double>& x, int ddof) {
    if (ddof < 0 || static_cast<size_t>(ddof) >= x.size())
        throw std::invalid_argument("invalid degrees of freedom");
    if (x.size() == 1) return 0.;
    const double n = static_cast<double>(x.size());
    const double pop = lib == 0 ? boost::math::statistics::variance(x)
        : gsl_stats_variance(x.data(),1,x.size()) * ((n-1)/n);
    return pop * (n / (n-ddof));
}
// SciPy tail convention: 0 two-sided, 1 less, 2 greater. These calls use
// public Student-t CDFs from the same library as the test statistic.
double student_p(int lib, double t, double df, int alt) {
    if (lib == 0) {
        const boost::math::students_t_distribution<double> dist(df);
        if (alt == 1) return boost::math::cdf(dist,t);
        if (alt == 2) return boost::math::cdf(boost::math::complement(dist,t));
        return 2 * boost::math::cdf(boost::math::complement(dist,std::abs(t)));
    }
    if (alt == 1) return gsl_cdf_tdist_P(t,df);
    if (alt == 2) return gsl_cdf_tdist_Q(t,df);
    return 2 * gsl_cdf_tdist_Q(std::abs(t),df);
}
std::vector<double> filter(const char* name, const std::vector<double>& x,
                           const double* p, int np) {
    const bool extrema = same(name,"tmin") || same(name,"tmax");
    if (extrema && np != 2) throw std::invalid_argument("tmin/tmax need limit,inclusive");
    if (!extrema && np != 5) throw std::invalid_argument("truncated stats need low,high,two inclusive flags,ddof");
    const double low = extrema && same(name,"tmin") ? p[0] :
                       extrema ? -INFINITY : p[0];
    const double high = extrema && same(name,"tmax") ? p[0] :
                        extrema ? INFINITY : p[1];
    const bool incl_low = extrema ? flag(p[1]) : flag(p[2]);
    const bool incl_high = extrema ? flag(p[1]) : flag(p[3]);
    if (low > high) throw std::invalid_argument("lower truncation bound exceeds upper");
    std::vector<double> kept;
    kept.reserve(x.size());
    for (double value : x)
        if ((value > low || (incl_low && value == low)) &&
            (value < high || (incl_high && value == high)))
            kept.push_back(value);
    if (kept.empty()) throw std::invalid_argument("empty truncated sample");
    return kept;
}
void test_evaluate(int lib, const char* name, const std::vector<double>& x,
                   const std::vector<double>& y, const double* p, int np,
                   double* result) {
    const bool single = same(name,"ttest_1samp");
    const bool paired = same(name,"ttest_rel");
    const bool independent = same(name,"ttest_ind");
    if (np != (paired ? 1 : 2))
        throw std::invalid_argument("incorrect t-test parameter arity");
    if (x.size() < 2 || (!single && (y.size() < 2 ||
        (paired && y.size() != x.size()))))
        throw std::invalid_argument("insufficient t-test samples");
    const int alt = alternative(p[paired ? 0 : 1]);
    double t = 0, df = 0;
    if (single || paired) {
        std::vector<double> diff;
        if (paired) {
            diff.reserve(x.size());
            for (size_t i = 0; i < x.size(); ++i) diff.push_back(x[i] - y[i]);
        }
        const auto& sample = paired ? diff : x;
        const double target = single ? p[0] : 0;
        df = static_cast<double>(sample.size() - 1);
        if (lib == 0) {
            // Invoke the actual public Boost t-test (which returns t and a
            // two-sided p-value); convert its tail with the public t CDF.
            const auto stats = paired ?
                boost::math::statistics::paired_samples_t_test(x,y) :
                boost::math::statistics::one_sample_t_test(sample,target);
            t = stats.first;
        } else {
            t = (mean(lib,sample)-target) /
                std::sqrt(variance(lib,sample,1)/sample.size());
        }
    } else if (independent) {
        const bool pooled = flag(p[0]);
        const double n1 = static_cast<double>(x.size());
        const double n2 = static_cast<double>(y.size());
        const double s1 = variance(lib,x,1);
        const double s2 = variance(lib,y,1);
        const double diff = mean(lib,x)-mean(lib,y);
        if (pooled) {
            df = n1+n2-2;
            if (lib == 1) t = gsl_stats_ttest(x.data(),1,x.size(),y.data(),1,y.size());
            else {
                const double sp = ((n1-1)*s1+(n2-1)*s2)/df;
                t = diff / std::sqrt(sp*(1/n1+1/n2));
            }
        } else {
            const double a = s1/n1, b = s2/n2;
            t = diff/std::sqrt(a+b);
            df = ((a+b)*(a+b)) /
                 (a*a/(n1-1)+b*b/(n2-1));
        }
    }
    result[0] = t;
    result[1] = student_p(lib,t,df,alt);
    result[2] = df;
}
void regress(int lib, const std::vector<double>& x, const std::vector<double>& y,
             const double* p, int np, double* result) {
    if (np != 1 || x.size() != y.size() || x.size() < 3)
        throw std::invalid_argument("linregress requires matching sample lengths >=3 and alternative");
    const int alt = alternative(p[0]);
    double intercept, slope;
    if (lib == 0) {
        const auto fit = boost::math::statistics::simple_ordinary_least_squares(x,y);
        intercept = fit.first;
        slope = fit.second;
    } else {
        double cov00, cov01, cov11, sumsq;
        const int status = gsl_fit_linear(x.data(),1,y.data(),1,x.size(),
            &intercept,&slope,&cov00,&cov01,&cov11,&sumsq);
        if (status != GSL_SUCCESS)
            throw std::runtime_error(gsl_strerror(status));
    }
    const double vx = variance(lib,x,0), vy = variance(lib,y,0);
    if (!(vx > 0)) throw std::invalid_argument("all x values identical");
    if (!(vy > 0)) throw std::invalid_argument("constant y: Sci correlation undefined");
    double r = lib == 0 ?
        boost::math::statistics::correlation_coefficient(x,y) :
        gsl_stats_correlation(x.data(),1,y.data(),1,x.size());
    // A mathematically exact fit may leave r = 1 - ulp in GSL; at this
    // boundary Sci's slope/intercept standard errors are exactly zero.
    bool exact_fit = true;
    for (size_t i = 0; i < x.size(); ++i)
        if (std::fma(slope,x[i],intercept) != y[i]) { exact_fit = false; break; }
    if (exact_fit && slope != 0) r = std::copysign(1.,slope);
    const double clipped_r = std::clamp(r,-1.,1.);
    const double df = static_cast<double>(x.size()-2);
    const double slope_stderr = std::sqrt(std::max(0.,1-clipped_r*clipped_r)*vy/vx/df);
    const double intercept_stderr = slope_stderr*std::sqrt(vx+std::pow(mean(lib,x),2));
    const double t = std::abs(clipped_r) == 1
        ? std::copysign(INFINITY,clipped_r)
        : clipped_r*std::sqrt(df/((1-clipped_r)*(1+clipped_r)));
    double pvalue = std::isinf(t) ?
        (alt == 0 ? 0 : alt == 1 ? (t > 0 ? 1 : 0) : (t > 0 ? 0 : 1)) :
        student_p(lib,t,df,alt);
    result[0] = slope;
    result[1] = intercept;
    result[2] = r;
    result[3] = pvalue;
    result[4] = slope_stderr;
    result[5] = intercept_stderr;
}
double gsl_multivariate(const char* name, const std::vector<double>& x,
                        const std::vector<double>& y, int np) {
    if (np != 0 || x.empty() || x.size() != y.size())
        throw std::invalid_argument("multivariate density needs matching vectors");
    if (same(name,"dirichlet_pdf") || same(name,"dirichlet_logpdf")) {
        double total = 0;
        for (size_t i = 0; i < x.size(); ++i) {
            if (!(x[i] > 0 && x[i] < 1 && y[i] > 0))
                throw std::invalid_argument("GSL Dirichlet interior domain requires positive components");
            total += x[i];
        }
        if (std::abs(total-1.) > 1e-12)
            throw std::invalid_argument("Dirichlet components do not sum to one");
        return same(name,"dirichlet_pdf") ?
            gsl_ran_dirichlet_pdf(x.size(),y.data(),x.data()) :
            gsl_ran_dirichlet_lnpdf(x.size(),y.data(),x.data());
    }
    double total = 0;
    std::vector<unsigned int> counts;
    counts.reserve(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        if (!(x[i] >= 0 && x[i] <= 1048576 && std::floor(x[i]) == x[i]) ||
            y[i] < 0 || y[i] > 1)
            throw std::invalid_argument("invalid GSL multinomial count or probability");
        counts.push_back(static_cast<unsigned int>(x[i]));
        total += y[i];
    }
    if (std::abs(total-1.) > 1e-12)
        throw std::invalid_argument("multinomial probabilities do not sum to one");
    return same(name,"multinomial_pmf") ?
        gsl_ran_multinomial_pdf(x.size(),y.data(),counts.data()) :
        gsl_ran_multinomial_lnpdf(x.size(),y.data(),counts.data());
}
double gsl_gaussian(const char* name, std::vector<double>& x,
                    std::vector<double>& y, const double* p, int np) {
    const size_t dimension = x.size();
    if (y.size() != dimension || dimension == 0 ||
        dimension > 1024 || np != static_cast<int>(dimension*dimension))
        throw std::invalid_argument("Gaussian requires square covariance matrix");
    std::vector<double> covariance(p,p+np), work_values(dimension);
    auto matrix = gsl_matrix_view_array(covariance.data(),dimension,dimension);
    auto observation = gsl_vector_view_array(x.data(),dimension);
    auto location = gsl_vector_view_array(y.data(),dimension);
    auto work = gsl_vector_view_array(work_values.data(),dimension);
    const int cholesky = gsl_linalg_cholesky_decomp1(&matrix.matrix);
    if (cholesky != GSL_SUCCESS) throw std::runtime_error(gsl_strerror(cholesky));
    double result = 0;
    const int status = same(name,"multivariate_normal_pdf") ?
        gsl_ran_multivariate_gaussian_pdf(&observation.vector,&location.vector,
            &matrix.matrix,&result,&work.vector) :
        gsl_ran_multivariate_gaussian_log_pdf(&observation.vector,&location.vector,
            &matrix.matrix,&result,&work.vector);
    if (status != GSL_SUCCESS) throw std::runtime_error(gsl_strerror(status));
    return result;
}
}

extern "C" {
const char* stats_extra_error() { return extra_error.c_str(); }
int stats_extra_supported(int lib, const char* name) { return available(lib,name) ? 1 : 0; }
int stats_extra_output_count(const char* name) {
    return available(0,name) || available(1,name) ? count(name) : 0;
}
int stats_extra_evaluate(int lib, const char* name,
                         const double* x, int nx, const double* y, int ny,
                         const double* p, int np, double* out, int capacity) {
    extra_error.clear();
    gsl_set_error_handler_off();
    try {
        if (!available(lib,name)) throw std::invalid_argument("unsupported direct comparison");
        if (!x || nx < 1 || nx > 1048576 || !out ||
            np < 0 || (np && !p) || ny < 0 || (ny && !y) ||
            capacity < count(name)) throw std::invalid_argument("invalid ABI inputs");
        std::vector<double> a(x,x+nx), b;
        if (ny) b.assign(y,y+ny);
        for (double v : a) if (!std::isfinite(v)) throw std::invalid_argument("nonfinite data");
        for (double v : b) if (!std::isfinite(v)) throw std::invalid_argument("nonfinite data");
        for (int i=0;i<np;++i) if (std::isnan(p[i])) throw std::invalid_argument("NaN parameter");
        if (gaussian(name)) {
            out[0] = gsl_gaussian(name,a,b,p,np);
        } else if (multivariate(name)) {
            out[0] = gsl_multivariate(name,a,b,np);
        } else if (same(name,"linregress")) {
            regress(lib,a,b,p,np,out);
        } else if (is_test(name)) {
            test_evaluate(lib,name,a,b,p,np,out);
        } else {
            if (ny) throw std::invalid_argument("unexpected secondary sample");
            if (same(name,"trim_mean")) {
                if (np != 1 || p[0] < 0 || p[0] >= .5 || !std::isfinite(p[0]))
                    throw std::invalid_argument("invalid trim proportion");
                std::sort(a.begin(),a.end());
                out[0] = gsl_stats_trmean_from_sorted_data(p[0],a.data(),1,a.size());
            } else {
                auto data = filter(name,a,p,np);
                if (same(name,"tmin")) out[0] = gsl_stats_min(data.data(),1,data.size());
                else if (same(name,"tmax")) out[0] = gsl_stats_max(data.data(),1,data.size());
                else if (same(name,"tmean")) out[0] = mean(lib,data);
                else {
                    const int correction = flag(p[4]); // saved fixtures: ddof 0/1
                    double var = variance(lib,data,correction);
                    if (same(name,"tvar")) out[0] = var;
                    else if (same(name,"tstd")) out[0] = std::sqrt(var);
                    else if (same(name,"tsem")) out[0] = std::sqrt(var/data.size());
                }
            }
        }
        for (int i=0;i<count(name);++i)
            if (!std::isfinite(out[i])) throw std::runtime_error("nonfinite statistic");
        return 0;
    } catch (const std::exception& ex) { extra_error = ex.what(); return -1; }
      catch (...) { extra_error = "unknown native exception"; return -1; }
}
}
