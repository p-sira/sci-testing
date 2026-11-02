// SciPy-compatible real sample reductions using PUBLIC GSL 2.7/2.8
// gsl_stats_mean, gsl_stats_variance and gsl_stats_tss.
//
// These are GSL-derived formulas, not purported GSL direct equivalents:
// gstd = exp(sqrt(variance(log(x), ddof))),
// gzscore = (log(x)-mean(log(x)))/std(log(x), ddof),
// pmean = mean(x**p)**(1/p), with p=0 the geometric mean,
// moment = mean((x-mean(x))**order).
// Only scalar one-dimensional, unweighted samples are currently supported;
// saved recipes use x=[...], y=[], g=[], p=[one scalar].

#include <gsl/gsl_statistics_double.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
thread_local std::string last_error;

bool named(const char* given, const char* expected) {
    return given && std::strcmp(given,expected)==0;
}

bool supported(int lib, const char* name) {
    return lib == 1 && (named(name,"gstd") || named(name,"gzscore") ||
                        named(name,"pmean") || named(name,"moment"));
}

int output_count(const char* name, int nx) {
    if (nx <= 0 || nx > 100000 || !supported(1,name)) return 0;
    return named(name,"gzscore") ? nx : 1;
}

void require(bool valid, const char* message) {
    if (!valid) throw std::invalid_argument(message);
}

int integer(double value, int min_value, int max_value, const char* label) {
    if (!std::isfinite(value) || value != std::trunc(value) ||
        value < min_value || value > max_value)
        throw std::invalid_argument(std::string("invalid ")+label);
    return static_cast<int>(value);
}

double mean(const std::vector<double>& values) {
    return gsl_stats_mean(values.data(),1,values.size());
}

double popvar(const std::vector<double>& values) {
    if (values.size() == 1) return 0.0;
    // GSL variance uses Bessel correction (n-1). Recover population var.
    return gsl_stats_variance(values.data(),1,values.size()) *
        ((values.size()-1.0)/values.size());
}

double gmean_or_gstd(const std::vector<double>& x, int ddof) {
    for (double value : x)
        require(value > 0,"geometric standard deviation requires positive observations");
    require(static_cast<size_t>(ddof) < x.size(),"ddof must be less than sample size");
    std::vector<double> transformed;
    transformed.reserve(x.size());
    for (double value:x) transformed.push_back(std::log(value));
    const double variance = popvar(transformed) *
        (x.size()/static_cast<double>(x.size()-ddof));
    require(variance >= 0,"negative logarithmic variance");
    return std::exp(std::sqrt(variance));
}

void geometric_scores(const std::vector<double>& x, int ddof,
                      double* out) {
    for (double value:x)
        require(value > 0,"geometric z-score requires positive observations");
    require(static_cast<size_t>(ddof) < x.size(),"ddof must be less than sample size");
    std::vector<double> transformed;
    transformed.reserve(x.size());
    for (double value:x) transformed.push_back(std::log(value));
    const double center = mean(transformed);
    const double variance = popvar(transformed) *
        (x.size()/static_cast<double>(x.size()-ddof));
    require(variance > 0,"geometric z-score undefined for a constant sample");
    const double scale = std::sqrt(variance);
    for (size_t i=0;i<x.size();++i) out[i] = (transformed[i]-center)/scale;
}

double power_mean(const std::vector<double>& x, double power) {
    require(std::isfinite(power) && power >= -32 && power <= 32,
            "pmean supports finite powers between -32 and 32");
    // SciPy's real pmean returns NaN for negative values even for integer
    // powers. Explicitly reject these as outside the shared real domain.
    for (double v:x)
        require(v>=0,"power mean requires nonnegative real samples");
    if (std::any_of(x.begin(),x.end(),[](double a){return a==0;})) {
        if (power <= 0) return 0.; // SciPy convention for p<=0 with zeros.
    }
    if (power == 0) {
        std::vector<double> logs;
        logs.reserve(x.size());
        for (double v:x) logs.push_back(std::log(v));
        return std::exp(mean(logs));
    }
    if (power == 1) return mean(x);
    // Anchor guarantees each transformed term is within [0,1] for both
    // positive and negative powers; prevents overflow and cancellation.
    const double anchor = power < 0 ?
        *std::min_element(x.begin(),x.end()) :
        *std::max_element(x.begin(),x.end());
    if (anchor == 0.0) return 0.0;
    std::vector<double> normalized;
    normalized.reserve(x.size());
    for (double v:x) normalized.push_back(std::pow(v/anchor,power));
    const double averaged = mean(normalized);
    require(averaged >= 0,"negative transformed power mean");
    return anchor*std::pow(averaged,1.0/power);
}

double central_moment(const std::vector<double>& x, int order) {
    if (order == 0) return 1.;
    if (order == 1 || x.size() == 1) return 0.;
    if (order == 2) {
        if (x.size()==1) return 0.;
        return gsl_stats_tss(x.data(),1,x.size())/x.size();
    }
    const double center = mean(x);
    std::vector<double> powered;
    powered.reserve(x.size());
    for (double v:x) powered.push_back(std::pow(v-center,order));
    return mean(powered);
}
}

extern "C" {
int gsl_stats_aggregates_supported(int library,const char* name) {
    return supported(library,name) ? 1 : 0;
}

int gsl_stats_aggregates_output_count(const char* name,int nx) {
    return output_count(name,nx);
}

const char* gsl_stats_aggregates_error() { return last_error.c_str(); }

int gsl_stats_aggregates_evaluate(int library,const char* name,
                                 const double* x,int nx,
                                 const double* y,int ny,
                                 const double* p,int np,
                                 const double* g,int ng,
                                 double* out,int capacity) {
    last_error.clear();
    try {
        require(supported(library,name),"unsupported library or statistical function");
        const int expected = output_count(name,nx);
        require(expected>0 && capacity>=expected && x && out,
                "invalid sample size or insufficient output capacity");
        require(ny==0 && ng==0 && !y && !g && np==1 && p,
                "expected unweighted x sample, empty y and g, and one scalar parameter");
        for (int i=0;i<nx;++i)
            require(std::isfinite(x[i]),"nonfinite sample observation");
        require(std::isfinite(p[0]),"nonfinite parameter");
        const std::vector<double> values(x,x+nx);
        std::vector<double> results(expected,0.0);
        if (named(name,"gstd")) {
            const int ddof = integer(p[0],0,nx-1,"geometric standard deviation ddof");
            results[0] = gmean_or_gstd(values,ddof);
        } else if (named(name,"gzscore")) {
            const int ddof = integer(p[0],0,nx-1,"geometric z-score ddof");
            geometric_scores(values,ddof,results.data());
        } else if (named(name,"pmean")) {
            results[0] = power_mean(values,p[0]);
        } else {
            const int order = integer(p[0],0,32,"central moment order");
            results[0] = central_moment(values,order);
        }
        for (double v:results)
            require(std::isfinite(v),"nonfinite GSL-derived statistical output");
        std::copy(results.begin(),results.end(),out);
        return 0;
    } catch(const std::exception& exc) {
        last_error=exc.what();
    } catch(...) {
        last_error="unknown GSL aggregate error";
    }
    return -1;
}
}
