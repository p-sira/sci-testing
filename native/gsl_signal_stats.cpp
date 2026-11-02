// GSL-derived SciPy signal detrend, median and Wiener filters.
//
// GSL has no native "scipy.signal.detrend", "medfilt", "medfilt2d" or
// "wiener" functions. This translation unit composes PUBLIC GSL 2.7+
// linear-fit, statistics and sorting APIs, preserving saved SciPy output
// ordering and zero-padding semantics. Provider GSL only (library=1).
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_sort_double.h>
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
constexpr int max_rows=1024;
constexpr int max_cols=8192;
constexpr int max_points=65536;
constexpr int max_kernel=127;

bool named(const char* name,const char* target) {
    return name && std::strcmp(name,target)==0;
}
bool supported(const char* name) {
    return named(name,"detrend") || named(name,"medfilt") ||
           named(name,"medfilt2d") || named(name,"wiener");
}
void require(bool ok,const char* message) {
    if (!ok)throw std::invalid_argument(message);
}
int output_count(const char* name,int rows,int cols) {
    if (!supported(name) || rows<=0 || rows>max_rows ||
        cols<=0 || cols>max_cols ||
        static_cast<long long>(rows)*cols>max_points ||
        (!named(name,"medfilt2d") && rows!=1))return 0;
    return rows*cols;
}
void detrend(const double* data,int n,int linear,std::vector<double>& result) {
    require(linear==0 || linear==1,"detrend mode must be constant(0) or linear(1)");
    // A one-point constant or linear trend is the point itself. A linear
    // fit with one observation is singular; SciPy returns [0] instead.
    if (n==1) {
        result[0]=0.;
        return;
    }
    if (!linear) {
        const double mean=gsl_stats_mean(data,1,static_cast<size_t>(n));
        require(std::isfinite(mean),"GSL detrend mean failed");
        for (int i=0;i<n;++i) result[i]=data[i]-mean;
        return;
    }
    std::vector<double> positions(n);
    for (int i=0;i<n;++i)
        positions[i]=static_cast<double>(i+1)/n;
    double c0=0.,c1=0.,cov00=0.,cov01=0.,cov11=0.,sum_sq=0.;
    const int status=gsl_fit_linear(
        positions.data(),1,data,1,static_cast<size_t>(n),
        &c0,&c1,&cov00,&cov01,&cov11,&sum_sq);
    require(status==GSL_SUCCESS && std::isfinite(c0) && std::isfinite(c1),
            "GSL linear detrend fit failed");
    for (int i=0;i<n;++i)
        result[i]=data[i]-(c0+c1*positions[i]);
}

void validate_kernel(int height,int width) {
    require(height>0 && width>0 && height<=max_kernel &&
            width<=max_kernel && (height&1) && (width&1),
            "median/Wiener kernel dimensions must be positive bounded odd integers");
}
std::vector<double> window_at(const double* data,int rows,int cols,
                              int center_row,int center_col,int height,int width) {
    std::vector<double> window;
    window.reserve(static_cast<size_t>(height)*width);
    const int rh=height/2;
    const int cw=width/2;
    for (int dr=-rh;dr<=rh;++dr)
        for (int dc=-cw;dc<=cw;++dc) {
            const int r=center_row+dr;
            const int c=center_col+dc;
            window.push_back(r<0 || c<0 || r>=rows || c>=cols
                                 ? 0. : data[r*cols+c]);
        }
    return window;
}
void median_filter(const double* data,int rows,int cols,
                   int height,int width,std::vector<double>& result) {
    validate_kernel(height,width);
    // Zero pad, including for windows larger than the input dimension.
    for (int r=0;r<rows;++r)
        for (int c=0;c<cols;++c) {
            auto window=window_at(data,rows,cols,r,c,height,width);
            gsl_sort(window.data(),1,window.size());
            result[r*cols+c]=gsl_stats_median_from_sorted_data(
                window.data(),1,window.size());
        }
}
void wiener(const double* data,int n,int width,int has_noise,
            double noise,std::vector<double>& result) {
    require(has_noise==0 || has_noise==1,
            "Wiener noise mode must be estimated(0) or specified(1)");
    validate_kernel(1,width);
    require(!has_noise || (std::isfinite(noise) && noise>=0.),
            "Wiener noise must be finite and nonnegative");
    std::vector<double> means(n),variances(n);
    for (int i=0;i<n;++i) {
        auto window=window_at(data,1,n,0,i,1,width);
        means[i]=gsl_stats_mean(window.data(),1,window.size());
        // Public GSL variance_with_fixed_mean uses population denominator N,
        // as in SciPy's local E[x^2]-E[x]^2 with zeros outside the signal.
        variances[i]=gsl_stats_variance_with_fixed_mean(
            window.data(),1,window.size(),means[i]);
    }
    if (!has_noise) {
        noise=gsl_stats_mean(variances.data(),1,variances.size());
    }
    require(std::isfinite(noise),"GSL Wiener estimated noise is nonfinite");
    for (int i=0;i<n;++i) {
        const double variance=variances[i];
        const double mean=means[i];
        // Match SciPy's out=where(lVar<noise,localMean,original - ...)
        // (exactly variance==noise and noise==0 yield original/NaN based on
        // floating arithmetic; such degenerate results are rejected below).
        result[i]=variance<noise ? mean :
                  mean+(data[i]-mean)*(1.-noise/variance);
    }
}
}

extern "C" {
const char* gsl_signal_stats_error() {return last_error.c_str();}
int gsl_signal_stats_supported(int library,const char* name) {
    return library==1 && supported(name) ? 1 : 0;
}
int gsl_signal_stats_output_count(const char* name,int rows,int cols) {
    return output_count(name,rows,cols);
}
int gsl_signal_stats_evaluate(int library,const char* name,
                              const double* data,int rows,int cols,
                              int kernel_rows,int kernel_cols,
                              int mode,double noise,double* out,int capacity) {
    last_error.clear();
    try {
        require(gsl_signal_stats_supported(library,name)==1,
                "unsupported GSL signal function or provider");
        const int count=output_count(name,rows,cols);
        require(count>0 && data && out && capacity>=count,
                "invalid GSL signal dimensions, data or output capacity");
        for (int i=0;i<count;++i)
            require(std::isfinite(data[i]),"GSL signal requires finite data");
        std::vector<double> result(count);
        if (named(name,"detrend")) {
            require(kernel_rows==1 && kernel_cols==1,
                    "detrend does not accept kernel dimensions");
            detrend(data,cols,mode,result);
        } else if (named(name,"medfilt") || named(name,"medfilt2d")) {
            require(mode==0,"median filter does not accept a mode");
            require(named(name,"medfilt2d") || kernel_rows==1,
                    "1D median kernel must have one row");
            median_filter(data,rows,cols,kernel_rows,kernel_cols,result);
        } else if (named(name,"wiener")) {
            require(kernel_rows==1,"1D Wiener kernel must have one row");
            wiener(data,cols,kernel_cols,mode,noise,result);
        }
        for (const auto value:result)
            require(std::isfinite(value),"nonfinite GSL-derived signal result");
        std::copy(result.begin(),result.end(),out);
        return 0;
    } catch (const std::exception& error) {last_error=error.what();}
      catch (...) {last_error="unknown GSL signal statistics failure";}
    return -1;
}
}
