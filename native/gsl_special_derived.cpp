// Independently callable public GSL (2.7/2.8) special-function adapters.
// Library index 1 = GSL. Composed results are labeled GSL-DERIVED in the
// documentation; no SciPy code, private GSL symbols or Boost APIs are used.
#include <gsl/gsl_errno.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
thread_local std::string last_error;
bool is(const char* name, const char* expected) {
    return name && std::strcmp(name, expected) == 0;
}
int logical_arity(const char* name) {
    if (is(name,"ai_zeros") || is(name,"bi_zeros") ||
        is(name,"erfinv") || is(name,"erfcinv") ||
        is(name,"erfi") || is(name,"erfcx") ||
        is(name,"log_ndtr")) return 1;
    if (is(name,"gdtria") || is(name,"elliprg")) return 3;
    if (is(name,"ln_multinomial")) return 2; // n + counts vector
    return 0;
}
void require(bool pass, const char* description) {
    if (!pass) throw std::domain_error(description);
}
void check(int status, const char* what) {
    if (status != GSL_SUCCESS)
        throw std::runtime_error(std::string(what)+": "+gsl_strerror(status));
}
double finite(double value) {
    require(std::isfinite(value), "nonfinite GSL-derived result");
    return value;
}
int exact_count(double value, int upper, const char* name) {
    require(std::isfinite(value) && value>=0 && value<=upper &&
            std::floor(value)==value,name);
    return static_cast<int>(value);
}
int output_count(const char* name,const double* args,int count) {
    if (!args || !name || !logical_arity(name)) return 0;
    if (is(name,"ln_multinomial")) return count>=1 && count<=100001 ? 1 : 0;
    if (count!=logical_arity(name)) return 0;
    if (is(name,"ai_zeros") || is(name,"bi_zeros")) {
        if (!std::isfinite(args[0]) || args[0]<1 ||
            args[0]>256 || std::trunc(args[0])!=args[0]) return 0;
        return 4*static_cast<int>(args[0]);
    }
    return 1;
}
double airy(bool ai,bool derivative,double x) {
    gsl_sf_result r{};
    int status=ai ?
        (derivative ? gsl_sf_airy_Ai_deriv_e(x,GSL_PREC_DOUBLE,&r)
                    : gsl_sf_airy_Ai_e(x,GSL_PREC_DOUBLE,&r)) :
        (derivative ? gsl_sf_airy_Bi_deriv_e(x,GSL_PREC_DOUBLE,&r)
                    : gsl_sf_airy_Bi_e(x,GSL_PREC_DOUBLE,&r));
    check(status, "GSL Airy evaluation");
    return finite(r.val);
}
double airy_zero(bool ai,bool derivative,int serial) {
    gsl_sf_result r{};
    int status=ai ?
        (derivative ? gsl_sf_airy_zero_Ai_deriv_e(serial,&r)
                    : gsl_sf_airy_zero_Ai_e(serial,&r)) :
        (derivative ? gsl_sf_airy_zero_Bi_deriv_e(serial,&r)
                    : gsl_sf_airy_zero_Bi_e(serial,&r));
    check(status,"GSL Airy zero");
    return finite(r.val);
}
void airy_zeros(bool ai,int nt,double* output) {
    for (int k=0;k<nt;++k) {
        const double a=airy_zero(ai,false,k+1);
        const double ap=airy_zero(ai,true,k+1);
        require(a<0 && ap<0,"unexpected Airy root outside negative real axis");
        output[k]=a;
        output[nt+k]=ap;
        output[2*nt+k]=airy(ai,false,ap);
        output[3*nt+k]=airy(ai,true,a);
        // Public GSL accuracy is considerably better than this tolerance;
        // avoid claiming roots if an older GSL returns unusable values.
        require(std::abs(airy(ai,false,a))<=1e-7 &&
                std::abs(airy(ai,true,ap))<=1e-7,
                "GSL Airy zero fails function/derivative residual");
        if (k) require(output[k]<output[k-1] && output[nt+k]<output[nt+k-1],
                       "GSL Airy zero ordering mismatch");
    }
}
double scalar(const char* name,const double* a,int count) {
    for (int i=0;i<count;++i) require(std::isfinite(a[i]),"nonfinite special input");
    constexpr double sqrt2=1.4142135623730950488;
    constexpr double sqrtpi=1.7724538509055160273;
    if (is(name,"erfinv")) {
        require(a[0]>-1 && a[0]<1,"erfinv requires -1 < p < 1 for finite output");
        if (a[0]==0) return 0.;
        // Use the smaller normal tail to avoid cancellation near +1.
        const double p=(1-std::abs(a[0]))/2;
        return finite(-std::copysign(1.0,a[0])*
            gsl_cdf_ugaussian_Pinv(p)/sqrt2);
    }
    if (is(name,"erfcinv")) {
        require(a[0]>0 && a[0]<2,"erfcinv requires 0 < p < 2 for finite output");
        if (a[0]==1) return 0.;
        const double p=a[0]<1 ? a[0]/2 : (2-a[0])/2;
        return finite((a[0]<1 ? -1 : 1)*gsl_cdf_ugaussian_Pinv(p)/sqrt2);
    }
    if (is(name,"erfi")) {
        require(std::abs(a[0])<=26,"erfi exceeds stable finite GSL-derived range");
        gsl_sf_result d{};check(gsl_sf_dawson_e(a[0],&d),"GSL Dawson");
        return finite(2./sqrtpi*std::exp(a[0]*a[0])*d.val);
    }
    if (is(name,"erfcx")) {
        require(a[0]>=-26 && a[0]<=26,"erfcx outside bounded finite range");
        gsl_sf_result r{};
        check(gsl_sf_log_erfc_e(a[0],&r),"GSL log(erfc)");
        return finite(std::exp(a[0]*a[0]+r.val));
    }
    if (is(name,"log_ndtr")) {
        gsl_sf_result r{};
        if (a[0]>0) {
            // log(1-Q(x)) avoids the log(2) cancellation for x >> 0.
            check(gsl_sf_erf_Q_e(a[0],&r),"GSL Gaussian survival");
            return finite(std::log1p(-r.val));
        }
        check(gsl_sf_log_erfc_e(-a[0]/sqrt2,&r),"GSL log(erfc)");
        return finite(-std::log(2.)+r.val);
    }
    if (is(name,"gdtria")) {
        // SciPy gdtria(p, shape, x): inverse rate parameter 'a', not x.
        require(a[0]>0 && a[0]<1 && a[1]>0 && a[2]>0,
                "gdtria requires 0<p<1, shape>0, x>0");
        return finite(gsl_cdf_gamma_Pinv(a[0],a[1],1.)/a[2]);
    }
    if (is(name,"elliprg")) {
        // Carlson RG homogeneity; sort to use positive largest z in
        // RG(x,y,z) = (z RF - (x-z)(y-z) RD/3 + sqrt(x*y/z))/2.
        // Scaling arguments to [0,1] bounds over/underflow in intermediates.
        require(a[0]>=0 && a[1]>=0 && a[2]>=0,"elliprg requires nonnegative inputs");
        double x=a[0], y=a[1], z=a[2];
        if (x>y) std::swap(x,y);
        if (y>z) std::swap(y,z);
        if (x>y) std::swap(x,y);
        if (z==0) return 0.;
        if (y==0) return finite(std::sqrt(z)/2);
        const double scale=z;
        x/=scale;y/=scale;z=1.;
        gsl_sf_result rf{},rd{};
        check(gsl_sf_ellint_RF_e(x,y,z,GSL_PREC_DOUBLE,&rf),"GSL Carlson RF");
        check(gsl_sf_ellint_RD_e(x,y,z,GSL_PREC_DOUBLE,&rd),"GSL Carlson RD");
        return finite(0.5*std::sqrt(scale)*
            (z*rf.val-(x-z)*(y-z)*rd.val/3+std::sqrt(x*y/z)));
    }
    if (is(name,"ln_multinomial")) {
        // Native shape is flattened [n, c0, ..., c_(k-1)]; empty
        // count vector is valid exactly for n=0.
        const int n=exact_count(a[0],1000000,"invalid multinomial total");
        long long total=0;
        double log_denominator=0.;
        for (int i=1;i<count;++i) {
            const int k=exact_count(a[i],1000000,"invalid multinomial category count");
            total+=k;
            gsl_sf_result term{};
            check(gsl_sf_lnfact_e(k,&term),"GSL log factorial category");
            log_denominator+=term.val;
        }
        require(total==n,"multinomial counts must sum to n");
        gsl_sf_result numerator{};
        check(gsl_sf_lnfact_e(n,&numerator),"GSL log factorial total");
        return finite(numerator.val-log_denominator);
    }
    throw std::invalid_argument("unknown GSL-derived special scalar");
}
} // namespace

extern "C" {
const char* gsl_special_derived_error() {return last_error.c_str();}
int gsl_special_derived_supported(int library,const char* name) {
    return library==1 ? logical_arity(name) : 0;
}
int gsl_special_derived_output_count(const char* name,const double* args,int count) {
    return output_count(name,args,count);
}
int gsl_special_derived_evaluate(int library,const char* name,const double* args,
                                 int count,double* output,int capacity) {
    last_error.clear();
    try {
        require(gsl_special_derived_supported(library,name)>0,
                "GSL-derived function is unavailable for this library");
        const int needed=output_count(name,args,count);
        require(needed>0 && output && capacity>=needed,
                "invalid GSL-derived output shape or capacity");
        gsl_set_error_handler_off();
        if (is(name,"ai_zeros") || is(name,"bi_zeros"))
            airy_zeros(is(name,"ai_zeros"),static_cast<int>(args[0]),output);
        else output[0]=scalar(name,args,count);
        return 0;
    } catch(const std::exception& e) {last_error=e.what();}
      catch(...) {last_error="unknown GSL-derived special failure";}
    return -1;
}
}
