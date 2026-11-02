// More saved Sci 0.3 reference oracles derived from public GSL primitives.
// Provider index 1 means GSL; the C ABI is independent of other native files.
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_statistics_double.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
thread_local std::string last_error;
constexpr double pi = 3.141592653589793238462643383279502884;

bool is(const char* name, const char* label) {
    return name && std::strcmp(name, label) == 0;
}
void require(bool pass, const char* message) {
    if (!pass) throw std::invalid_argument(message);
}
void check(int status) {
    if (status != GSL_SUCCESS) throw std::runtime_error(gsl_strerror(status));
}
double checked_finite(double v) {
    require(std::isfinite(v), "nonfinite GSL-derived result");
    return v;
}
int integer(double value, int lower, int upper) {
    require(std::isfinite(value) && value >= lower && value <= upper &&
            std::floor(value) == value, "invalid integral parameter");
    return static_cast<int>(value);
}
bool stats(const char* name) {
    return is(name,"stats.entropy") || is(name,"stats.rankdata") ||
           is(name,"stats.tiecorrect") || is(name,"stats.mode") ||
           is(name,"stats.wasserstein_distance") ||
           is(name,"stats.energy_distance") || is(name,"stats.boxcox_llf") ||
           is(name,"stats.zmap") || is(name,"stats.expectile");
}
bool special(const char* name) {
    return is(name,"special.exp2") || is(name,"special.exp10") ||
           is(name,"special.cosm1") || is(name,"special.sindg") ||
           is(name,"special.cosdg") || is(name,"special.tandg") ||
           is(name,"special.cotdg") || is(name,"special.entr") ||
           is(name,"special.xlogy") || is(name,"special.xlog1py") ||
           is(name,"special.rel_entr") || is(name,"special.kl_div");
}
int output_count(int provider,const char* name,int nx,int ny,int np) {
    if (provider != 1 || !name || nx < 0 || nx > 100000 ||
        ny < 0 || ny > 100000 || np < 0 || np > 100000) return 0;
    if (special(name)) {
        if (ny || (is(name,"special.xlogy") || is(name,"special.xlog1py") ||
                   is(name,"special.rel_entr") || is(name,"special.kl_div")
                   ? (nx != 2 || np) : (nx != 1 || np))) return 0;
        return 1;
    }
    if (!stats(name)) return 0;
    if (is(name,"stats.mode")) return nx && !ny && !np ? 2 : 0;
    if (is(name,"stats.rankdata")) return nx && !ny && np==1 ? nx : 0;
    if (is(name,"stats.zmap")) return nx && ny>=2 && np==1 ? nx : 0;
    if (is(name,"stats.entropy")) return nx && (ny==0 || ny==nx) &&
        np==1 ? 1 : 0;
    if (is(name,"stats.wasserstein_distance") ||
        is(name,"stats.energy_distance")) return nx && ny && !np ? 1 : 0;
    if (is(name,"stats.boxcox_llf")) return nx>=2 && !ny && np==1 ? 1 : 0;
    if (is(name,"stats.expectile")) return nx && !ny && np==1 ? 1 : 0;
    if (is(name,"stats.tiecorrect")) return nx && !ny && !np ? 1 : 0;
    return 0;
}
std::vector<double> sorted(const double* x, int n) {
    std::vector<double> values(x,x+n);
    gsl_sort(values.data(),1,values.size());
    return values;
}
double log_gsl(double value) {
    gsl_sf_result result{};
    check(gsl_sf_log_e(value,&result));
    return result.val;
}
double log1p_gsl(double value) {
    gsl_sf_result result{};
    check(gsl_sf_log_1plusx_e(value,&result));
    return result.val;
}
double exp_gsl(double value) {
    gsl_sf_result result{};
    check(gsl_sf_exp_e(value,&result));
    return result.val;
}
double trig_gsl(bool cosine,double value) {
    gsl_sf_result result{};
    check(cosine ? gsl_sf_cos_e(value,&result) : gsl_sf_sin_e(value,&result));
    return result.val;
}
void scalar(const char* name,const double* x,double* out) {
    const double a=x[0];
    if (is(name,"special.exp2")) out[0]=exp_gsl(a*std::log(2.));
    else if (is(name,"special.exp10")) out[0]=exp_gsl(a*std::log(10.));
    else if (is(name,"special.cosm1")) {
        const double s=trig_gsl(false,a/2.);
        out[0]=-2.*s*s;
    } else if (is(name,"special.sindg")) out[0]=trig_gsl(false,a*pi/180.);
    else if (is(name,"special.cosdg")) out[0]=trig_gsl(true,a*pi/180.);
    else if (is(name,"special.tandg") || is(name,"special.cotdg")) {
        const double rad=a*pi/180.;
        const double s=trig_gsl(false,rad),c=trig_gsl(true,rad);
        require(is(name,"special.cotdg") ? s!=0. : c!=0.,
                "nonfinite trigonometric ratio");
        out[0]=is(name,"special.cotdg") ? c/s : s/c;
    } else if (is(name,"special.entr")) {
        require(a>=0.,"entr requires a nonnegative input");
        out[0]=a==0. ? 0. : -a*log_gsl(a);
    } else if (is(name,"special.xlogy") || is(name,"special.xlog1py")) {
        if (a==0.) out[0]=0.;
        else out[0]=a*(is(name,"special.xlogy")
                            ? log_gsl(x[1]):log1p_gsl(x[1]));
    } else if (is(name,"special.rel_entr") || is(name,"special.kl_div")) {
        const double b=x[1];
        require(a>=0. && b>=0.,"relative entropy requires nonnegative inputs");
        if (a==0.) out[0]=is(name,"special.kl_div") ? b : 0.;
        else {
            require(b>0.,"positive x requires positive y");
            // The logarithm of x/y loses accuracy near unity. GSL log1p
            // avoids this cancellation on the saved near-equal profiles.
            const double t=(a-b)/b;
            const double term=a*log1p_gsl(t);
            out[0]=is(name,"special.kl_div") ? term-a+b : term;
        }
    } else throw std::invalid_argument("unsupported GSL special scalar");
}
void ranks(const double* x,int n,int method,double* out) {
    std::vector<size_t> index(n);
    gsl_sort_index(index.data(),x,1,n);
    if (method==4) {
        // GSL correctly orders values, but makes no promise about the
        // original order of equal values. Ordinal ranks require that order.
        for (int start=0;start<n;) {
            int end=start+1;
            while(end<n && x[index[end]]==x[index[start]])++end;
            std::sort(index.begin()+start,index.begin()+end);
            start=end;
        }
        for (int i=0;i<n;++i) out[index[i]]=i+1.;
        return;
    }
    int dense=0;
    for (int start=0;start<n;) {
        int end=start+1;
        while (end<n && x[index[end]]==x[index[start]]) ++end;
        ++dense;
        const double rank=method==0 ? (start+end+1.)/2. :
                          method==1 ? start+1. : method==2 ? end : dense;
        for (int j=start;j<end;++j) out[index[j]]=rank;
        start=end;
    }
}
double tiecorrect(const double* x,int n) {
    if (n<2) return 1.;
    const auto v=sorted(x,n);
    double excess=0.;
    for (int start=0;start<n;) {
        int end=start+1;
        while(end<n && v[end]==v[start])++end;
        const double t=end-start;
        excess+=t*t*t-t;
        start=end;
    }
    const double nn=n;
    return 1.-excess/(nn*nn*nn-nn);
}
double entropy(const double* x,int n,const double* y,int ny,double base) {
    require(base==0. || (base>0. && base!=1.), "entropy base must be 0, positive and not 1");
    const double sx=gsl_stats_mean(x,1,n)*n;
    require(sx>0. && std::isfinite(sx),"entropy requires positive total mass");
    double sy=1.;
    if (ny) {
        sy=gsl_stats_mean(y,1,ny)*ny;
        require(sy>0. && std::isfinite(sy),"cross entropy requires positive comparison mass");
    }
    double total=0.;
    for (int i=0;i<n;++i) {
        require(x[i]>=0. && (!ny || y[i]>=0.),"entropy requires nonnegative masses");
        if (x[i]==0.) continue;
        const double px=x[i]/sx;
        if (!ny) total-=px*log_gsl(px);
        else {
            const double qx=y[i]/sy;
            require(qx>0., "positive probability compared with zero");
            total+=px*log_gsl(px/qx);
        }
    }
    return base==0. ? total : total/log_gsl(base);
}
double distance(const double* x,int nx,const double* y,int ny,bool energy) {
    const auto sx=sorted(x,nx),sy=sorted(y,ny);
    int ix=0,iy=0;
    double last=std::min(sx.front(),sy.front()),integral=0.;
    // GSL supplies the sorted empirical samples; this is the exact
    // integral of |F-G| (Wasserstein) or 2*(F-G)^2 (energy squared).
    while(ix<nx || iy<ny) {
        const double next=std::min(ix<nx ? sx[ix] : std::numeric_limits<double>::infinity(),
                                   iy<ny ? sy[iy] : std::numeric_limits<double>::infinity());
        const double delta=next-last;
        const double d=static_cast<double>(ix)/nx-static_cast<double>(iy)/ny;
        integral+=delta*(energy ? d*d : std::abs(d));
        while(ix<nx && sx[ix]==next)++ix;
        while(iy<ny && sy[iy]==next)++iy;
        last=next;
    }
    return energy ? std::sqrt(2.*integral) : integral;
}
double boxcox_llf(const double* x,int n,double lambda) {
    std::vector<double> transformed(n);
    double sumlog=0.;
    for (int i=0;i<n;++i) {
        require(x[i]>0.,"Box-Cox log likelihood requires positive samples");
        const double lx=log_gsl(x[i]);
        sumlog+=lx;
        transformed[i]=lambda==0. ? lx : std::expm1(lambda*lx)/lambda;
    }
    const double variation=gsl_stats_tss(transformed.data(),1,n)/n;
    require(variation>0. && std::isfinite(variation),"Box-Cox requires nonconstant sample");
    return (lambda-1.)*sumlog-.5*n*log_gsl(variation);
}
double expectile(const double* x,int n,double alpha) {
    require(alpha>=0. && alpha<=1.,"expectile alpha must be in [0,1]");
    const auto range=std::minmax_element(x,x+n);
    if (alpha==0.) return *range.first;
    if (alpha==1.) return *range.second;
    if (*range.first==*range.second) return *range.first;
    double lo=*range.first,hi=*range.second;
    std::vector<double> weights(n);
    for(int k=0;k<90;++k) {
        const double estimate=lo+(hi-lo)/2.;
        for(int i=0;i<n;++i) weights[i]=x[i]<estimate ? 1.-alpha : alpha;
        const double expected=gsl_stats_wmean(weights.data(),1,x,1,n);
        if (expected > estimate) lo=estimate;
        else hi=estimate;
    }
    return lo+(hi-lo)/2.;
}
void evaluate_stats(const char* name,const double* x,int nx,
                    const double* y,int ny,const double* p,double* out) {
    if (is(name,"stats.rankdata")) {
        ranks(x,nx,integer(p[0],0,4),out);
    } else if (is(name,"stats.tiecorrect")) out[0]=tiecorrect(x,nx);
    else if (is(name,"stats.mode")) {
        auto sorted_values=sorted(x,nx);
        int best_count=0;
        for(int i=0;i<nx;) {
            int j=i+1;while(j<nx && sorted_values[j]==sorted_values[i])++j;
            if (j-i>best_count) {out[0]=sorted_values[i];best_count=j-i;}
            i=j;
        }
        out[1]=best_count;
    } else if (is(name,"stats.entropy")) out[0]=entropy(x,nx,y,ny,p[0]);
    else if (is(name,"stats.wasserstein_distance"))
        out[0]=distance(x,nx,y,ny,false);
    else if (is(name,"stats.energy_distance"))
        out[0]=distance(x,nx,y,ny,true);
    else if (is(name,"stats.boxcox_llf")) out[0]=boxcox_llf(x,nx,p[0]);
    else if (is(name,"stats.expectile")) out[0]=expectile(x,nx,p[0]);
    else if (is(name,"stats.zmap")) {
        const int ddof=integer(p[0],0,ny-1);
        const double mean=gsl_stats_mean(y,1,ny);
        const double variance=gsl_stats_tss_m(y,1,ny,mean)/(ny-ddof);
        require(variance>0.,"zmap requires nonconstant comparison samples");
        const double scale=std::sqrt(variance);
        for(int i=0;i<nx;++i)out[i]=(x[i]-mean)/scale;
    } else throw std::invalid_argument("unsupported GSL stats operation");
}
}

extern "C" {
const char* extra_native_stats_error() {return last_error.c_str();}
int extra_native_stats_supported(int provider,const char* name) {
    return provider==1 && (special(name) || stats(name)) ? 1 : 0;
}
int extra_native_stats_output_count(int provider,const char* name,int nx,int ny,int np) {
    return output_count(provider,name,nx,ny,np);
}
int extra_native_stats_evaluate(int provider,const char* name,
        const double* x,int nx,const double* y,int ny,const double* p,int np,
        double* output,int capacity) {
    last_error.clear();
    try {
        const int needed=output_count(provider,name,nx,ny,np);
        require(needed>0 && capacity>=needed && output && (nx==0 || x) &&
                (ny==0 || y) && (np==0 || p), "invalid GSL output shape/buffers");
        for(int i=0;i<nx;++i) require(std::isfinite(x[i]),"nonfinite primary input");
        for(int i=0;i<ny;++i) require(std::isfinite(y[i]),"nonfinite secondary input");
        for(int i=0;i<np;++i) require(std::isfinite(p[i]),"nonfinite parameters");
        gsl_set_error_handler_off();
        if (special(name)) scalar(name,x,output);
        else evaluate_stats(name,x,nx,y,ny,p,output);
        for(int i=0;i<needed;++i) checked_finite(output[i]);
        return 0;
    } catch(const std::exception& e) {last_error=e.what();}
      catch(...) {last_error="unknown GSL stats/special failure";}
    return -1;
}
}
