// Public GSL interpolation competitors for the saved Sci scalar/regular-grid
// cases. GSL polynomial interpolation applies ONLY to distinct abscissae:
// confluent (Hermite) Krogh nodes are not represented by this GSL API.
//
// Encoding for barycentric/krogh/PCHIP:
//   [n, nq, xa[n], ya[n], queries[nq]]
// Encoding for interpn:
//   [ndim, lengths[ndim], axes[sum(lengths)], values[product(lengths)],
//    nq, points[nq*ndim], method(0=linear,1=nearest),
//    bounds_error(0/1), has_fill(0/1), fill]
// For all operations, results are Sci's flattened values, no extra metadata.
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
thread_local std::string last_error;
constexpr int max_n = 64;
constexpr int max_queries = 4096;
constexpr int max_grid_values = 32768;
bool same(const char* a, const char* b) { return a && std::strcmp(a,b)==0; }

bool supports(int library, const char* name) {
    return library==1 && (same(name,"interpolate.barycentric_interpolate") ||
            same(name,"interpolate.krogh_interpolate") ||
            same(name,"interpolate.pchip_interpolate") ||
            same(name,"interpolate.interpn"));
}

struct Reader {
    const double* data;
    int size, index=0;
    Reader(const double* d,int n):data(d),size(n) {
        if (!d || n<0 || n>200000)
            throw std::invalid_argument("invalid GSL interpolation input buffer");
    }
    double real() {
        if (index>=size) throw std::invalid_argument("truncated interpolation input");
        const double x=data[index++];
        if (!std::isfinite(x)) throw std::invalid_argument("nonfinite interpolation input");
        return x;
    }
    int integer(int upper) {
        double v=real();
        if (v<0 || v>upper || v!=std::trunc(v))
            throw std::invalid_argument("invalid interpolation count or option");
        return static_cast<int>(v);
    }
    std::vector<double> values(int n) {
        std::vector<double> a;
        a.reserve(n);
        for(int i=0;i<n;++i) a.push_back(real());
        return a;
    }
    void done() {
        if(index!=size) throw std::invalid_argument("unexpected extra interpolation input");
    }
};

void check_order(const std::vector<double>& x) {
    for(size_t i=1;i<x.size();++i)
        if (!(x[i]>x[i-1]))
            throw std::invalid_argument("interpolation nodes must be strictly ascending");
}

struct SplineDeleter {
    void operator()(gsl_spline* a) const { if(a) gsl_spline_free(a); }
};
struct AccelDeleter {
    void operator()(gsl_interp_accel* a) const { if(a) gsl_interp_accel_free(a); }
};

std::vector<double> interpolate_polynomial(const char* name, Reader& r) {
    const int n=r.integer(max_n), nq=r.integer(max_queries);
    if(n<2 || nq<1) throw std::invalid_argument("invalid polynomial node/query size");
    const auto xa=r.values(n),ya=r.values(n),q=r.values(nq);
    r.done();
    check_order(xa);
    const bool steffen=same(name,"interpolate.pchip_interpolate");
    if(steffen) {
        // Steffen's monotone slope selection differs from SciPy's PCHIP.
        // Only allow exact query knots or affine data, where the two agree.
        if(n<3) throw std::invalid_argument("Steffen requires at least 3 knots");
        const double slope=(ya[1]-ya[0])/(xa[1]-xa[0]);
        bool affine=true;
        for(int i=2;i<n;++i) {
            const double next=(ya[i]-ya[i-1])/(xa[i]-xa[i-1]);
            if (next!=slope) affine=false;
        }
        if(!affine) {
            for(double x:q)
                if(!std::binary_search(xa.begin(),xa.end(),x))
                    throw std::invalid_argument("GSL Steffen is not SciPy PCHIP");
        }
    }
    const gsl_interp_type* type=steffen?gsl_interp_steffen:gsl_interp_polynomial;
    std::unique_ptr<gsl_spline,SplineDeleter> spline(gsl_spline_alloc(type,n));
    std::unique_ptr<gsl_interp_accel,AccelDeleter> accel(gsl_interp_accel_alloc());
    if(!spline || !accel || gsl_spline_init(spline.get(),xa.data(),ya.data(),n)!=GSL_SUCCESS)
        throw std::runtime_error("GSL spline initialization failed");
    std::vector<double> result;
    result.reserve(nq);
    for(double x:q) {
        double value=0;
        int status=GSL_SUCCESS;
        if(x>=xa.front() && x<=xa.back()) {
            status=gsl_spline_eval_e(spline.get(),x,accel.get(),&value);
        } else if(!steffen) {
            // GSL's public polynomial interpolation type is global; its
            // public eval callback admits extrapolation. gsl_spline_eval_e
            // enforces a *different* (bounded) spline domain, so dispatch
            // directly to this documented public interpolation type callback.
            status=spline->interp->type->eval(
                spline->interp->state,xa.data(),ya.data(),n,x,accel.get(),&value);
        } else {
            throw std::invalid_argument("Steffen extrapolation unsupported");
        }
        if(status!=GSL_SUCCESS || !std::isfinite(value))
            throw std::runtime_error("GSL polynomial/Steffen evaluation failed");
        result.push_back(value);
    }
    return result;
}

std::vector<double> interpolate_grid(Reader& r) {
    const int ndim=r.integer(4);
    if(ndim<1) throw std::invalid_argument("regular grid requires >=1 axis");
    std::vector<int> lens(ndim);
    int total=1;
    for(int d=0;d<ndim;++d) {
        lens[d]=r.integer(max_n);
        if(lens[d]<1 || total>max_grid_values/lens[d])
            throw std::invalid_argument("invalid regular grid axis length");
        total*=lens[d];
    }
    std::vector<std::vector<double>> axes;
    std::vector<bool> descending;
    for(int d=0;d<ndim;++d) {
        auto axis=r.values(lens[d]);
        if(lens[d]>1) {
            bool asc=axis[1]>axis[0],desc=axis[1]<axis[0];
            if(!asc && !desc)
                throw std::invalid_argument("repeated regular grid coordinate");
            if(desc) std::reverse(axis.begin(),axis.end());
            check_order(axis);
            descending.push_back(desc);
        } else descending.push_back(false);
        axes.push_back(std::move(axis));
    }
    auto data=r.values(total);
    const int nq=r.integer(max_queries);
    if(nq<1) throw std::invalid_argument("empty regular grid query set");
    std::vector<double> queries=r.values(nq*ndim);
    const int method=r.integer(1), bounds=r.integer(1), has_fill=r.integer(1);
    const double fill=r.real();
    r.done();
    std::vector<int> strides(ndim,1);
    for(int d=ndim-2;d>=0;--d) strides[d]=strides[d+1]*lens[d+1];
    std::vector<double> result;
    result.reserve(nq);
    for(int k=0;k<nq;++k) {
        bool outside=false;
        std::vector<int> indices(ndim);
        std::vector<double> weights(ndim);
        for(int d=0;d<ndim;++d) {
            const auto& ax=axes[d];
            const double x=queries[k*ndim+d];
            if(x<ax.front() || x>ax.back()) outside=true;
            if(lens[d]==1) {
                indices[d]=0;
                weights[d]=0;
            } else {
                // This is the public GSL bracketing routine. On the two
                // exteriors it returns the endpoint interval, making
                // SciPy-compatible linear extrapolation possible.
                const size_t index=gsl_interp_bsearch(ax.data(),x,0,ax.size()-1);
                indices[d]=static_cast<int>(index);
                weights[d]=(x-ax[index])/(ax[index+1]-ax[index]);
            }
        }
        if(outside && bounds)
            throw std::invalid_argument("regular grid query outside bounds_error=True");
        if(outside && has_fill) {
            result.push_back(fill);
            continue;
        }
        if(method==1) {
            int index=0;
            for(int d=0;d<ndim;++d) {
                int offset=indices[d]+(weights[d]>0.5 && lens[d]>1?1:0);
                if(descending[d]) offset=lens[d]-1-offset;
                index+=offset*strides[d];
            }
            result.push_back(data[index]);
        } else {
            // Multilinear tensor interpolation with the GSL interval
            // bracketing. Does not confuse GSL 1D spline with Sci interpn.
            double y=0;
            for(int corners=0;corners<(1<<ndim);++corners) {
                double weight=1;
                int index=0;
                for(int d=0;d<ndim;++d) {
                    const bool hi=((corners>>d)&1) && lens[d]>1;
                    if(lens[d]==1 && ((corners>>d)&1)) {weight=0;break;}
                    weight*=hi?weights[d]:(1-weights[d]);
                    int offset=indices[d]+(hi?1:0);
                    if(descending[d]) offset=lens[d]-1-offset;
                    index+=offset*strides[d];
                }
                if(weight!=0) y+=weight*data[index];
            }
            if(!std::isfinite(y))
                throw std::runtime_error("nonfinite regular grid result");
            result.push_back(y);
        }
    }
    return result;
}

std::vector<double> compute(const char* name, const double* data,int length) {
    if(!supports(1,name))
        throw std::invalid_argument("unsupported GSL interpolation operation");
    gsl_set_error_handler_off();
    Reader reader(data,length);
    return same(name,"interpolate.interpn")?
        interpolate_grid(reader):interpolate_polynomial(name,reader);
}
} // namespace

extern "C" {
int gsl_interpolation_supported(int library,const char* name) {
    return supports(library,name)?2:0;
}
int gsl_interpolation_output_count(const char* name,const double* data,int count) {
    try {return static_cast<int>(compute(name,data,count).size());}
    catch (...) {return 0;}
}
int gsl_interpolation_evaluate(int library,const char* name,const double* data,
                                int count,double* output,int capacity) {
    last_error.clear();
    try {
        if(!supports(library,name))
            throw std::invalid_argument("unsupported GSL interpolation library/function");
        auto result=compute(name,data,count);
        if(!output || capacity<static_cast<int>(result.size()))
            throw std::invalid_argument("insufficient GSL interpolation output capacity");
        std::copy(result.begin(),result.end(),output);
        return 0;
    } catch(const std::exception& e) {last_error=e.what();}
      catch(...) {last_error="unknown GSL interpolation failure";}
    return -1;
}
const char* gsl_interpolation_error() {return last_error.c_str();}
}
