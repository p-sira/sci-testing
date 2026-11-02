// Public Boost.Math 1.83-compatible vector distance primitives.
// Matrix/condensed results use SciPy row-major or upper-triangle traversal.
// Weighted cases pre-transform inputs and then call Boost norms.
// Library=0 (Boost.Math) only. No unsupported metric fallback.
#include <boost/math/tools/norms.hpp>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
thread_local std::string last_error;
constexpr int MAX_ROWS=128, MAX_DIM=512, MAX_POWER=32;
bool is(const char* actual,const char* expected) {
    return actual && std::strcmp(actual,expected)==0;
}
void require(bool ok,const char* message) {
    if (!ok)throw std::domain_error(message);
}
bool known(const char* name) {
    return is(name,"euclidean") || is(name,"sqeuclidean") ||
           is(name,"cityblock") || is(name,"chebyshev") ||
           is(name,"distance_matrix") || is(name,"pdist") ||
           is(name,"cdist") || is(name,"minkowski") ||
           is(name,"minkowski_distance") || is(name,"minkowski_distance_p") ||
           is(name,"hamming");
}
bool paired(const char* name) {
    return is(name,"distance_matrix") || is(name,"pdist") || is(name,"cdist");
}
int needed(const char* name,int urows,int vrows,int dim) {
    if (!known(name) || urows<1 || urows>MAX_ROWS ||
        vrows<0 || vrows>MAX_ROWS || dim<0 || dim>MAX_DIM)return 0;
    if (is(name,"pdist"))return vrows==0 && dim>0 ? urows*(urows-1)/2 : 0;
    if (is(name,"cdist") || is(name,"distance_matrix"))
        return vrows>=1 && dim>0 ? urows*vrows : 0;
    return urows==1 && vrows==1 ? 1 : 0;
}
enum class Metric {euclidean,sqeuclidean,cityblock,chebyshev,
                   minkowski,minkowski_p,hamming};
Metric which(const char* name,const char* metric) {
    const char* check=paired(name) ? metric : name;
    if (is(check,"euclidean") || is(check,"distance_matrix"))return Metric::euclidean;
    if (is(check,"sqeuclidean"))return Metric::sqeuclidean;
    if (is(check,"cityblock"))return Metric::cityblock;
    if (is(check,"chebyshev"))return Metric::chebyshev;
    if (is(check,"minkowski") || is(check,"minkowski_distance"))return Metric::minkowski;
    if (is(check,"minkowski_distance_p"))return Metric::minkowski_p;
    if (is(check,"hamming"))return Metric::hamming;
    throw std::domain_error("metric has no public Boost.Math norm equivalent");
}
double calculate(Metric metric,const double* u,const double* v,
                 int dim,const double* w,int wcount,double p) {
    using boost::math::tools::l1_distance;
    using boost::math::tools::l2_distance;
    using boost::math::tools::sup_distance;
    using boost::math::tools::lp_distance;
    using boost::math::tools::hamming_distance;
    require(dim>=0 && dim<=MAX_DIM, "invalid spatial feature dimension");
    require(wcount==0 || wcount==dim, "invalid spatial weight dimension");
    std::vector<double> a,b;
    a.reserve(dim);b.reserve(dim);
    double sum_w=0.;
    if (metric==Metric::minkowski || metric==Metric::minkowski_p)
        require((std::isfinite(p) && p>=1. && p<=MAX_POWER &&
                 p==std::floor(p)) || (std::isinf(p) && p>0.),
                "Boost lp_distance requires integer p>=1, or +infinity");
    for (int i=0;i<dim;++i) {
        require(std::isfinite(u[i]) && std::isfinite(v[i]),
                "spatial coordinates must be finite");
        const double wi=wcount ? w[i] : 1.;
        require(std::isfinite(wi) && wi>=0.,
                "spatial weights must be finite and nonnegative");
        if (metric==Metric::hamming) {
            a.push_back(u[i]==v[i] ? 0. : wi);
            b.push_back(0.);sum_w+=wi;
            continue;
        }
        if ((metric==Metric::chebyshev ||
            ((metric==Metric::minkowski || metric==Metric::minkowski_p) &&
             std::isinf(p))) && wi==0.)continue;
        if (wcount) {
            const double factor =
                metric==Metric::sqeuclidean ? std::sqrt(wi) :
                metric==Metric::cityblock ? wi :
                (metric==Metric::minkowski || metric==Metric::minkowski_p) ?
                   std::pow(wi,1./p) : 1.;
            if (metric==Metric::chebyshev) {
                a.push_back(u[i]); b.push_back(v[i]);
            } else if (metric==Metric::euclidean) {
                throw std::domain_error("weighted Euclidean not mapped");
            } else {
                a.push_back((u[i]-v[i])*factor);
                b.push_back(0.);
            }
        } else {
            a.push_back(u[i]);b.push_back(v[i]);
        }
    }
    double result=0.;
    switch(metric) {
        case Metric::euclidean:result=l2_distance(a,b);break;
        case Metric::sqeuclidean: {
            const double d=l2_distance(a,b);result=d*d;break;
        }
        case Metric::cityblock:result=l1_distance(a,b);break;
        case Metric::chebyshev:
            result=a.empty() ? 0. : sup_distance(a,b);break;
        case Metric::minkowski:
        case Metric::minkowski_p:
            if (a.empty())result=0.;
            else if (std::isinf(p))result=sup_distance(a,b);
            else result=lp_distance(a,b,static_cast<unsigned>(p));
            if (metric==Metric::minkowski_p && std::isfinite(p))
                result=std::pow(result,p);
            break;
        case Metric::hamming:
            require(dim==0 || sum_w>0.,"Hamming weights must have positive total");
            if (dim==0)result=std::numeric_limits<double>::quiet_NaN();
            else if (wcount)result=l1_distance(a,b)/sum_w;
            else result=double(hamming_distance(u,u+dim,v))/dim;
            break;
    }
    require(std::isfinite(result),"Boost.Math spatial metric returned nonfinite output");
    return result;
}
std::vector<double> compute(const char* name,const double* u,int urows,
                            const double* v,int vrows,int dim,
                            const double* w,int wcount,const char* metric,double p) {
    const int count=needed(name,urows,vrows,dim);
    require(count>0,"invalid spatial output shape or unsupported function");
    require(u && (v || is(name,"pdist")) &&
            (wcount==0 || w) && wcount>=0,"missing spatial input buffer");
    if (is(name,"distance_matrix")) {
        require(wcount==0,"distance_matrix has no weights");
        require((std::isfinite(p) && (p==1. || p==2.)) ||
                (std::isinf(p) && p>0.),
                "distance_matrix supports Boost L1/L2/sup powers only");
    }
    Metric kind=is(name,"distance_matrix") ? Metric::euclidean : which(name,metric);
    if (is(name,"distance_matrix"))
        kind=p==1. ? Metric::cityblock :
             p==2. ? Metric::euclidean : Metric::chebyshev;
    if (paired(name))
        require(wcount==0,"pairwise Boost metrics do not accept weights");
    if (kind==Metric::minkowski || kind==Metric::minkowski_p)
        require((std::isfinite(p) && p>=1. && p<=MAX_POWER &&
                 p==std::floor(p)) || (std::isinf(p) && p>0.),
                "fractional/invalid Minkowski p has no Boost lp_distance equivalent");
    std::vector<double> result;
    result.reserve(count);
    if (is(name,"pdist")) {
        for (int i=0;i<urows;++i)
            for (int j=i+1;j<urows;++j)
                result.push_back(calculate(kind,u+i*dim,u+j*dim,dim,nullptr,0,p));
    } else if (is(name,"cdist") || is(name,"distance_matrix")) {
        for (int i=0;i<urows;++i)
            for (int j=0;j<vrows;++j)
                result.push_back(calculate(kind,u+i*dim,v+j*dim,dim,nullptr,0,p));
    } else {
        result.push_back(calculate(kind,u,v,dim,w,wcount,p));
    }
    require(static_cast<int>(result.size())==count,
            "unexpected Boost spatial result length");
    return result;
}
}
extern "C" {
const char* boost_spatial_error() {return last_error.c_str();}
int boost_spatial_supported(int library,const char* name) {
    return library==0 && known(name) ? 1 : 0;
}
int boost_spatial_output_count(const char* name,int urows,int vrows,int dim) {
    return needed(name,urows,vrows,dim);
}
int boost_spatial_evaluate(int library,const char* name,
                            const double* u,int urows,const double* v,
                            int vrows,int dim,const double* w,int wcount,
                            const char* metric,double p,double* out,int capacity) {
    last_error.clear();
    try {
        require(boost_spatial_supported(library,name)==1,
                "unsupported library or spatial function");
        const int count=needed(name,urows,vrows,dim);
        require(count>0 && out && capacity>=count,
                "invalid spatial shape/output capacity");
        const auto result=compute(name,u,urows,v,vrows,dim,w,wcount,metric,p);
        std::copy(result.begin(),result.end(),out);
        return 0;
    } catch (const std::exception& ex) {last_error=ex.what();}
      catch (...) {last_error="unknown Boost spatial failure";}
    return -1;
}
}
