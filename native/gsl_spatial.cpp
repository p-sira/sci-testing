// Spatial distances DERIVED from PUBLIC GSL 2.7 BLAS/statistical primitives.
// GSL has no named SciPy spatial.distance equivalents. These compositions
// deliberately preserve SciPy's weighting, zero-dimensional and edge rules.
// 1 = GSL provider index, other library values are unsupported.
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
thread_local std::string last_error;
constexpr int max_dim = 65536;

bool named(const char* x, const char* expected) {
    return x && std::strcmp(x, expected) == 0;
}
bool supported(const char* name) {
    return named(name,"euclidean") || named(name,"sqeuclidean") ||
           named(name,"cityblock") || named(name,"chebyshev") ||
           named(name,"cosine") || named(name,"correlation");
}
void require(bool pass, const char* reason) {
    if (!pass) throw std::domain_error(reason);
}
void gsl_check(int status, const char* call) {
    if (status != GSL_SUCCESS)
        throw std::runtime_error(std::string(call)+": "+gsl_strerror(status));
}
double dot(const std::vector<double>& a, const std::vector<double>& b) {
    require(a.size()==b.size(),"GSL BLAS dot size mismatch");
    if (a.empty()) return 0.;
    gsl_vector_const_view av = gsl_vector_const_view_array(a.data(),a.size());
    gsl_vector_const_view bv = gsl_vector_const_view_array(b.data(),b.size());
    double result=0.;
    gsl_check(gsl_blas_ddot(&av.vector,&bv.vector,&result),"GSL blas_ddot");
    return result;
}
double norm(const std::vector<double>& a) {
    if (a.empty()) return 0.;
    gsl_vector_const_view av=gsl_vector_const_view_array(a.data(),a.size());
    return gsl_blas_dnrm2(&av.vector);
}
double absolute_sum(const std::vector<double>& a) {
    if (a.empty()) return 0.;
    gsl_vector_const_view av=gsl_vector_const_view_array(a.data(),a.size());
    return gsl_blas_dasum(&av.vector);
}
double maximum_absolute(const std::vector<double>& a) {
    if (a.empty()) return 0.;
    gsl_vector_const_view av=gsl_vector_const_view_array(a.data(),a.size());
    return std::abs(a.at(gsl_blas_idamax(&av.vector)));
}

double calculate(const char* name, const double* u, const double* v,
                 int dim, const double* weights, int wcount) {
    require(supported(name) && dim>=0 && dim<=max_dim,
            "unknown GSL-derived spatial distance or invalid dimension");
    require(wcount==0 || (wcount==dim && weights),
            "spatial weights must match feature dimension");
    require(dim==0 || (u && v),"spatial coordinate arrays are missing");
    std::vector<double> diff;
    std::vector<double> first;
    std::vector<double> second;
    diff.reserve(dim);
    first.reserve(dim);
    second.reserve(dim);

    double total_weight=0.;
    for (int i=0;i<dim;++i) {
        const double a=u[i],b=v[i];
        const double wi=wcount ? weights[i] : 1.;
        require(std::isfinite(a) && std::isfinite(b) &&
                std::isfinite(wi) && wi>=0.,
                "spatial input or weight must be finite and weight nonnegative");
        total_weight+=wi;
        if (named(name,"chebyshev")) {
            if (wi>0.) diff.push_back(a-b); // SciPy uses only weight mask.
        } else if (named(name,"cityblock")) {
            diff.push_back((a-b)*wi);
        } else if (named(name,"sqeuclidean") ||
                   named(name,"euclidean")) {
            diff.push_back((a-b)*std::sqrt(wi));
        } else {
            first.push_back(a);
            second.push_back(b);
        }
    }
    if (named(name,"euclidean")) return norm(diff);   // gsl_blas_dnrm2
    if (named(name,"sqeuclidean")) return dot(diff,diff); // gsl_blas_ddot
    if (named(name,"cityblock")) return absolute_sum(diff); // gsl_blas_dasum
    if (named(name,"chebyshev")) {
        require(dim>0,"chebyshev undefined for empty feature vectors");
        return maximum_absolute(diff); // idamax; all-zero weights => zero
    }

    require(dim>0,"cosine/correlation undefined for empty vectors");
    require(total_weight>0.,"cosine/correlation require positive weight total");

    if (named(name,"correlation")) {
        if (!wcount) {
            // SciPy uses one minus sample correlation, no Bessel correction.
            const double c=gsl_stats_correlation(first.data(),1,second.data(),1,
                                                 static_cast<std::size_t>(dim));
            require(std::isfinite(c),"correlation undefined for constant vectors");
            // SciPy clamps only upper positive similarity for cosine style.
            return 1.-std::clamp(c,-1.,1.);
        }
        // SciPy's weighted correlation normalizes weights, centers using
        // weighted means, then computes the cosine of sqrt(w)-scaled
        // centered vectors. The normalization cancels algebraically.
        std::vector<double> w(weights,weights+dim);
        const double mu=gsl_stats_wmean(w.data(),1,first.data(),1,dim);
        const double mv=gsl_stats_wmean(w.data(),1,second.data(),1,dim);
        for (int i=0;i<dim;++i) {
            const double root=std::sqrt(w[i]);
            first[i]=(first[i]-mu)*root;
            second[i]=(second[i]-mv)*root;
        }
    } else if (wcount) {
        for (int i=0;i<dim;++i) {
            const double root=std::sqrt(weights[i]);
            first[i]*=root;
            second[i]*=root;
        }
    }

    // SciPy's cosine intentionally forms sqrt(dot(u,u)*dot(v,v));
    // that intermediate may overflow for very large finite vectors even
    // when a stable product-of-norms would not. Preserve saved SciPy
    // behavior, including its large-vector case returning distance 1.
    // Weighted correlation uses centered sqrt(weight) vectors and the
    // stable BLAS norm path, as its normalization differs.
    double similarity;
    if (named(name,"cosine")) {
        const double uu=dot(first,first);
        const double vv=dot(second,second);
        require(uu>0 && vv>0,"cosine undefined for zero-norm vectors");
        similarity=dot(first,second)/std::sqrt(uu*vv);
    } else {
        const double nu=norm(first), nv=norm(second);
        require(nu>0. && nv>0.,"correlation undefined for zero-norm vectors");
        similarity=dot(first,second)/(nu*nv);
    }
    require(std::isfinite(similarity),"nonfinite cosine/correlation similarity");
    // SciPy cosine has clip(1 - similarity, 0, 2).
    return std::clamp(1.-similarity,0.,2.);
}
} // namespace

extern "C" {
const char* gsl_spatial_error() {return last_error.c_str();}
int gsl_spatial_supported(int library,const char* name) {
    return library==1 && supported(name) ? 1 : 0;
}
int gsl_spatial_output_count(const char* name,int dim) {
    return supported(name) && dim>=0 && dim<=max_dim ? 1 : 0;
}
int gsl_spatial_evaluate(int library,const char* name,
                         const double* u,const double* v,int dim,
                         const double* w,int wcount,double* out,int capacity) {
    last_error.clear();
    try {
        require(gsl_spatial_supported(library,name)==1,
                "spatial distance does not support requested GSL provider");
        require(gsl_spatial_output_count(name,dim)==1 && out && capacity>=1,
                "invalid spatial result capacity or feature count");
        gsl_set_error_handler_off();
        const double value=calculate(name,u,v,dim,w,wcount);
        require(std::isfinite(value),"nonfinite GSL spatial output");
        out[0]=value;
        return 0;
    } catch(const std::exception& e) {last_error=e.what();}
      catch(...) {last_error="unknown GSL spatial evaluation failure";}
    return -1;
}
}
