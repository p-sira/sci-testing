// Public Boost.Math 1.83+ interpolation adapters for saved SciPy value queries.
//
// * pchip_interpolate: boost::math::interpolators::pchip with explicit
//   SciPy/PCHIP one-sided endpoint derivatives, not Boost's default secants.
// * barycentric_interpolate: boost::math::interpolators::barycentric_rational
//   with order=n-1, which is mathematically polynomial Lagrange interpolation.
// * krogh_interpolate: the SAME polynomial representation only for pairwise
//   DISTINCT nodes. This is not a general Krogh/Hermite confluent interpolator:
//   repeated nodes/derivative tuples are intentionally rejected.
//
// All outputs are query values only. Boost's PCHIP refuses extrapolation;
// this adapter explicitly rejects such requests rather than fabricating
// SciPy's extrapolating contract outside the proven intersection.
#include <boost/math/interpolators/pchip.hpp>
#include <boost/math/interpolators/barycentric_rational.hpp>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
thread_local std::string last_error;
// The saved polynomial degrees are at most nine; avoid untested, poorly
// conditioned high-degree extrapolation rather than implying broad fidelity.
constexpr int MAX_NODES = 16;
constexpr int MAX_QUERIES = 2048;
bool named(const char* name,const char* match) {
    return name && std::strcmp(name,match)==0;
}
bool supported(const char* name) {
    return named(name,"pchip_interpolate") ||
           named(name,"barycentric_interpolate") ||
           named(name,"krogh_interpolate");
}

// SciPy PchipInterpolator._edge_case: endpoint one-sided shape-preserving
// slope, distinct from public Boost.Math pchip's default boundary secant.
double pchip_endpoint(double h0,double h1,double slope0,double slope1) {
    double d=((2*h0+h1)*slope0-h0*slope1)/(h0+h1);
    if(slope0==0. || (d>0)!=(slope0>0))return 0.;
    if((slope1==0. || (slope0>0)!=(slope1>0)) &&
       std::abs(d)>3*std::abs(slope0))
        return 3*slope0;
    return d;
}
void validate(const char* name,const double* nodes,const double* y,
              int n,const double* query,int nq,double* output,int capacity) {
    if(!supported(name))throw std::invalid_argument("unknown public Boost interpolation API");
    if(n<2 || n>MAX_NODES || (named(name,"pchip_interpolate") && n<4))
        throw std::domain_error("invalid interpolation node count");
    if(nq<1 || nq>MAX_QUERIES)
        throw std::domain_error("invalid interpolation query count");
    if(!nodes || !y || !query || !output || capacity<nq)
        throw std::invalid_argument("null input or insufficient interpolation output capacity");
    for(int i=0;i<n;++i) {
        if(!std::isfinite(nodes[i]) || !std::isfinite(y[i]))
            throw std::domain_error("interpolation nodes and values must be finite real scalars");
        if(i && !(nodes[i]>nodes[i-1]))
            throw std::domain_error("nodes must be strictly increasing and distinct; confluent Krogh is unsupported");
    }
    for(int i=0;i<nq;++i) {
        if(!std::isfinite(query[i]))
            throw std::domain_error("interpolation queries must be finite real scalars");
        if(named(name,"pchip_interpolate") &&
           (query[i]<nodes[0] || query[i]>nodes[n-1]))
            throw std::domain_error("Boost pchip cannot extrapolate; query outside nodes");
    }
}
}

extern "C" {
const char* boost_interpolation_error() {return last_error.c_str();}
int boost_interpolation_supported(const char* name) {return supported(name)?1:0;}
int boost_interpolation_evaluate(const char* name,
        const double* nodes,const double* y,int n,
        const double* query,int nq,double* output,int capacity) {
    last_error.clear();
    try {
        validate(name,nodes,y,n,query,nq,output,capacity);
        std::vector<double> result(static_cast<size_t>(nq));
        if(named(name,"pchip_interpolate")) {
            const double h0=nodes[1]-nodes[0], h1=nodes[2]-nodes[1];
            const double slope0=(y[1]-y[0])/h0, slope1=(y[2]-y[1])/h1;
            const double hend=nodes[n-1]-nodes[n-2],hprev=nodes[n-2]-nodes[n-3];
            const double send=(y[n-1]-y[n-2])/hend;
            const double sprev=(y[n-2]-y[n-3])/hprev;
            const double left=pchip_endpoint(h0,h1,slope0,slope1);
            const double right=pchip_endpoint(hend,hprev,send,sprev);
            boost::math::interpolators::pchip<std::vector<double>> interp(
                std::vector<double>(nodes,nodes+n),std::vector<double>(y,y+n),
                left,right);
            for(int i=0;i<nq;++i)result[i]=interp(query[i]);
        } else {
            // n-1 barycentric rational order is EXACT polynomial interpolation,
            // so this is a valid public Boost basis for distinct-node Krogh.
            boost::math::interpolators::barycentric_rational<double> interp(
                nodes,y,static_cast<size_t>(n),static_cast<size_t>(n-1));
            for(int i=0;i<nq;++i)result[i]=interp(query[i]);
        }
        for(double value:result)
            if(!std::isfinite(value))
                throw std::runtime_error("nonfinite interpolation output");
        std::copy(result.begin(),result.end(),output);
        return 0;
    } catch(const std::exception& ex) {
        last_error=ex.what();
        return -1;
    } catch(...) {
        last_error="unknown public Boost interpolation failure";
        return -1;
    }
}
}
