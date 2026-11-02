// Derived statistics comparisons from public Boost.Math 1.83 APIs.
// All statistical distributions/moments are evaluated by Boost.Math.
// SciPy-style parameterization, degrees of freedom, and tail assembly are
// explicit here so provider coverage is not conflated with a direct API.
#include <boost/math/statistics/univariate_statistics.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
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
bool named(const char* a, const char* b) { return a && std::strcmp(a,b)==0; }
int output_count(const char* name) {
    if (named(name,"gmean") || named(name,"hmean") || named(name,"gstd")) return 1;
    if (named(name,"jarque_bera") || named(name,"chisquare") ||
        named(name,"power_divergence") || named(name,"mode") ||
        named(name,"binomtest") || named(name,"fisher_exact")) return 2;
    return 0;
}
void ensure(bool condition, const char* error) {
    if (!condition) throw std::invalid_argument(error);
}
unsigned integer(double v, unsigned limit) {
    ensure(std::isfinite(v) && v >= 0 && v <= limit &&
           std::floor(v) == v, "invalid nonnegative integer");
    return static_cast<unsigned>(v);
}
double positive_sum(const std::vector<double>& v) {
    double sum = std::accumulate(v.begin(),v.end(),0.);
    ensure(std::isfinite(sum), "nonfinite sum");
    return sum;
}
std::vector<double> vector(const double* input, int count) {
    ensure(count >= 0 && count <= 100000 && (input || !count),
           "invalid vector length");
    std::vector<double> result;
    if (count) result.assign(input,input+count);
    for (double x:result) ensure(std::isfinite(x),"nonfinite sample");
    return result;
}
double chi_tail(double statistic, double degrees) {
    ensure(std::isfinite(statistic) && statistic >= 0 && degrees > 0,
           "invalid chi-squared statistic/degrees");
    return boost::math::cdf(boost::math::complement(
        boost::math::chi_squared_distribution<double>(degrees),statistic));
}
// Binomial probabilities are evaluated by Boost's public discrete distribution.
double binomial_test(const double* p) {
    const unsigned k=integer(p[0],10000), n=integer(p[1],10000);
    const double probability=p[2];
    const unsigned alternative=integer(p[3],2);
    ensure(k<=n && probability>=0 && probability<=1,
           "invalid binomial test inputs");
    if (n==0) return 1.;
    if (probability==0 || probability==1) {
        const bool possible = probability==0 ? k==0 : k==n;
        if (alternative==0) return possible ? 1. : 0.;
        if (alternative==1) return probability==0 ? 1. : (k==n ? 1. : 0.);
        return probability==1 ? 1. : (k==0 ? 1. : 0.);
    }
    const boost::math::binomial_distribution<double> dist(n,probability);
    if (alternative==1) return boost::math::cdf(dist,k);
    if (alternative==2)
        return k==0 ? 1. : boost::math::cdf(boost::math::complement(dist,k-1));
    const double observed=boost::math::pdf(dist,k);
    double tail=0.;
    for (unsigned i=0;i<=n;++i) {
        const double current=boost::math::pdf(dist,i);
        if (current<=observed*(1.+1e-7))
            tail+=current;
    }
    return std::min(1.,tail);
}
void fisher(const std::vector<double>& x,const double* p,double* out) {
    ensure(x.size()==4, "Fisher test requires 2x2 data");
    const unsigned a=integer(x[0],10000),b=integer(x[1],10000);
    const unsigned c=integer(x[2],10000),d=integer(x[3],10000);
    const unsigned alternative=integer(p[0],2);
    const unsigned total=a+b+c+d;
    ensure(total>0 && total<=20000,"Fisher table outside validated range");
    const unsigned rows=a+b, cols=a+c;
    ensure(rows>0 && cols>0 && rows<total && cols<total,
           "degenerate Fisher table has undefined odds ratio");
    out[0]=static_cast<double>(a)*d/(static_cast<double>(b)*c);
    if (b==0 || c==0)
        throw std::domain_error("nonfinite Fisher odds ratio");
    boost::math::hypergeometric_distribution<double> dist(cols,rows,total);
    const unsigned lower=std::max(0,static_cast<int>(rows+cols)-static_cast<int>(total));
    const unsigned upper=std::min(rows,cols);
    if (alternative==1) {
        out[1]=boost::math::cdf(dist,a);
    } else if (alternative==2) {
        out[1]=a==lower ? 1. : boost::math::cdf(boost::math::complement(dist,a-1));
    } else {
        const double observed=boost::math::pdf(dist,a);
        double sum=0.;
        for (unsigned j=lower;j<=upper;++j)
            if (boost::math::pdf(dist,j)<=observed*(1.+1e-14))
                sum+=boost::math::pdf(dist,j);
        out[1]=std::min(1.,sum);
    }
}
void evaluate(const char* name,const std::vector<double>& x,
              const std::vector<double>& y,const std::vector<double>& p,
              double* out) {
    namespace stats=boost::math::statistics;
    if (named(name,"binomtest")) {
        ensure(x.empty() && y.empty() && p.size()==4, "invalid binomtest shape");
        const unsigned k=integer(p[0],10000),n=integer(p[1],10000);
        ensure(k<=n,"binomial successes exceed trials");
        out[0]=n ? static_cast<double>(k)/n : 0.;
        out[1]=binomial_test(p.data());
        return;
    }
    ensure(!x.empty(), "requires a nonempty sample vector");
    if (named(name,"fisher_exact")) {
        ensure(p.size()==1, "invalid Fisher parameter shape");
        fisher(x,p.data(),out);
        return;
    }
    if (named(name,"mode")) {
        ensure(p.empty() && y.empty(), "mode does not accept parameters");
        // Public Boost mode returns all tied modes, unlike SciPy's smallest
        // matching value. Sorting before the public call makes its order stable.
        auto sorted=x;
        std::sort(sorted.begin(),sorted.end());
        const auto modes=stats::mode(sorted);
        ensure(!modes.empty(),"empty Boost mode result");
        out[0]=*std::min_element(modes.begin(),modes.end());
        out[1]=static_cast<double>(std::count(sorted.begin(),sorted.end(),out[0]));
        return;
    }
    if (named(name,"gmean") || named(name,"hmean") || named(name,"gstd")) {
        ensure(y.empty(),"unexpected second input");
        for (double v:x) ensure(v>0,"geometric/harmonic inputs must be positive");
        if (named(name,"hmean")) {
            std::vector<double> reciprocal;
            reciprocal.reserve(x.size());
            for (double v:x) reciprocal.push_back(1./v);
            out[0]=1./stats::mean(reciprocal);
            return;
        }
        std::vector<double> logs;
        logs.reserve(x.size());
        for (double v:x) logs.push_back(std::log(v));
        if (named(name,"gmean")) {
            ensure(p.empty(),"unexpected gmean parameters");
            out[0]=std::exp(stats::mean(logs));
            return;
        }
        ensure(p.size()==1, "invalid gstd parameters");
        const unsigned ddof=integer(p[0],1000);
        ensure(x.size()>ddof,"gstd degrees of freedom");
        out[0]=std::exp(std::sqrt(stats::variance(logs)*
                                  x.size()/(x.size()-ddof)));
        return;
    }
    if (named(name,"jarque_bera")) {
        ensure(p.empty() && y.empty() && x.size()>=3,"invalid Jarque-Bera input");
        const double skew=stats::skewness(x);
        const double excess=stats::kurtosis(x)-3.;
        out[0]=(x.size()/6.)*(skew*skew+excess*excess/4.);
        out[1]=chi_tail(out[0],2.);
        return;
    }
    if (named(name,"chisquare") || named(name,"power_divergence")) {
        ensure(y.size()==x.size() && x.size()>=2,
               "expected frequencies must match observed shape");
        const bool power=named(name,"power_divergence");
        ensure(p.size()==(power ? 2u : 1u),"invalid divergence parameters");
        unsigned ddof=integer(p[0],1000);
        ensure(x.size()>ddof+1,"chi-squared degrees of freedom");
        const double lambda=power ? p[1] : 1.;
        for (std::size_t i=0;i<x.size();++i)
            ensure(x[i]>=0 && y[i]>0,"invalid observed/expected frequencies");
        const double sum_x=positive_sum(x),sum_y=positive_sum(y);
        ensure(std::abs(sum_x-sum_y) <= 1e-8*std::max(sum_x,sum_y),
               "observed and expected totals differ");
        double total=0.;
        for (std::size_t i=0;i<x.size();++i) {
            const double a=x[i],b=y[i];
            if (lambda==1.) {
                total+=(a-b)*(a-b)/b;
            } else if (lambda==0.) {
                if (a>0.) total+=2*a*std::log(a/b);
            } else if (lambda==-1.) {
                ensure(a>0,"divergence infinite for zero observed frequency");
                total+=2*b*std::log(b/a);
            } else if (lambda==-0.5) {
                total+=4*std::pow(std::sqrt(a)-std::sqrt(b),2);
            } else {
                ensure(a>0 || lambda>0, "divergence undefined at zero frequency");
                total+=2*a*(std::pow(a/b,lambda)-1)/(lambda*(lambda+1));
            }
        }
        out[0]=std::max(0.,total);
        out[1]=chi_tail(out[0],x.size()-1-ddof);
        return;
    }
    throw std::invalid_argument("unknown statistics function");
}
}

extern "C" {
int boost_stats_more_supported(int library,const char* name) {
    return library==0 ? output_count(name) : 0;
}
int boost_stats_more_output_count(const char* name) { return output_count(name); }
const char* boost_stats_more_error() { return last_error.c_str(); }
int boost_stats_more_evaluate(int library,const char* name,
                              const double* x,int nx,const double* y,int ny,
                              const double* p,int np,double* out,int capacity) {
    last_error.clear();
    try {
        int n=boost_stats_more_supported(library,name);
        ensure(n && out && capacity>=n,"unsupported provider or output capacity");
        auto xs=vector(x,nx),ys=vector(y,ny),ps=vector(p,np);
        double values[2]{};
        evaluate(name,xs,ys,ps,values);
        for (int i=0;i<n;++i) ensure(std::isfinite(values[i]),"nonfinite result");
        std::copy_n(values,n,out);
        return 0;
    } catch (const std::exception& exception) {
        last_error=exception.what();
    } catch (...) { last_error="unknown Boost statistics error"; }
    return -1;
}
}
