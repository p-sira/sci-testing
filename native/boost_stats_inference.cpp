// Public Boost.Math statistics and probability distributions composed into
// Sci's saved inference and circular-statistics signatures.
// Supports the Boost 1.83 surface: mean(), variance(), students_t,
// fisher_f_distribution, chi_squared_distribution.
#include <boost/math/statistics/univariate_statistics.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
thread_local std::string error_text;
constexpr double pi = 3.14159265358979323846264338327950288;
bool eq(const char* a, const char* b) {
    return a && std::strcmp(a, b) == 0;
}
int output_count(const char* name) {
    if (eq(name, "ttest_ind_from_stats")) return 3;
    if (eq(name, "f_oneway") || eq(name, "bartlett")) return 2;
    if (eq(name, "circmean") || eq(name, "circvar") ||
        eq(name, "circstd")) return 1;
    return 0;
}
void require(bool valid, const char* description) {
    if (!valid) throw std::invalid_argument(description);
}
bool flag(double value) {
    require(value == 0. || value == 1., "expected boolean 0 or 1");
    return value == 1.;
}
int alternative(double value) {
    require(value == 0. || value == 1. || value == 2.,
            "alternative must be 0(two-sided), 1(less), or 2(greater)");
    return static_cast<int>(value);
}
std::vector<double> collect(const double* values, int n) {
    require(n >= 0 && n <= 100000 && (values || n == 0),
            "invalid vector pointer or length");
    std::vector<double> result;
    if (n > 0) result.assign(values, values + n);
    for (double value : result)
        require(std::isfinite(value), "nonfinite sample");
    return result;
}
double mean(const std::vector<double>& values) {
    return boost::math::statistics::mean(values);
}
double sample_variance(const std::vector<double>& values) {
    require(values.size() >= 2, "at least two values per group required");
    return boost::math::statistics::variance(values) *
        (static_cast<double>(values.size()) / (values.size()-1));
}
double student_tail(double statistic, double df, int alt) {
    require(df > 0 && std::isfinite(df) && std::isfinite(statistic),
            "nonfinite t statistic or degrees of freedom");
    const boost::math::students_t_distribution<double> t(df);
    if (alt == 1) return boost::math::cdf(t, statistic);
    if (alt == 2)
        return boost::math::cdf(boost::math::complement(t,statistic));
    return 2 * boost::math::cdf(
        boost::math::complement(t,std::abs(statistic)));
}
void ttest_from_stats(const std::vector<double>& p,double* output) {
    require(p.size() == 8, "ttest_ind_from_stats requires 8 parameters");
    const auto& v = p;
    const double n1=v[2], n2=v[5];
    require(n1 >= 2 && n2 >= 2 && n1 == std::floor(n1) &&
            n2 == std::floor(n2) && n1 <= 1e9 && n2 <= 1e9 &&
            v[1] >= 0 && v[4] >= 0,
            "invalid sample sizes or standard deviations");
    const bool equal = flag(v[6]);
    const int alt=alternative(v[7]);
    const double delta=v[0]-v[3];
    double df=0, variance=0;
    if (equal) {
        df = n1+n2-2;
        const double pooled=((n1-1)*v[1]*v[1]+(n2-1)*v[4]*v[4])/df;
        variance=pooled*(1/n1+1/n2);
    } else {
        const double a=v[1]*v[1]/n1,b=v[4]*v[4]/n2;
        variance=a+b;
        df=(a+b)*(a+b)/(a*a/(n1-1)+b*b/(n2-1));
    }
    require(variance > 0 && std::isfinite(variance),
            "zero or nonfinite t-test variance");
    output[0] = delta/std::sqrt(variance);
    output[1] = student_tail(output[0],df,alt);
    output[2] = df;
}
void f_oneway(const std::vector<std::vector<double>>& groups,
              double* output) {
    const auto k=groups.size();
    require(k>=2, "ANOVA needs at least two groups");
    double total=0, n=0;
    std::vector<double> group_mean;
    group_mean.reserve(k);
    for (const auto& group:groups) {
        require(group.size()>=2, "ANOVA groups require two or more samples");
        const double m=mean(group);
        group_mean.push_back(m);
        total+=static_cast<double>(group.size())*m;
        n+=group.size();
    }
    const double grand_mean=total/n;
    double between=0, within=0;
    for (std::size_t i=0;i<k;++i) {
        const auto& group=groups[i];
        const double difference=group_mean[i]-grand_mean;
        between+=group.size()*difference*difference;
        within+=(group.size()-1)*sample_variance(group);
    }
    require(within>0 && std::isfinite(within),
            "ANOVA undefined with zero within-group variance");
    output[0]=(between/(k-1))/(within/(n-k));
    require(output[0]>=0 && std::isfinite(output[0]),
            "ANOVA nonfinite statistic");
    const boost::math::fisher_f_distribution<double> dist(k-1,n-k);
    output[1]=boost::math::cdf(boost::math::complement(dist,output[0]));
}
void bartlett(const std::vector<std::vector<double>>& groups,
              double* output) {
    const auto k=groups.size();
    require(k>=2, "Bartlett test needs at least two groups");
    double degrees=0,weighted=0,weighted_log=0,reciprocal=0;
    for (const auto& group:groups) {
        const double df=group.size()-1.;
        const double s2=sample_variance(group);
        require(s2>0 && std::isfinite(s2),
                "Bartlett requires positive sample variance");
        degrees+=df;
        weighted+=df*s2;
        weighted_log+=df*std::log(s2);
        reciprocal+=1/df;
    }
    const double pooled=weighted/degrees;
    const double correction=1+(reciprocal-1/degrees)/(3*(k-1));
    const double statistic=(degrees*std::log(pooled)-weighted_log)/correction;
    require(std::isfinite(statistic) && statistic>=-1e-11,
            "Bartlett nonfinite/negative statistic");
    output[0]=std::max(0.,statistic);
    boost::math::chi_squared_distribution<double> chi(k-1.);
    output[1]=boost::math::cdf(boost::math::complement(chi,output[0]));
}
void circular(const char* name, const std::vector<double>& x,
              const std::vector<double>& p, double* output) {
    const bool is_std=eq(name,"circstd");
    require(p.size() == (is_std ? 3u : 2u) && !x.empty(),
            "circular stats require sample, low, high, optional normalize");
    const double low=p[0], high=p[1],span=high-low;
    require(span>0 && std::isfinite(span), "invalid circular angle range");
    bool normalize=is_std ? flag(p[2]) : false;
    std::vector<double> sine,cosine;
    sine.reserve(x.size());
    cosine.reserve(x.size());
    const double factor=2*pi/span;
    for (double value:x) {
        const double angle=(value-low)*factor;
        sine.push_back(std::sin(angle));
        cosine.push_back(std::cos(angle));
    }
    const double s=mean(sine),c=mean(cosine);
    if (eq(name,"circmean")) {
        double angle=std::atan2(s,c);
        if (angle<0) angle+=2*pi;
        output[0]=low+angle/factor;
    } else {
        const double R=std::clamp(std::hypot(s,c),0.,1.);
        if (eq(name,"circvar")) output[0]=1-R;
        else output[0]=std::sqrt(-2*std::log(R))*
              (normalize ? 1. : 1/factor);
    }
}
std::vector<std::vector<double>> split_groups(
    const std::vector<double>& flat,const int* sizes,int n_groups) {
    require(n_groups>=2 && n_groups<=10000 && sizes,
            "group lengths not supplied");
    std::vector<std::vector<double>> groups;
    std::size_t position=0;
    for (int i=0;i<n_groups;++i) {
        require(sizes[i]>=2 && static_cast<std::size_t>(sizes[i])<=flat.size()-position,
                "invalid group size");
        const std::size_t end=position+sizes[i];
        groups.emplace_back(flat.begin()+position,flat.begin()+end);
        position=end;
    }
    require(position==flat.size(),"group lengths differ from flattened samples");
    return groups;
}
void evaluate(const char* name,const std::vector<double>& x,
              const std::vector<double>& y,const std::vector<double>& p,
              const std::vector<double>& g,const int* sizes,int ng,
              double* out) {
    require(y.empty(), "unexpected y vector");
    if (eq(name,"ttest_ind_from_stats")) {
        require(x.empty() && g.empty() && ng==0, "ttest accepts only p");
        ttest_from_stats(p,out);
    } else if (eq(name,"f_oneway") || eq(name,"bartlett")) {
        require(x.empty() && p.empty(), "group tests accept only g");
        const auto groups=split_groups(g,sizes,ng);
        if (eq(name,"f_oneway")) f_oneway(groups,out);
        else bartlett(groups,out);
    } else {
        require(g.empty() && ng==0, "circular stats accept only x and p");
        circular(name,x,p,out);
    }
}
}  // namespace

extern "C" {
int boost_stats_inference_supported(int library,const char* name) {
    return library==0 ? output_count(name) : 0;
}
int boost_stats_inference_output_count(const char* name) {
    return output_count(name);
}
const char* boost_stats_inference_error() { return error_text.c_str(); }
int boost_stats_inference_evaluate(
    int library,const char* name,const double* x,int nx,const double* y,int ny,
    const double* p,int np,const double* g,int ng,const int* group_sizes,
    int n_groups,double* out,int capacity) {
    error_text.clear();
    try {
        const int count=boost_stats_inference_supported(library,name);
        require(count && out && capacity>=count,"unsupported provider or invalid output capacity");
        const auto xs=collect(x,nx),ys=collect(y,ny);
        const auto ps=collect(p,np),gs=collect(g,ng);
        double result[3]{};
        evaluate(name,xs,ys,ps,gs,group_sizes,n_groups,result);
        for (int i=0;i<count;++i)
            require(std::isfinite(result[i]),"nonfinite statistics result");
        std::copy_n(result,count,out);
        return 0;
    } catch (const std::exception& error) {
        error_text=error.what();
    } catch (...) { error_text="unexpected Boost inference exception"; }
    return -1;
}
}
