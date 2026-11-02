// Public Boost.Math Bessel functions and real positive-root primitives.
// The separate ABI preserves SciPy vector/tuple shapes, including the
// globally merged metadata of jnjnp_zeros.
//
// Intended for Boost >=1.83: cyl_bessel_j, cyl_neumann,
// cyl_bessel_j_zero, cyl_neumann_zero. No private Boost detail APIs.
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/airy.hpp>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace {
thread_local std::string last_error;
constexpr int max_count = 256;
constexpr int max_joint = 128;  // Sci Rust contract; SciPy permits 1200
constexpr int max_order = 256;

int arity(const char* name) {
    if (!name) return 0;
    if (!std::strcmp(name, "jnp_zeros") ||
        !std::strcmp(name, "ynp_zeros") ||
        !std::strcmp(name, "jnyn_zeros")) return 2;
    if (!std::strcmp(name, "y0_zeros") ||
        !std::strcmp(name, "y1_zeros") ||
        !std::strcmp(name, "y1p_zeros") ||
        !std::strcmp(name, "ai_zeros") ||
        !std::strcmp(name, "bi_zeros") ||
        !std::strcmp(name, "jnjnp_zeros")) return 1;
    return 0;  // Existing jn_zeros/yn_zeros handled in special_extra.cpp.
}

int lanes(const char* name) {
    if (!name) return 0;
    if (!std::strcmp(name, "jnyn_zeros") ||
        !std::strcmp(name, "ai_zeros") ||
        !std::strcmp(name, "bi_zeros") ||
        !std::strcmp(name, "jnjnp_zeros")) return 4;
    if (!std::strcmp(name, "y0_zeros") ||
        !std::strcmp(name, "y1_zeros") ||
        !std::strcmp(name, "y1p_zeros")) return 2;
    return arity(name) ? 1 : 0;
}

void require(bool condition, const char* message) {
    if (!condition) throw std::domain_error(message);
}

int input_integer(double input, int lo, int hi, const char* message) {
    require(std::isfinite(input) && input >= lo && input <= hi &&
            input == std::floor(input), message);
    return static_cast<int>(input);
}

int requested_count(const char* name, const double* args, int count) {
    require(args && arity(name) == count, "unsupported Bessel zero API or arity");
    if (count == 2) input_integer(args[0], 0, max_order, "invalid Bessel order");
    return input_integer(args[count-1], 1,
                         !std::strcmp(name, "jnjnp_zeros") ? max_joint : max_count,
                         "invalid Bessel zero count");
}

double j(int n, double x) { return boost::math::cyl_bessel_j(double(n), x); }
double y(int n, double x) { return boost::math::cyl_neumann(double(n), x); }
double jp(int n, double x) {
    return n == 0 ? -j(1,x) : n / x * j(n,x) - j(n+1,x);
}
double yp(int n, double x) {
    return n == 0 ? -y(1,x) : n / x * y(n,x) - y(n+1,x);
}

// Brackets use neighboring zeros of the *undifferentiated* function:
// J'_n first root in (n, j[n,1]) for n>=1, then in successive J_n
// root intervals; Y'_n roots lie between successive Y_n roots. For n=0,
// J'_0=-J_1 and Y'_0=-Y_1 exactly. This avoids phase-grid skipping.
template <class F>
double bisect(F function, double left, double right) {
    double fleft = function(left), fright = function(right);
    require(std::isfinite(fleft) && std::isfinite(fright),
            "nonfinite Bessel derivative at root bracket");
    require((fleft < 0) != (fright < 0),
            "Bessel derivative root was not bracketed");
    for (int i=0; i<100; ++i) {
        double mid = left + (right-left)*.5;
        if (mid == left || mid == right) break;
        double value = function(mid);
        require(std::isfinite(value), "nonfinite Bessel derivative during inversion");
        if (value == 0) return mid;
        if ((value < 0) == (fleft < 0)) {
            left = mid;
            fleft = value;
        } else right = mid;
    }
    double root = left + (right-left)*.5;
    require(std::abs(function(root)) < 5e-10,
            "Bessel derivative root failed residual check");
    return root;
}

double root(int n, int serial, int kind) {
    require(n >= 0 && n <= max_order && serial >= 1 && serial <= max_count+1,
            "invalid Bessel root rank/order");
    if (kind == 0) return boost::math::cyl_bessel_j_zero(double(n), serial);
    if (kind == 2) return boost::math::cyl_neumann_zero(double(n), serial);
    if (kind == 1) {
        if (n == 0) return root(1,serial,0);
        double hi = root(n,serial,0);
        double lo = serial == 1 ? double(n) : root(n,serial-1,0);
        return bisect([&](double x){return jp(n,x);},lo,hi);
    }
    if (kind == 3) {
        if (n == 0) return root(1,serial,2);
        double lo = root(n,serial,2);
        double hi = root(n,serial+1,2);
        return bisect([&](double x){return yp(n,x);},lo,hi);
    }
    throw std::invalid_argument("invalid root kind");
}

struct Entry {
    double x;
    int n, rank, derivative;
};

double airy_root(bool first_kind, int serial) {
    return first_kind ? boost::math::airy_ai_zero<double>(serial) :
                        boost::math::airy_bi_zero<double>(serial);
}

double airy_function(bool first_kind, double x) {
    return first_kind ? boost::math::airy_ai(x) : boost::math::airy_bi(x);
}

double airy_prime(bool first_kind, double x) {
    return first_kind ? boost::math::airy_ai_prime(x) :
                        boost::math::airy_bi_prime(x);
}

// Airy Ai' has its first root between Ai(1) and 0; successive roots
// between Ai(k) and Ai(k-1). Bi' has its first root strictly between
// Bi(1) and Bi(2), with the same interlacing for subsequent ranks.
double airy_derivative_root(bool first_kind, int serial) {
    const double left = first_kind ? airy_root(true,serial) :
                                     airy_root(false,serial+1);
    const double right = first_kind ?
        (serial == 1 ? 0. : airy_root(true,serial-1)) :
        airy_root(false,serial);
    return bisect([&](double x){return airy_prime(first_kind,x);}, left,right);
}

std::vector<Entry> joint_zeros(int nt) {
    std::vector<Entry> entries{{0.,0,0,1}};
    // Increase a bound until >=nt qualifying roots have been found.
    // The first positive J_n/J'_n root is >n for n>0, so all higher
    // degrees beyond cutoff cannot add an earlier positive root.
    double cutoff = 5 + 3*std::sqrt(double(nt));
    for (int expansion=0; expansion<7; ++expansion) {
        entries.resize(1);
        int degree_limit = std::min(max_order, int(std::ceil(cutoff)));
        for (int n=0; n<=degree_limit; ++n) {
            for (int kind=0; kind<=1; ++kind) {
                for (int serial=1; serial<=nt; ++serial) {
                    double x = root(n,serial,kind);
                    if (x > cutoff) break;
                    entries.push_back({x,n,serial,kind});
                }
            }
        }
        if (int(entries.size()) >= nt) break;
        cutoff *= 1.7;
    }
    require(int(entries.size()) >= nt, "insufficient joint Bessel roots");
    std::sort(entries.begin(),entries.end(),[](const Entry& a,const Entry& b){
        return a.x < b.x;
    });
    // Exact identities J'_0=-J_1 must have both members: J_1 first,
    // J'_0 second. Group nearby floating-point ties without using a
    // non-transitive tolerance in the std::sort comparator.
    for (std::size_t begin=0; begin<entries.size();) {
        std::size_t end=begin+1;
        while (end<entries.size() && std::abs(entries[end].x-entries[begin].x)<1e-11)
            ++end;
        std::stable_sort(entries.begin()+begin,entries.begin()+end,
                         [](const Entry& a,const Entry& b){
                             return std::tie(a.derivative,a.n,a.rank) <
                                    std::tie(b.derivative,b.n,b.rank);
                         });
        begin=end;
    }
    entries.resize(nt);
    return entries;
}

std::vector<double> calculate(const char* name, const double* args, int count) {
    int nt = requested_count(name,args,count);
    int order = count == 2 ? static_cast<int>(args[0]) : 0;
    int components = lanes(name);
    std::vector<double> result(static_cast<std::size_t>(nt)*components);
    if (!std::strcmp(name,"jnjnp_zeros")) {
        auto entries = joint_zeros(nt);
        for (int k=0; k<nt; ++k) {
            result[k] = entries[k].x;
            result[nt+k] = entries[k].n;
            result[2*nt+k] = entries[k].rank;
            result[3*nt+k] = entries[k].derivative;
        }
    } else if (!std::strcmp(name,"ai_zeros") || !std::strcmp(name,"bi_zeros")) {
        bool first_kind=!std::strcmp(name,"ai_zeros");
        for (int k=1; k<=nt; ++k) {
            double z = airy_root(first_kind,k);
            double zp = airy_derivative_root(first_kind,k);
            result[k-1] = z;
            result[nt+k-1] = zp;
            result[2*nt+k-1] = airy_function(first_kind,zp);
            result[3*nt+k-1] = airy_prime(first_kind,z);
        }
    } else if (!std::strcmp(name,"jnyn_zeros")) {
        for (int kind=0; kind<4; ++kind)
            for (int k=1; k<=nt; ++k) result[kind*nt+k-1] = root(order,k,kind);
    } else if (!std::strcmp(name,"jnp_zeros") ||
               !std::strcmp(name,"ynp_zeros")) {
        int kind = !std::strcmp(name,"jnp_zeros") ? 1 : 3;
        for (int k=1; k<=nt; ++k) result[k-1] = root(order,k,kind);
    } else {
        int kind = !std::strcmp(name,"y0_zeros") ? 2 :
                   !std::strcmp(name,"y1_zeros") ? 2 : 3;
        int degree = !std::strcmp(name,"y0_zeros") ? 0 : 1;
        for (int k=1; k<=nt; ++k) {
            double x = root(degree,k,kind);
            result[k-1] = x;
            // SciPy's ACTUAL y0_zeros paired values are Y1 (not Y0'=-Y1),
            // despite its documentation calling them derivatives. This
            // matches the Sci Rust facade and the saved binary64 oracle.
            // y1_zeros -> Y1'=Y0-Y1/x; y1p_zeros -> Y1.
            result[nt+k-1] = degree == 0 ? y(1,x) :
                             kind == 2 ? yp(1,x) : y(1,x);
        }
    }
    for (double v:result) require(std::isfinite(v),"nonfinite Bessel zero/value");
    return result;
}

// Caching by exact integer name/order/count prevents the 2,001-row generated
// fixture sweeps from recomputing the same vector thousands of times. Bound
// memory; this is thread-local and does not share mutable state between jobs.
std::vector<double> cached(const char* name, const double* args, int count) {
    int nt=requested_count(name,args,count);
    int n=count==2 ? static_cast<int>(args[0]) : 0;
    std::string key=std::string(name)+":"+std::to_string(n)+":"+std::to_string(nt);
    thread_local std::unordered_map<std::string,std::vector<double>> cache;
    auto found=cache.find(key);
    if (found!=cache.end()) return found->second;
    auto values=calculate(name,args,count);
    if (cache.size()>=256) cache.clear();
    cache.emplace(std::move(key),values);
    return values;
}
}  // namespace

extern "C" {
const char* bessel_zeros_error() {return last_error.c_str();}
int bessel_zeros_supported(int library,const char* name) {
    return library == 0 ? arity(name) : 0;  // Boost.Math only.
}
int bessel_zeros_output_count(const char* name,const double* args,int count) {
    try {return requested_count(name,args,count)*lanes(name);}
    catch (...) {return 0;}
}
int bessel_zeros_evaluate(int library,const char* name,const double* args,
                          int count,double* output,int capacity) {
    last_error.clear();
    try {
        require(bessel_zeros_supported(library,name)==count && count>0,
                "unsupported library or Bessel zero API");
        int needed=bessel_zeros_output_count(name,args,count);
        require(needed>0 && output && capacity>=needed,"invalid Bessel vector output capacity");
        auto values=cached(name,args,count);
        std::copy(values.begin(),values.end(),output);
        return 0;
    } catch(const std::exception& e) {last_error=e.what();}
      catch(...) {last_error="unknown Bessel zero failure";}
    return -1;
}
}
