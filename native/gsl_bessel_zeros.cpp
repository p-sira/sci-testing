// GSL-DERIVED real integer-order Bessel zeros; no public GSL Y-zero API exists.
// Uses public gsl_sf_bessel_{Jn,Yn,zero_Jnu}_e plus gsl_root_fsolver_brent.
// Full SciPy/Sci contracts: scalar vector or column-major flattened tuples.
// This translation unit deliberately provides GSL (library=1) only.
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_roots.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace {
thread_local std::string last_error;
constexpr int MAX_N = 64;
constexpr int MAX_COUNT = 128;
constexpr int MAX_JOINT = 128;

void require(bool ok, const char* why) {
    if (!ok) throw std::domain_error(why);
}

bool named(const char* actual, const char* expected) {
    return actual && std::strcmp(actual, expected) == 0;
}

int arity(const char* name) {
    if (named(name, "jnp_zeros") || named(name, "yn_zeros") ||
        named(name, "ynp_zeros") || named(name, "jnyn_zeros")) return 2;
    if (named(name, "jnjnp_zeros") || named(name, "y0_zeros") ||
        named(name, "y1_zeros") || named(name, "y1p_zeros")) return 1;
    return 0;
}

int lanes(const char* name) {
    if (named(name, "jnyn_zeros") || named(name, "jnjnp_zeros")) return 4;
    if (named(name, "y0_zeros") || named(name, "y1_zeros") ||
        named(name, "y1p_zeros")) return 2;
    return arity(name) ? 1 : 0;
}

int checked_int(double x, int min, int max, const char* description) {
    require(std::isfinite(x) && x == std::floor(x) && x >= min && x <= max,
            description);
    return static_cast<int>(x);
}

int count_checked(const char* name, const double* args, int count) {
    require(name && args && arity(name) == count, "unknown Bessel zeros API or wrong arity");
    if (count == 2)
        checked_int(args[0], 0, MAX_N, "Bessel order must be an integer in [0,64]");
    return checked_int(args[count-1], 1,
                       named(name, "jnjnp_zeros") ? MAX_JOINT : MAX_COUNT,
                       "Bessel zero count must be an integer in [1,128]");
}

double eval_jy(int family, int n, double x) {
    require(n >= 0 && n <= MAX_N + 1 && std::isfinite(x) && x > 0.,
            "nonpositive or invalid real Bessel argument");
    gsl_sf_result result{};
    // IMPORTANT: gsl_sf_bessel_Ynu_e(0,J0-zero) sometimes reports an
    // incorrect sign even with status GSL_SUCCESS; use Yn_e at all orders.
    const int rc = family == 0 ? gsl_sf_bessel_Jn_e(n, x, &result)
                               : gsl_sf_bessel_Yn_e(n, x, &result);
    require(rc == GSL_SUCCESS && std::isfinite(result.val),
            "GSL integer Bessel evaluation failed or returned nonfinite output");
    return result.val;
}

double zero_j(int n, int serial) {
    require(n >= 0 && n <= MAX_N + 1 && serial >= 1 &&
            serial <= MAX_COUNT + 2, "invalid GSL J zero order/rank");
    gsl_sf_result result{};
    const int rc = gsl_sf_bessel_zero_Jnu_e(static_cast<double>(n),
                                             static_cast<unsigned>(serial), &result);
    require(rc == GSL_SUCCESS && std::isfinite(result.val) && result.val > 0.,
            "GSL J_n zero evaluation failed");
    double x = result.val;
    // GSL's zero approximation can be inaccurate for higher orders.
    // Newton polish uses public integer-order GSL values and their recurrence.
    for (int i = 0; i < 5; ++i) {
        const double value = eval_jy(0, n, x);
        const double derivative = (n == 0 ? -eval_jy(0, 1, x)
                                          : double(n)/x*value - eval_jy(0, n+1, x));
        require(std::isfinite(derivative) && std::abs(derivative) > 1e-15,
                "GSL J_n zero had vanishing derivative");
        const double step = value / derivative;
        require(std::isfinite(step) && std::abs(step) < 0.1 * x,
                "GSL J_n zero polishing was not bracket-safe");
        const double next = x - step;
        if (next == x) break;
        require(next > 0., "GSL J_n Newton refinement escaped positive axis");
        x = next;
    }
    require(std::abs(eval_jy(0, n, x)) <= 5e-10,
            "GSL J_n polished zero failed residual check");
    return x;
}

// kind: J=0, J'=1, Y=2, Y'=3.
struct Callback {
    int kind;
    int n;
    bool failed = false;
};

double function_value(int kind, int n, double x) {
    if (kind == 0 || kind == 2) return eval_jy(kind == 0 ? 0 : 1, n, x);
    const int family = kind == 1 ? 0 : 1;
    if (n == 0) return -eval_jy(family, 1, x);
    return double(n)/x * eval_jy(family, n, x) - eval_jy(family, n+1, x);
}

double gsl_callback(double x, void* ptr) {
    auto& context = *static_cast<Callback*>(ptr);
    try {
        return function_value(context.kind, context.n, x);
    } catch (...) {
        context.failed = true;
        return std::numeric_limits<double>::quiet_NaN();
    }
}

double brent(int kind, int n, double a, double b) {
    require(std::isfinite(a) && std::isfinite(b) && a > 0. && b > a,
            "invalid Bessel root bracket");
    const double left_val = function_value(kind, n, a);
    const double right_val = function_value(kind, n, b);
    if (left_val == 0.) return a;
    if (right_val == 0.) return b;
    require((left_val < 0.) != (right_val < 0.),
            "Bessel root bracket endpoints have same sign");
    Callback context{kind, n, false};
    gsl_function fn{&gsl_callback, &context};
    std::unique_ptr<gsl_root_fsolver, decltype(&gsl_root_fsolver_free)> solver(
        gsl_root_fsolver_alloc(gsl_root_fsolver_brent), &gsl_root_fsolver_free);
    require(static_cast<bool>(solver), "cannot allocate GSL Brent solver");
    require(gsl_root_fsolver_set(solver.get(), &fn, a, b) == GSL_SUCCESS &&
            !context.failed, "GSL Brent root setup failed");
    double root = 0.;
    bool converged = false;
    for (int iter = 0; iter < 100; ++iter) {
        const int status = gsl_root_fsolver_iterate(solver.get());
        require(status == GSL_SUCCESS && !context.failed,
                "GSL Brent root iteration failed");
        root = gsl_root_fsolver_root(solver.get());
        const double lo = gsl_root_fsolver_x_lower(solver.get());
        const double hi = gsl_root_fsolver_x_upper(solver.get());
        require(std::isfinite(root) && std::isfinite(lo) && std::isfinite(hi),
                "GSL Brent produced nonfinite root");
        if (gsl_root_test_interval(lo, hi, 5e-14, 4e-15) == GSL_SUCCESS) {
            converged = true;
            break;
        }
    }
    require(converged && root > 0., "GSL Brent failed to converge");
    require(std::abs(function_value(kind, n, root)) < 5e-9,
            "GSL Bessel root failed residual check");
    return root;
}

// Cache root computations across overlapping scalar/vector fixture calls.
// The counts are small, but 2,001 generated cases repeat few integer ranks.
double root(int kind, int n, int serial) {
    require(kind >= 0 && kind <= 3 && n >= 0 && n <= MAX_N &&
            serial > 0 && serial <= MAX_COUNT+1, "invalid Bessel root request");
    thread_local std::unordered_map<std::uint64_t, double> cache;
    const auto key = (std::uint64_t(kind) << 32) |
                     (std::uint64_t(n) << 16) | std::uint64_t(serial);
    if (auto found = cache.find(key); found != cache.end()) return found->second;
    double value = 0.;
    if (kind == 0) {
        value = zero_j(n, serial);
    } else if (kind == 1) {
        value = n == 0 ? root(0, 1, serial) :
            brent(kind, n, serial == 1 ? static_cast<double>(n)
                                      : root(0, n, serial-1),
                  root(0, n, serial));
    } else if (kind == 2) {
        const double lo = serial == 1 ? (n == 0 ? 0.1 : double(n))
                                       : root(0, n, serial-1);
        value = brent(kind, n, lo, root(0, n, serial));
    } else {
        value = n == 0 ? root(2, 1, serial) :
            brent(kind, n, root(2, n, serial), root(2, n, serial+1));
    }
    require(std::isfinite(value) && value > 0., "nonfinite positive Bessel root");
    if (cache.size() > 16000) cache.clear();
    cache.emplace(key, value);
    return value;
}

struct Entry { double x; int n, rank, derivative; };

std::vector<Entry> joint(int count) {
    // J'_0(0)=0 is included by the SciPy contract, with rank=0.
    std::vector<Entry> entries{{0.,0,0,1}};
    double cutoff = 5. + 3.*std::sqrt(double(count));
    for (int expansion = 0; expansion < 7; ++expansion) {
        entries.resize(1);
        const int upper_n = std::min(MAX_N, int(std::ceil(cutoff)));
        for (int n = 0; n <= upper_n; ++n) {
            for (int kind = 0; kind <= 1; ++kind) {
                for (int k = 1; k <= count; ++k) {
                    const double x = root(kind, n, k);
                    if (x > cutoff) break;
                    entries.push_back({x,n,k,kind});
                }
            }
        }
        if (int(entries.size()) >= count) break;
        cutoff *= 1.7;
    }
    require(int(entries.size()) >= count, "not enough GSL joint J/J' roots");
    std::sort(entries.begin(), entries.end(),
              [](const Entry& a, const Entry& b) {return a.x < b.x;});
    // Tie at J1=J0' must be in SciPy's J-first then J'-second ordering.
    for (std::size_t begin=0; begin<entries.size();) {
        std::size_t end=begin+1;
        while (end<entries.size() && std::abs(entries[end].x-entries[begin].x)<1e-11)
            ++end;
        std::stable_sort(entries.begin()+begin, entries.begin()+end,
                         [](const Entry& a, const Entry& b) {
                             return std::tie(a.derivative,a.n,a.rank) <
                                    std::tie(b.derivative,b.n,b.rank);
                         });
        begin=end;
    }
    entries.resize(count);
    return entries;
}

std::vector<double> evaluate(const char* name, const double* args, int ar) {
    const int count = count_checked(name,args,ar);
    const int n = ar == 2 ? static_cast<int>(args[0]) : 0;
    std::vector<double> out(static_cast<std::size_t>(count) * lanes(name));
    if (named(name, "jnjnp_zeros")) {
        const auto merged = joint(count);
        for (int k=0;k<count;++k) {
            out[k]=merged[k].x;
            out[count+k]=merged[k].n;
            out[2*count+k]=merged[k].rank;
            out[3*count+k]=merged[k].derivative;
        }
    } else if (named(name, "jnyn_zeros")) {
        for (int kind=0;kind<4;++kind)
            for (int k=1;k<=count;++k) out[kind*count+k-1]=root(kind,n,k);
    } else if (named(name, "jnp_zeros") || named(name, "yn_zeros") ||
               named(name, "ynp_zeros")) {
        const int kind = named(name, "jnp_zeros") ? 1 :
                         named(name, "yn_zeros") ? 2 : 3;
        for (int k=1;k<=count;++k) out[k-1]=root(kind,n,k);
    } else {
        const int kind = named(name, "y1p_zeros") ? 3 : 2;
        const int degree = named(name, "y0_zeros") ? 0 : 1;
        for (int k=1;k<=count;++k) {
            const double x=root(kind,degree,k);
            out[k-1]=x;
            // Actual SciPy's y0_zeros paired output is +Y1, NOT -Y1.
            out[count+k-1] = degree == 0 || kind == 3 ?
                function_value(2,1,x) : function_value(3,1,x);
        }
    }
    for (const double x : out)
        require(std::isfinite(x), "nonfinite GSL Bessel zero output");
    return out;
}

std::vector<double> cached(const char* name, const double* args, int ar) {
    const int count = count_checked(name,args,ar);
    const int n = ar == 2 ? static_cast<int>(args[0]) : 0;
    const auto key = std::string(name) + ':' + std::to_string(n) +
                     ':' + std::to_string(count);
    thread_local std::unordered_map<std::string,std::vector<double>> cache;
    if (auto found=cache.find(key);found!=cache.end())return found->second;
    auto value=evaluate(name,args,ar);
    if (cache.size() >= 256)cache.clear();
    cache.emplace(key,value);
    return value;
}
} // namespace

extern "C" {
const char* gsl_bessel_zeros_error() {return last_error.c_str();}
int gsl_bessel_zeros_supported(int library, const char* name) {
    return library == 1 ? arity(name) : 0;
}
int gsl_bessel_zeros_output_count(const char* name,const double* args,int ar) {
    try {return count_checked(name,args,ar)*lanes(name);}
    catch (...) {return 0;}
}
int gsl_bessel_zeros_evaluate(int library,const char* name,const double* args,
                              int ar,double* output,int capacity) {
    last_error.clear();
    // GSL otherwise aborts the process on a nonconvergent root bracket.
    static const bool no_gsl_abort = [] {gsl_set_error_handler_off();return true;}();
    (void)no_gsl_abort;
    try {
        require(gsl_bessel_zeros_supported(library,name) == ar && ar > 0,
                "unsupported GSL Bessel zeros library/API or arity");
        const int required=gsl_bessel_zeros_output_count(name,args,ar);
        require(required > 0 && output && capacity >= required,
                "invalid GSL Bessel zero input/capacity");
        const auto value=cached(name,args,ar);
        std::copy(value.begin(),value.end(),output); // failure never partially writes
        return 0;
    } catch (const std::exception& e) {last_error=e.what();}
      catch (...) {last_error="unknown GSL Bessel zero failure";}
    return -1;
}
}
