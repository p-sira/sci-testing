// GSL-derived nonlinear minimization of the SAVED explicit 2-D callback.
//
// Objective: f(x,y) = 3 + (x-1)^2 + (y+2)^2.
// Gradient:  (2(x-1), 2(y+2)).
//
// These are actual public GSL 2.7 multimin methods: Nelder-Mead simplex2
// (derivative-free) and vector_bfgs2 (analytic gradient). They are independent
// minimizers, NOT shortcuts that simply return the analytic minimizer.
// SciPy and GSL iterate differently; the saved output contract is ONLY
// [optimized x, optimized y, objective], with no trajectory or status fields.
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
thread_local std::string last_error;
constexpr int library_gsl = 1;
constexpr int input_count = 2;
constexpr int output_count = 3;
constexpr int maximum_iterations = 512;
constexpr double tolerance = 1e-8;
constexpr double maximum_start = 1000.;

bool named(const char* value,const char* other) {
    return value && std::strcmp(value,other)==0;
}
bool known(const char* name) {
    return named(name,"minimize") || named(name,"try_minimize") ||
           named(name,"minimize_with_gradient") ||
           named(name,"try_minimize_with_gradient");
}
bool gradient(const char* name) {
    return named(name,"minimize_with_gradient") ||
           named(name,"try_minimize_with_gradient");
}
void require(bool yes,const std::string& message) {
    if (!yes)throw std::domain_error(message);
}
void gsl_checked(int status,const char* operation) {
    if (status!=GSL_SUCCESS)
        throw std::runtime_error(std::string(operation)+": "+gsl_strerror(status));
}

double objective(const gsl_vector* point,void*) {
    const double dx=gsl_vector_get(point,0)-1.;
    const double dy=gsl_vector_get(point,1)+2.;
    // The additive constant is restored in the public output. GSL simplex2
    // ranks objective values near the minimizer: adding 3 inside the solver
    // would round all sufficiently small residuals to exactly 3 and stall
    // contraction. Minimizing f-3 is mathematically equivalent, but avoids
    // this preventable binary64 loss of ordering information.
    return dx*dx+dy*dy;
}
void objective_gradient(const gsl_vector* point,void*,gsl_vector* output) {
    gsl_vector_set(output,0,2.*(gsl_vector_get(point,0)-1.));
    gsl_vector_set(output,1,2.*(gsl_vector_get(point,1)+2.));
}
void objective_and_gradient(const gsl_vector* point,void* context,
                            double* f,gsl_vector* output) {
    *f=objective(point,context);
    objective_gradient(point,context,output);
}

using vector_ptr=std::unique_ptr<gsl_vector,decltype(&gsl_vector_free)>;
using simplex_ptr=std::unique_ptr<gsl_multimin_fminimizer,
                                  decltype(&gsl_multimin_fminimizer_free)>;
using bfgs_ptr=std::unique_ptr<gsl_multimin_fdfminimizer,
                               decltype(&gsl_multimin_fdfminimizer_free)>;
vector_ptr make_vector() {
    vector_ptr ptr(gsl_vector_alloc(2),&gsl_vector_free);
    require(static_cast<bool>(ptr),"GSL vector allocation failed");
    return ptr;
}

std::vector<double> minimize(const double* start,bool with_gradient) {
    vector_ptr location=make_vector();
    gsl_vector_set(location.get(),0,start[0]);
    gsl_vector_set(location.get(),1,start[1]);

    double x=0.,y=0.;
    bool converged=false;
    if (with_gradient) {
        gsl_multimin_function_fdf callback{
            &objective,&objective_gradient,&objective_and_gradient,2,nullptr};
        bfgs_ptr minimizer(gsl_multimin_fdfminimizer_alloc(
            gsl_multimin_fdfminimizer_vector_bfgs2,2),
            &gsl_multimin_fdfminimizer_free);
        require(static_cast<bool>(minimizer),"GSL BFGS2 allocation failed");
        // Initial line-search step and line-search tolerance, not the
        // convergence accuracy; convergence is based on gradient norm.
        gsl_checked(gsl_multimin_fdfminimizer_set(
            minimizer.get(),&callback,location.get(),0.1,0.1),
            "GSL BFGS2 initialization");
        for (int i=0;i<=maximum_iterations;++i) {
            const int status=gsl_multimin_test_gradient(
                gsl_multimin_fdfminimizer_gradient(minimizer.get()),tolerance);
            if (status==GSL_SUCCESS) {
                converged=true;
                break;
            }
            require(status==GSL_CONTINUE,"GSL BFGS2 gradient test failed");
            if (i==maximum_iterations)break;
            gsl_checked(gsl_multimin_fdfminimizer_iterate(minimizer.get()),
                        "GSL BFGS2 iteration");
        }
        if (converged) {
            const gsl_vector* solution=gsl_multimin_fdfminimizer_x(minimizer.get());
            x=gsl_vector_get(solution,0);
            y=gsl_vector_get(solution,1);
        }
    } else {
        gsl_multimin_function callback{&objective,2,nullptr};
        simplex_ptr minimizer(gsl_multimin_fminimizer_alloc(
            gsl_multimin_fminimizer_nmsimplex2,2),
            &gsl_multimin_fminimizer_free);
        require(static_cast<bool>(minimizer),"GSL simplex2 allocation failed");
        vector_ptr steps=make_vector();
        // Nondegenerate positive simplex sides; scale with the local
        // starting magnitudes without assuming the known exact minimizer.
        for (int i=0;i<2;++i)
            gsl_vector_set(steps.get(),i,0.1*std::max(1.,std::abs(start[i])));
        gsl_checked(gsl_multimin_fminimizer_set(
            minimizer.get(),&callback,location.get(),steps.get()),
            "GSL simplex2 initialization");
        for (int i=0;i<=maximum_iterations;++i) {
            const double size=gsl_multimin_fminimizer_size(minimizer.get());
            require(std::isfinite(size),"GSL simplex2 size nonfinite");
            // The public solver's cached minimum() and x() need not be
            // treated as a single atomic (position,cost) observation.
            // Re-evaluate the objective at the CURRENT public position;
            // this avoids accepting/rejecting based on a stale cached
            // simplex cost if GSL has updated its geometry in between.
            const gsl_vector* current=gsl_multimin_fminimizer_x(minimizer.get());
            require(current!=nullptr,"GSL simplex2 position unavailable");
            const double residual=objective(current,nullptr);
            require(std::isfinite(residual) && residual>=0.,
                    "GSL simplex2 returned invalid residual");
            const int status=gsl_multimin_test_size(size,tolerance);
            // For this strictly convex saved quadratic, an objective
            // residual <= tolerance^2 directly bounds Euclidean distance
            // to the unique optimum by tolerance. This independent
            // convergence check permits genuine high-accuracy solutions
            // even if simplex geometry stalls at machine precision.
            if (status==GSL_SUCCESS || residual<=tolerance*tolerance) {
                converged=true;
                break;
            }
            require(status==GSL_CONTINUE,"GSL simplex2 convergence test failed");
            if (i==maximum_iterations)break;
            gsl_checked(gsl_multimin_fminimizer_iterate(minimizer.get()),
                        "GSL simplex2 iteration");
        }
        if (converged) {
            const gsl_vector* solution=gsl_multimin_fminimizer_x(minimizer.get());
            x=gsl_vector_get(solution,0);
            y=gsl_vector_get(solution,1);
        }
    }
    require(converged,"GSL multimin did not converge in 512 iterations");
    require(std::isfinite(x) && std::isfinite(y),
            "GSL multimin produced a nonfinite solution");
    const double dx=x-1.,dy=y+2.;
    // Convergence is independently checked against the known strictly
    // convex quadratic optimum, not inferred from a successful GSL status.
    require(std::abs(dx)<=1e-6 && std::abs(dy)<=1e-6,
            "GSL multimin converged away from analytic minimizer");
    // The SAVED output's cost is f(REPORTED x,y). A solver's internally
    // cached minimum is not an output oracle: even a finite optimizer
    // cache can describe an earlier simplex vertex rather than x(). Compute
    // the objective from the final coordinates after independent validation.
    const double f=3.+dx*dx+dy*dy;
    require(std::isfinite(f),"nonfinite independently evaluated objective");
    return {x,y,f};
}
}

extern "C" {
const char* gsl_optimize_error() {return last_error.c_str();}
int gsl_optimize_supported(int library,const char* name) {
    return library==library_gsl && known(name) ? input_count : 0;
}
int gsl_optimize_output_count(const char* name,const double* values,int count) {
    if (!known(name) || !values || count!=input_count)return 0;
    if (!std::isfinite(values[0]) || !std::isfinite(values[1]) ||
        std::abs(values[0])>maximum_start ||
        std::abs(values[1])>maximum_start)return 0;
    return output_count;
}
int gsl_optimize_evaluate(int library,const char* name,
                          const double* inputs,int count,
                          double* output,int capacity) {
    last_error.clear();
    // GSL normally invokes an aborting error handler. Unreliable status
    // paths must always return an error instead of terminating the process.
    static const bool safe=[] {gsl_set_error_handler_off();return true;}();
    (void)safe;
    try {
        require(gsl_optimize_supported(library,name)==input_count,
                "unsupported GSL optimize provider/function");
        require(gsl_optimize_output_count(name,inputs,count)==output_count &&
                output && capacity>=output_count,
                "GSL optimize requires finite 2D start in [-1000,1000] and 3 outputs");
        const auto result=minimize(inputs,gradient(name));
        std::copy(result.begin(),result.end(),output); // no partial writes
        return 0;
    } catch (const std::exception& error) {last_error=error.what();}
      catch (...) {last_error="unknown GSL optimize failure";}
    return -1;
}
}
