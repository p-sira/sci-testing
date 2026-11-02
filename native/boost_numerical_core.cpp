// Public Boost.Math numerical oracles for saved Sci callbacks.
// Input by name: fixed_quad=[a,b]; derivative=[x,step];
// gradient=[x,y,step]; jacobian/hessian=[x,y].
// The saved step is validated but Boost determines its own step (or uses
// automatic differentiation). This is an independent mathematical oracle,
// not a promise of identical Sci finite-difference methodology.
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/differentiation/autodiff.hpp>
#include <boost/math/differentiation/finite_difference.hpp>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <string>
#include <tuple>

namespace {
thread_local std::string last_error;
bool named(const char* name, const char* match) {
    return name && std::strcmp(name, match) == 0;
}
int input_arity(const char* name) {
    if (named(name,"fixed_quad") || named(name,"derivative") ||
        named(name,"jacobian") || named(name,"hessian")) return 2;
    if (named(name,"gradient")) return 3;
    return 0;
}
int output_arity(const char* name) {
    if (named(name,"gradient")) return 2;
    if (named(name,"jacobian") || named(name,"hessian")) return 4;
    return input_arity(name) ? 1 : 0;
}
void require(bool yes, const char* reason) {
    if (!yes) throw std::invalid_argument(reason);
}
void validate(const char* name, const double* args, int count) {
    require(args && input_arity(name) == count && count>0,
            "unsupported numerical function or input arity");
    for (int i=0; i<count; ++i)
        require(std::isfinite(args[i]), "Boost numerical inputs must be finite");
    if (named(name,"derivative"))
        require(args[1]>0,"derivative step must be positive");
    if (named(name,"gradient"))
        require(args[2]>0,"gradient step must be positive");
}
double sine(double x) { return std::sin(x); }
void calculate(const char* name, const double* args, double* out) {
    if (named(name,"fixed_quad")) {
        // Sci fixture: fixed_quad(exp, a, b, n=16), not adaptive quadrature.
        out[0]=boost::math::quadrature::gauss<double,16>::integrate(
            [](double x){return std::exp(x);},args[0],args[1]);
    } else if (named(name,"derivative")) {
        // Public sixth-order Boost finite difference: f(x)=sin(x).
        out[0]=boost::math::differentiation::
            finite_difference_derivative<decltype(&sine),double,6>(&sine,args[0]);
    } else {
        // Public forward-mode autodiff, two independent variables.
        auto variables=boost::math::differentiation::
            make_ftuple<double,2,2>(args[0],args[1]);
        auto x=std::get<0>(variables);
        auto y=std::get<1>(variables);
        if (named(name,"gradient")) {
            auto f=x*x+3*x*y+2*y*y;
            out[0]=f.derivative(1,0);
            out[1]=f.derivative(0,1);
        } else if (named(name,"jacobian")) {
            auto f=x*y;
            auto g=x+y;
            // Row-major: [f_x,f_y,g_x,g_y].
            out[0]=f.derivative(1,0);
            out[1]=f.derivative(0,1);
            out[2]=g.derivative(1,0);
            out[3]=g.derivative(0,1);
        } else if (named(name,"hessian")) {
            auto f=x*x+3*x*y+2*y*y;
            out[0]=f.derivative(2,0);
            out[1]=f.derivative(1,1);
            out[2]=f.derivative(1,1);
            out[3]=f.derivative(0,2);
        }
    }
    for (int i=0;i<output_arity(name);++i)
        require(std::isfinite(out[i]),"Boost numerical output must be finite");
}
}

extern "C" {
const char* boost_numerical_core_error() {return last_error.c_str();}
int boost_numerical_core_supported(int library,const char* name) {
    return library==0 ? input_arity(name) : 0;
}
int boost_numerical_core_output_count(
        const char* name,const double* args,int count) {
    try {validate(name,args,count);return output_arity(name);}
    catch (...) {return 0;}
}
int boost_numerical_core_evaluate(
        int library,const char* name,const double* args,int count,
        double* out,int capacity) {
    last_error.clear();
    try {
        require(library==0,"Boost numerical core supports Boost.Math only");
        validate(name,args,count);
        require(out && capacity>=output_arity(name),
                "invalid Boost numerical output pointer/capacity");
        double result[4]={};
        calculate(name,args,result);
        std::copy(result,result+output_arity(name),out);
        return 0;
    } catch(const std::exception& e) {last_error=e.what();}
      catch(...) {last_error="unknown Boost numerical core failure";}
    return -1;
}
}
