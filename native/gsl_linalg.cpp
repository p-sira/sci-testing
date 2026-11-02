// GSL 2.7-era public LU, BLAS, Cholesky, symmetric eigensystem and SVD.
// This independent bridge mirrors only the validated Sci saved signatures.
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
thread_local std::string last_error;
std::once_flag handler_once;
using Mat=std::unique_ptr<gsl_matrix,decltype(&gsl_matrix_free)>;
using Vec=std::unique_ptr<gsl_vector,decltype(&gsl_vector_free)>;
using Perm=std::unique_ptr<gsl_permutation,decltype(&gsl_permutation_free)>;
using EigWs=std::unique_ptr<gsl_eigen_symm_workspace,decltype(&gsl_eigen_symm_free)>;
bool same(const char* a,const char* b) { return a && std::strcmp(a,b)==0; }
void need(bool valid,const char* message) {
    if (!valid) throw std::invalid_argument(message);
}
void status(int code,const char* where) {
    if (code!=GSL_SUCCESS)
        throw std::runtime_error(std::string(where)+": "+gsl_strerror(code));
}
void finite(const double* values,int n) {
    need(n>=0 && n<=1000000 && (values || !n), "invalid array pointer or length");
    for (int i=0;i<n;++i) need(std::isfinite(values[i]),"nonfinite input");
}
Mat matrix(const double* values,int rows,int cols) {
    need(rows>0 && cols>0 && rows<=1024 && cols<=1024 &&
         static_cast<long long>(rows)*cols<=1000000, "invalid matrix dimensions");
    finite(values,rows*cols);
    Mat mat(gsl_matrix_alloc(rows,cols),gsl_matrix_free);
    need(bool(mat),"GSL matrix allocation failed");
    for (int i=0;i<rows;++i)
        for (int j=0;j<cols;++j) gsl_matrix_set(mat.get(),i,j,values[i*cols+j]);
    return mat;
}
Vec vector(const double* values,int n) {
    need(n>0 && n<=1000000,"invalid vector size");
    finite(values,n);
    Vec v(gsl_vector_alloc(n),gsl_vector_free);
    need(bool(v),"GSL vector allocation failed");
    for (int i=0;i<n;++i) gsl_vector_set(v.get(),i,values[i]);
    return v;
}
void output_matrix(const gsl_matrix* m,double* out) {
    for (size_t i=0;i<m->size1;++i)
        for (size_t j=0;j<m->size2;++j)
            out[i*m->size2+j]=gsl_matrix_get(m,i,j);
}
void output_vector(const gsl_vector* v,double* out) {
    for (size_t i=0;i<v->size;++i) out[i]=gsl_vector_get(v,i);
}
bool is_square(const char* name) {
    return same(name,"det") || same(name,"solve") || same(name,"inv") ||
           same(name,"cholesky") || same(name,"cho_solve") ||
           same(name,"eigvalsh");
}
int count(const char* name,int ar,int ac,int br,int bc,int nu,int nv) {
    if (!name || ar<=0 || ac<=0 || ar>1024 || ac>1024) return 0;
    if (is_square(name) && ar!=ac) return 0;
    if (same(name,"det") || same(name,"dot")) return same(name,"dot") ?
        ((nu>0 && nu==nv) ? 1 : 0) : 1;
    if (same(name,"solve") || same(name,"cho_solve"))
        return br==ar && bc==1 ? ar : 0;
    if (same(name,"inv") || same(name,"cholesky"))
        return ar*ac;
    if (same(name,"eigvalsh")) return ar;
    if (same(name,"svdvals")) return std::min(ar,ac);
    if (same(name,"matvec")) return nu==ac ? ar : 0;
    if (same(name,"matmul")) return br==ac && bc>0 ? ar*bc : 0;
    return 0;
}
void evaluate(const char* name,const double* a,int ar,int ac,
              const double* b,int br,int bc,const double* u,int nu,
              const double* v,int nv,int lower,double* out) {
    if (same(name,"dot")) {
        finite(u,nu); finite(v,nv);
        Vec x=vector(u,nu),y=vector(v,nv);
        status(gsl_blas_ddot(x.get(),y.get(),out),"gsl_blas_ddot");
        return;
    }
    Mat A=matrix(a,ar,ac);
    if (same(name,"matvec")) {
        Vec x=vector(u,nu);
        Vec result(gsl_vector_alloc(ar),gsl_vector_free);
        need(bool(result),"GSL result vector allocation failed");
        gsl_vector_set_zero(result.get());
        status(gsl_blas_dgemv(CblasNoTrans,1.,A.get(),x.get(),0.,
                              result.get()),"gsl_blas_dgemv");
        output_vector(result.get(),out);
        return;
    }
    if (same(name,"matmul")) {
        Mat B=matrix(b,br,bc);
        Mat product(gsl_matrix_alloc(ar,bc),gsl_matrix_free);
        need(bool(product),"GSL product matrix allocation failed");
        gsl_matrix_set_zero(product.get());
        status(gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,A.get(),
                              B.get(),0.,product.get()),"gsl_blas_dgemm");
        output_matrix(product.get(),out);
        return;
    }
    if (same(name,"det") || same(name,"solve") || same(name,"inv")) {
        Perm perm(gsl_permutation_alloc(ar),gsl_permutation_free);
        need(bool(perm),"GSL permutation allocation failed");
        int sign=1;
        status(gsl_linalg_LU_decomp(A.get(),perm.get(),&sign),"gsl_linalg_LU_decomp");
        if (same(name,"det")) { out[0]=gsl_linalg_LU_det(A.get(),sign); return; }
        // Avoid GSL's built-in error handler on singular matrices.
        for (int i=0;i<ar;++i)
            need(gsl_matrix_get(A.get(),i,i)!=0., "singular LU decomposition");
        if (same(name,"inv")) {
            Mat inverse(gsl_matrix_alloc(ar,ar),gsl_matrix_free);
            need(bool(inverse),"GSL inverse allocation failed");
            status(gsl_linalg_LU_invert(A.get(),perm.get(),inverse.get()),
                   "gsl_linalg_LU_invert");
            output_matrix(inverse.get(),out);
        } else {
            Vec rhs=vector(b,br), solution(gsl_vector_alloc(ar),gsl_vector_free);
            need(bool(solution),"GSL solution allocation failed");
            status(gsl_linalg_LU_solve(A.get(),perm.get(),rhs.get(),
                                       solution.get()),"gsl_linalg_LU_solve");
            output_vector(solution.get(),out);
        }
        return;
    }
    if (same(name,"cholesky") || same(name,"cho_solve")) {
        need(lower==0 || lower==1,"lower flag must be zero or one");
        status(gsl_linalg_cholesky_decomp(A.get()),"gsl_linalg_cholesky_decomp");
        if (same(name,"cholesky")) {
            // GSL returns L; SciPy expects L for lower=True and L^T otherwise.
            for (int i=0;i<ar;++i)
                for (int j=0;j<ar;++j)
                    out[i*ar+j]=(lower ? j<=i : i<=j) ?
                         (lower ? gsl_matrix_get(A.get(),i,j) :
                                  gsl_matrix_get(A.get(),j,i)) : 0.;
        } else {
            Vec rhs=vector(b,br),solution(gsl_vector_alloc(ar),gsl_vector_free);
            need(bool(solution),"GSL solution allocation failed");
            status(gsl_linalg_cholesky_solve(A.get(),rhs.get(),solution.get()),
                   "gsl_linalg_cholesky_solve");
            output_vector(solution.get(),out);
        }
        return;
    }
    if (same(name,"eigvalsh")) {
        need(lower==0 || lower==1,"lower flag must be zero or one");
        // GSL eigen_symm reads the lower half. Respect SciPy's lower flag
        // even for saved intentionally asymmetric matrices.
        for (int i=0;i<ar;++i)
            for (int j=0;j<i;++j) {
                const double element=lower ? gsl_matrix_get(A.get(),i,j)
                                           : gsl_matrix_get(A.get(),j,i);
                gsl_matrix_set(A.get(),i,j,element);
                gsl_matrix_set(A.get(),j,i,element);
            }
        EigWs workspace(gsl_eigen_symm_alloc(ar),gsl_eigen_symm_free);
        Vec values(gsl_vector_alloc(ar),gsl_vector_free);
        need(bool(workspace) && bool(values),"GSL eigensystem allocation failed");
        status(gsl_eigen_symm(A.get(),values.get(),workspace.get()),
               "gsl_eigen_symm");
        output_vector(values.get(),out);
        std::sort(out,out+ar);
        return;
    }
    if (same(name,"svdvals")) {
        // gsl_linalg_SV_decomp requires rows >= columns; transpose wide A.
        if (ar<ac) {
            Mat transposed(gsl_matrix_alloc(ac,ar),gsl_matrix_free);
            need(bool(transposed),"GSL transpose allocation failed");
            gsl_matrix_transpose_memcpy(transposed.get(),A.get());
            A=std::move(transposed);
        }
        const int columns=std::min(ar,ac);
        Mat V(gsl_matrix_alloc(columns,columns),gsl_matrix_free);
        Vec S(gsl_vector_alloc(columns),gsl_vector_free);
        Vec work(gsl_vector_alloc(columns),gsl_vector_free);
        need(bool(V) && bool(S) && bool(work),"GSL SVD allocation failed");
        status(gsl_linalg_SV_decomp(A.get(),V.get(),S.get(),work.get()),
               "gsl_linalg_SV_decomp");
        output_vector(S.get(),out);
        std::sort(out,out+columns,std::greater<double>());
        return;
    }
    throw std::invalid_argument("unsupported GSL linalg API");
}
}

extern "C" {
const char* gsl_linalg_error() { return last_error.c_str(); }
int gsl_linalg_supported(int library,const char* name) {
    if (library!=1 || !name) return 0;
    return same(name,"det") || same(name,"solve") || same(name,"inv") ||
           same(name,"dot") || same(name,"matvec") || same(name,"matmul") ||
           same(name,"cholesky") || same(name,"cho_solve") ||
           same(name,"eigvalsh") || same(name,"svdvals");
}
int gsl_linalg_output_count(const char* name,int ar,int ac,
                            int br,int bc,int nu,int nv) {
    return count(name,ar,ac,br,bc,nu,nv);
}
int gsl_linalg_evaluate(int library,const char* name,const double* a,int ar,int ac,
                         const double* b,int br,int bc,const double* u,int nu,
                         const double* v,int nv,int lower,double* out,int capacity) {
    last_error.clear();
    try {
        std::call_once(handler_once,[]{gsl_set_error_handler_off();});
        need(gsl_linalg_supported(library,name),"unsupported GSL linalg API/provider");
        const int size=count(name,ar,ac,br,bc,nu,nv);
        need(size>0 && size<=1000000 && out && capacity>=size,
             "invalid GSL linalg dimensions or output capacity");
        finite(a,ar*ac);
        finite(b,br*bc);
        finite(u,nu);
        finite(v,nv);
        evaluate(name,a,ar,ac,b,br,bc,u,nu,v,nv,lower,out);
        for (int i=0;i<size;++i)
            need(std::isfinite(out[i]),"GSL linalg produced nonfinite output");
        return 0;
    } catch (const std::exception& ex) {
        last_error=ex.what();
    } catch (...) { last_error="unknown GSL linalg exception"; }
    return -1;
}
}
