// Public GSL 2.7 sparse COO/CSR and sparse BLAS comparison adapter.
//
// Wire format: [rows, cols, nnz, (row,col,value)*nnz,
//               (brows,bcols,bnnz,(row,col,value)*bnnz for add/matmul),
//               (x[cols] for matvec)]
// C output: [rows, cols, dense C-order values...]. The Python wrapper
// subsequently adds Sci's finite-value (class,value) tags.
#include <cstddef> // GSL 2.8 gsl_spmatrix.h uses size_t without including stddef.h.
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr int maxdim=128;
constexpr int max_entries=32768;
constexpr int max_input=160000;
thread_local std::string last_error;
std::once_flag disable_handler_once;
using Sp=std::unique_ptr<gsl_spmatrix,decltype(&gsl_spmatrix_free)>;
using Dense=std::unique_ptr<gsl_matrix,decltype(&gsl_matrix_free)>;
using Vector=std::unique_ptr<gsl_vector,decltype(&gsl_vector_free)>;

bool eq(const char* a,const char* b) {
    return a && std::strcmp(a,b)==0;
}
void require(bool condition,const char* message) {
    if (!condition) throw std::invalid_argument(message);
}
void check(int status,const char* procedure) {
    if (status!=GSL_SUCCESS)
        throw std::runtime_error(std::string(procedure)+": "+gsl_strerror(status));
}
bool supported(int provider,const char* name) {
    return provider==1 && (
        eq(name,"sparse.matvec") || eq(name,"sparse.matrix_transpose") ||
        eq(name,"sparse.transpose") || eq(name,"sparse.from_dense") ||
        eq(name,"sparse.to_dense") || eq(name,"sparse.add") ||
        eq(name,"sparse.matmul"));
}
class Reader {
public:
    Reader(const double* values,int size) : data(values), length(size) {
        require(size>=0 && size<=max_input && (values || size==0),
                "invalid input pointer or size");
    }
    double real() {
        require(pos<length,"truncated GSL sparse wire input");
        double value=data[pos++];
        require(std::isfinite(value),"nonfinite sparse matrix input");
        return value;
    }
    int integer(int bound) {
        const double v=real();
        require(v>=0 && v<=bound && std::floor(v)==v,
                "invalid sparse index or dimension");
        return static_cast<int>(v);
    }
    void done() { require(pos==length,"trailing sparse wire input"); }
private:
    const double* data;
    const int length;
    int pos=0;
};
struct Matrix {
    int rows;
    int cols;
    Sp sparse;
};
Sp alloc(int rows,int cols,int count,int type) {
    require(rows>0 && cols>0 && rows<=maxdim && cols<=maxdim &&
            rows*cols<=max_entries && count>=0 && count<=max_entries,
            "invalid sparse matrix dimensions or capacity");
    Sp result(gsl_spmatrix_alloc_nzmax(rows,cols,std::max(1,count),type),
              gsl_spmatrix_free);
    require(bool(result),"GSL sparse allocation failed");
    return result;
}
Dense dense(int rows,int cols) {
    require(rows>0 && cols>0 && rows*cols<=max_entries,
            "invalid dense conversion shape");
    Dense output(gsl_matrix_calloc(rows,cols),gsl_matrix_free);
    require(bool(output),"GSL dense matrix allocation failed");
    return output;
}
Sp column_major(const Matrix& source) {
    // gsl_spmatrix_compress requires a COO source; converting an already
    // compressed CSR matrix directly is invalid and can corrupt its indices.
    auto values=dense(source.rows,source.cols);
    check(gsl_spmatrix_sp2d(values.get(),source.sparse.get()),
          "gsl_spmatrix_sp2d(CSR to dense)");
    auto coo=alloc(source.rows,source.cols,source.rows*source.cols,
                   GSL_SPMATRIX_COO);
    check(gsl_spmatrix_d2sp(coo.get(),values.get()),
          "gsl_spmatrix_d2sp(dense to COO)");
    Sp csc(gsl_spmatrix_compress(coo.get(),GSL_SPMATRIX_CSC),
           gsl_spmatrix_free);
    require(bool(csc),"GSL COO->CSC compression failed");
    return csc;
}
Matrix parse(Reader& reader) {
    const int rows=reader.integer(maxdim);
    const int cols=reader.integer(maxdim);
    const int triples=reader.integer(max_entries);
    auto coo=alloc(rows,cols,triples,GSL_SPMATRIX_COO);
    for (int k=0;k<triples;++k) {
        const int i=reader.integer(maxdim);
        const int j=reader.integer(maxdim);
        const double v=reader.real();
        require(i<rows && j<cols,"out-of-bounds sparse COO index");
        // Unlike gsl_spmatrix_set's replacement semantics, SciPy's
        // coo_matrix(...).tocsr() coalesces repeated coordinates by summing.
        // Explicit zeros remain numerically zero after CSR conversion.
        const double old=gsl_spmatrix_get(coo.get(),i,j);
        require(std::isfinite(old+v),"nonfinite duplicate sum");
        check(gsl_spmatrix_set(coo.get(),i,j,old+v),
              "gsl_spmatrix_set(coalesced COO)");
    }
    Sp csr(gsl_spmatrix_compress(coo.get(),GSL_SPMATRIX_CSR),
           gsl_spmatrix_free);
    require(bool(csr),"GSL COO to CSR compression failed");
    return {rows,cols,std::move(csr)};
}
void write_dense(const gsl_matrix* source,std::vector<double>& values) {
    values.reserve(2+source->size1*source->size2);
    values.push_back(static_cast<double>(source->size1));
    values.push_back(static_cast<double>(source->size2));
    for (size_t i=0;i<source->size1;++i)
        for (size_t j=0;j<source->size2;++j)
            values.push_back(gsl_matrix_get(source,i,j));
}
std::vector<double> evaluate(const char* name,const double* input,int input_len) {
    Reader reader(input,input_len);
    Matrix a=parse(reader);
    if (eq(name,"sparse.matvec")) {
        Vector x(gsl_vector_alloc(a.cols),gsl_vector_free);
        Vector y(gsl_vector_calloc(a.rows),gsl_vector_free);
        require(bool(x) && bool(y),"GSL sparse matvec allocation failed");
        for (int i=0;i<a.cols;++i) gsl_vector_set(x.get(),i,reader.real());
        reader.done();
        check(gsl_spblas_dgemv(CblasNoTrans,1.,a.sparse.get(),
                               x.get(),0.,y.get()),"gsl_spblas_dgemv");
        std::vector<double> result{1.,static_cast<double>(a.rows)};
        for (int i=0;i<a.rows;++i) result.push_back(gsl_vector_get(y.get(),i));
        return result;
    }
    if (eq(name,"sparse.transpose") || eq(name,"sparse.matrix_transpose")) {
        reader.done();
        auto dest=alloc(a.cols,a.rows,a.sparse->nz,GSL_SPMATRIX_CSR);
        check(gsl_spmatrix_transpose_memcpy(dest.get(),a.sparse.get()),
              "gsl_spmatrix_transpose_memcpy");
        auto output=dense(a.cols,a.rows);
        check(gsl_spmatrix_sp2d(output.get(),dest.get()),"gsl_spmatrix_sp2d");
        std::vector<double> values;
        write_dense(output.get(),values);
        return values;
    }
    if (eq(name,"sparse.from_dense")) {
        reader.done();
        auto matrix=dense(a.rows,a.cols);
        check(gsl_spmatrix_sp2d(matrix.get(),a.sparse.get()),"gsl_spmatrix_sp2d");
        auto coo=alloc(a.rows,a.cols,a.rows*a.cols,GSL_SPMATRIX_COO);
        check(gsl_spmatrix_d2sp(coo.get(),matrix.get()),"gsl_spmatrix_d2sp");
        Sp converted(gsl_spmatrix_compress(coo.get(),GSL_SPMATRIX_CSR),
                     gsl_spmatrix_free);
        require(bool(converted),"GSL dense->sparse CSR compression failed");
        auto output=dense(a.rows,a.cols);
        check(gsl_spmatrix_sp2d(output.get(),converted.get()),"gsl_spmatrix_sp2d");
        std::vector<double> values;
        write_dense(output.get(),values);
        return values;
    }
    if (eq(name,"sparse.to_dense")) {
        reader.done();
        auto output=dense(a.rows,a.cols);
        check(gsl_spmatrix_sp2d(output.get(),a.sparse.get()),"gsl_spmatrix_sp2d");
        std::vector<double> values;
        write_dense(output.get(),values);
        return values;
    }
    Matrix b=parse(reader);
    reader.done();
    const bool multiplication=eq(name,"sparse.matmul");
    require(multiplication ? a.cols==b.rows :
            (a.rows==b.rows && a.cols==b.cols),
            "incompatible sparse operand shapes");
    const int result_cols=multiplication ? b.cols : a.cols;
    // gsl_spblas_dgemm is the public CSparse/CCS algorithm: its 2.7/2.8
    // contract expressly requires all THREE matrices to be compressed-column
    // (CSC/CCS), unlike gsl_spblas_dgemv which accepts CSR.
    const int storage=multiplication ? GSL_SPMATRIX_CSC : GSL_SPMATRIX_CSR;
    auto result=alloc(a.rows,result_cols,a.rows*result_cols,storage);
    if (multiplication) {
        Sp column_a=column_major(a);
        Sp column_b=column_major(b);
        require(bool(column_a) && bool(column_b),"GSL CSR->CSC conversion failed");
        check(gsl_spblas_dgemm(1.,column_a.get(),column_b.get(),result.get()),
              "gsl_spblas_dgemm");
    } else {
        check(gsl_spmatrix_add(result.get(),a.sparse.get(),b.sparse.get()),
              "gsl_spmatrix_add");
    }
    auto output=dense(a.rows,result_cols);
    check(gsl_spmatrix_sp2d(output.get(),result.get()),"gsl_spmatrix_sp2d");
    std::vector<double> values;
    write_dense(output.get(),values);
    return values;
}
int output_size(const char* name,const double* input,int n) {
    Reader reader(input,n);
    const int rows=reader.integer(maxdim),cols=reader.integer(maxdim);
    // The rest of the encoded COO is checked during evaluation.
    if (rows<=0 || cols<=0 || rows*cols>max_entries) return 0;
    if (eq(name,"sparse.matvec")) return 2+rows;
    if (eq(name,"sparse.matrix_transpose") || eq(name,"sparse.transpose"))
        return 2+rows*cols;
    if (eq(name,"sparse.from_dense") || eq(name,"sparse.to_dense") ||
        eq(name,"sparse.add")) return 2+rows*cols;
    // Dimensions of B are only known after reading all entries of A.
    if (eq(name,"sparse.matmul")) {
        int count=reader.integer(max_entries);
        for (int k=0;k<count;++k) {
            reader.integer(maxdim);reader.integer(maxdim);reader.real();
        }
        const int brows=reader.integer(maxdim),bcols=reader.integer(maxdim);
        if (brows!=cols || bcols<=0 || rows*bcols>max_entries) return 0;
        return 2+rows*bcols;
    }
    return 0;
}
} // namespace

extern "C" {
int gsl_sparse_supported(int provider,const char* name) {
    return supported(provider,name) ? 1 : 0;
}
const char* gsl_sparse_error() { return last_error.c_str(); }
int gsl_sparse_output_count(const char* name,const double* input,int n) {
    try {
        return supported(1,name) ? output_size(name,input,n) : 0;
    } catch (...) { return 0; }
}
int gsl_sparse_evaluate(int provider,const char* name,const double* input,
                        int n,double* output,int capacity) {
    last_error.clear();
    try {
        std::call_once(disable_handler_once,[]{gsl_set_error_handler_off();});
        require(supported(provider,name),"unsupported GSL sparse operation/provider");
        require(output,"null GSL sparse output pointer");
        auto values=evaluate(name,input,n);
        require(static_cast<int>(values.size())<=capacity,
                "insufficient GSL sparse output capacity");
        for (double value:values)
            require(std::isfinite(value),"nonfinite GSL sparse result");
        std::copy(values.begin(),values.end(),output);
        return 0;
    } catch (const std::exception& ex) {
        last_error=ex.what();
    } catch (...) {
        last_error="unknown GSL sparse adapter exception";
    }
    return -1;
}
}
