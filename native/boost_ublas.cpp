// Public Boost.Numeric.uBLAS dense/sparse linear algebra for Sci saved cases.
//
// Input conventions (all double, integer dimensions are exact):
// Dense:
//   dot:         [n, u[n], v[n]]
//   transpose:   [rows, cols, a[rows*cols]]
//   matvec:      [rows, cols, a[rows*cols], u[cols]]
//   matmul:      [rows, shared, cols, a[rows*shared], b[shared*cols]]
//   det, solve:  [n, n, a[n*n], (b[n] for solve)]
// Sparse:
//   [rows, cols, nnz, (row,col,value)*nnz,
//      (brows,bcols,bnnz,(row,col,value)*bnnz for add/multiply/matmul),
//      (x[cols] for matvec)]
//
// Results: dense row-major data (scalar for dot/det), sparse [rows,cols,
// row-major data]. Python applies the Sci sparse finite-value tagging exactly.
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace ub = boost::numeric::ublas;
using Matrix = ub::matrix<double, ub::row_major>;
using Vector = ub::vector<double>;
using Sparse = ub::compressed_matrix<double, ub::row_major>;
using ExtendedMatrix = ub::matrix<long double, ub::row_major>;
using ExtendedVector = ub::vector<long double>;
thread_local std::string last_error;
constexpr int bound = 128;
constexpr int matrix_bound = 32768;

bool named(const char* name, const char* expected) {
    return name && std::strcmp(name, expected) == 0;
}

bool supported(int library, const char* name) {
    if (library != 0) return false;
    return named(name, "linalg.det") || named(name, "linalg.solve") ||
        named(name, "linalg.dot") || named(name, "linalg.matvec") ||
        named(name, "linalg.matmul") || named(name, "linalg.transpose") ||
        named(name, "sparse.add") || named(name, "sparse.multiply") ||
        named(name, "sparse.matvec") || named(name, "sparse.matmul") ||
        named(name, "sparse.transpose") || named(name, "sparse.from_dense") ||
        named(name, "sparse.to_dense");
}

class Reader {
public:
    Reader(const double* values, int size) : data(values), length(size) {
        if (!values || size < 0) throw std::invalid_argument("invalid input buffer");
    }
    double real() {
        if (index >= length) throw std::invalid_argument("truncated uBLAS input");
        const double v = data[index++];
        if (!std::isfinite(v)) throw std::invalid_argument("uBLAS input must be finite");
        return v;
    }
    int number(int upper = bound) {
        const double v = real();
        if (v < 0 || v > upper || std::trunc(v) != v)
            throw std::invalid_argument("invalid uBLAS integer dimension/index");
        return static_cast<int>(v);
    }
    void finished() {
        if (index != length) throw std::invalid_argument("unexpected extra uBLAS input");
    }
private:
    const double* data;
    int length;
    int index = 0;
};

Matrix dense(Reader& r, int rows, int cols) {
    if (rows * cols > matrix_bound)
        throw std::invalid_argument("uBLAS matrix too large");
    Matrix a(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            a(i, j) = r.real();
    return a;
}

Vector vector(Reader& r, int size) {
    Vector v(size);
    for (int i = 0; i < size; ++i) v(i) = r.real();
    return v;
}

Sparse sparse(Reader& r, int& rows, int& cols) {
    rows = r.number();
    cols = r.number();
    const int count = r.number(matrix_bound);
    if (rows * cols > matrix_bound)
        throw std::invalid_argument("uBLAS sparse matrix too large");
    Sparse a(rows, cols, std::min(count, rows*cols));
    for (int k = 0; k < count; ++k) {
        const int i = r.number(), j = r.number();
        const double value = r.real();
        if (i >= rows || j >= cols)
            throw std::invalid_argument("invalid sparse COO coordinate");
        // Sci COO -> CSR sums duplicate coordinates, including explicit zeros.
        a(i, j) += value;
    }
    return a;
}

std::vector<double> flattened(const Matrix& matrix, bool tagged_shape = false) {
    std::vector<double> output;
    output.reserve(matrix.size1()*matrix.size2() + (tagged_shape ? 2 : 0));
    if (tagged_shape) {
        output.push_back(static_cast<double>(matrix.size1()));
        output.push_back(static_cast<double>(matrix.size2()));
    }
    for (std::size_t i = 0; i < matrix.size1(); ++i)
        for (std::size_t j = 0; j < matrix.size2(); ++j)
            output.push_back(matrix(i, j));
    return output;
}

std::vector<double> flattened_sparse(const Sparse& matrix) {
    Matrix dense_result(matrix.size1(), matrix.size2());
    dense_result = matrix;
    return flattened(dense_result, true);
}

std::vector<double> evaluate_dense(const char* name, Reader& input) {
    if (named(name, "linalg.dot")) {
        const int n = input.number();
        const Vector a = vector(input, n), b = vector(input, n);
        input.finished();
        return {ub::inner_prod(a, b)};
    }
    const int rows = input.number(), cols = input.number();
    const Matrix a = dense(input, rows, cols);
    if (named(name, "linalg.transpose")) {
        input.finished();
        return flattened(Matrix(ub::trans(a)));
    }
    if (named(name, "linalg.matvec")) {
        const Vector v = vector(input, cols);
        input.finished();
        const Vector result = ub::prod(a, v);
        return std::vector<double>(result.begin(), result.end());
    }
    if (named(name, "linalg.matmul")) {
        const int result_cols = input.number();
        const Matrix b = dense(input, cols, result_cols);
        input.finished();
        return flattened(Matrix(ub::prod(a, b)));
    }
    if (rows != cols)
        throw std::invalid_argument("uBLAS determinant/solve needs a square matrix");
    const Vector rhs = named(name, "linalg.solve") ? vector(input, rows) : Vector();
    input.finished();
    // This independently tests uBLAS pivoted LU while retaining a wider
    // accumulator for poorly conditioned saved matrices (e.g. Hilbert n=8).
    // Inputs are still the exact saved binary64 values.
    ExtendedMatrix decomposition(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            decomposition(i, j) = static_cast<long double>(a(i, j));
    ub::permutation_matrix<std::size_t> permutation(rows);
    const std::size_t status = ub::lu_factorize(decomposition, permutation);
    if (status != 0)
        throw std::runtime_error("singular uBLAS LU decomposition");
    if (named(name, "linalg.solve")) {
        ExtendedVector answer(rows);
        for (int i = 0; i < rows; ++i)
            answer(i) = static_cast<long double>(rhs(i));
        ub::lu_substitute(decomposition, permutation, answer);
        std::vector<double> result;
        result.reserve(rows);
        for (int i = 0; i < rows; ++i)
            result.push_back(static_cast<double>(answer(i)));
        return result;
    }
    long double determinant = 1.0L;
    for (int i = 0; i < rows; ++i) {
        determinant *= decomposition(i, i);
        if (permutation(i) != static_cast<std::size_t>(i)) determinant = -determinant;
    }
    return {static_cast<double>(determinant)};
}

std::vector<double> evaluate_sparse(const char* name, Reader& input) {
    int rows, cols;
    const Sparse a = sparse(input, rows, cols);
    if (named(name, "sparse.from_dense")) {
        input.finished();
        // Explicitly create a dense intermediate then convert to a new uBLAS
        // compressed sparse representation for the Sci dense->sparse contract.
        Matrix original(rows, cols);
        original = a;
        Sparse converted(rows, cols);
        converted = original;
        return flattened_sparse(converted);
    }
    if (named(name, "sparse.to_dense")) {
        input.finished();
        Matrix dense_result(rows, cols);
        dense_result = a;
        return flattened(dense_result, true);
    }
    if (named(name, "sparse.transpose")) {
        input.finished();
        Sparse result(ub::trans(a));
        return flattened_sparse(result);
    }
    if (named(name, "sparse.matvec")) {
        const Vector x = vector(input, cols);
        input.finished();
        const Vector result = ub::prod(a, x);
        std::vector<double> data = {1., static_cast<double>(rows)};
        data.insert(data.end(), result.begin(), result.end());
        return data;
    }
    int brows, bcols;
    const Sparse b = sparse(input, brows, bcols);
    input.finished();
    if (named(name, "sparse.matmul")) {
        if (cols != brows) throw std::invalid_argument("sparse matmul shape mismatch");
        Sparse product(ub::prod(a, b));
        return flattened_sparse(product);
    }
    if (rows != brows || cols != bcols)
        throw std::invalid_argument("sparse pair shape mismatch");
    if (named(name, "sparse.add")) {
        Sparse result(a + b);
        return flattened_sparse(result);
    }
    Sparse result(ub::element_prod(a, b));
    return flattened_sparse(result);
}

std::vector<double> compute(const char* name, const double* data, int length) {
    if (!supported(0, name))
        throw std::invalid_argument("unsupported Boost.Numeric.uBLAS function");
    Reader input(data, length);
    auto result = std::strncmp(name, "linalg.", 7) == 0 ?
        evaluate_dense(name, input) : evaluate_sparse(name, input);
    for (double value : result)
        if (!std::isfinite(value))
            throw std::runtime_error("nonfinite uBLAS result");
    return result;
}
} // namespace

extern "C" {
int boost_ublas_supported(int library, const char* name) {
    return supported(library, name) ? 2 : 0;
}

const char* boost_ublas_error() {
    return last_error.c_str();
}

int boost_ublas_output_count(const char* name, const double* data, int length) {
    try {
        return static_cast<int>(compute(name, data, length).size());
    } catch (...) {
        return 0;
    }
}

int boost_ublas_evaluate(int library, const char* name, const double* data,
                         int length, double* output, int capacity) {
    last_error.clear();
    try {
        if (!supported(library, name))
            throw std::invalid_argument("unsupported Boost.Numeric.uBLAS library or function");
        const auto values = compute(name, data, length);
        if (!output || capacity < static_cast<int>(values.size()))
            throw std::invalid_argument("insufficient Boost.uBLAS output capacity");
        std::copy(values.begin(), values.end(), output);
        return 0;
    } catch (const std::exception& error) {
        last_error = error.what();
    } catch (...) {
        last_error = "unknown Boost.uBLAS evaluation failure";
    }
    return -1;
}
}
