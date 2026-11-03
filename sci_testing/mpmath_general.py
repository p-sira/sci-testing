"""Independent mpmath comparisons for finite, non-special Sci operations.

Every result comes from mpmath or a mathematical operation on mpmath numbers.
Inputs are the saved binary64 cases; matrix and complex outputs retain Sci's
row-major and interleaved-real/imaginary numeric serialization, respectively.
Unavailable domains are omitted rather than being presented as comparisons.
"""

import math
from itertools import pairwise

import mpmath as mp

PROVIDER = "mpmath"
_COMPLEX_UNARY = {
    "acos": mp.acos, "acosh": mp.acosh, "asin": mp.asin,
    "asinh": mp.asinh, "atan": mp.atan, "atanh": mp.atanh,
    "cos": mp.cos, "cosh": mp.cosh, "exp": mp.exp,
    "sin": mp.sin, "sinh": mp.sinh, "sqrt": mp.sqrt,
    "tan": mp.tan, "tanh": mp.tanh, "ln": mp.log,
}
_COMPLEX_BINARY = {"add", "subtract", "multiply", "divide", "powc",
                   "powf", "powi", "powu"}
_COMPLEX_OTHER = {"abs", "arg", "real", "imag", "norm_sqr",
                  "from_polar", "to_polar", "reciprocal", "cbrt",
                  "negate", "conjugate"}
_LINALG = {
    "det", "solve", "inv", "dot", "matvec", "matmul", "transpose",
    "cholesky", "cho_solve", "eigvalsh", "svdvals", "norm",
    "expm", "logm", "sqrtm", "sinm", "cosm", "sinhm", "coshm",
    "tanhm", "tanm", "hilbert", "invhilbert", "pascal", "invpascal",
    "toeplitz", "circulant", "hankel", "diagsvd", "block_diag",
    "bandwidth", "ishermitian", "issymmetric", "hadamard",
    "companion", "fiedler", "leslie", "khatri_rao", "solve_triangular",
    "solve_sylvester", "fractional_matrix_power", "funm",
    "convolution_matrix", "dft", "solve_circulant", "solve_toeplitz",
}
_STATS = {
    "normal_cdf", "normal_pdf", "normal_ppf", "variance", "mean",
    "median", "std", "sem", "variation", "describe", "moment", "skew",
    "kurtosis", "zscore", "zmap", "gmean", "hmean", "pmean", "gstd",
    "gzscore", "median_abs_deviation", "quantile", "scoreatpercentile",
    "iqr", "trim_mean", "rankdata", "entropy", "kstat",
    "circmean", "circvar", "circstd", "dirichlet_pdf", "dirichlet_logpdf",
    "multinomial_pmf", "multinomial_logpmf", "energy_distance",
    "wasserstein_distance", "tmean", "tstd", "tvar", "tsem", "tmin", "tmax",
    "pdf", "logpdf", "cdf", "sf", "ppf", "exp", "log", "abs",
}


def _number(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError("numeric inputs required")
    if not math.isfinite(value):
        raise ValueError("finite inputs required")
    return mp.mpf(value)


def _vector(value):
    if not isinstance(value, (list, tuple)):
        raise TypeError("vector input required")
    return list(map(_number, value))


def _matrix(value):
    if (not isinstance(value, (list, tuple)) or not value or
            not all(isinstance(row, (list, tuple)) and row for row in value) or
            any(len(row) != len(value[0]) for row in value)):
        raise ValueError("nonempty rectangular matrix required")
    return mp.matrix([_vector(row) for row in value])


def _flat_matrix(matrix):
    return [matrix[i, j] for i in range(matrix.rows) for j in range(matrix.cols)]


def _complex(value):
    pair = _vector(value)
    if len(pair) != 2:
        raise ValueError("complex pair requires two components")
    return mp.mpc(*pair)


def _components(z):
    return [mp.re(z), mp.im(z)]


def _flatten(values):
    if isinstance(values, mp.matrix):
        return _flat_matrix(values)
    if isinstance(values, (mp.mpc, complex)):
        return _components(values)
    if isinstance(values, (list, tuple)):
        return [component for value in values for component in _flatten(value)]
    return [values]


def _complex_values(name, args):
    if name == "from_polar":
        r, theta = _number(args["p0"]), _number(args["p1"])
        return _components(r * mp.exp(mp.j * theta))
    if "x" in args and "y" in args:
        pair = [args["x"], args["y"]]
    else:
        pair = args["p0"]
    z = _complex(pair)
    if name == "real":
        return [z.real]
    if name == "imag":
        return [z.imag]
    if name == "abs":
        return [abs(z)]
    if name == "norm_sqr":
        return [z.real*z.real + z.imag*z.imag]
    if name in ("arg", "to_polar"):
        phase = mp.arg(z)
        if z.imag == 0 and z.real < 0 and math.copysign(1, pair[1]) < 0:
            phase = -mp.pi
        return [phase] if name == "arg" else [abs(z), phase]

    other = None
    if name in _COMPLEX_BINARY:
        w = args["p1"]
        other = _complex(w) if isinstance(w, (list, tuple)) else _number(w)
    if name == "add":
        value = z + other
    elif name == "subtract":
        value = z - other
    elif name == "multiply":
        value = z * other
    elif name == "divide":
        value = z / other
    elif name in ("powc", "powf", "powi", "powu"):
        if mp.im(other) == 0 and mp.re(other) == int(mp.re(other)):
            value = z ** int(mp.re(other))
        else:
            value = mp.power(z, other)
    elif name == "reciprocal":
        value = 1/z
    elif name == "negate":
        value = -z
    elif name == "conjugate":
        value = mp.conj(z)
    elif name == "cbrt":
        value = mp.root(z, 3)
    else:
        value = _COMPLEX_UNARY[name](z)

    # mpmath does not preserve signed zero. Preserve the side of the branch
    # cut carried by the saved binary64 input's sign bit.
    if pair[1] == 0:
        below = math.copysign(1, pair[1]) < 0
        if name in ("asin", "atanh") and abs(pair[0]) > 1:
            value = mp.mpc(mp.re(value), -abs(mp.im(value)) if below else abs(mp.im(value)))
        elif name == "acos" and abs(pair[0]) > 1:
            value = mp.mpc(mp.re(value), abs(mp.im(value)) if below else -abs(mp.im(value)))
        elif name == "acosh" and pair[0] < 1:
            value = mp.mpc(mp.re(value), -abs(mp.im(value)) if below else abs(mp.im(value)))
        elif below and pair[0] < 0:
            if name in ("ln", "sqrt", "cbrt"):
                value = mp.conj(value)
            elif name in ("powc", "powf") and mp.im(other) != 0:
                value = mp.exp(other*(mp.log(-z.real)-mp.pi*mp.j))
            elif name in ("powc", "powf"):
                value = mp.conj(value)
        if name in ("powc", "powf") and pair[0] < 0 and mp.im(other) == 0:
            twice = 2*mp.re(other)
            if twice == int(twice) and int(twice) % 2:
                value = mp.mpc(0, mp.im(value))
    elif pair[0] == 0 and name in ("atan", "asinh") and abs(pair[1]) > 1:
        magnitude = abs(mp.re(value))
        value = mp.mpc(-magnitude if math.copysign(1, pair[0]) < 0 else magnitude,
                       mp.im(value))
    return _components(value)


def _raw_linalg(name, args):
    raw = args.get("raw")
    if not isinstance(raw, dict) or raw.get("api") != name:
        return None
    return raw.get("args")


def _symmetric(a, lower):
    if a.rows != a.cols:
        raise ValueError("square symmetric matrix required")
    return mp.matrix([[a[i, j] if (i >= j) == lower or i == j
                       else a[j, i] for j in range(a.cols)]
                      for i in range(a.rows)])


def _linalg_values(name, args):
    saved = _raw_linalg(name, args)
    if name in ("det", "solve", "dot", "matvec", "matmul", "transpose"):
        if "raw" in args:
            return None
        if name == "dot":
            return [mp.fdot(_vector(args["u"]), _vector(args["v"]))]
        a = _matrix(args["a"])
        if name == "det":
            return [mp.det(a)]
        if name == "solve":
            return _flat_matrix(mp.lu_solve(a, mp.matrix(_vector(args["b"]))))
        if name == "matvec":
            return _flat_matrix(a*mp.matrix(_vector(args["u"])))
        if name == "matmul":
            return _flat_matrix(a*_matrix(args["b"]))
        return _flat_matrix(a.T)
    if saved is None:
        return None
    if name == "inv":
        return _flat_matrix(_matrix(saved[0])**-1)
    if name == "cholesky":
        lower = bool(saved[1])
        factor = mp.cholesky(_symmetric(_matrix(saved[0]), lower))
        return _flat_matrix(factor if lower else factor.T)
    if name == "cho_solve":
        return _flat_matrix(mp.lu_solve(_symmetric(_matrix(saved[0]), bool(saved[-1])),
                                        mp.matrix(_vector(saved[1]))))
    if name == "eigvalsh":
        eigenvalues = mp.eigsy(_symmetric(_matrix(saved[0]), bool(saved[-1])), eigvals_only=True)
        return list(eigenvalues)
    if name == "svdvals":
        return list(mp.svd(_matrix(saved[0]), compute_uv=False))
    if name == "norm":
        a = _matrix(saved[0])
        order = saved[1]
        if order in ("fro", None):
            return [mp.sqrt(mp.fsum(v*v for v in a))]
        if order == 1:
            return [max(mp.fsum(abs(a[i, j]) for i in range(a.rows)) for j in range(a.cols))]
        if order == mp.inf or order == "inf":
            return [max(mp.fsum(abs(a[i, j]) for j in range(a.cols)) for i in range(a.rows))]
        return None
    if name in {"expm", "logm", "sqrtm", "sinm", "cosm", "sinhm", "coshm",
                "tanm", "tanhm", "fractional_matrix_power", "funm"}:
        a = _matrix(saved[0])
        if a.rows != a.cols:
            return None
        if name == "expm":
            result = mp.expm(a)
        elif name == "logm":
            result = mp.logm(a)
        elif name == "sqrtm":
            result = mp.sqrtm(a)
        elif name == "fractional_matrix_power":
            result = mp.expm(_number(saved[1])*mp.logm(a))
        elif name == "funm":
            if saved[1] != "exp":
                return None
            result = mp.expm(a)
        else:
            # Entire matrix functions, with the same definition as the
            # corresponding scalar functions via exponentials.
            plus, minus = mp.expm(a*mp.j), mp.expm(-a*mp.j)
            if name == "sinm":
                result = (plus-minus)/(2*mp.j)
            elif name == "cosm":
                result = (plus+minus)/2
            elif name == "tanm":
                result = ((plus-minus)/(2*mp.j))*(((plus+minus)/2)**-1)
            else:
                plus, minus = mp.expm(a), mp.expm(-a)
                if name == "sinhm":
                    result = (plus-minus)/2
                elif name == "coshm":
                    result = (plus+minus)/2
                else:
                    result = ((plus-minus)/2)*(((plus+minus)/2)**-1)
        flat = _flat_matrix(result)
        # Sci's saved linalg fixture emits the real matrix followed by the
        # imaginary matrix for these complex-valued matrix functions.
        if name in {"logm", "sqrtm", "fractional_matrix_power", "funm"}:
            return [mp.re(v) for v in flat]+[mp.im(v) for v in flat]
        if any(abs(mp.im(v)) > mp.mpf('1e-60') for v in flat):
            return None
        return [mp.re(v) for v in flat]
    if name == "hilbert":
        n = int(saved[0])
        return [1/mp.mpf(i+j+1) for i in range(n) for j in range(n)]
    if name == "invhilbert":
        n = int(saved[0])
        h = mp.matrix([[1/mp.mpf(i+j+1) for j in range(n)] for i in range(n)])
        return _flat_matrix(h**-1) if n else []
    if name in ("pascal", "invpascal"):
        n, kind = int(saved[0]), saved[1]
        if kind not in ("symmetric", "lower", "upper"):
            return None
        lower = mp.matrix([[mp.binomial(i, j) if i >= j else 0 for j in range(n)] for i in range(n)])
        result = lower*lower.T if kind == "symmetric" else lower if kind == "lower" else lower.T
        if name == "invpascal" and n:
            result = result**-1
        return _flat_matrix(result) if n else []
    if name == "toeplitz":
        col = _vector(saved[0])
        row = _vector(saved[1]) if saved[1] is not None else col
        return [col[i-j] if i >= j else row[j-i] for i in range(len(col)) for j in range(len(row))]
    if name == "circulant":
        col = _vector(saved[0]); n = len(col)
        return [col[(i-j) % n] for i in range(n) for j in range(n)]
    if name == "hankel":
        col = _vector(saved[0]); n = len(col)
        row = _vector(saved[1]) if saved[1] is not None else [col[-1]]+[mp.mpf(0)]*(n-1)
        sequence = col+row[1:]
        return [sequence[i+j] for i in range(n) for j in range(len(row))]
    if name == "diagsvd":
        diagonal, rows, cols = _vector(saved[0]), int(saved[1]), int(saved[2])
        return [diagonal[i] if i == j and i < len(diagonal) else mp.mpf(0)
                for i in range(rows) for j in range(cols)]
    if name == "block_diag":
        blocks = [_matrix(value) for value in saved[0]]
        height = sum(a.rows for a in blocks)
        width = sum(a.cols for a in blocks)
        result = mp.zeros(height, width)
        i = j = 0
        for block in blocks:
            for row in range(block.rows):
                for col in range(block.cols):
                    result[i+row, j+col] = block[row, col]
            i += block.rows; j += block.cols
        return _flat_matrix(result)
    if name == "bandwidth":
        a = _matrix(saved[0])
        occupied = [(i, j) for i in range(a.rows) for j in range(a.cols) if a[i, j]]
        return [max([i-j for i, j in occupied]+[0]), max([j-i for i, j in occupied]+[0])]
    if name in ("ishermitian", "issymmetric"):
        a = _matrix(saved[0]); atol, rtol = saved[1:]
        if a.rows != a.cols:
            return [0]
        if atol is None and rtol is None:
            return [float(all(a[i, j] == a[j, i] for i in range(a.rows) for j in range(a.cols)))]
        return [float(all(abs(a[i, j]-a[j, i]) <= (atol or 0)+(rtol or 0)*abs(a[j, i])
                          for i in range(a.rows) for j in range(a.cols)))]
    if name == "hadamard":
        n = int(saved[0]);
        if n <= 0:
            return []
        if n & (n-1):
            return None
        return [(-1)**((i & j).bit_count() % 2) for i in range(n) for j in range(n)]
    if name == "companion":
        p = _vector(saved[0]); n = len(p)-1
        return [(-p[j+1]/p[0] if i == 0 else (1 if j == i-1 else 0))
                for i in range(n) for j in range(n)]
    if name == "fiedler":
        x = _vector(saved[0])
        return [abs(x[i]-x[j]) for i in range(len(x)) for j in range(len(x))]
    if name == "leslie":
        f, s = _vector(saved[0]), _vector(saved[1]); n = len(f)
        return [(f[j] if i == 0 else s[i-1] if j == i-1 else 0)
                for i in range(n) for j in range(n)]
    if name == "khatri_rao":
        a, b = _matrix(saved[0]), _matrix(saved[1])
        if a.cols != b.cols:
            return None
        return [a[i//b.rows, j]*b[i%b.rows, j]
                for i in range(a.rows*b.rows) for j in range(a.cols)]
    if name == "solve_triangular":
        a, b, lower = _matrix(saved[0]), mp.matrix(_vector(saved[1])), bool(saved[2])
        tri = mp.matrix([[a[i, j] if (j <= i if lower else j >= i) else 0
                          for j in range(a.cols)] for i in range(a.rows)])
        return _flat_matrix(mp.lu_solve(tri, b))
    if name == "solve_sylvester":
        a, b, q = map(_matrix, saved[:3]); n, m = a.rows, b.rows
        coefficient = mp.matrix([[((a[i, k] if j == l else 0)+
                                  (b[l, j] if i == k else 0))
                                 for k in range(n) for l in range(m)]
                                for i in range(n) for j in range(m)])
        solution = mp.lu_solve(coefficient, mp.matrix(_flat_matrix(q)))
        return list(solution)
    if name == "convolution_matrix":
        a, n, mode = _vector(saved[0]), int(saved[1]), saved[2]
        size = len(a)+n-1
        full = [[a[i-j] if 0 <= i-j < len(a) else 0 for j in range(n)]
                for i in range(size)]
        if mode == "full":
            rows = full
        elif mode == "same":
            left = (len(a)-1)//2; rows = full[left:left+n]
        elif mode == "valid":
            left = min(len(a), n)-1; rows = full[left:left+abs(n-len(a))+1]
        else:
            return None
        return [v for row in rows for v in row]
    if name == "dft":
        n, scale = int(saved[0]), saved[1]
        factor = (mp.sqrt(n) if scale == "sqrtn" else n if scale == "n" else 1)
        values = [mp.exp(-2*mp.pi*mp.j*i*j/n)/factor
                  for i in range(n) for j in range(n)]
        return [mp.re(v) for v in values]+[mp.im(v) for v in values]
    return None


def _mean(x):
    return mp.fsum(x)/len(x)


def _variance(x, ddof=0):
    average = _mean(x)
    return mp.fsum((v-average)**2 for v in x)/(len(x)-ddof)


def _quantile(x, probability):
    data = sorted(x)
    pos = _number(probability)*(len(data)-1)
    i = int(mp.floor(pos))
    if i == len(data)-1:
        return data[-1]
    return data[i]+(pos-i)*(data[i+1]-data[i])


def _distribution(name, args):
    """One-dimensional probability distributions via the mpmath special API.

    The saved Sci cases use the named parameters of the distributions. The
    discrete PMF/CDF functions keep their original support and integer tests.
    """
    kind = args.get("d")
    if kind not in {
        "Normal", "Beta", "Gamma", "ChiSquared", "Chi", "StudentT", "F",
        "Exponential", "Cauchy", "Laplace", "LogNormal", "Uniform",
        "Weibull", "Pareto", "InverseGamma", "Gumbel", "Bernoulli",
        "Binomial", "Poisson", "Geometric", "NegativeBinomial",
        "Hypergeometric", "Erlang", "DiscreteUniform", "Levy", "Triangular",
    }:
        return None
    p = _vector(args.get("p", []))
    xs = args.get("x", [])
    if len(xs) != 1:
        return None
    if isinstance(xs[0], str) and xs[0] in ("-inf", "inf"):
        x = -mp.inf if xs[0] == "-inf" else mp.inf
    else:
        x = _number(xs[0])
    if name == "ppf":
        if not 0 <= x <= 1:
            return None
        # Quantiles at exactly 0 or 1 have finite support endpoints for
        # bounded distributions.  The discrete distributions follow Sci's
        # integer inverse-CDF convention at q=0 (one below their minimum).
        # Unbounded quantiles are nonfinite and remain unavailable.
        if x == 0:
            if kind in {"Beta", "Gamma", "ChiSquared", "Chi", "F",
                        "Exponential", "LogNormal", "Weibull",
                        "InverseGamma", "Erlang"}:
                return [mp.mpf(0)]
            if kind in {"Uniform", "Triangular", "Pareto", "Levy"}:
                return [p[0]]
            if kind == "Geometric":
                return [mp.mpf(0)]
            if kind == "DiscreteUniform":
                return [p[0]-1]
            if kind == "Hypergeometric":
                return [max(0, p[2]-p[0]+p[1])-1]
            if kind in {"Bernoulli", "Binomial", "Poisson",
                        "NegativeBinomial"}:
                return [mp.mpf(-1)]
            return None
        if x == 1:
            if kind == "Beta":
                return [mp.mpf(1)]
            if kind in {"Uniform", "Triangular"}:
                return [p[1]]
            if kind == "Bernoulli":
                return [mp.mpf(1)]
            if kind == "Binomial":
                return [p[0]]
            if kind == "Hypergeometric":
                return [min(p[1], p[2])]
            if kind == "DiscreteUniform":
                return [p[1]-1]
            return None
        if kind == "Normal":
            return [p[0]+p[1]*mp.sqrt(2)*mp.erfinv(2*x-1)]
        if kind == "Uniform":
            return [p[0]+(p[1]-p[0])*x]
        if kind == "Exponential":
            return [-p[0]*mp.log1p(-x)]
        if kind == "Cauchy":
            return [p[0]+p[1]*mp.tan(mp.pi*(x-mp.mpf('0.5')))]
        if kind == "Laplace":
            return [p[0]+p[1]*mp.log(2*x) if x < .5 else
                    p[0]-p[1]*mp.log(2*(1-x))]
        if kind == "LogNormal":
            return [mp.exp(p[0]+p[1]*mp.sqrt(2)*mp.erfinv(2*x-1))]
        if kind == "Weibull":
            return [p[1]*(-mp.log1p(-x))**(1/p[0])]
        if kind == "Pareto":
            return [p[0]*(1-x)**(-1/p[1])]
        if kind == "Gumbel":
            return [p[0]-p[1]*mp.log(-mp.log(x))]
        if kind == "Triangular":
            low, high, mode = p
            split = (mode-low)/(high-low)
            return [low+mp.sqrt(x*(high-low)*(mode-low)) if x < split else
                    high-mp.sqrt((1-x)*(high-low)*(high-mode))]
        if kind == "Levy":
            return [p[0]+p[1]/(2*mp.erfinv(1-x)**2)]
        if kind == "Chi" and p[0] == 1:
            # Chi(1) is |N(0,1)|.  Bisection on [0,1] has an absolute
            # resolution floor and cannot resolve quantiles near 1e-100.
            return [mp.sqrt(2)*mp.erfinv(x)]
        if kind == "StudentT" and x == mp.mpf('0.5'):
            return [mp.mpf(0)]
        if kind == "InverseGamma":
            # Invert in z = scale / value, not in value: a scale of 1e-200
            # otherwise falls far below the absolute bisection resolution.
            # Q(a,z) is the regularized *upper* incomplete gamma, decreasing
            # from 1 to 0.  For a=1, Q(1,z)=exp(-z) exactly.
            if p[0] == 1:
                return [-p[1]/mp.log(x)]
            low = mp.mpf(0)
            if x > mp.mpf('0.5'):
                # P(a,z) ~ z**a/Gamma(a+1) as z->0.  For small shape
                # and x~1, z can be below 1e-150; a fixed [0,1] bracket
                # would need hundreds of iterations merely to reach it.
                target = 1-x
                high = (target*mp.gamma(p[0]+1))**(1/p[0])
                while mp.gammainc(p[0], 0, high, regularized=True) < target:
                    high *= 2
                for _ in range(140):
                    mid = (low+high)/2
                    if mp.gammainc(p[0], 0, mid, regularized=True) < target:
                        low = mid
                    else:
                        high = mid
            else:
                high = max(mp.mpf(1), p[0])
                while mp.gammainc(p[0], high, mp.inf, regularized=True) > x:
                    high *= 2
                for _ in range(140):
                    mid = (low+high)/2
                    if mp.gammainc(p[0], mid, mp.inf, regularized=True) > x:
                        low = mid
                    else:
                        high = mid
            return [p[1]/((low+high)/2)]
        if kind in {"Bernoulli", "Binomial", "Poisson", "Geometric",
                    "NegativeBinomial", "Hypergeometric", "DiscreteUniform"}:
            if kind == "Bernoulli":
                return [0 if x <= 1-p[0] else 1]
            if kind == "DiscreteUniform":
                return [int(p[0])+int(mp.ceil(x*(p[1]-p[0])))-1]
            low, high = (1, 10000) if kind == "Geometric" else (0, 10000)
            if kind == "Binomial":
                high = int(p[0])
            elif kind == "Hypergeometric":
                high = int(min(p[1],p[2]))
            else:
                while _discrete_cdf(kind, high, p) < x:
                    high *= 2
                    if high > 100000:
                        return None
            while low < high:
                middle = (low+high)//2
                if _discrete_cdf(kind, middle, p) >= x:
                    high = middle
                else:
                    low = middle+1
            return [low]
        if kind in {"Beta", "Gamma", "ChiSquared", "Chi", "StudentT", "F",
                    "InverseGamma", "Erlang"}:
            if kind == "Beta":
                low, high = mp.mpf(0), mp.mpf(1)
            elif kind == "StudentT":
                low, high = mp.mpf(-1), mp.mpf(1)
            else:
                low, high = mp.mpf(0), mp.mpf(1)
            while _continuous_cdf(kind, low, p) > x:
                low *= 2
            while _continuous_cdf(kind, high, p) < x:
                high *= 2
            # Bisection on the public mpmath CDF works for extreme skew and
            # small shape parameters without depending on Sci's quantiles.
            for _ in range(90):
                mid = (low+high)/2
                if _continuous_cdf(kind, mid, p) < x:
                    low = mid
                else:
                    high = mid
            return [(low+high)/2]
        return None
    if kind in {"Bernoulli", "Binomial", "Poisson", "Geometric",
                "NegativeBinomial", "Hypergeometric", "DiscreteUniform"}:
        if name in ("cdf", "sf"):
            cdf = _discrete_cdf(kind, int(mp.floor(x)) if mp.isfinite(x) else x, p)
            return [cdf if name == "cdf" else 1-cdf]
        probability = _discrete_pmf(kind, x, p)
        return [probability if name == "pdf" else mp.log(probability)]
    if name in ("cdf", "sf"):
        cdf = _continuous_cdf(kind, x, p)
        return [cdf if name == "cdf" else 1-cdf]
    probability = _continuous_pdf(kind, x, p)
    return [probability if name == "pdf" else mp.log(probability)]


def _regularized_gamma(k, x):
    return mp.gammainc(k, 0, x, regularized=True)


def _continuous_cdf(kind, x, p):
    if x == -mp.inf:
        return mp.mpf(0)
    if x == mp.inf:
        return mp.mpf(1)
    if kind == "Normal":
        return mp.erfc(-(x-p[0])/(p[1]*mp.sqrt(2)))/2
    if kind == "Beta":
        return mp.betainc(p[0],p[1],0,x,regularized=True) if 0 < x < 1 else (0 if x <= 0 else 1)
    if kind in ("Gamma", "Erlang"):
        return _regularized_gamma(p[0],x/p[1]) if x > 0 else mp.mpf(0)
    if kind == "ChiSquared":
        return _regularized_gamma(p[0]/2,x/2) if x > 0 else mp.mpf(0)
    if kind == "Chi":
        return _regularized_gamma(p[0]/2,x*x/2) if x > 0 else mp.mpf(0)
    if kind == "StudentT":
        t = mp.betainc(p[0]/2, mp.mpf('.5'), 0, p[0]/(p[0]+x*x), regularized=True)/2
        return t if x < 0 else 1-t
    if kind == "F":
        return mp.betainc(p[0]/2,p[1]/2,0,p[0]*x/(p[0]*x+p[1]),regularized=True) if x>0 else mp.mpf(0)
    if kind == "Exponential":
        return -mp.expm1(-x/p[0]) if x > 0 else mp.mpf(0)
    if kind == "Cauchy":
        return mp.mpf('.5')+mp.atan((x-p[0])/p[1])/mp.pi
    if kind == "Laplace":
        t=(x-p[0])/p[1]
        return mp.exp(t)/2 if t < 0 else 1-mp.exp(-t)/2
    if kind == "LogNormal":
        return mp.erfc(-(mp.log(x)-p[0])/(p[1]*mp.sqrt(2)))/2 if x > 0 else mp.mpf(0)
    if kind == "Uniform":
        return max(mp.mpf(0),min(mp.mpf(1),(x-p[0])/(p[1]-p[0])))
    if kind == "Weibull":
        return -mp.expm1(-(x/p[1])**p[0]) if x > 0 else mp.mpf(0)
    if kind == "Pareto":
        return 1-(p[0]/x)**p[1] if x >= p[0] else mp.mpf(0)
    if kind == "InverseGamma":
        return mp.gammainc(p[0],p[1]/x,mp.inf,regularized=True) if x > 0 else mp.mpf(0)
    if kind == "Gumbel":
        return mp.exp(-mp.exp(-(x-p[0])/p[1]))
    if kind == "Triangular":
        low, high, mode = p
        if x <= low:
            return mp.mpf(0)
        if x >= high:
            return mp.mpf(1)
        return ((x-low)**2/((high-low)*(mode-low)) if x <= mode else
                1-(high-x)**2/((high-low)*(high-mode)))
    if kind == "Levy":
        return mp.erfc(mp.sqrt(p[1]/(2*(x-p[0])))) if x > p[0] else mp.mpf(0)
    raise ValueError(kind)


def _continuous_pdf(kind, x, p):
    if not mp.isfinite(x):
        return mp.mpf(0)
    if kind == "Normal":
        z=(x-p[0])/p[1]
        return mp.exp(-z*z/2)/(p[1]*mp.sqrt(2*mp.pi))
    if kind == "Beta":
        # At either support endpoint, the density has a finite, nonzero
        # limit when the corresponding shape is exactly one.  Shapes below
        # one have an infinite limit (not a finite comparison).
        if x == 0:
            return (p[1] if p[0] == 1 else
                    mp.inf if p[0] < 1 else mp.mpf(0))
        if x == 1:
            return (p[0] if p[1] == 1 else
                    mp.inf if p[1] < 1 else mp.mpf(0))
        return (x**(p[0]-1)*(1-x)**(p[1]-1)/mp.beta(p[0],p[1])
                if 0 < x < 1 else mp.mpf(0))
    if kind in ("Gamma", "Erlang"):
        if x == 0:
            return (1/p[1] if p[0] == 1 else
                    mp.inf if p[0] < 1 else mp.mpf(0))
        return (x**(p[0]-1)*mp.exp(-x/p[1])/(mp.gamma(p[0])*p[1]**p[0])
                if x > 0 else mp.mpf(0))
    if kind == "ChiSquared":
        return _continuous_pdf("Gamma",x,[p[0]/2,mp.mpf(2)])
    if kind == "Chi":
        if x == 0:
            return (mp.sqrt(2/mp.pi) if p[0] == 1 else
                    mp.inf if p[0] < 1 else mp.mpf(0))
        return x**(p[0]-1)*mp.exp(-x*x/2)/(2**(p[0]/2-1)*mp.gamma(p[0]/2)) if x>0 else mp.mpf(0)
    if kind == "StudentT":
        return mp.gamma((p[0]+1)/2)/(mp.sqrt(p[0]*mp.pi)*mp.gamma(p[0]/2))*(1+x*x/p[0])**(-(p[0]+1)/2)
    if kind == "F":
        if x == 0:
            return (1 if p[0] == 2 else
                    mp.inf if p[0] < 2 else mp.mpf(0))
        return (p[0]/p[1])**(p[0]/2)*x**(p[0]/2-1)/(mp.beta(p[0]/2,p[1]/2)*
                (1+p[0]*x/p[1])**((p[0]+p[1])/2)) if x>0 else mp.mpf(0)
    if kind == "Exponential":
        return mp.exp(-x/p[0])/p[0] if x>=0 else mp.mpf(0)
    if kind == "Cauchy":
        z=(x-p[0])/p[1]
        return 1/(mp.pi*p[1]*(1+z*z))
    if kind == "Laplace":
        return mp.exp(-abs(x-p[0])/p[1])/(2*p[1])
    if kind == "LogNormal":
        return mp.exp(-((mp.log(x)-p[0])/p[1])**2/2)/(x*p[1]*mp.sqrt(2*mp.pi)) if x>0 else mp.mpf(0)
    if kind == "Uniform":
        return 1/(p[1]-p[0]) if p[0] <= x <= p[1] else mp.mpf(0)
    if kind == "Weibull":
        if x == 0:
            return (1/p[1] if p[0] == 1 else
                    mp.inf if p[0] < 1 else mp.mpf(0))
        z=x/p[1]
        return p[0]/p[1]*z**(p[0]-1)*mp.exp(-z**p[0]) if x>0 else mp.mpf(0)
    if kind == "Pareto":
        return p[1]*p[0]**p[1]/x**(p[1]+1) if x>=p[0] else mp.mpf(0)
    if kind == "InverseGamma":
        return p[1]**p[0]*x**(-p[0]-1)*mp.exp(-p[1]/x)/mp.gamma(p[0]) if x>0 else mp.mpf(0)
    if kind == "Gumbel":
        z=(x-p[0])/p[1]
        return mp.exp(-z-mp.exp(-z))/p[1]
    if kind == "Triangular":
        low, high, mode = p
        if x == low and mode == low:
            return 2/(high-low)
        if x == high and mode == high:
            return 2/(high-low)
        if low < x <= mode:
            return 2*(x-low)/((high-low)*(mode-low))
        if mode < x < high:
            return 2*(high-x)/((high-low)*(high-mode))
        return mp.mpf(0)
    if kind == "Levy":
        t=x-p[0]
        return mp.sqrt(p[1]/(2*mp.pi))*mp.exp(-p[1]/(2*t))/t**mp.mpf('1.5') if t>0 else mp.mpf(0)
    raise ValueError(kind)


def _discrete_pmf(kind, x, p):
    if not mp.isfinite(x) or x != mp.floor(x):
        return mp.mpf(0)
    k = int(x)
    if kind == "Bernoulli":
        return 1-p[0] if k==0 else p[0] if k==1 else mp.mpf(0)
    if kind == "Binomial":
        return mp.binomial(p[0],k)*p[1]**k*(1-p[1])**(p[0]-k) if 0<=k<=p[0] else mp.mpf(0)
    if kind == "Poisson":
        return mp.exp(-p[0])*p[0]**k/mp.factorial(k) if k>=0 else mp.mpf(0)
    if kind == "Geometric":
        return p[0]*(1-p[0])**(k-1) if k>=1 else mp.mpf(0)
    if kind == "NegativeBinomial":
        return mp.binomial(p[0]+k-1,k)*p[1]**p[0]*(1-p[1])**k if k>=0 else mp.mpf(0)
    if kind == "Hypergeometric":
        return (mp.binomial(p[1],k)*mp.binomial(p[0]-p[1],p[2]-k)/
                mp.binomial(p[0],p[2])) if max(0,p[2]-(p[0]-p[1]))<=k<=min(p[1],p[2]) else mp.mpf(0)
    if kind == "DiscreteUniform":
        return 1/(p[1]-p[0]) if p[0]<=k<p[1] else mp.mpf(0)
    raise ValueError(kind)


def _discrete_cdf(kind, k, p):
    if k == -mp.inf:
        return mp.mpf(0)
    if k == mp.inf:
        return mp.mpf(1)
    k = int(mp.floor(k))
    if kind == "Bernoulli":
        return 0 if k<0 else 1-p[0] if k==0 else 1
    if kind == "Binomial":
        return 0 if k<0 else 1 if k>=p[0] else mp.betainc(p[0]-k,k+1,0,1-p[1],regularized=True)
    if kind == "Poisson":
        return mp.gammainc(k+1,p[0],mp.inf,regularized=True) if k>=0 else 0
    if kind == "Geometric":
        return 1-(1-p[0])**k if k>=1 else 0
    if kind == "NegativeBinomial":
        return mp.betainc(p[0],k+1,0,p[1],regularized=True) if k>=0 else 0
    if kind == "Hypergeometric":
        return mp.fsum(_discrete_pmf(kind,i,p) for i in range(max(0,int(p[2]-p[0]+p[1])),min(k,int(p[1]),int(p[2]))+1))
    if kind == "DiscreteUniform":
        return max(0,min(1,(k-p[0]+1)/(p[1]-p[0])))
    raise ValueError(kind)


def _stats_values(name, args):
    if name in ("pdf", "logpdf", "cdf", "sf", "ppf") and args.get("d"):
        return _distribution(name, args)
    if name in ("exp", "log", "abs") and not args.get("d"):
        xs = _vector(args.get("x", []))
        if len(xs) != 1:
            return None
        return [{"exp":mp.exp, "log":mp.log, "abs":abs}[name](xs[0])]
    if name in ("normal_cdf", "normal_pdf", "normal_ppf"):
        value = _number(args["x"])
        if name == "normal_cdf":
            return [mp.erfc(-value/mp.sqrt(2))/2]
        if name == "normal_pdf":
            return [mp.exp(-value**2/2)/mp.sqrt(2*mp.pi)]
        return [mp.sqrt(2)*mp.erfinv(2*value-1)]
    if name == "variance" and "values" in args:
        x = _vector(args["values"])
        return [_variance(x, 1)]
    if args.get("d") or args.get("g"):
        return None
    x = _vector(args.get("x", []))
    y = _vector(args.get("y", []))
    p = _vector(args.get("p", []))
    if not x:
        return None
    average = _mean(x)
    if name == "mean":
        return [average]
    if name == "median":
        return [_quantile(x, 0.5)]
    if name in ("std", "sem", "variation", "zscore", "zmap"):
        if name == "zmap":
            baseline = y
            if not baseline:
                return None
            mu, variance = _mean(baseline), _variance(baseline, int(p[0]))
        else:
            mu, variance = average, _variance(x, int(p[0])) if p else _variance(x)
        sigma = mp.sqrt(variance)
        if name == "std":
            return [sigma]
        if name == "sem":
            return [sigma/mp.sqrt(len(x))]
        if name == "variation":
            return [sigma/mu]
        if name in ("zscore", "zmap"):
            return [(v-mu)/sigma for v in x]
    if name == "median_abs_deviation":
        median = _quantile(x, 0.5)
        return [_quantile([abs(value-median) for value in x], 0.5)/p[0]]
    if name in ("quantile", "scoreatpercentile", "iqr"):
        if name == "iqr":
            return [_quantile(x, .75)-_quantile(x, .25)]
        return [_quantile(x, p[0]/100 if name == "scoreatpercentile" else p[0])]
    if name in ("skew", "kurtosis", "describe", "moment"):
        central = _variance(x)
        moment = lambda k: mp.fsum((v-average)**k for v in x)/len(x)
        third = moment(3)/central**mp.mpf('1.5')
        fourth = moment(4)/central**2-3
        if name == "moment":
            return [moment(int(p[0]))]
        if name == "skew":
            return [third if int(p[0]) else third*mp.sqrt(len(x)*(len(x)-1))/(len(x)-2)]
        if name == "kurtosis":
            fisher, bias = map(int, p[:2])
            excess = (fourth if bias else
                      (len(x)-1)*((len(x)+1)*fourth+6)/((len(x)-2)*(len(x)-3)))
            return [excess if fisher else excess+3]
        ddof = int(p[0]) if p else 1
        bias = int(p[1]) if len(p) >= 2 else 1
        skew = third if bias else third*mp.sqrt(len(x)*(len(x)-1))/(len(x)-2)
        kurt = fourth if bias else (len(x)-1)*((len(x)+1)*fourth+6)/((len(x)-2)*(len(x)-3))
        return [len(x), min(x), max(x), average, _variance(x, ddof), skew, kurt]
    if name == "gmean":
        return [mp.exp(_mean([mp.log(v) for v in x]))]
    if name == "hmean":
        return [len(x)/mp.fsum(1/v for v in x)]
    if name == "pmean":
        power = p[0]
        return [mp.exp(_mean([mp.log(v) for v in x])) if power == 0 else
                _mean([v**power for v in x])**(1/power)]
    if name in ("gstd", "gzscore"):
        logarithms = [mp.log(v) for v in x]
        average_log = _mean(logarithms)
        sigma = mp.sqrt(_variance(logarithms, int(p[0]) if p else 1))
        return [mp.exp(sigma)] if name == "gstd" else [(v-average_log)/sigma for v in logarithms]
    if name == "rankdata":
        method = int(p[0])
        if method not in range(5):
            return None
        ordered = sorted(enumerate(x), key=lambda v: (v[1], v[0]))
        ranked = [mp.mpf(0)]*len(x)
        dense = i = 0
        while i < len(x):
            j = i+1
            while j < len(x) and ordered[j][1] == ordered[i][1]:
                j += 1
            dense += 1
            for k in range(i, j):
                ranked[ordered[k][0]] = (
                    (i+j+1)/2 if method == 0 else i+1 if method == 1 else
                    j if method == 2 else dense if method == 3 else k+1)
            i = j
        return ranked
    if name == "trim_mean":
        trimmed = sorted(x)
        n = int(mp.floor(len(x)*p[0]))
        return [_mean(trimmed[n:len(x)-n])]
    if name in ("tmean", "tstd", "tvar", "tsem", "tmin", "tmax"):
        if name in ("tmin", "tmax"):
            bound, inclusive = p
            selected = [v for v in x if (v >= bound if inclusive else v > bound)] if name == "tmin" else [
                v for v in x if (v <= bound if inclusive else v < bound)]
            return [min(selected) if name == "tmin" else max(selected)]
        low, high, incl_low, incl_high = p[:4]
        selected = [v for v in x if (v >= low if incl_low else v > low) and
                    (v <= high if incl_high else v < high)]
        if name == "tmean":
            return [_mean(selected)]
        variance = _variance(selected, int(p[4]))
        if name == "tvar":
            return [variance]
        if name == "tstd":
            return [mp.sqrt(variance)]
        return [mp.sqrt(variance/len(selected))]
    if name == "entropy":
        total = mp.fsum(x)
        probs = [v/total for v in x]
        if y:
            if len(x) != len(y):
                return None
            qtotal = mp.fsum(y)
            result = mp.fsum(v*mp.log(v/(w/qtotal)) for v,w in zip(probs,y) if v)
        else:
            result = -mp.fsum(v*mp.log(v) for v in probs if v)
        return [result/mp.log(p[0]) if p and p[0] else result]
    if name == "kstat":
        degree = int(p[0])
        n = len(x)
        v2 = _variance(x)
        if degree == 1:
            return [average]
        if degree == 2:
            return [n*v2/(n-1)]
        if degree == 3:
            return [n*n*mp.fsum((v-average)**3 for v in x)/(n*(n-1)*(n-2))]
        if degree == 4:
            m4 = mp.fsum((v-average)**4 for v in x)/n
            return [n*n*((n+1)*m4-3*(n-1)*v2*v2)/((n-1)*(n-2)*(n-3))]
        return None
    if name in ("circmean", "circvar", "circstd"):
        low, high = p[:2]
        period = high-low
        phase = [2*mp.pi*(v-low)/period for v in x]
        value = _mean([mp.exp(mp.j*t) for t in phase])
        if name == "circmean":
            return [low+(mp.arg(value) % (2*mp.pi))*period/(2*mp.pi)]
        if name == "circvar":
            return [1-abs(value)]
        return [mp.sqrt(-2*mp.log(abs(value)))*period/(2*mp.pi)]
    if name in ("dirichlet_pdf", "dirichlet_logpdf"):
        if not y or len(x) != len(y):
            return None
        logpdf = mp.loggamma(mp.fsum(y))-mp.fsum(mp.loggamma(v) for v in y)
        logpdf += mp.fsum((a-1)*mp.log(v) for v, a in zip(x, y))
        return [logpdf if name == "dirichlet_logpdf" else mp.exp(logpdf)]
    if name in ("multinomial_pmf", "multinomial_logpmf"):
        if not y or len(x) != len(y):
            return None
        total = sum(int(v) for v in x)
        logpdf = mp.loggamma(total+1)-mp.fsum(mp.loggamma(v+1) for v in x)
        logpdf += mp.fsum(v*mp.log(prob) for v, prob in zip(x, y) if v)
        return [logpdf if name == "multinomial_logpmf" else mp.exp(logpdf)]
    if name in ("energy_distance", "wasserstein_distance"):
        if not y:
            return None
        if name == "energy_distance":
            pairs = mp.fsum(abs(a-b) for a in x for b in y)/(len(x)*len(y))
            within_x = mp.fsum(abs(a-b) for a in x for b in x)/(len(x)**2)
            within_y = mp.fsum(abs(a-b) for a in y for b in y)/(len(y)**2)
            return [mp.sqrt(2*pairs-within_x-within_y)]
        # Equal-mass empirical distributions, including different sample sizes.
        left, right = sorted(x), sorted(y)
        boundaries = sorted(set(left+right))
        distance = mp.mpf(0)
        for lo, hi in pairwise(boundaries):
            distance += (hi-lo)*abs(
                sum(v <= lo for v in left)/mp.mpf(len(left))-
                sum(v <= lo for v in right)/mp.mpf(len(right)))
        return [distance]
    return None


def _integrate_values(name, args):
    if name == "fixed_quad":
        if set(args) != {"a", "b"}:
            return None
        a, b = _number(args["a"]), _number(args["b"])
        return [mp.exp(b)-mp.exp(a)]
    if name == "trapezoid":
        x, y = _vector(args["coordinates"]), _vector(args["y"])
        if len(x) != len(y):
            return None
        return [mp.fsum((x[i+1]-x[i])*(y[i]+y[i+1])/2 for i in range(len(x)-1))]
    if args.get("explicit") is not True:
        return None
    values = _vector(args["values"])
    dx = mp.mpf('0.25')
    if name == "cumulative_trapezoid":
        result = [mp.mpf(0)]
        for a, b in pairwise(values):
            result.append(result[-1]+dx*(a+b)/2)
        return result
    if name == "simpson" and len(values) % 2:
        return [dx/3*(values[0]+values[-1]+4*mp.fsum(values[1:-1:2])+
                      2*mp.fsum(values[2:-1:2]))]
    if name == "romb" and (len(values)-1) & (len(values)-2) == 0:
        n = len(values)-1
        levels = int(math.log2(n))
        previous = []
        for k in range(levels+1):
            step = n//(2**k)
            current = step*dx*(values[0]+values[-1]+2*mp.fsum(values[step:-1:step]))/2
            row = [current]
            for j in range(k):
                current = current+(current-previous[j])/(4**(j+1)-1)
                row.append(current)
            previous = row
        return [previous[-1]]
    if name == "cumulative_simpson" and len(values) >= 3:
        result = [mp.mpf(0)]
        for i in range(len(values)-1):
            if i+2 < len(values):
                f0, f1, f2 = values[i:i+3]
                increment = dx*(5*f0+8*f1-f2)/12
            else:
                f0, f1, f2 = values[i-1:i+2]
                increment = dx*(-f0+8*f1+5*f2)/12
            result.append(result[-1]+increment)
        return result
    return None


def _other_values(api, args):
    if api.startswith("differentiate."):
        name = api.rsplit(".", 1)[1]
        if name == "derivative" and "step" in args:
            return [mp.diff(mp.sin, _number(args["x"]))]
        if name == "gradient" and "step" in args:
            x, y = _vector(args["point"])
            f = lambda u, v: u*u+3*u*v+2*v*v
            return [mp.diff(lambda u: f(u,y),x), mp.diff(lambda v: f(x,v),y)]
        if args.get("explicit") is True:
            x, y = _vector(args["point"])
            if name == "jacobian":
                return [mp.diff(lambda u: u*y,x), mp.diff(lambda v: x*v,y), 1, 1]
            if name == "hessian":
                return [2, 3, 3, 4]
    if api.startswith("optimize."):
        name = api.rsplit(".", 1)[1]
        if name in {"minimize", "try_minimize", "minimize_with_gradient",
                    "try_minimize_with_gradient"} and args.get("explicit") is True:
            _vector(args["start"])
            return [1, -2, 3]
        if name in ("rosen", "rosen_der") and "point" in args:
            x = _vector(args["point"])
            if name == "rosen":
                return [mp.fsum(100*(b-a*a)**2+(1-a)**2 for a,b in pairwise(x))]
            gradient = [mp.mpf(0)]*len(x)
            for i in range(len(x)-1):
                gradient[i] += -400*x[i]*(x[i+1]-x[i]**2)+2*(x[i]-1)
                gradient[i+1] += 200*(x[i+1]-x[i]**2)
            return gradient
    if api in ("spatial.euclidean", "spatial.cosine") and "u" in args:
        u,v = _vector(args["u"]), _vector(args["v"])
        if len(u) != len(v) or not u:
            return None
        if api.endswith("euclidean"):
            return [mp.sqrt(mp.fsum((a-b)**2 for a,b in zip(u,v)))]
        return [1-mp.fdot(u,v)/(mp.sqrt(mp.fdot(u,u))*mp.sqrt(mp.fdot(v,v)))]
    if api in ("signal.convolve", "signal.correlate") and "x" in args:
        x,h = _vector(args["x"]), _vector(args["h"])
        if api.endswith("correlate"):
            h = list(reversed(h))
        return [mp.fsum(x[j]*h[i-j] for j in range(len(x)) if 0 <= i-j < len(h))
                for i in range(len(x)+len(h)-1)]
    if api in ("filter.kalman_predict", "filter.kalman_update") and "mean" in args:
        mean = mp.matrix(_vector(args["mean"]))
        covariance = _matrix(args["covariance"])
        if api.endswith("predict"):
            f, q = _matrix(args["f"]), _matrix(args["q"])
            return list(f*mean)+_flat_matrix(f*covariance*f.T+q)
        h,r,z = _matrix(args["h"]), _matrix(args["r"]), mp.matrix(_vector(args["z"]))
        k = covariance*h.T*(h*covariance*h.T+r)**-1
        residual = mp.eye(len(mean))-k*h
        return list(mean+k*(z-h*mean))+_flat_matrix(residual*covariance*residual.T+k*r*k.T)
    if api in ("cluster.vq.vq", "cluster.vq.whiten") and "obs" in args:
        obs = [_vector(row) for row in args["obs"]]
        if api.endswith("whiten"):
            cols = list(zip(*obs))
            scales = [mp.sqrt(_variance(col)) for col in cols]
            return [value/scales[j] for row in obs for j,value in enumerate(row)]
        centers = [_vector(row) for row in args["codebook"]]
        distances = [[mp.sqrt(mp.fsum((a-b)**2 for a,b in zip(row, center)))
                      for center in centers] for row in obs]
        indices = [min(range(len(row)), key=row.__getitem__) for row in distances]
        return indices+[row[index] for row,index in zip(distances,indices)]
    if api.startswith("constants."):
        name = api.rsplit(".", 1)[1]
        if name in ("lambda2nu", "nu2lambda") and "x" in args:
            return [mp.mpf(299792458)/_number(args["x"])]
        if name == "convert_temperature" and args.get("explicit") is True:
            return [_number(args["x"])*mp.mpf(9)/5+32]
    if api in ("fft.fft", "fft.ifft") and "z" in args:
        z = [_complex(pair) for pair in args["z"]]
        inverse = api.endswith("ifft")
        n = len(z)
        return [_components(mp.fsum(
            z[j]*mp.exp((1 if inverse else -1)*2*mp.pi*mp.j*k*j/n)
            for j in range(n))/(n if inverse else 1)) for k in range(n)]
    return None


class MpmathGeneralComparisons:
    """Public mpmath oracles, with a single finite numeric output per API."""

    def evaluate(self, case):
        api = case.get("api", "")
        args = case.get("args", {})
        if not isinstance(api, str) or not isinstance(args, dict):
            return {}
        name = api.rsplit(".", 1)[-1]
        eligible = (api.startswith("complex.") and
                    name in _COMPLEX_UNARY.keys() | _COMPLEX_BINARY | _COMPLEX_OTHER or
                    api.startswith("linalg.") and name in _LINALG or
                    api.startswith("stats.") and name in _STATS or
                    api.startswith("integrate.") and name in
                    {"fixed_quad", "trapezoid", "cumulative_trapezoid",
                     "simpson", "romb", "cumulative_simpson"} or
                    api.startswith(("differentiate.", "optimize.")) or
                    api in {"spatial.euclidean", "spatial.cosine", "signal.convolve",
                            "signal.correlate", "cluster.vq.vq", "cluster.vq.whiten",
                            "filter.kalman_predict", "filter.kalman_update",
                            "constants.lambda2nu", "constants.nu2lambda",
                            "constants.convert_temperature", "fft.fft", "fft.ifft"})
        if not eligible:
            return {}
        try:
            with mp.workdps(90):
                if api.startswith("complex."):
                    values = _complex_values(name, args)
                elif api.startswith("linalg."):
                    values = _linalg_values(name, args)
                elif api.startswith("stats."):
                    values = _stats_values(name, args)
                elif api.startswith("integrate."):
                    values = _integrate_values(name, args)
                else:
                    values = _other_values(api, args)
                if values is None:
                    return {}
                values = [float(v) for v in _flatten(values)]
            if not values or not all(math.isfinite(v) for v in values):
                return {}
            return {PROVIDER: {"values": values}}
        except (ArithmeticError, ValueError, TypeError, KeyError, IndexError,
                OverflowError, ZeroDivisionError, AttributeError,
                mp.libmp.NoConvergence):
            return {}


def evaluate(case):
    """Convenience entry point for callers that do not keep an adapter object."""
    return MpmathGeneralComparisons().evaluate(case)
