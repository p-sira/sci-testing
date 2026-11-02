//! Independent numerical comparisons via the public RustFFT, nalgebra and
//! statrs APIs. The numeric opcodes are part of the Python adapter contract.

use nalgebra::{DMatrix, DVector};
use rustfft::{FftPlanner, num_complex::Complex};
use statrs::distribution::{
    Bernoulli, Beta, Binomial, Cauchy, Chi, ChiSquared, Continuous, ContinuousCDF, Discrete,
    DiscreteCDF, Erlang, Exp, FisherSnedecor, Gamma, Geometric, Gumbel, InverseGamma, Laplace,
    LogNormal, NegativeBinomial, Normal, Pareto, Poisson, StudentsT, Uniform, Weibull,
};
use std::panic::{AssertUnwindSafe, catch_unwind};

fn fft_axis(data: &mut [Complex<f64>], shape: &[usize], axis: usize, inverse: bool) {
    let length = shape[axis];
    let stride: usize = shape[axis + 1..].iter().product();
    let outer = data.len() / (length * stride);
    let mut planner = FftPlanner::<f64>::new();
    let plan = if inverse {
        planner.plan_fft_inverse(length)
    } else {
        planner.plan_fft_forward(length)
    };
    let mut lane = vec![Complex::new(0.0, 0.0); length];
    for block in 0..outer {
        for offset in 0..stride {
            let base = block * length * stride + offset;
            for k in 0..length {
                lane[k] = data[base + k * stride];
            }
            plan.process(&mut lane);
            let factor = if inverse { 1.0 / length as f64 } else { 1.0 };
            for k in 0..length {
                data[base + k * stride] = lane[k] * factor;
            }
        }
    }
}

fn dense_fft(
    op: u32,
    input: &[f64],
    ndim: usize,
    axes_len: usize,
    metadata: &[f64],
) -> Option<Vec<f64>> {
    // Shape and axis metadata are integer-valued, finite, compact and explicit.
    if ndim == 0 || ndim > 4 || axes_len > ndim || metadata.len() != ndim + axes_len {
        return None;
    }
    let integers = metadata
        .iter()
        .map(|&x| {
            if x >= 0.0 && x <= 131_072.0 && x.fract() == 0.0 {
                Some(x as usize)
            } else {
                None
            }
        })
        .collect::<Option<Vec<_>>>()?;
    let shape = &integers[..ndim];
    let axes = &integers[ndim..];
    if shape.contains(&0) || axes.iter().any(|&a| a >= ndim) {
        return None;
    }
    let mut used = [false; 4];
    for &axis in axes {
        if used[axis] {
            return None;
        }
        used[axis] = true;
    }
    let n: usize = shape
        .iter()
        .try_fold(1usize, |acc, s| acc.checked_mul(*s))?;
    if n > 65_536 {
        return None;
    }
    let forward_real = matches!(op, 16 | 22 | 24);
    let inverse_real = matches!(op, 17 | 23 | 25);
    if (forward_real || inverse_real) && axes.is_empty() {
        return None;
    }
    let real_axis = if forward_real || inverse_real {
        Some(*axes.last()?)
    } else {
        None
    };
    let packed_shape = if let Some(axis) = real_axis {
        let mut dims = shape.to_vec();
        dims[axis] = shape[axis] / 2 + 1;
        dims
    } else {
        shape.to_vec()
    };
    let packed_n: usize = packed_shape.iter().product();
    let input_n = if inverse_real { packed_n } else { n };
    let input_real = forward_real;
    if input.len() != input_n * if input_real { 1 } else { 2 } {
        return None;
    }
    let mut buffer = vec![Complex::new(0.0, 0.0); n];
    let mut strides = vec![1; ndim];
    for axis in (0..ndim - 1).rev() {
        strides[axis] = strides[axis + 1] * shape[axis + 1];
    }
    if inverse_real {
        // Inverse on the ordinary complex axes comes before Hermitian
        // reconstruction, matching the public irfftn contract.
        let mut packed = input
            .chunks_exact(2)
            .map(|v| Complex::new(v[0], v[1]))
            .collect::<Vec<_>>();
        let axis = real_axis?;
        for &a in axes {
            if a != axis {
                fft_axis(&mut packed, &packed_shape, a, true);
            }
        }
        let stride = strides[axis];
        for index in 0..n {
            let frequency = (index / stride) % shape[axis];
            let paired = frequency > shape[axis] / 2;
            let src_frequency = if paired {
                shape[axis] - frequency
            } else {
                frequency
            };
            let low = index % stride;
            let high = index / (stride * shape[axis]);
            let src = high * (shape[axis] / 2 + 1) * stride + src_frequency * stride + low;
            buffer[index] = if paired {
                packed[src].conj()
            } else {
                packed[src]
            };
        }
        fft_axis(&mut buffer, shape, axis, true);
        return Some(buffer.iter().flat_map(|x| [x.re, 0.0]).collect());
    }
    for (i, chunk) in input
        .chunks_exact(if input_real { 1 } else { 2 })
        .enumerate()
    {
        buffer[i] = Complex::new(chunk[0], if input_real { 0.0 } else { chunk[1] });
    }
    let inverse = matches!(op, 19 | 21);
    for &axis in axes {
        fft_axis(&mut buffer, shape, axis, inverse);
    }
    if forward_real {
        let axis = real_axis?;
        let stride = strides[axis];
        let mut result = Vec::with_capacity(packed_n * 2);
        for high in 0..n / (shape[axis] * stride) {
            for frequency in 0..=shape[axis] / 2 {
                for low in 0..stride {
                    let z = buffer[high * shape[axis] * stride + frequency * stride + low];
                    result.extend([z.re, z.im]);
                }
            }
        }
        return Some(result);
    }
    Some(buffer.iter().flat_map(|z| [z.re, z.im]).collect())
}

fn distribution(op: u32, observations: &[f64], id: usize, params: &[f64]) -> Option<Vec<f64>> {
    if params.iter().any(|v| !v.is_finite()) {
        return None;
    }
    // Distribution opcodes: 0=PDF, 1=log PDF, 2=CDF, 3=SF,
    // 4=PPF. Constructors are public statrs 0.19.1 APIs.
    macro_rules! evaluate {
        ($value:expr) => {{
            let dist = $value.ok()?;
            observations
                .iter()
                .map(|&x| match op {
                    0 => Some(dist.pdf(x)),
                    1 => Some(dist.ln_pdf(x)),
                    2 => Some(dist.cdf(x)),
                    3 => Some(dist.sf(x)),
                    4 if (0.0..=1.0).contains(&x) => Some(dist.inverse_cdf(x)),
                    _ => None,
                })
                .collect::<Option<Vec<_>>>()
        }};
    }
    macro_rules! discrete {
        ($value:expr) => {{
            let dist = $value.ok()?;
            observations
                .iter()
                .map(|&x| {
                    if op == 4 {
                        return (0.0 < x && x < 1.0).then(|| dist.inverse_cdf(x) as f64);
                    }
                    if x < 0.0 {
                        return match op {
                            0 => Some(0.0),
                            2 => Some(0.0),
                            3 => Some(1.0),
                            _ => None,
                        };
                    }
                    if x > u64::MAX as f64 {
                        return None;
                    }
                    let k = x as u64;
                    match op {
                        0 => Some(if x.fract() == 0.0 { dist.pmf(k) } else { 0.0 }),
                        1 if x.fract() == 0.0 => Some(dist.ln_pmf(k)),
                        2 => Some(dist.cdf(k)),
                        3 => Some(dist.sf(k)),
                        _ => None,
                    }
                })
                .collect::<Option<Vec<_>>>()
        }};
    }
    match (id, params) {
        (1, [mean, sd]) => evaluate!(Normal::new(*mean, *sd)),
        (2, [a, b]) => evaluate!(Beta::new(*a, *b)),
        (3, [shape, scale]) if *scale > 0.0 => evaluate!(Gamma::new(*shape, 1.0 / scale)),
        (4, [df]) => evaluate!(ChiSquared::new(*df)),
        (5, [df]) if *df >= 1.0 && df.fract() == 0.0 && *df <= 65536.0 => {
            evaluate!(Chi::new(*df as u64))
        }
        (6, [df]) => evaluate!(StudentsT::new(0.0, 1.0, *df)),
        (7, [df1, df2]) => evaluate!(FisherSnedecor::new(*df1, *df2)),
        (8, [scale]) if *scale > 0.0 => evaluate!(Exp::new(1.0 / scale)),
        (9, [location, scale]) => evaluate!(Cauchy::new(*location, *scale)),
        (10, [location, scale]) => evaluate!(Laplace::new(*location, *scale)),
        (11, [location, scale]) => evaluate!(LogNormal::new(*location, *scale)),
        (12, [low, high]) => evaluate!(Uniform::new(*low, *high)),
        (13, [shape, scale]) => evaluate!(Weibull::new(*shape, *scale)),
        (14, [scale, shape]) => evaluate!(Pareto::new(*scale, *shape)),
        (15, [location, scale]) => evaluate!(Gumbel::new(*location, *scale)),
        (16, [shape, scale])
            if *scale > 0.0
                && (op != 4
                    || observations
                        .iter()
                        .all(|&q| (1e-10..=1.0 - 1e-10).contains(&q))) =>
        {
            evaluate!(InverseGamma::new(*shape, *scale))
        }
        (17, [shape, scale])
            if *shape >= 1.0 && shape.fract() == 0.0 && *shape <= 65536.0 && *scale > 0.0 =>
        {
            evaluate!(Erlang::new(*shape as u64, 1.0 / scale))
        }
        (18, [p]) => discrete!(Bernoulli::new(*p)),
        (19, [n, p]) if *n >= 0.0 && n.fract() == 0.0 && *n <= 65536.0 => {
            discrete!(Binomial::new(*p, *n as u64))
        }
        (20, [rate]) => discrete!(Poisson::new(*rate)),
        (21, [p]) if op != 4 || *p < 1.0 => discrete!(Geometric::new(*p)),
        (22, [r, p]) => discrete!(NegativeBinomial::new(*r, *p)),
        _ => None,
    }
}

fn compute(op: u32, input: &[f64], rows: usize, cols: usize, extra: &[f64]) -> Option<Vec<f64>> {
    if input.is_empty()
        || input.len() > 131_072
        || extra.len() > 131_072
        || !input.iter().chain(extra).all(|x| x.is_finite())
    {
        return None;
    }
    match op {
        // FFT inputs and outputs are interleaved [real, imaginary] pairs.
        1 | 2 => {
            if input.len() % 2 != 0 || !extra.is_empty() {
                return None;
            }
            let n = input.len() / 2;
            let mut data: Vec<Complex<f64>> = input
                .chunks_exact(2)
                .map(|pair| Complex::new(pair[0], pair[1]))
                .collect();
            let mut planner = FftPlanner::<f64>::new();
            let fft = if op == 1 {
                planner.plan_fft_forward(n)
            } else {
                planner.plan_fft_inverse(n)
            };
            fft.process(&mut data);
            let scale = if op == 2 { 1.0 / n as f64 } else { 1.0 };
            Some(
                data.iter()
                    .flat_map(|z| [z.re * scale, z.im * scale])
                    .collect(),
            )
        }
        // The standard normal distribution matches the saved standard
        // normal Sci contract (zero mean, unit standard deviation).
        3 | 4 | 15 => {
            if input.len() != 1 || !extra.is_empty() {
                return None;
            }
            let normal = Normal::new(0.0, 1.0).ok()?;
            if op == 4 && !(0.0 < input[0] && input[0] < 1.0) {
                return None;
            }
            Some(vec![if op == 3 {
                normal.cdf(input[0])
            } else if op == 4 {
                normal.inverse_cdf(input[0])
            } else {
                normal.pdf(input[0])
            }])
        }
        // General dense matrix primitives; row-major data and outputs.
        5..=14 => {
            if op == 8 {
                return (input.len() == extra.len() && !extra.is_empty())
                    .then(|| vec![input.iter().zip(extra).map(|(x, y)| x * y).sum()]);
            }
            if rows == 0 || cols == 0 || rows > 64 || cols > 64 || input.len() != rows * cols {
                return None;
            }
            let a = DMatrix::from_row_slice(rows, cols, input);
            match op {
                5 if rows == cols && extra.is_empty() => Some(vec![a.determinant()]),
                6 if rows == cols && extra.len() == rows => {
                    let b = DVector::from_row_slice(extra);
                    Some(a.lu().solve(&b)?.iter().copied().collect())
                }
                7 if rows == cols && extra.is_empty() => {
                    let inverse = a.try_inverse()?;
                    Some(matrix_rows(&inverse))
                }
                9 if extra.len() == cols => {
                    let vector = DVector::from_row_slice(extra);
                    Some((a * vector).iter().copied().collect())
                }
                10 if extra.len() % cols == 0 && !extra.is_empty() => {
                    let b_cols = extra.len() / cols;
                    if b_cols > 64 {
                        return None;
                    }
                    let b = DMatrix::from_row_slice(cols, b_cols, extra);
                    Some(matrix_rows(&(a * b)))
                }
                11 if rows == cols && extra.len() == 1 => {
                    let lower = extra[0] == 1.0;
                    let factor = a.cholesky()?.l();
                    Some(matrix_rows(&if lower {
                        factor
                    } else {
                        factor.transpose()
                    }))
                }
                12 if rows == cols && extra.len() == rows + 1 => {
                    let vector = DVector::from_row_slice(&extra[..rows]);
                    let solution = a.cholesky()?.solve(&vector);
                    Some(solution.iter().copied().collect())
                }
                13 if rows == cols && extra.len() == 1 => {
                    // SciPy's lower flag selects the authoritative triangle.
                    let lower = extra[0] == 1.0;
                    let mut symmetric = a;
                    for i in 0..rows {
                        for j in (i + 1)..cols {
                            if lower {
                                symmetric[(i, j)] = symmetric[(j, i)];
                            } else {
                                symmetric[(j, i)] = symmetric[(i, j)];
                            }
                        }
                    }
                    let mut eigenvalues: Vec<f64> = symmetric
                        .symmetric_eigen()
                        .eigenvalues
                        .iter()
                        .copied()
                        .collect();
                    eigenvalues.sort_by(f64::total_cmp);
                    Some(eigenvalues)
                }
                14 if extra.is_empty() => {
                    let mut singular_values: Vec<f64> = a
                        .svd(false, false)
                        .singular_values
                        .iter()
                        .copied()
                        .collect();
                    singular_values.sort_by(|a, b| b.total_cmp(a));
                    Some(singular_values)
                }
                _ => None,
            }
        }
        16..=25 => dense_fft(op, input, rows, cols, extra),
        30 => distribution(cols as u32, input, rows, extra),
        _ => None,
    }
}

fn matrix_rows(matrix: &DMatrix<f64>) -> Vec<f64> {
    (0..matrix.nrows())
        .flat_map(|i| (0..matrix.ncols()).map(move |j| matrix[(i, j)]))
        .collect()
}

/// # Safety
/// All input and output pointers must be valid for the specified f64 lengths,
/// which are validated by the calling Python ctypes adapter.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rust_oracle_evaluate(
    op: u32,
    input: *const f64,
    input_len: usize,
    rows: usize,
    cols: usize,
    extra: *const f64,
    extra_len: usize,
    output: *mut f64,
    output_len: usize,
) -> i32 {
    if input.is_null() || output.is_null() || (extra_len != 0 && extra.is_null()) {
        return 1;
    }
    // Refuse pathological lengths before creating any slices.
    if input_len > 131_072 || extra_len > 131_072 || output_len > 131_072 {
        return 1;
    }
    let result = catch_unwind(AssertUnwindSafe(|| {
        // SAFETY: The Python FFI provides arrays with the stated lengths.
        let inputs = unsafe { std::slice::from_raw_parts(input, input_len) };
        let extras = if extra_len == 0 {
            &[]
        } else {
            unsafe { std::slice::from_raw_parts(extra, extra_len) }
        };
        compute(op, inputs, rows, cols, extras)
    }));
    match result {
        Ok(Some(values)) if values.len() == output_len && values.iter().all(|v| v.is_finite()) => {
            // SAFETY: Python allocated output_len writable f64 elements.
            unsafe { std::ptr::copy_nonoverlapping(values.as_ptr(), output, output_len) };
            0
        }
        Ok(_) => 1,  // Unsupported numerical domain or shape.
        Err(_) => 2, // Unexpected internal failure; report to caller.
    }
}
