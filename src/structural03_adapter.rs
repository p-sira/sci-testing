//! Sci 0.3 structural/spectral comparisons through the public Sci facade.
//!
//! The Sci facade does not currently expose `sci::core`; core CSR tests call
//! the public `sci_core::csr` exports directly through a dedicated dependency.

use sci::{SparseMatrix, signal, sparse};
use sci_core::csr;
use serde_json::Value;
use std::error::Error;

type Result<T> = std::result::Result<T, Box<dyn Error>>;

fn integer(value: &Value) -> Result<usize> {
    Ok(usize::try_from(
        value.as_u64().ok_or("expected nonnegative integer")?,
    )?)
}

fn float(value: &Value) -> Result<f64> {
    value.as_f64().ok_or_else(|| "expected float".into())
}

fn ints(value: &Value) -> Result<Vec<usize>> {
    value
        .as_array()
        .ok_or("expected integer array")?
        .iter()
        .map(integer)
        .collect()
}

fn floats(value: &Value) -> Result<Vec<f64>> {
    value
        .as_array()
        .ok_or("expected float array")?
        .iter()
        .map(float)
        .collect()
}

fn triples(value: &Value) -> Result<Vec<(usize, usize, f64)>> {
    value
        .as_array()
        .ok_or("expected coordinate array")?
        .iter()
        .map(|entry| {
            let fields = entry.as_array().ok_or("expected COO triple")?;
            if fields.len() != 3 {
                return Err("expected three COO fields".into());
            }
            Ok((
                integer(&fields[0])?,
                integer(&fields[1])?,
                float(&fields[2])?,
            ))
        })
        .collect()
}

fn storage(a: &Value) -> Result<SparseMatrix> {
    Ok(SparseMatrix {
        rows: integer(&a["rows"])?,
        cols: integer(&a["cols"])?,
        indptr: ints(&a["indptr"])?,
        indices: ints(&a["indices"])?,
        data: floats(&a["data"])?,
    })
}

fn float_pair(out: &mut Vec<f64>, value: f64) {
    out.push(value);
    out.push(f64::from(u8::from(value == 0. && value.is_sign_negative())));
}

fn matrix_storage(matrix: &SparseMatrix) -> Vec<f64> {
    let mut out = vec![
        matrix.rows as f64,
        matrix.cols as f64,
        matrix.indptr.len() as f64,
    ];
    out.extend(matrix.indptr.iter().map(|&v| v as f64));
    out.push(matrix.indices.len() as f64);
    out.extend(matrix.indices.iter().map(|&v| v as f64));
    out.push(matrix.data.len() as f64);
    for &value in &matrix.data {
        float_pair(&mut out, value);
    }
    out
}

fn index_storage<T>(indptr: &[T], indices: &[T]) -> Vec<f64>
where
    T: Copy + Into<i64>,
{
    let mut out = vec![indptr.len() as f64];
    out.extend(indptr.iter().map(|&v| v.into() as f64));
    out.push(indices.len() as f64);
    out.extend(indices.iter().map(|&v| v.into() as f64));
    out
}

fn evaluate_success(api: &str, a: &Value) -> Result<Vec<f64>> {
    if api == "signal.utilities.closest_stft_dual_window" {
        let win = floats(&a["win"])?;
        let desired = floats(&a["desired"])?;
        let hop = integer(&a["hop"])?;
        let scaled = a["scaled"].as_bool().ok_or("missing scaled flag")?;
        let (dual, alpha) =
            signal::utilities::closest_stft_dual_window(&win, hop, &desired, scaled)?;
        let mut out = vec![dual.len() as f64, alpha];
        out.extend(dual);
        return Ok(out);
    }

    if api == "core.csr.from_coo" {
        let rows = integer(&a["rows"])?;
        let cols = integer(&a["cols"])?;
        let matrix = csr::from_coo(rows, cols, &triples(&a["coordinates"])?)?;
        return Ok(matrix_storage(&matrix));
    }

    let matrix = storage(a)?;
    Ok(match api {
        "core.csr.validate" => {
            csr::validate(&matrix)?;
            vec![1.]
        }
        "core.csr.to_dense" => {
            let dense = csr::to_dense(&matrix)?;
            let mut out = vec![matrix.rows as f64, matrix.cols as f64];
            for row in dense {
                for value in row {
                    float_pair(&mut out, value);
                }
            }
            out
        }
        "sparse.array.csr_indices_i32" => {
            let indices = sparse::array::csr_indices_i32(&matrix)?;
            index_storage(&indices.indptr, &indices.indices)
        }
        "sparse.array.csr_indices_i64" => {
            let indices = sparse::array::csr_indices_i64(&matrix)?;
            index_storage(&indices.indptr, &indices.indices)
        }
        _ => return Err(format!("unregistered structural03 API: {api}").into()),
    })
}

/// Evaluate the six Sci 0.3 structural exports; invalid fixtures verify errors.
pub fn evaluate(api: &str, a: &Value) -> Result<Vec<f64>> {
    match a["operation"]
        .as_str()
        .ok_or("missing structural03 operation")?
    {
        "success" => evaluate_success(api, a),
        "reject" => Ok(vec![f64::from(u8::from(evaluate_success(api, a).is_err()))]),
        value => Err(format!("unknown structural03 operation: {value}").into()),
    }
}
