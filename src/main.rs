//! Evaluates fixed cases through the public Sci facade; references live in Python.
use sci::{
    cluster,
    complex::{self, Complex},
    constants, differentiate, fft, filter, integrate, interpolate, linalg, optimize, signal,
    spatial, special, stats,
};
use serde_json::{Value, json};
use std::{
    error::Error,
    io::{self, BufRead},
};

type Result<T> = std::result::Result<T, Box<dyn Error>>;
fn vector(v: &Value) -> Result<Vec<f64>> {
    Ok(serde_json::from_value(v.clone())?)
}
fn matrix(v: &Value) -> Result<Vec<Vec<f64>>> {
    Ok(serde_json::from_value(v.clone())?)
}
fn scalar(v: &Value) -> Result<f64> {
    v.as_f64().ok_or_else(|| "expected number".into())
}
fn flatten(z: impl IntoIterator<Item = Complex>) -> Vec<f64> {
    z.into_iter().flat_map(|z| [z.re, z.im]).collect()
}
fn gaussian(g: filter::Gaussian) -> Vec<f64> {
    g.mean
        .into_iter()
        .chain(g.covariance.into_iter().flatten())
        .collect()
}
fn evaluate(case: &Value) -> Result<Vec<f64>> {
    let a = &case["args"];
    let x = || scalar(&a["x"]);
    let v = |key: &str| vector(&a[key]);
    let m = |key: &str| matrix(&a[key]);
    let api = case["api"].as_str().ok_or("missing API")?;
    Ok(match api {
        "cluster.vq" => {
            let (codes, distances) = cluster::vq::vq(&m("obs")?, &m("codebook")?)?;
            codes
                .into_iter()
                .map(|c| c as f64)
                .chain(distances)
                .collect()
        }
        "cluster.whiten" => cluster::vq::whiten(&m("obs")?)?
            .into_iter()
            .flatten()
            .collect(),
        "complex.exp" | "complex.sqrt" => {
            let z = Complex {
                re: x()?,
                im: scalar(&a["y"])?,
            };
            flatten([if api == "complex.exp" {
                complex::exp(z)
            } else {
                complex::sqrt(z)
            }])
        }
        "constants.lambda2nu" => vec![constants::lambda2nu(x()?)],
        "constants.nu2lambda" => vec![constants::nu2lambda(x()?)],
        "differentiate.derivative" => vec![differentiate::derivative(
            f64::sin,
            x()?,
            scalar(&a["step"])?,
        )?],
        "differentiate.gradient" => differentiate::gradient(
            |p| p[0] * p[0] + 3.0 * p[0] * p[1] + 2.0 * p[1] * p[1],
            &v("point")?,
            scalar(&a["step"])?,
        )?,
        "fft.fft" | "fft.ifft" => {
            let z: Vec<_> = m("z")?
                .into_iter()
                .map(|r| Complex { re: r[0], im: r[1] })
                .collect();
            flatten(if api == "fft.fft" {
                fft::fft(&z)
            } else {
                fft::ifft(&z)
            })
        }
        "filter.kalman_predict" | "filter.kalman_update" => {
            let g = filter::Gaussian {
                mean: v("mean")?,
                covariance: m("covariance")?,
            };
            gaussian(if api == "filter.kalman_predict" {
                filter::kalman_predict(&g, &m("f")?, &m("q")?)?
            } else {
                filter::kalman_update(&g, &v("z")?, &m("h")?, &m("r")?)?.estimate
            })
        }
        "integrate.trapezoid" => vec![integrate::trapezoid(
            &v("y")?,
            Some(&v("coordinates")?),
            1.0,
        )?],
        "integrate.fixed_quad" => vec![integrate::fixed_quad(
            f64::exp,
            scalar(&a["a"])?,
            scalar(&a["b"])?,
            16,
        )?],
        "interpolate.barycentric_interpolate" => {
            interpolate::barycentric_interpolate(&v("nodes")?, &v("values")?, &v("query")?)?
        }
        "interpolate.krogh_interpolate" => {
            interpolate::krogh_interpolate(&v("nodes")?, &v("values")?, &v("query")?)?
        }
        "linalg.solve" => linalg::solve(&m("a")?, &v("b")?)?,
        "linalg.det" => vec![linalg::det(&m("a")?)?],
        "optimize.rosen" => vec![optimize::rosen(&v("point")?)?],
        "optimize.rosen_der" => optimize::rosen_der(&v("point")?)?,
        "signal.convolve" => signal::convolve(&v("x")?, &v("h")?, sci::Mode::Full)?,
        "signal.correlate" => signal::correlate(&v("x")?, &v("h")?, sci::Mode::Full)?,
        "spatial.euclidean" => vec![spatial::euclidean(&v("u")?, &v("v")?, None)?],
        "spatial.cosine" => vec![spatial::cosine(&v("u")?, &v("v")?, None)?],
        "special.erf" => vec![special::erf(x()?)],
        "special.gamma" => vec![special::gamma(x()?)],
        "stats.normal_cdf" => vec![stats::normal_cdf(x()?, 0.0, 1.0)?],
        "stats.variance" => vec![stats::variance(&v("values")?, 1)?],
        _ => return Err(format!("unsupported API: {api}").into()),
    })
}
fn main() -> Result<()> {
    for line in io::stdin().lock().lines() {
        let case: Value = serde_json::from_str(&line?)?;
        let result = match evaluate(&case) {
            Ok(values) if values.iter().all(|v| v.is_finite()) => {
                json!({"id":case["id"],"values":values})
            }
            Ok(_) => json!({"id":case["id"],"error":"nonfinite output"}),
            Err(error) => json!({"id":case["id"],"error":error.to_string()}),
        };
        println!("{result}");
    }
    Ok(())
}
