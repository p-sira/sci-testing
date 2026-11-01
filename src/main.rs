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

#[path = "../sci/tests/common/mod.rs"]
mod common;
mod extra;
mod extracted;
mod fixture_adapters;
mod generated;
mod interpolate_adapter;
mod io_adapter;
mod linalg_adapter;
mod numeric;
mod sparse_adapter;
mod special03_adapter;
mod special_complex03_adapter;
mod structural03_adapter;

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
    if a["special03_case"].as_bool() == Some(true) {
        return special03_adapter::evaluate(case["api"].as_str().ok_or("missing API")?, a);
    }
    if a["structural03_case"].as_bool() == Some(true) {
        return structural03_adapter::evaluate(case["api"].as_str().ok_or("missing API")?, a);
    }
    if a["special_complex03_case"].as_bool() == Some(true) {
        return special_complex03_adapter::evaluate(case["api"].as_str().ok_or("missing API")?, a);
    }
    if a["io_case"].as_bool() == Some(true) {
        return io_adapter::evaluate(case["api"].as_str().ok_or("missing API")?, a);
    }
    if a["explicit"].as_bool() == Some(true) {
        return extra::evaluate(case["api"].as_str().ok_or("missing API")?, a);
    }
    if let Some(fields) = a["fields"].as_array() {
        let fields = fields
            .iter()
            .map(|v| v.as_str().ok_or("expected string field"))
            .collect::<std::result::Result<Vec<_>, _>>()?;
        let api = case["api"].as_str().ok_or("missing API")?;
        let values = match a["suite"].as_str().ok_or("missing fixture suite")? {
            "interpolation" => interpolate_adapter::run(&fields),
            "linalg" => linalg_adapter::run(&fields),
            "sparse" => sparse_adapter::core::run(&fields),
            "sparse_parity" => sparse_adapter::parity::run(&fields),
            "fft" => extracted::fft::run(&fields),
            "periodic" => extracted::periodic::run(&fields),
            "cluster" => extracted::cluster::run(&fields),
            "spatial" => extracted::spatial::run(&fields),
            "optimize" => extracted::optimize::run(&fields),
            "differentiate" => extracted::differentiate::run(&fields),
            "filter" => extracted::filter::run(&fields),
            _ => return Err(format!("unregistered fixture suite: {api}").into()),
        };
        return Ok(if api.starts_with("sparse.") {
            values
                .into_iter()
                .flat_map(|v| {
                    if v.is_finite() {
                        [0., v]
                    } else if v.is_nan() {
                        [2., 0.]
                    } else {
                        [v.signum(), 0.]
                    }
                })
                .collect()
        } else {
            values
        });
    }
    if a.get("arguments").is_some() || a.get("p").is_some() {
        return fixture_adapters::evaluate(case["api"].as_str().ok_or("missing API")?, a);
    }
    if !a["p0"].is_null() {
        return generated::evaluate(case["api"].as_str().ok_or("missing API")?, a);
    }
    let x = || scalar(&a["x"]);
    let v = |key: &str| vector(&a[key]);
    let m = |key: &str| matrix(&a[key]);
    let api = case["api"].as_str().ok_or("missing API")?;
    Ok(match api {
        "cluster.vq.vq" => {
            let (codes, distances) = cluster::vq::vq(&m("obs")?, &m("codebook")?)?;
            codes
                .into_iter()
                .map(|c| c as f64)
                .chain(distances)
                .collect()
        }
        "cluster.vq.whiten" => cluster::vq::whiten(&m("obs")?)?
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
                fft::fourier::fft(&z)
            } else {
                fft::fourier::ifft(&z)
            })
        }
        "filter.kalman_predict" | "filter.kalman_update" => {
            let g = filter::Gaussian {
                mean: v("mean")?,
                covariance: m("covariance")?,
            };
            gaussian(if api == "filter.kalman_predict" {
                filter::kalman::kalman_predict(&g, &m("f")?, &m("q")?)?
            } else {
                filter::kalman::kalman_update(&g, &v("z")?, &m("h")?, &m("r")?)?.estimate
            })
        }
        "integrate.trapezoid" => vec![integrate::sampled::trapezoid(
            &v("y")?,
            Some(&v("coordinates")?),
            1.0,
        )?],
        "integrate.fixed_quad" => vec![integrate::quadrature::fixed_quad(
            f64::exp,
            scalar(&a["a"])?,
            scalar(&a["b"])?,
            16,
        )?],
        "interpolate.barycentric_interpolate" => interpolate::polynomial::barycentric_interpolate(
            &v("nodes")?,
            &v("values")?,
            &v("query")?,
        )?,
        "interpolate.krogh_interpolate" => {
            interpolate::polynomial::krogh_interpolate(&v("nodes")?, &v("values")?, &v("query")?)?
        }
        "linalg.solve" => linalg::solvers::solve(&m("a")?, &v("b")?)?,
        "linalg.det" => vec![linalg::basic::det(&m("a")?)?],
        "optimize.rosen" => vec![optimize::derivatives::rosen(&v("point")?)?],
        "optimize.rosen_der" => optimize::derivatives::rosen_der(&v("point")?)?,
        "signal.convolve" => signal::convolution::convolve(&v("x")?, &v("h")?, sci::Mode::Full)?,
        "signal.correlate" => signal::convolution::correlate(&v("x")?, &v("h")?, sci::Mode::Full)?,
        "spatial.euclidean" => vec![spatial::distance::euclidean(&v("u")?, &v("v")?, None)?],
        "spatial.cosine" => vec![spatial::distance::cosine(&v("u")?, &v("v")?, None)?],
        "special.erf" => vec![special::elementary::erf(x()?)],
        "special.gamma" => vec![special::gamma::gamma(x()?)],
        "stats.normal_cdf" => vec![stats::distributions::normal_cdf(x()?, 0.0, 1.0)?],
        "stats.variance" => vec![stats::descriptive::variance(&v("values")?, 1)?],
        _ => return Err(format!("unsupported API: {api}").into()),
    })
}
fn main() -> Result<()> {
    std::panic::set_hook(Box::new(|_| {}));
    for line in io::stdin().lock().lines() {
        let case: Value = serde_json::from_str(&line?)?;
        let evaluation = std::panic::catch_unwind(|| evaluate(&case)).unwrap_or_else(|payload| {
            let message = payload
                .downcast_ref::<String>()
                .map(String::as_str)
                .or_else(|| payload.downcast_ref::<&str>().copied())
                .unwrap_or("adapter panic");
            Err(message.to_owned().into())
        });
        let result = match evaluation {
            Ok(values) if values.iter().all(|v| v.is_finite()) => {
                json!({"id":case["id"],"values":values})
            }
            Ok(values) => {
                let values: Vec<Value> = values
                    .into_iter()
                    .map(|v| {
                        if v.is_finite() {
                            json!(v)
                        } else {
                            json!(v.to_string())
                        }
                    })
                    .collect();
                json!({"id":case["id"],"error":"nonfinite output","values":values})
            }
            Err(error) => json!({"id":case["id"],"error":error.to_string()}),
        };
        println!("{result}");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn json_preserves_adjacent_floats() {
        for text in ["-0.9999999999999999", "-1.9999999999999998"] {
            let expected: f64 = text.parse().unwrap();
            let parsed: Value = serde_json::from_str(text).unwrap();
            assert_eq!(scalar(&parsed).unwrap().to_bits(), expected.to_bits());
        }
    }

    #[test]
    fn gamma_near_poles_survives_json_transport() {
        // Adjacent binary64 inputs must reach Sci unchanged. Parsing them as
        // exact negative integers instead turns finite gamma values into NaNs.
        for text in ["-0.9999999999999999", "-1.9999999999999998"] {
            let case: Value = serde_json::from_str(&format!(
                r#"{{"api":"special.gamma","args":{{"x":{text}}}}}"#
            ))
            .unwrap();
            let expected: f64 = special::gamma::gamma(text.parse().unwrap());
            assert!(expected.is_finite());
            assert_eq!(evaluate(&case).unwrap(), vec![expected]);
        }
    }

    #[test]
    fn namespaced_dispatch_preserves_independent_values() {
        // Exact polynomial, linear algebra and sampled integral identities.
        for (api, args, expected) in [
            (
                "linalg.solve",
                json!({"a": [[2., 0.], [0., 4.]], "b": [6., 8.]}),
                vec![3., 2.],
            ),
            (
                "signal.convolve",
                json!({"x": [1., 2.], "h": [1., -1.]}),
                vec![1., 1., -2.],
            ),
            (
                "integrate.trapezoid",
                json!({"y": [0., 1., 4.], "coordinates": [0., 1., 2.]}),
                vec![3.],
            ),
            ("stats.variance", json!({"values": [1., 3., 5.]}), vec![4.]),
            (
                "special.assoc_laguerre",
                json!({"p0": 0.5, "p1": 1, "p2": 0.0}),
                vec![0.5],
            ),
        ] {
            let actual = evaluate(&json!({"api": api, "args": args})).unwrap();
            assert_eq!(actual, expected, "{api}");
        }
    }

    #[test]
    fn complex_special_namespace_preserves_real_gamma_and_complex_zeta() {
        let real = json!({"api": "special.gamma", "args": {"x": 5.0}});
        assert_eq!(evaluate(&real).unwrap(), vec![24.0]);

        let complex_gamma = json!({
            "api": "special.complex.gamma",
            "args": {"special_complex03_case": true, "z": [5.0, 0.0]}
        });
        let result = evaluate(&complex_gamma).unwrap();
        assert!((result[0] - 24.0).abs() < 1e-11);
        assert!(result[1].abs() < 1e-11);

        let complex_zeta = json!({
            "api": "special.complex.zeta",
            "args": {"special_complex03_case": true, "z": [2.0, 0.0]}
        });
        let result = evaluate(&complex_zeta).unwrap();
        assert!((result[0] - std::f64::consts::PI.powi(2) / 6.0).abs() < 1e-11);
        assert!(result[1].abs() < 1e-11);

        // A marker cannot accidentally redirect an unrelated real-valued API.
        let misrouted = json!({
            "api": "special.gamma",
            "args": {"special_complex03_case": true, "z": [5.0, 0.0]}
        });
        assert!(evaluate(&misrouted).is_err());
    }

    #[test]
    fn new_special_and_structural_markers_reach_distinct_adapters() {
        let special = json!({
            "api": "special.gamma.gamma_ratio",
            "args": {"special03_case": true, "a": 6.0, "b": 5.0}
        });
        assert_eq!(evaluate(&special).unwrap(), vec![5.0]);

        let structural = json!({
            "api": "core.csr.from_coo",
            "args": {
                "structural03_case": true,
                "operation": "success",
                "rows": 1,
                "cols": 1,
                "coordinates": [[0, 0, 2.5]]
            }
        });
        // rows, columns, indptr count and values, index count and values, data count,
        // data with an explicit signed-zero bit in the structural transport.
        assert_eq!(
            evaluate(&structural).unwrap(),
            vec![1.0, 1.0, 2.0, 0.0, 1.0, 1.0, 0.0, 1.0, 2.5, 0.0]
        );

        let rejected = json!({
            "api": "core.csr.validate",
            "args": {
                "structural03_case": true,
                "operation": "reject",
                "rows": 1,
                "cols": 1,
                "indptr": [0, 1],
                "indices": [2],
                "data": [1.0]
            }
        });
        assert_eq!(evaluate(&rejected).unwrap(), vec![1.0]);
    }
}
