//! Explicit recipes for callback APIs and exact metadata contracts.
use super::{Result, matrix, scalar, vector};
use crate::numeric::Numbers;
use sci::{constants, differentiate, integrate, optimize, stats};
use serde_json::Value;

fn text(value: &str) -> Vec<f64> {
    std::iter::once(value.len() as f64)
        .chain(value.bytes().map(f64::from))
        .collect()
}

pub fn evaluate(api: &str, a: &Value) -> Result<Vec<f64>> {
    let v = |key: &str| vector(&a[key]);
    let x = || scalar(&a["x"]);
    Ok(match api {
        "linalg.dot" => vec![sci::linalg::basic::dot(&v("u")?, &v("v")?)?],
        "linalg.matvec" => sci::linalg::basic::matvec(&matrix(&a["a"])?, &v("u")?)?,
        "linalg.matmul" => {
            sci::linalg::basic::matmul(&matrix(&a["a"])?, &matrix(&a["b"])?)?.numbers()?
        }
        "linalg.transpose" => sci::linalg::basic::transpose(&matrix(&a["a"])?)?.numbers()?,
        "stats.normal_pdf" => vec![stats::distributions::normal_pdf(x()?, 0., 1.)?],
        "stats.normal_ppf" => vec![stats::distributions::normal_ppf(x()?, 0., 1.)?],
        "constants.value" => vec![constants::value(a["key"].as_str().ok_or("missing key")?)?],
        "constants.precision" => vec![constants::precision(
            a["key"].as_str().ok_or("missing key")?,
        )?],
        "constants.unit" => text(constants::unit(a["key"].as_str().ok_or("missing key")?)?),
        "constants.find" => constants::find(a["key"].as_str())
            .into_iter()
            .flat_map(text)
            .collect(),
        "constants.convert_temperature" => vec![constants::convert_temperature(
            x()?,
            constants::TemperatureScale::Celsius,
            constants::TemperatureScale::Fahrenheit,
        )],
        "differentiate.jacobian" => {
            differentiate::jacobian(|p| vec![p[0] * p[1], p[0] + p[1]], &v("point")?, 1e-4)?
                .numbers()?
        }
        "differentiate.hessian" => differentiate::hessian(
            |p| p[0] * p[0] + 3. * p[0] * p[1] + 2. * p[1] * p[1],
            &v("point")?,
            1e-4,
        )?
        .numbers()?,
        "integrate.simpson" => vec![integrate::sampled::simpson(&v("values")?, None, 0.25)?],
        "integrate.cumulative_trapezoid" => {
            integrate::sampled::cumulative_trapezoid(&v("values")?, None, 0.25, Some(0.))?
        }
        "integrate.cumulative_simpson" => {
            integrate::sampled::cumulative_simpson(&v("values")?, None, 0.25, Some(0.))?
        }
        "integrate.romb" => vec![integrate::sampled::romb(&v("values")?, 0.25)?],
        "integrate.newton_cotes" => {
            integrate::quadrature::newton_cotes(&v("values")?)?.numbers()?
        }
        "interpolate.pchip_interpolate" => sci::interpolate::polynomial::pchip_interpolate(
            &v("nodes")?,
            &v("values")?,
            &v("query")?,
        )?,
        "interpolate.lagrange" => {
            sci::interpolate::polynomial::lagrange(&v("nodes")?, &v("values")?)?
        }
        "interpolate.pade" => {
            sci::interpolate::polynomial::pade(&v("values")?, 2, Some(2))?.numbers()?
        }
        "cluster.vq.kmeans" => {
            sci::cluster::vq::kmeans(&matrix(&a["obs"])?, &matrix(&a["initial"])?, 1e-8, 32)?
                .numbers()?
        }
        "cluster.vq.kmeans2" => {
            sci::cluster::vq::kmeans2(&matrix(&a["obs"])?, &matrix(&a["initial"])?, 8)?.numbers()?
        }
        name if name.starts_with("optimize.") => {
            let start = v("start")?;
            let objective = |p: &[f64]| 3. + (p[0] - 1.).powi(2) + (p[1] + 2.).powi(2);
            let gradient = |p: &[f64]| vec![2. * (p[0] - 1.), 2. * (p[1] + 2.)];
            let options = optimize::minimize::Options {
                max_iterations: 512,
                tolerance: 1e-8,
                learning_rate: 0.25,
            };
            let result = match name {
                "optimize.minimize" => optimize::minimize::minimize(objective, &start, options),
                "optimize.try_minimize" => {
                    optimize::minimize::try_minimize(|p| Ok(objective(p)), &start, options)
                }
                "optimize.minimize_with_gradient" => {
                    optimize::minimize::minimize_with_gradient(objective, gradient, &start, options)
                }
                "optimize.try_minimize_with_gradient" => {
                    optimize::minimize::try_minimize_with_gradient(
                        |p| Ok(objective(p)),
                        |p| Ok(gradient(p)),
                        &start,
                        options,
                    )
                }
                _ => return Err(format!("unknown callback API: {api}").into()),
            }?;
            if result.status != optimize::minimize::Status::Converged {
                return Err("optimizer did not converge within its work budget".into());
            }
            result.point.into_iter().chain([result.cost]).collect()
        }
        _ => return Err(format!("unregistered explicit recipe: {api}").into()),
    })
}
