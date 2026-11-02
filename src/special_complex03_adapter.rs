//! Direct public Sci 0.3 complex-special API comparisons.
//!
//! Each input complex number is `[real, imaginary]`. Scalar complex results
//! are flattened to two f64 values, Airy arrays to eight, and Fresnel's
//! `(S, C)` pair to four. Distinct namespace and marker avoid the real facade.

use super::Result;
use crate::numeric::{Argument, Numbers};
use sci::{complex::Complex, special::complex as c};
use serde_json::Value;

fn z(a: &Value, key: &str) -> Result<Complex> {
    Complex::read(&a[key])
}

fn v(a: &Value) -> Result<f64> {
    f64::read(&a["v"])
}

/// Evaluate exactly the named public complex overload, preserving result order.
pub fn evaluate(api: &str, a: &Value) -> Result<Vec<f64>> {
    if a["special_complex03_case"].as_bool() != Some(true) {
        return Err("missing special_complex03_case marker".into());
    }
    match api {
        "special.complex.airy" => c::airy(z(a, "z")?).numbers(),
        "special.complex.airye" => c::airye(z(a, "z")?).numbers(),
        "special.complex.dawsn" => c::dawsn(z(a, "z")?).numbers(),
        "special.complex.digamma" => c::digamma(z(a, "z")?).numbers(),
        "special.complex.erf" => c::erf(z(a, "z")?).numbers(),
        "special.complex.erfc" => c::erfc(z(a, "z")?).numbers(),
        "special.complex.erfcx" => c::erfcx(z(a, "z")?).numbers(),
        "special.complex.erfi" => c::erfi(z(a, "z")?).numbers(),
        "special.complex.exp1" => c::exp1(z(a, "z")?).numbers(),
        "special.complex.expi" => c::expi(z(a, "z")?).numbers(),
        "special.complex.fresnel" => c::fresnel(z(a, "z")?).numbers(),
        "special.complex.gamma" => c::gamma(z(a, "z")?).numbers(),
        "special.complex.hankel1" => c::hankel1(v(a)?, z(a, "z")?).numbers(),
        "special.complex.hankel1e" => c::hankel1e(v(a)?, z(a, "z")?).numbers(),
        "special.complex.hankel2" => c::hankel2(v(a)?, z(a, "z")?).numbers(),
        "special.complex.hankel2e" => c::hankel2e(v(a)?, z(a, "z")?).numbers(),
        "special.complex.iv" => c::iv(v(a)?, z(a, "z")?).numbers(),
        "special.complex.ive" => c::ive(v(a)?, z(a, "z")?).numbers(),
        "special.complex.jv" => c::jv(v(a)?, z(a, "z")?).numbers(),
        "special.complex.jve" => c::jve(v(a)?, z(a, "z")?).numbers(),
        "special.complex.kv" => c::kv(v(a)?, z(a, "z")?).numbers(),
        "special.complex.kve" => c::kve(v(a)?, z(a, "z")?).numbers(),
        "special.complex.lambertw" => c::lambertw(z(a, "z")?, i32::read(&a["k"])?).numbers(),
        "special.complex.loggamma" => c::loggamma(z(a, "z")?).numbers(),
        "special.complex.rgamma" => c::rgamma(z(a, "z")?).numbers(),
        "special.complex.spence" => c::spence(z(a, "z")?).numbers(),
        "special.complex.wofz" => c::wofz(z(a, "z")?).numbers(),
        "special.complex.xlog1py" => c::xlog1py(z(a, "x")?, z(a, "y")?).numbers(),
        "special.complex.xlogy" => c::xlogy(z(a, "x")?, z(a, "y")?).numbers(),
        "special.complex.yv" => c::yv(v(a)?, z(a, "z")?).numbers(),
        "special.complex.yve" => c::yve(v(a)?, z(a, "z")?).numbers(),
        "special.complex.zeta" => c::zeta(z(a, "z")?).numbers(),
        _ => Err(format!("unknown public complex-special API: {api}").into()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn complex_shape_order_and_namespaced_gamma() {
        let args = json!({"special_complex03_case": true, "z": [1.0, 0.0]});
        let airy = evaluate("special.complex.airy", &args).unwrap();
        assert_eq!(airy.len(), 8);
        assert!(airy[0] > 0.0 && airy[2] < 0.0 && airy[4] > 0.0 && airy[6] > 0.0);
        let fresnel = evaluate("special.complex.fresnel", &args).unwrap();
        assert_eq!(fresnel.len(), 4);
        assert!(fresnel[0] > 0.0 && fresnel[2] > 0.0);
        let gamma = evaluate("special.complex.gamma", &args).unwrap();
        assert_eq!(gamma, [1.0, 0.0]);
    }

    #[test]
    fn typed_parameters_and_marker_are_checked() {
        assert!(evaluate("special.complex.gamma", &json!({"z":[1.0,0.0]})).is_err());
        assert!(
            evaluate(
                "special.complex.erf",
                &json!({"special_complex03_case":true,"z":[1.0]})
            )
            .is_err()
        );
        assert!(
            evaluate(
                "special.complex.lambertw",
                &json!({"special_complex03_case":true,"z":[1.0,0.0],"k":0.5})
            )
            .is_err()
        );
        assert!(
            evaluate(
                "special.complex.iv",
                &json!({"special_complex03_case":true,"z":[1.0,0.0]})
            )
            .is_err()
        );
        assert!(
            evaluate(
                "special.gamma",
                &json!({"special_complex03_case":true,"z":[1.0,0.0]})
            )
            .is_err()
        );
    }
}
