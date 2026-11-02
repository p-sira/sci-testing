//! Typed adapters for nine public Sci 0.3 special-function exports.
//!
//! Comparisons use their own independent Python mathematical references.
//! Complex inputs are JSON `[re, im]`; complex outputs are `[re, im]`.

use super::Result;
use crate::numeric::{Argument, Numbers};
use sci::{complex::Complex, special};
use serde_json::Value;

fn number(args: &Value, name: &str) -> Result<f64> {
    f64::read(&args[name])
}

fn index(args: &Value, name: &str) -> Result<u32> {
    u32::read(&args[name])
}

fn complex(args: &Value, name: &str) -> Result<Complex> {
    Complex::read(&args[name])
}

/// Return each public Sci function's value in the harness' numeric ordering.
pub fn evaluate(api: &str, a: &Value) -> Result<Vec<f64>> {
    if a["special03_case"].as_bool() != Some(true) {
        return Err("missing special03_case marker".into());
    }
    match api {
        "special.bessel.iv_ratio" => {
            special::bessel::iv_ratio(number(a, "v")?, number(a, "x")?).numbers()
        }
        "special.bessel.iv_ratio_c" => {
            special::bessel::iv_ratio_c(number(a, "v")?, number(a, "x")?).numbers()
        }
        "special.bessel.iv_ratioinv" => {
            special::bessel::iv_ratioinv(number(a, "v")?, number(a, "r")?).numbers()
        }
        "special.combinatorics.normalized_gen_harmonic" => {
            special::combinatorics::normalized_gen_harmonic(
                index(a, "j")?,
                index(a, "k")?,
                index(a, "n")?,
                number(a, "a")?,
            )
            .numbers()
        }
        "special.complex.hurwitz_zeta" => {
            special::complex::hurwitz_zeta(complex(a, "s")?, complex(a, "q")?).numbers()
        }
        "special.elementary.log1pmx" => special::elementary::log1pmx(number(a, "x")?).numbers(),
        "special.gamma.digammainv" => special::gamma::digammainv(number(a, "y")?).numbers(),
        "special.gamma.gamma_ratio" => {
            special::gamma::gamma_ratio(number(a, "a")?, number(a, "b")?).numbers()
        }
        "special.integrals.scaled_exp1" => {
            special::integrals::scaled_exp1(number(a, "x")?).numbers()
        }
        _ => Err(format!("unregistered Sci 0.3 special API: {api}").into()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn typed_adapter_checks_functions_and_shapes() {
        assert_eq!(
            evaluate(
                "special.bessel.iv_ratio",
                &json!({"special03_case":true, "v":0.5, "x":0.0})
            )
            .unwrap(),
            [0.]
        );
        assert_eq!(
            evaluate(
                "special.bessel.iv_ratio_c",
                &json!({"special03_case":true, "v":0.5, "x":0.0})
            )
            .unwrap(),
            [1.]
        );
        assert_eq!(
            evaluate(
                "special.bessel.iv_ratioinv",
                &json!({"special03_case":true, "v":1.0, "r":0.0})
            )
            .unwrap(),
            [0.]
        );
        assert_eq!(
            evaluate(
                "special.combinatorics.normalized_gen_harmonic",
                &json!({"special03_case":true, "j":2, "k":3, "n":4, "a":0.0})
            )
            .unwrap(),
            [0.5]
        );
        let z = evaluate(
            "special.complex.hurwitz_zeta",
            &json!({"special03_case":true, "s":[2.0,0.0], "q":[1.0,0.0]}),
        )
        .unwrap();
        assert_eq!(z.len(), 2);
        assert!((z[0] - std::f64::consts::PI.powi(2) / 6.).abs() < 1e-12);
        assert!(z[1].abs() < 1e-12);
        assert_eq!(
            evaluate(
                "special.elementary.log1pmx",
                &json!({"special03_case":true, "x":0.0})
            )
            .unwrap(),
            [0.]
        );
        assert!(
            (evaluate(
                "special.gamma.digammainv",
                &json!({"special03_case":true, "y":-std::f64::consts::EULER_GAMMA})
            )
            .unwrap()[0]
                - 1.)
                .abs()
                < 1e-13
        );
        let ratio = evaluate(
            "special.gamma.gamma_ratio",
            &json!({"special03_case":true, "a":201.0, "b":200.0}),
        )
        .unwrap()[0];
        assert!((ratio - 200.).abs() < 1e-10);
        assert_eq!(
            evaluate(
                "special.integrals.scaled_exp1",
                &json!({"special03_case":true, "x":0.0})
            )
            .unwrap(),
            [0.]
        );
    }

    #[test]
    fn typed_adapter_rejects_bad_domains_and_shapes() {
        assert!(evaluate("special.bessel.iv_ratio", &json!({"v":1.0,"x":1.0})).is_err());
        assert!(
            evaluate(
                "special.bessel.iv_ratio",
                &json!({"special03_case":true,"v":0.25,"x":1.0})
            )
            .unwrap()[0]
                .is_nan()
        );
        assert!(
            evaluate(
                "special.bessel.iv_ratioinv",
                &json!({"special03_case":true,"v":1.0,"r":1.0})
            )
            .unwrap()[0]
                .is_infinite()
        );
        assert!(
            evaluate(
                "special.elementary.log1pmx",
                &json!({"special03_case":true,"x":-2.0})
            )
            .unwrap()[0]
                .is_nan()
        );
        assert!(
            evaluate(
                "special.gamma.gamma_ratio",
                &json!({"special03_case":true,"a":-1.0,"b":2.0})
            )
            .unwrap()[0]
                .is_nan()
        );
        assert!(
            evaluate(
                "special.integrals.scaled_exp1",
                &json!({"special03_case":true,"x":-1.0})
            )
            .unwrap()[0]
                .is_nan()
        );
        assert!(
            evaluate(
                "special.combinatorics.normalized_gen_harmonic",
                &json!({"special03_case":true,"j":0,"k":1,"n":5,"a":1.0})
            )
            .is_err()
        );
        assert!(
            evaluate(
                "special.combinatorics.normalized_gen_harmonic",
                &json!({"special03_case":true,"j":1.5,"k":2,"n":5,"a":1.0})
            )
            .is_err()
        );
        assert!(
            evaluate(
                "special.complex.hurwitz_zeta",
                &json!({"special03_case":true,"s":[1.0,0.0],"q":[1.0,0.0]})
            )
            .is_err()
        );
        assert!(
            evaluate(
                "special.complex.hurwitz_zeta",
                &json!({"special03_case":true,"s":[2.0],"q":[1.0,0.0]})
            )
            .is_err()
        );
    }
}
