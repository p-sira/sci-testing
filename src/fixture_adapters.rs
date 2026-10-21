//! Reuse the pinned Sci fixture dispatchers for structured functional APIs.
//! Their input and output conventions are independently defined by the paired
//! SciPy generators. Keep all numerical calls on the public Sci facade.
use super::Result;
use serde_json::Value;

#[allow(dead_code, unused_imports)]
mod statistics {
    include!(concat!(env!("OUT_DIR"), "/statistics.rs"));
    pub fn run(fields: &[&str]) -> sci::Result<Vec<f64>> {
        evaluate(fields)
    }
}
#[allow(dead_code, unused_imports)]
mod signal {
    include!(concat!(env!("OUT_DIR"), "/signal.rs"));
    pub fn run(name: &str, args: &[&str]) -> Vec<Vec<f64>> {
        evaluate(name, args)
    }
}

fn encoded(value: &Value) -> String {
    match value {
        Value::Null => "none".into(),
        Value::String(text) => text.clone(),
        Value::Array(values) => {
            let separator = if values.first().is_some_and(Value::is_array) {
                ":"
            } else {
                ";"
            };
            values
                .iter()
                .map(encoded)
                .collect::<Vec<_>>()
                .join(separator)
        }
        _ => value.to_string(),
    }
}

pub fn evaluate(api: &str, args: &Value) -> Result<Vec<f64>> {
    if let Some(public_name) = api.strip_prefix("stats.") {
        let name = args["dispatch_api"].as_str().unwrap_or(public_name);
        let groups = args["g"]
            .as_array()
            .ok_or("expected groups")?
            .iter()
            .map(encoded)
            .collect::<Vec<_>>()
            .join("/");
        let fields = [
            "0".into(),
            name.into(),
            encoded(&args["p"]),
            encoded(&args["x"]),
            encoded(&args["y"]),
            groups,
            encoded(&args["d"]),
        ];
        let mut values = statistics::run(&fields.iter().map(String::as_str).collect::<Vec<_>>())?;
        if name == "weightedtau" {
            values.truncate(1);
        }
        return Ok(values);
    }
    if let Some(name) = api.strip_prefix("signal.") {
        let mut fields = args["arguments"]
            .as_array()
            .ok_or("expected arguments")?
            .iter()
            .map(encoded)
            .collect::<Vec<_>>();
        if matches!(
            name,
            "lp2lp_zpk"
                | "lp2hp_zpk"
                | "lp2bp_zpk"
                | "lp2bs_zpk"
                | "bilinear_zpk"
                | "zpk2tf"
                | "freqs_zpk"
                | "freqz_zpk"
                | "zpk2ss"
                | "zpk2sos"
        ) {
            fields[0] = args["arguments"][0]
                .as_array()
                .ok_or("expected ZPK")?
                .iter()
                .map(encoded)
                .collect::<Vec<_>>()
                .join("|");
        }
        return Ok(
            signal::run(name, &fields.iter().map(String::as_str).collect::<Vec<_>>())
                .into_iter()
                .flatten()
                .collect(),
        );
    }
    Err(format!("unsupported structured API: {api}").into())
}
