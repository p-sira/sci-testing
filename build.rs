//! Reuse fixture modules with one shared parser and absolute fixture paths.
use std::{env, fs, path::PathBuf};

fn main() {
    let manifest = PathBuf::from(env::var_os("CARGO_MANIFEST_DIR").unwrap());
    let tests = manifest.join("sci/tests").canonicalize().unwrap();
    let output = PathBuf::from(env::var_os("OUT_DIR").unwrap());
    let escaped = tests
        .to_string_lossy()
        .replace('\\', "\\\\")
        .replace('"', "\\\"");
    for name in ["statistics", "signal", "filter"] {
        let path = tests.join(format!("{name}.rs"));
        println!("cargo:rerun-if-changed={}", path.display());
        let source = fs::read_to_string(path)
            .unwrap()
            .replacen("mod common;", "use crate::common;", 1)
            .replace(
                "include_str!(\"fixtures/",
                &format!("include_str!(\"{escaped}/fixtures/"),
            );
        fs::write(output.join(format!("{name}.rs")), source).unwrap();
    }
}
