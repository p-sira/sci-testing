# sci-testing

Independent accuracy comparisons for the public Sci API. Keep this Git repository
next to a checkout named `sci`; the Rust path dependency is `../sci`. Generated
Markdown is published under `sci/docs/accuracy`, while all comparison code, inputs,
lockfiles and raw measurements remain here. No remote is configured yet.

```sh
uv sync --frozen
uv run --frozen python -m unittest -v
uv run --frozen python report.py --output ../sci/docs/accuracy
cargo fmt --all -- --check
cargo clippy --locked --all-targets -- -D warnings
```

Requires Python 3.13, Rust, a C++17 compiler and libclang. Default Sci features
exercise all 14 modules. For shared build artifacts, set `CARGO_TARGET_DIR` to
an existing Cargo target directory. The runner is built with Cargo's dev profile;
custom compiler/profile/target environment overrides are refused. Cargo.lock
pins the comparison dependency graph separately from the Sci workspace lockfile.
`uv.lock` pins Python dependencies and their transitive dependencies.

`cases.json` contains only fixed inputs, without expected values. `src/main.rs`
reads one case per line from stdin and writes one result per line, recording Sci
errors and nonfinite results explicitly. `report.py` calls public comparison
APIs, computes references independently at 100 and 160 decimal digits using
mpmath and analytic definitions, checks precision stability, validates output
lengths and codebook indices, and generates tables of symmetric relative error.

The finite nonzero metric is `abs(actual-reference)/min(abs(actual),abs(reference))`,
scaled by `2^52`. Effective zeros (magnitudes below `2^-1022`) agree with each
other; a single effective zero contributes `2^52` epsilon. Opposite signs use the
formula. Maximum absolute error retains subnormal magnitudes. Nonfinite samples
are failures. The report states coverage, normalization and library adaptations;
the 135 cases cover 28 representative APIs, not every Sci export or domain.

Raw results are in `results/measurements.json`; references are unrounded decimal
strings and library outputs round-trip as binary64 JSON numbers. Provenance in
`results/environment.json` records the Sci source commit, source cleanliness,
compiler/platform versions, NumPy/SciPy build configuration and a source digest.
The digest includes this README, both lockfiles, manifests, cases, runner,
generator and tests; it excludes generated files. Numerical Sci sources must be
clean; documentation changes are permitted. Resolve a digest by checking the
recorded file list in `report.py`. A second identical run reproduces all pages
and raw measurements without timestamps.

Existing Sci fixtures and API tests remain in the Sci repository. This repository
owns the additional cross-library accuracy harness only.
