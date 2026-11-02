# core: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-core)

Measured binary64 (`f64`) results using **Sci-core as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers 3053 cases and 3 public API names;
coverage is limited to the operations and inputs below.

## Agreement with Sci

Each library cell is **maximum / mean difference from Sci in ε** over cases
where both evaluations succeed. `Cases` counts input cases; `Sci failed N/T`
identifies missing baseline evaluations. `failed N/T` counts that library's own
failures; `Sci unavailable N/T` counts its successful evaluations excluded
because Sci failed. A dash means unavailable or no paired finite results.
Means include all compared scalar components, including exact zeros.
Full outputs, references, warnings and errors remain in the per-function CSVs.

| Function | Cases | SciPy |
| --- | ---: | --- |
| `csr.from_coo` | 1014 | 0.0 / 0.0 |
| `csr.to_dense` | 1020 | 0.0 / 0.0 |
| `csr.validate` | 1019 | 0.0 / 0.0 |


## Conventions and limitations

CSR constructors, validators and conversion routines are exercised against independent SciPy sparse and explicit storage invariants. Canonical pointer arrays, index ordering, zeros and rejection behavior are retained in the per-case measurements.

Inputs are generated deterministically into `sci-testing/cases/generated/<module>/<function>.csv` (optionally `.gz`); raw outputs, high-precision
references and the environment are stored in that repository's `results/`.
NaN, infinity, signed-zero contracts, invalid domains and work-limit enforcement
are outside this finite-input accuracy sample. Existing Sci API tests cover
additional behavior. The [coverage audit](../../results/coverage.csv) checks every public function
against the compiler-resolved inventory. These samples do not establish
worst-case bounds or full domain coverage. SciPy aliases, compatibility namespaces and Cython
variants are not counted as independent implementations.

## Reproduction

- Sci source: `c781e8259b175ebeb3c7d27a522ef3f3a7e24998` (numerical sources clean: true).
- Testing source SHA-256: `a9a4095493112b6b04c27596db78ecaef74ca844a93c23be4da6563ce7cbfcbe`.
- NumPy 2.5.4; SciPy 1.18.1; FilterPy 1.4.5; mpmath 1.3.0.
- Native comparison versions: {'Boost.Math': '1_92', 'GSL': '2.8'}.
- rustc 1.99.1 (58cf8f9c3 2026-10-13); Python 3.13.3; dev; opt-level=0; debug=line-tables-only; no custom RUSTFLAGS.
- Linux-7.2.11-arch1-1-x86_64-with-glibc2.44; AMD Ryzen 5 4600H with Radeon Graphics.

See the [shared metric, provenance and reproduction instructions](README.md).

Sci-only refresh: independent references and other library outputs are retained from the saved snapshot. Their original environment is recorded in `results/environment.json` under `reference_environment`.

Sci results were refreshed by namespace. Per-function Sci revisions are recorded in `results/environment.json` under `sci_measurement_commits`.
