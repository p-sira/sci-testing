# sci-testing

Library comparisons using Sci as the baseline for every public function in the sibling
[Sci](../sci) checkout. Reports, inputs, recipes and measurements live here;
start with the [comparison tables](docs/accuracy/README.md).

```sh
uv sync --frozen
uv run --frozen python -m sci_testing.report --full
uv run --frozen python -m sci_testing.report --check
```

After the initial full run, refresh only Sci outputs with:

```sh
uv run --frozen python -m sci_testing.report
uv run --frozen python -m sci_testing.report --module special
uv run --frozen python -m sci_testing.report --module sparse.linalg
uv run --frozen python -m sci_testing.report --module linalg --module interpolate
```

The default reuses saved independent references and all other library results,
including their errors and warnings. `--module` selects a module or nested
namespace and can be repeated; unselected measurements remain unchanged.
`--sci-only` explicitly selects the default mode. Use `--full` to recompute all
references and library results, including when inputs change or no snapshot exists.
`--module` applies only to Sci-only refreshes.

Before measuring an updated Sci checkout, update `sci-revision.txt` to its commit
and refresh `Cargo.lock` if its dependencies changed. Sci-only runs require the
same generated inputs and intact saved measurements. They rebuild the Rust runner,
refresh coverage and reports, and avoid evaluating comparison libraries or
independent references. Original reference provenance is retained under
`reference_environment` in `results/environment.json`; partial refreshes also
record per-function Sci revisions under `sci_measurement_commits`.

The full suite has about 2.26 million cases for 1,109 distinct public functions.
Each recipe adds 1,000 endpoint-inclusive NumPy linspace cases and 1,000 uniform
random cases from seed `20261020`, with independent streams per API. Original
boundary profiles are retained. Constrained integer domains, fixed shapes and
metadata queries can repeat bounded inputs; `results/coverage.csv` reports unique
inputs separately. The audit requires every compiler-inventoried function to
have a successful public-API comparison and both dense sample groups.
Distribution density sweeps use 1%–99% quantiles to preserve continuous support
and integer arguments for probability masses.
The runner evaluates identical deterministic inputs once and expands the measured
result to each repeated record. Error-policy APIs still execute every call.

## Build and validation

Requires Python 3.13, Rust, a C++17 compiler, libclang, Boost.Math headers and GSL
development files. The Linux workflow installs `libboost-math-dev` and `libgsl-dev`.
`uv.lock` pins NumPy 2.5.4, SciPy 1.18.1, mpmath 1.3.0 and FilterPy 1.4.5;
`Cargo.lock` independently pins the Rust dependency graph.

```sh
uv run --frozen python -m unittest -v
cargo fmt --all -- --check
cargo test --locked
cargo clippy --locked --all-targets -- -D warnings
RUSTDOCFLAGS='-D warnings' cargo doc --no-deps
```

Keep the repositories next to each other: the Rust dependency is `../sci`.
Sci must match `sci-revision.txt`; modified numerical sources or fixture providers
are refused. Documentation edits are permitted. The runner uses Cargo's dev
profile without custom compiler/profile/target overrides. Set `CARGO_TARGET_DIR`
to reuse an existing build directory.

## Files and recipes

- `sci_testing/`: Python reporting, reference evaluation, CSV handling and generators.
- `tests/`: Python regression tests, discovered by `python -m unittest`.
- `src/`: Rust public-API runner and generated fixture adapters.
- `native/`: C++ bridge for Boost.Math and GSL comparisons.
- `cases/base/<module>/<function>.csv`: original hand-picked inputs.
- `cases/recipes/<module>/<function>.csv`: imported valid fixture profiles.
- `cases/generated/<module>/<function>.csv`: complete generated input cases.
- `results/measurements/<module>/<function>.csv[.gz]`: independent references, actual
  library outputs, exceptions and warnings, keyed by case ID and library.
- `results/coverage.csv`: per-function counts, uniqueness and successful coverage.
- `results/environment.json`: versions, compiler settings and snapshot digests.
- `docs/accuracy/`: compact tables and methodology, generated here exclusively.

Nested namespaces stay explicit, for example `signal.windows/hann.csv` and
`cluster.vq/kmeans.csv`. Inputs contain no expected results. Common numerical
cells use `s:` scalars, `i:` exact integers, `v:` space-separated vectors and `m:`
semicolon-separated matrix rows. `b:`, `t:` and `n:` represent booleans, text and
None. Heterogeneous structures use compact JSON cells inside CSV. Binary64
numbers retain their round-trip decimal representations, signed zero and
subnormals; large integer parameters retain their exact integer values.
Measurement reference vectors are written once on the first library row for a
case; later rows for that case leave the reference cell empty. Warning columns
preserve multiple messages, including embedded newlines.
Measurement CSVs above 16 MiB use deterministic gzip compression. The loader
handles both forms without changing stored precision or numerical values.

`sci_testing/generate_cases.py` generates inputs only. Python family decorators
and typed Rust macro registrations reduce repetitive adapters.
`sci_testing/prepare_recipes.py` imports inputs from the pinned Sci fixtures.
The Rust code generators reuse existing public-API fixture dispatchers, retaining
their structural invariants and observable ordering. Tests check that generated
Rust is current. Run generators as modules from the repository root:

```sh
uv run --frozen python -m sci_testing.api_registry
uv run --frozen python -m sci_testing.generate_adapters
uv run --frozen python -m sci_testing.generate_linalg_adapter
uv run --frozen python -m sci_testing.generate_interpolate_adapter
uv run --frozen python -m sci_testing.generate_sparse_adapter
```

`sci_testing/inventory.py` records compiler-resolved exports, including macro
expansions and aliases. Ordinary runs read its checked-in CSV and require stable
Rust only.

## References and failures

Original module recipes use independent mpmath references at 100 and 160 decimal
digits. Added recipes use pinned public SciPy/stdlib calls, independent identities
or quadrature; binary64 references measure agreement with that reference library.
Sci and its private backends never generate expectations. Nonunique factors use
projectors, reconstructions and ordered eigenvalues. Module pages document the
remaining adaptations and reference accuracy limits.

The separate C++ bridge calls Boost.Math and GSL where matching scalar APIs exist.
It converts elliptic parameter m to modulus sqrt(m), contains native exceptions,
disables GSL's abort handler and retains GSL error statuses. Versions and native
compiler flags are recorded. See the [Boost.Math reference](https://www.boost.org/doc/libs/latest/libs/math/doc/html/index.html)
and [GSL special-function conventions](https://www.gnu.org/software/gsl/doc/html/specfunc.html).

Tables use the corresponding **Sci-* crate as the reference**, with columns
**Function**, **Cases**, then the relevant libraries (SciPy, Boost.Math, GSL and
others). Cells show maximum / mean symmetric relative difference from Sci in
epsilon over cases where both evaluations succeed. Agreement with Sci does not
establish numerical accuracy; independent accuracy references remain in the raw
CSV `reference` field. `failed N/T` counts each library’s own failures.
`Sci failed N/T` in Cases counts baseline failures; `Sci unavailable N/T` counts
successful library evaluations excluded because Sci failed. Full errors and outputs
remain in CSV. Metadata and discrete ordering are also covered by the underlying
Sci fixture assertions. These samples do not establish whole-domain accuracy.

`python -m sci_testing.report --check` verifies inputs, source and result digests, the complete
coverage audit and all generated pages without compiling or writing files.
Moving or editing harness sources invalidates the recorded source digest; rerun
the report against the pinned Sci checkout to create a new snapshot. Use `--full`
when changing reference recipes or comparison adapters; the default preserves
their saved results.
Sources, manifests, lockfiles, recipes and generated inputs are hashed; derived
results and reports are excluded from the source digest. Each run records its
own platform. Repeated runs on the same platform reproduce the snapshot without
timestamps. GitHub Actions regenerates and uploads reports and data on pushes
and pull requests when this repository is hosted on GitHub.

## Sparse comparisons

The sparse suite covers all 98 compiler-inventoried exports (93 distinct
functions), including `sparse.csgraph` and `sparse.linalg`. It imports 253 input
profiles from Sci's two sparse suites, adds two uniform-scale review cases, and generates 1,000 linspace and 1,000 seeded
random cases per export. `sci_testing/generate_sparse_adapter.py` reuses public Rust dispatch
and structural invariants; `sci_testing/sparse_cases.py` evaluates independent pinned SciPy
recipes without running their file writers. NPZ comparisons use memory archives.

Matrix outputs include dimensions. Nonunique eigenvectors, permutations, flow
solutions and random streams retain the fixture invariants and compare stable
observables. Every sparse scalar uses a `(class, value)` pair: class 0 carries a
finite value, class +1/-1 denotes positive/negative infinity, and class 2 denotes
NaN; nonfinite payloads are zero. This preserves unreachable graph vertices in
the finite report transport. Means include these exact structural tags. Numerical
input sweeps preserve dimensions, sparsity and integer flow capacities. See the
[sparse comparison table](docs/accuracy/sparse.md) and
[PR #7 review](docs/pr7-review.md).
