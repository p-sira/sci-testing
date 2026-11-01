# Sci library comparisons

Historical and newly measured comparisons use each Sci module as the reference on
the same saved inputs. The compiler inventory lists 1142 mapped canonical
names across 17 modules and 2115 exported names; these are not
all independently verified against Sci 0.3. There are 0 currently
unmeasured new API names, 0 other untested exports and
53 unresolved legacy mappings.
The [coverage audit](../../results/coverage.csv) explicitly marks pending names
alongside measured sample counts, unique inputs, successes and failures.

Tables list **Function**, **Cases**, then relevant libraries (SciPy, Boost.Math,
GSL and others). Each library cell is **maximum / mean symmetric relative
difference from Sci in machine epsilon**, over cases where both evaluations
succeed. `failed N/T` counts that library's failures; `Sci failed N/T` in Cases
counts baseline failures, and `Sci unavailable N/T` marks otherwise successful
evaluations excluded from comparison. Unavailable comparisons use a dash. The per-function
CSV preserves every reference, result, warning and exception.

| Module | Recipe names | Cases |
| --- | ---: | ---: |
| [cluster](cluster.md) | 31 | 62217 |
| [complex](complex.md) | 34 | 70919 |
| [constants](constants.md) | 7 | 14252 |
| [core](core.md) | 3 | 3053 |
| [differentiate](differentiate.md) | 4 | 8277 |
| [fft](fft.md) | 45 | 90336 |
| [filter](filter.md) | 10 | 20226 |
| [integrate](integrate.md) | 7 | 14135 |
| [interpolate](interpolate.md) | 25 | 50224 |
| [io](io.md) | 10 | 413 |
| [linalg](linalg.md) | 109 | 218680 |
| [optimize](optimize.md) | 10 | 20152 |
| [signal](signal.md) | 179 | 358507 |
| [sparse](sparse.md) | 100 | 198293 |
| [spatial](spatial.md) | 33 | 66272 |
| [special](special.md) | 405 | 811236 |
| [stats](stats.md) | 169 | 342414 |


## Inputs and references

Historical dense recipes and new special-function sweeps use 1,000
endpoint-inclusive linspace inputs plus 1,000 independently seeded random inputs
per API (seed `20261020`), with preserved boundary profiles. New structural CSR
and spectral APIs use 1,000 seeded shape-varied inputs plus edge and rejection
cases. I/O APIs use independently generated Matrix Market, WAV and MATLAB
binary fixtures covering formats, data types and rejected inputs, with 8–136
cases per API. Integer-only functions, metadata queries and constrained shapes
use bounded discrete profiles, so some inputs repeat; the coverage CSV reports
unique inputs explicitly. Identical deterministic inputs share one measured
evaluation; all repeated case records remain in the CSV and totals.
Error-policy APIs execute every call. These samples do not establish domain-wide
accuracy bounds. Ordinary Sci API tests retain invalid-input, NaN, infinity,
signed-zero and work-limit contracts.

Tables use the measured binary64 output from the corresponding **Sci-* crate**
as their reference. Agreement with Sci does not establish numerical accuracy.
Independent accuracy references remain in the raw CSV `reference` field: original
module recipes use mpmath at 100 and 160 decimal digits; additional recipes use
pinned public SciPy/stdlib APIs, independent identities and quadrature. These
references are independent of the table baseline. Public Sci APIs are
called through typed Rust registrations and the pinned Sci fixture dispatchers;
no Sci or private-backend result generates an expected value.

Nonunique decompositions use projectors, reconstructions and ordered eigenvalues.
Callback minimizers compare converged points and costs with analytic minimizers.
Text metadata uses exact UTF-8 length/byte encodings. Module pages state further
ordering, normalization, branch and reference limitations. Boost.Math and GSL
are evaluated by a separate C++ bridge where matching scalar APIs exist. Native
exceptions and GSL error statuses remain failures, including underflow statuses.
Elliptic parameter m is converted to modulus sqrt(m); negative parameters may
therefore fall outside the native comparator's supported real-modulus domain.
Boost betaln is the logarithm of its beta evaluation, with that expression's
intermediate range limitations. Versions and compiler flags are recorded in
`results/environment.json`.

## Comparison metric

For finite nonzero library output v and measured Sci output r,

\[
E(v,r)=\frac{|v-r|}{\min(|v|,|r|)},\qquad E_\varepsilon=E/\varepsilon,
\qquad\varepsilon=2^{-52}.
\]

Values below `2^-1022` are effective zeros: two agree exactly; one contributes
`2^52` epsilon. Opposite signs use the stated formula. Means aggregate scalar
components, including exact zeros; vector-quantization code indices are checked
exactly and omitted from its error aggregate. Epsilon-scaled error is not an ULP
count. Near cancellation or zeros, inspect absolute differences in the raw data.
This adapts [Boost.Math's error-report convention](https://www.boost.org/doc/libs/latest/libs/math/doc/html/math_toolkit/relative_error.html).

## Reproduce

Initialize the `sci/` submodule at the revision recorded in `sci-revision.txt`.
Install a C++17 compiler, libclang, Boost.Math headers and GSL development files.
The Cargo and uv lockfiles pin Rust/Python dependencies separately from Sci.

```sh
git submodule update --init --recursive
uv sync --frozen
uv run --frozen python -m unittest -v
cargo test --locked
uv run --frozen python -m sci_testing.report --full
uv run --frozen python -m sci_testing.report --check
```

`python -m sci_testing.report` refreshes Sci outputs while retaining saved independent
references and other library results. Select modules with `--module special` or
nested namespaces with `--module sparse.linalg`; repeat `--module` to select several.
Unselected measurements stay fixed. Original reference provenance is retained under
`reference_environment` in `results/environment.json`; partial refreshes record
per-function Sci revisions under `sci_measurement_commits`. Commit the updated
`sci/` submodule pointer and `sci-revision.txt` and, if necessary, `Cargo.lock` before
measuring a newer Sci commit.
Use `--full` for the initial snapshot, changed inputs, or refreshing all references
and comparison libraries. `--check` validates source/input/result digests,
the complete coverage audit and every generated page without compiling or writing
files. It validates the recorded platform snapshot; it does not require identical
floating-point measurements on another platform. Set `CARGO_TARGET_DIR` to reuse
an existing build directory. The runner uses Cargo's dev profile without custom
compiler flags. Numerical sources and fixture providers must be clean; Sci
documentation edits are permitted.

The GitHub workflow regenerates these reports and uploads inputs and measurements.
An identical run on the same platform reproduces the snapshot without timestamps.
All reports and data belong to this repository; Sci links to them externally.

## Snapshot

- Sci commit: `3fbdbf3a73bcad57fa91e30262a4e9f9d09fd57a`.
- Harness source SHA-256: `9b806b00dd16330a3acc24ed8c2ff9399bd069fd0a2017e2d4988c55f8cd6e81`.
- Input SHA-256: `f85348e5bcc668c154cd406478d139513eb93340546e43deef2386fe52da7854`.
- 2349606 saved cases covering 1181 recipe API names; 1142 mapped canonical public functions in the inventory.
- Native versions: {'Boost.Math': '1_92', 'GSL': '2.8'}.

Sci-only refresh: independent references and other library outputs are retained from the saved snapshot. See `reference_environment` in `results/environment.json` for their original provenance.

Per-function Sci revisions for partial refreshes are recorded under `sci_measurement_commits` in `results/environment.json`.

New API references were evaluated separately from the historical snapshot. Their provenance is recorded in `reference_extensions` in `results/environment.json`.
