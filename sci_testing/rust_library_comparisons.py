"""Independent public RustFFT, nalgebra, and statrs numerical competitors.

The lazily built Rust cdylib lives in ``rust_oracles/``, outside Sci's Cargo
workspace, and links no Sci crate. Providers identify their actual public
library; output shapes match the saved Sci-testing cases. Unsupported cases
return {} without creating a misleading comparison value.
"""

import ctypes as ct
from functools import lru_cache
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

from . import ROOT


VERSIONS = {"RustFFT": "6.4.0", "nalgebra": "0.35.0", "statrs": "0.19.1"}
_RAW_FFT = {
    "rfft": 16, "irfft": 17, "fftn": 18, "ifftn": 19,
    "fft2": 20, "ifft2": 21, "rfftn": 22, "irfftn": 23,
    "rfft2": 24, "irfft2": 25,
}
_REAL_FORWARD = {16, 22, 24}
_REAL_INVERSE = {17, 23, 25}
_DISTRIBUTIONS = {
    "Normal": 1, "Beta": 2, "Gamma": 3, "ChiSquared": 4,
    "Chi": 5, "StudentT": 6, "F": 7, "Exponential": 8,
    "Cauchy": 9, "Laplace": 10, "LogNormal": 11, "Uniform": 12,
    "Weibull": 13, "Pareto": 14, "Gumbel": 15,
    "InverseGamma": 16, "Erlang": 17, "Bernoulli": 18,
    "Binomial": 19, "Poisson": 20, "Geometric": 21,
    "NegativeBinomial": 22,
}
_DISTRIBUTION_OPERATIONS = {
    "pdf": 0, "logpdf": 1, "cdf": 2, "sf": 3, "ppf": 4,
}
_NAMES = {
    "fft.fft": ("RustFFT", 1),
    "fft.ifft": ("RustFFT", 2),
    "stats.normal_cdf": ("statrs", 3),
    "stats.normal_ppf": ("statrs", 4),
    "stats.normal_pdf": ("statrs", 15),
    "linalg.det": ("nalgebra", 5),
    "linalg.solve": ("nalgebra", 6),
    "linalg.inv": ("nalgebra", 7),
    "linalg.dot": ("nalgebra", 8),
    "linalg.matvec": ("nalgebra", 9),
    "linalg.matmul": ("nalgebra", 10),
    "linalg.cholesky": ("nalgebra", 11),
    "linalg.cho_solve": ("nalgebra", 12),
    "linalg.eigvalsh": ("nalgebra", 13),
    "linalg.svdvals": ("nalgebra", 14),
}
_NAMES.update({f"fft.{name}": ("RustFFT", opcode)
               for name, opcode in _RAW_FFT.items()})
_NAMES.update({f"stats.{name}": ("statrs", 30)
               for name in _DISTRIBUTION_OPERATIONS})


@lru_cache(maxsize=1)
def _library():
    """Compile once on demand; Cargo locks the target across processes."""
    if shutil.which("cargo") is None:
        return None
    manifest = ROOT / "rust_oracles" / "Cargo.toml"
    target = ROOT / "target" / "rust-oracles"
    env = os.environ.copy()
    env["CARGO_TARGET_DIR"] = str(target)
    command = ("cargo", "build", "--release", "--locked", "--manifest-path", str(manifest))
    result = subprocess.run(command, cwd=ROOT, env=env, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f"Rust oracle build failed: {result.stderr.strip()}")
    if sys.platform == "darwin":
        binary = target / "release" / "libsci_testing_rust_oracles.dylib"
    elif sys.platform == "win32":
        binary = target / "release" / "sci_testing_rust_oracles.dll"
    else:
        binary = target / "release" / "libsci_testing_rust_oracles.so"
    library = ct.CDLL(str(binary))
    function = library.rust_oracle_evaluate
    function.argtypes = (
        ct.c_uint32, ct.POINTER(ct.c_double), ct.c_size_t,
        ct.c_size_t, ct.c_size_t,
        ct.POINTER(ct.c_double), ct.c_size_t,
        ct.POINTER(ct.c_double), ct.c_size_t,
    )
    function.restype = ct.c_int
    return library


def _vector(value):
    if not isinstance(value, (tuple, list)) or not value or len(value) > 131_072:
        raise ValueError("expected a nonempty finite vector")
    output = [float(item) for item in value]
    if not all(math.isfinite(item) for item in output):
        raise ValueError("nonfinite vector input")
    return output


def _matrix(value):
    if not isinstance(value, (tuple, list)) or not 1 <= len(value) <= 64:
        raise ValueError("expected a nonempty matrix")
    rows = len(value)
    if not isinstance(value[0], (tuple, list)) or not 1 <= len(value[0]) <= 64:
        raise ValueError("expected nonempty rectangular matrix rows")
    cols = len(value[0])
    if any(not isinstance(row, (list, tuple)) or len(row) != cols for row in value):
        raise ValueError("inconsistent matrix rows")
    return [item for row in value for item in _vector(row)], rows, cols


def _raw_fft(api, op, args):
    """Use only the saved raw FFT recipe, with its exact C-order shape."""
    raw = args.get("raw")
    if not isinstance(raw, dict) or raw.get("api") != api.rsplit(".", 1)[-1]:
        return None
    real = _vector(raw["re"])
    if op in _REAL_FORWARD:
        data = real
    else:
        imaginary = raw.get("im")
        imaginary = [0.0] * len(real) if imaginary is None else _vector(imaginary)
        if len(real) != len(imaginary):
            return None
        data = [value for pair in zip(real, imaginary) for value in pair]

    if op in (16, 17):
        if op == 16:
            shape = [len(real)]
        else:
            specified = raw.get("n")
            shape = [2 * (len(real) - 1) if specified is None else specified]
        axes = [0]
    else:
        shape = raw["shape"]
        axes = raw.get("axes", list(range(len(shape)))) if op in (18, 19, 22, 23) else [0, 1]
    if (not isinstance(shape, (tuple, list)) or not 1 <= len(shape) <= 4
            or not all(type(v) is int and 0 < v <= 65536 for v in shape)
            or not isinstance(axes, (tuple, list))
            or not all(type(v) is int and 0 <= v < len(shape) for v in axes)
            or len(set(axes)) != len(axes)):
        return None
    if (op in _REAL_FORWARD | _REAL_INVERSE) and not axes:
        return None
    size = math.prod(shape)
    if size > 65536:
        return None
    if op in _REAL_INVERSE:
        reduced_shape = list(shape)
        reduced_shape[axes[-1]] = reduced_shape[axes[-1]] // 2 + 1
        if len(real) != math.prod(reduced_shape):
            # Padding/truncation for 1-D irfft is part of the public contract.
            if op != 17:
                return None
            required = reduced_shape[0]
            pairs = [data[i:i + 2] for i in range(0, len(data), 2)]
            pairs = (pairs + [[0.0, 0.0]] * required)[:required]
            data = [v for pair in pairs for v in pair]
    elif len(real) != size:
        return None
    count = (math.prod(shape if op in _REAL_INVERSE else
                       [shape[i] // 2 + 1 if i == axes[-1] else dim
                        for i, dim in enumerate(shape)]
                       if op in _REAL_FORWARD else shape)) * 2
    return ("RustFFT", op, data, len(shape), len(axes),
            [*shape, *axes], count)


def _distribution(api, args):
    operation = api.rsplit(".", 1)[-1]
    name = args.get("d")
    params = args.get("p")
    samples = args.get("x")
    if (name not in _DISTRIBUTIONS or not isinstance(params, (tuple, list))
            or not isinstance(samples, (tuple, list)) or not samples):
        return None
    try:
        # Boundary points may have mathematically nonfinite outputs. They
        # remain unavailable rather than being smuggled in as finite values.
        observations = _vector(samples)
        parameters = [float(p) for p in params]
        if not all(math.isfinite(p) for p in parameters):
            return None
        return ("statrs", 30, observations, _DISTRIBUTIONS[name],
                _DISTRIBUTION_OPERATIONS[operation], parameters, len(observations))
    except (ValueError, TypeError, OverflowError):
        return None


def _encode(case):
    """Return (provider, opcode, input, rows, cols, extra, output size)."""
    api = case.get("api")
    if api not in _NAMES:
        return None
    args = case.get("args")
    if not isinstance(args, dict):
        return None
    provider, opcode = _NAMES[api]
    try:
        if opcode in _RAW_FFT.values():
            return _raw_fft(api, opcode, args)
        if opcode == 30:
            return _distribution(api, args)
        if opcode in (1, 2):
            z = args["z"]
            if not isinstance(z, (list, tuple)) or not 1 <= len(z) <= 65_536:
                return None
            data = []
            for pair in z:
                if not isinstance(pair, (list, tuple)) or len(pair) != 2:
                    return None
                data.extend(_vector(pair))
            return provider, opcode, data, 0, 0, [], len(data)
        if opcode in (3, 4, 15):
            value = float(args["x"])
            if not math.isfinite(value) or (opcode == 4 and not 0.0 < value < 1.0):
                return None
            return provider, opcode, [value], 0, 0, [], 1
        if opcode == 8:
            u, v = _vector(args["u"]), _vector(args["v"])
            if len(u) != len(v):
                return None
            return provider, opcode, u, 0, 0, v, 1
        if opcode in (7, 11, 12, 13, 14):
            raw = args["raw"]
            saved = raw["args"]
            if not isinstance(saved, list) or not saved:
                return None
            data, rows, cols = _matrix(saved[0])
        else:
            data, rows, cols = _matrix(args["a"])
        if opcode in (5, 6, 7, 11, 12, 13) and rows != cols:
            return None
        extra = []
        if opcode == 6:
            extra = _vector(args["b"])
            if len(extra) != rows:
                return None
            count = rows
        elif opcode == 9:
            extra = _vector(args["u"])
            if len(extra) != cols:
                return None
            count = rows
        elif opcode == 10:
            extra, b_rows, b_cols = _matrix(args["b"])
            if b_rows != cols:
                return None
            count = rows * b_cols
        elif opcode in (11, 13):
            if not isinstance(saved[-1], bool):
                return None
            extra = [float(saved[-1])]
            count = rows * cols if opcode == 11 else rows
        elif opcode == 12:
            if len(saved) != 3 or not isinstance(saved[2], bool):
                return None
            extra = [*_vector(saved[1]), float(saved[2])]
            if len(extra) != rows + 1:
                return None
            count = rows
        elif opcode == 14:
            count = min(rows, cols)
        else:
            count = 1 if opcode == 5 else rows * cols
        return provider, opcode, data, rows, cols, extra, count
    except (ValueError, KeyError, TypeError, IndexError, OverflowError):
        return None


def _doubles(values):
    return (ct.c_double * len(values))(*values) if values else None


def evaluate(case):
    """Evaluate a single supported saved case through its Rust library."""
    encoded = _encode(case)
    if encoded is None:
        return {}
    library = _library()
    if library is None:
        return {}
    provider, opcode, data, rows, cols, extra, count = encoded
    results = (ct.c_double * count)()
    status = library.rust_oracle_evaluate(
        opcode, _doubles(data), len(data), rows, cols,
        _doubles(extra), len(extra), results, count,
    )
    if status == 1:  # Unsupported mathematical domain or singular matrix.
        return {}
    if status != 0:
        return {provider: {"error": f"Rust oracle evaluation failed (status {status})"}}
    values = list(results)
    if not all(math.isfinite(item) for item in values):
        return {}
    return {provider: {"values": values}}
