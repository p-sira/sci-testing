"""Independent SciPy/NumPy oracles for six structural Sci 0.3 exports.

The CSR comparison includes canonical pointers, column order, stored values,
explicit zero signs and dense row-major values.  The STFT comparison checks
both the SciPy projection and the NOLA/dual-window equations.  Invalid storage
is deliberately tested against the stricter, documented Sci core contract:
SciPy's low-level CSR constructor does not enforce canonical storage.

No reference depends on Sci Rust.  All returned values are finite scalars;
invalid cases return the one-element vector [1.] when the call must reject.
"""

import hashlib
import math
import random
from functools import lru_cache

import numpy as np
import scipy
from scipy import signal, sparse

SCIPY_VERSION = '1.18.1'
NUMPY_VERSION = '2.5.4'


@lru_cache(maxsize=1)
def _versions():
    if (scipy.__version__, np.__version__) != (SCIPY_VERSION, NUMPY_VERSION):
        raise RuntimeError('structural03 references require scipy==1.18.1 and numpy==2.5.4')


APIS = (
    'core.csr.from_coo',
    'core.csr.to_dense',
    'core.csr.validate',
    'sparse.array.csr_indices_i32',
    'sparse.array.csr_indices_i64',
    'signal.utilities.closest_stft_dual_window',
)


def _case(api, name, operation='success', **arguments):
    return {
        'id': f'{api}/{name}', 'api': api,
        'group': 'Sci 0.3 structural and spectral; independent SciPy/NumPy reference',
        'args': dict(structural03_case=True, operation=operation, **arguments),
    }


def _rng(api, seed):
    digest = hashlib.sha256(f'structural03:{seed}:{api}'.encode()).digest()
    return random.Random(int.from_bytes(digest, 'big'))


def _coo(rows, cols, triples):
    return dict(rows=rows, cols=cols, coordinates=triples)


def _csr(rows, cols, triples):
    """Build valid CSR independently, retaining explicit zeros and duplicates.

    The finite dyadic samples make duplicate sums exactly representable at
    their chosen scales; SciPy is used again to check reference consistency.
    """
    ordered = sorted(triples, key=lambda entry: (entry[0], entry[1]))
    indices, data = [], []
    indptr = [0] * (rows + 1)
    last = None
    for r, c, v in ordered:
        if last == (r, c):
            data[-1] += v
        else:
            indices.append(c)
            data.append(v)
            indptr[r + 1] += 1
            last = r, c
    for i in range(rows):
        indptr[i + 1] += indptr[i]
    return dict(rows=rows, cols=cols, indptr=indptr, indices=indices, data=data)


def _random_coo(rng, i):
    # Empty and singleton cases are covered explicitly by _sparse_edges;
    # randomized sweeps exercise distinct nontrivial storage arrangements.
    shapes = [(2, 3), (3, 5), (5, 3), (8, 8), (13, 21), (32, 17),
              (1, 8), (8, 1), (6, 9), (19, 5), (9, 9), (64, 32)]
    rows, cols = shapes[i % len(shapes)]
    if not rows or not cols:
        return _coo(rows, cols, [])
    cap = min(3 * rows * cols, 100)
    count = rng.randrange(1, cap + 1)
    entries = []
    for _ in range(count):
        row, col = rng.randrange(rows), rng.randrange(cols)
        # Dyadic values, mixed signs/scales, and retained explicit zeros.
        value = rng.randrange(-128, 129) * 2.0 ** rng.choice([-12, -4, 0, 4])
        if rng.randrange(11) == 0:
            value = -0.0 if rng.randrange(2) else 0.0
        entries.append([row, col, value])
    return _coo(rows, cols, entries)


def _invalid_csr():
    base = _csr(3, 4, [[0, 1, 2.], [1, 3, -1.], [2, 0, -0.]])
    mutations = [
        ('short-indptr', dict(indptr=[0, 1, 3])),
        ('nonzero-start', dict(indptr=[1, 1, 2, 3])),
        ('bad-end', dict(indptr=[0, 1, 2, 2])),
        ('reverse-pointer', dict(indptr=[0, 2, 1, 3])),
        ('pointer-past-end', dict(indptr=[0, 4, 4, 4])),
        ('short-indices', dict(indices=[1, 3])),
        ('short-data', dict(data=[2., -1.])),
        ('column-out-of-range', dict(indices=[1, 4, 0])),
        ('unsorted-columns', dict(rows=1, cols=4, indptr=[0, 2],
                                  indices=[3, 1], data=[2., 3.])),
        ('duplicate-columns', dict(rows=1, cols=4, indptr=[0, 2],
                                   indices=[1, 1], data=[2., 3.])),
        ('axis-exceeds-limit', dict(rows=3, cols=1_000_001)),
    ]
    for name, patch in mutations:
        yield name, dict(base, **patch)


def _sparse_edges(api):
    from_coo = api == 'core.csr.from_coo'
    inputs = [
        ('empty-0x0', _coo(0, 0, [])),
        ('empty-0x3', _coo(0, 3, [])),
        ('empty-3x0', _coo(3, 0, [])),
        ('single-negative-zero', _coo(1, 1, [[0, 0, -0.]])),
        ('duplicates-cancel', _coo(2, 3, [[1, 2, -3.], [0, 0, 4.],
                                            [1, 2, 3.], [0, 0, -4.]])),
        ('unsorted-mixed', _coo(3, 4, [[2, 3, 4.], [0, 2, -1.],
                                          [1, 1, 0.], [0, 0, 7.]])),
        ('subnormal', _coo(1, 2, [[0, 0, 2.0 ** -1074],
                                   [0, 1, -2.0 ** -1074]])),
        ('large-finite', _coo(1, 2, [[0, 0, 1e150], [0, 1, -1e150]])),
    ]
    for name, source in inputs:
        if from_coo:
            yield _case(api, name, **source)
        else:
            values = _csr(source['rows'], source['cols'], source['coordinates'])
            operation = ('reject' if api == 'core.csr.to_dense'
                         and source['rows'] == 0 and source['cols'] > 0 else 'success')
            yield _case(api, name, operation, **values)
    if from_coo:
        invalid = [
            ('bad-row', _coo(2, 3, [[2, 0, 1.]])),
            ('bad-column', _coo(2, 3, [[0, 3, 1.]])),
            ('bad-empty-coordinate', _coo(0, 3, [[0, 1, 1.]])),
            ('large-row-axis', _coo(1_000_001, 1, [])),
            ('large-column-axis', _coo(1, 1_000_001, [])),
            ('duplicate-overflow', _coo(1, 1, [[0, 0, 1e308], [0, 0, 1e308]])),
        ]
        for name, args in invalid:
            yield _case(api, name, 'reject', **args)
    else:
        for name, args in _invalid_csr():
            yield _case(api, name, 'reject', **args)
        if api == 'core.csr.to_dense':
            yield _case(api, 'dense-work-limit', 'reject',
                        **_csr(1001, 1000, []))


def _signal_edges(api):
    for name, win, hop, desired, scaled in [
        ('rectangular-unscaled', [1.] * 8, 4, [1.] * 8, False),
        ('rectangular-scaled', [1.] * 8, 4, [1.] * 8, True),
        ('nonconstant-desired', [1.] * 8, 4, list(range(1, 9)), True),
        ('single-tap', [2.], 1, [3.], False),
        ('single-tap-scaled', [2.], 1, [3.], True),
        ('hop-one', [1., 2., 3., 4.], 1, [1., 2., 3., 4.], True),
        ('hop-window-size', [1., 2., 3., 4.], 4, [-1., 2., 3., 1.], False),
        ('small-amplitudes', [2e-6, 3e-6, 4e-6, 5e-6], 2,
         [1e-3, 2e-3, -1e-3, 3e-3], False),
    ]:
        yield _case(api, name, win=win, hop=hop, desired=desired, scaled=scaled)
    for name, win, hop, desired, scaled in [
        ('empty-window', [], 1, [], False),
        ('zero-hop', [1.] * 4, 0, [1.] * 4, False),
        ('hop-too-large', [1.] * 4, 5, [1.] * 4, True),
        ('desired-too-short', [1.] * 4, 2, [1.] * 3, True),
        ('nola-failure', [1., 0., 1., 0.], 2, [1.] * 4, False),
        ('all-zero-window', [0.] * 4, 2, [1.] * 4, True),
        ('scaled-zero-target', [1.] * 8, 4, [0.] * 8, True),
    ]:
        yield _case(api, name, 'reject', win=win, hop=hop,
                    desired=desired, scaled=scaled)


def cases(count, seed, apis=None):
    """Generate `count` seeded, shape-varied samples per requested export.

    Each API also receives a fixed edge/error set.  The API-specific RNG
    ensures filters and iteration order cannot change any input case bytes.
    """
    if count < 0:
        raise ValueError('count must be nonnegative')
    selected = set(APIS if apis is None else apis)
    unknown = selected.difference(APIS)
    if unknown:
        raise ValueError(f'unknown structural03 API(s): {sorted(unknown)}')
    for api in APIS:
        if api not in selected:
            continue
        if api.startswith('signal.'):
            yield from _signal_edges(api)
        else:
            yield from _sparse_edges(api)
        rng = _rng(api, seed)
        for i in range(count):
            if api.startswith('signal.'):
                n = (2, 3, 4, 5, 8, 11, 16, 17, 32, 63, 64)[i % 11]
                hop = rng.randrange(1, n + 1)
                shape = i % 6
                if shape == 0:
                    win = [1.] * n
                elif shape == 1:
                    win = [1. + .5 * math.sin(2 * math.pi * j / n) for j in range(n)]
                elif shape == 2:
                    win = [(-1.) ** j * (1. + j / n) for j in range(n)]
                else:
                    win = [rng.randrange(-64, 65) / 32 + .125 for _ in range(n)]
                # Ensure NOLA even for full-hop cases and cancellation.
                if shape > 2:
                    win = [w + (0.125 if w >= 0 else -0.125) for w in win]
                scaled = bool(i % 2)
                desired = ([w * (1. + rng.uniform(-.3, .3)) + rng.uniform(-.25, .25)
                            for w in win] if scaled else
                           [rng.randrange(-128, 129) / 64 + .25 for _ in range(n)])
                yield _case(api, f'sampled-{i:04d}', win=win, hop=hop,
                            desired=desired, scaled=scaled)
            else:
                source = _random_coo(rng, i)
                args = source if api == 'core.csr.from_coo' else _csr(
                    source['rows'], source['cols'], source['coordinates'])
                operation = ('reject' if api == 'core.csr.to_dense'
                             and source['rows'] == 0 and source['cols'] > 0 else 'success')
                yield _case(api, f'sampled-{i:04d}', operation, **args)


def _float_pair(value):
    value = float(value)
    if not math.isfinite(value):
        raise AssertionError('reference contains a nonfinite value')
    return [value, float(value == 0. and math.copysign(1., value) < 0.)]


def _canonical_matrix(a):
    rows, cols = a['rows'], a['cols']
    matrix = sparse.csr_array((np.asarray(a['data'], dtype=np.float64),
                               np.asarray(a['indices'], dtype=np.int64),
                               np.asarray(a['indptr'], dtype=np.int64)),
                              shape=(rows, cols))
    assert matrix.has_sorted_indices and matrix.has_canonical_format
    return matrix


def _valid_csr(a):
    rows, cols = a['rows'], a['cols']
    p, j, v = a['indptr'], a['indices'], a['data']
    if min(rows, cols) < 0 or max(rows, cols) > 1_000_000:
        return False
    if len(p) != rows + 1 or not p or p[0] != 0 or p[-1] != len(v):
        return False
    if len(v) != len(j) or len(v) > 1_000_000:
        return False
    if not all(math.isfinite(vv) for vv in v):
        return False
    for r in range(rows):
        lo, hi = p[r:r + 2]
        if lo > hi or lo < 0 or hi > len(v):
            return False
        if any(c < 0 or c >= cols for c in j[lo:hi]):
            return False
        if any(x >= y for x, y in zip(j[lo:hi], j[lo + 1:hi])):
            return False
    return True


def _reference_success(api, a):
    if api == 'core.csr.from_coo':
        rows, cols = a['rows'], a['cols']
        triples = a['coordinates']
        if not (0 <= rows <= 1_000_000 and 0 <= cols <= 1_000_000):
            raise ValueError('invalid axis')
        if any(r >= rows or c >= cols or r < 0 or c < 0 or not math.isfinite(v)
               for r, c, v in triples):
            raise ValueError('invalid COO coordinate')
        if not triples:
            mat = sparse.csr_array((rows, cols), dtype=np.float64)
        else:
            r, c, v = map(np.asarray, zip(*triples))
            mat = sparse.coo_array((v.astype(float), (r, c)), shape=(rows, cols)).tocsr()
        if not np.all(np.isfinite(mat.data)):
            raise ValueError('nonfinite duplicate accumulation')
        return [float(rows), float(cols), float(len(mat.indptr)),
                *map(float, mat.indptr), float(len(mat.indices)),
                *map(float, mat.indices), float(len(mat.data)),
                *(v for x in mat.data for v in _float_pair(x))]
    if not _valid_csr(a):
        raise ValueError('invalid canonical CSR')
    rows, cols = a['rows'], a['cols']
    if api == 'core.csr.validate':
        return [1.]
    if api in ('sparse.array.csr_indices_i32', 'sparse.array.csr_indices_i64'):
        # NumPy also verifies the explicit integer representation is exact.
        dtype = np.int32 if api.endswith('i32') else np.int64
        pointers = np.asarray(a['indptr'], dtype=dtype)
        indices = np.asarray(a['indices'], dtype=dtype)
        assert np.array_equal(pointers.astype(np.int64), a['indptr'])
        assert np.array_equal(indices.astype(np.int64), a['indices'])
        return [float(len(pointers)), *map(float, pointers),
                float(len(indices)), *map(float, indices)]
    if api == 'core.csr.to_dense':
        if (rows == 0 and cols != 0) or rows * cols > 1_000_000:
            raise ValueError('unrepresentable/oversized dense matrix')
        # SciPy checks numerical content; NumPy assignment preserves the
        # explicit -0 sign that scipy.csr_array.toarray() can discard.
        sci = _canonical_matrix(a)
        dense = np.zeros((rows, cols), dtype=float)
        for r in range(rows):
            for p in range(a['indptr'][r], a['indptr'][r + 1]):
                dense[r, a['indices'][p]] = a['data'][p]
        np.testing.assert_allclose(dense, sci.toarray(), rtol=0, atol=0)
        return [float(rows), float(cols),
                *(v for x in dense.flat for v in _float_pair(x))]
    raise ValueError(f'unsupported structural API: {api}')


def _signal_reference(a):
    win, desired = np.asarray(a['win'], dtype=np.float64), np.asarray(a['desired'], dtype=np.float64)
    hop, scaled = a['hop'], a['scaled']
    if (len(win) == 0 or not 1 <= hop <= len(win)
            or len(win) != len(desired)
            or not np.all(np.isfinite(win)) or not np.all(np.isfinite(desired))):
        raise ValueError('invalid dual inputs')
    energy = np.bincount(np.arange(len(win)) % hop, weights=win ** 2, minlength=hop)
    if np.any(energy <= np.finfo(float).eps):
        raise ValueError('NOLA failure')
    # The Rust implementation's scaled solver deliberately rejects a zero
    # target or unstable least-squares scaling.
    if scaled and not np.any(desired):
        raise ValueError('zero scaled target')
    dual, alpha = signal.closest_STFT_dual_window(win, hop, desired, scaled=scaled)
    assert dual.shape == win.shape and math.isfinite(float(alpha))
    assert np.all(np.isfinite(dual))
    residues = np.bincount(np.arange(len(win)) % hop,
                           weights=win * dual, minlength=hop)
    np.testing.assert_allclose(residues, np.ones(hop), rtol=3e-12, atol=3e-12)
    # For a fixed alpha, the correction to alpha*desired must be perpendicular
    # to the null space: it is constant times win within each hop residue.
    gap = dual - float(alpha) * desired
    for residue in range(hop):
        idx = np.arange(residue, len(win), hop)
        projection = np.dot(gap[idx], win[idx]) / energy[residue]
        np.testing.assert_allclose(gap[idx], projection * win[idx], atol=3e-11, rtol=3e-11)
    return [float(len(dual)), float(alpha), *map(float, dual)]


def reference(case):
    """Return the finite flattened independent reference for one case."""
    _versions()
    api, a = case['api'], case['args']
    if api not in APIS or a.get('structural03_case') is not True:
        raise ValueError('expected a structural03 case')
    if a['operation'] == 'reject':
        # Check the documented rejection contract independently. SciPy's
        # public signal API rejects these invalid-window cases too, while its
        # low-level CSR constructors can accept noncanonical storage.
        try:
            if api.startswith('signal.'):
                _signal_reference(a)
            else:
                _reference_success(api, a)
        except (ValueError, OverflowError):
            return [1.]
        raise AssertionError(f'invalid fixture unexpectedly valid: {case["id"]}')
    if a['operation'] != 'success':
        raise ValueError('unknown structural03 operation')
    result = _signal_reference(a) if api.startswith('signal.') else _reference_success(api, a)
    if not all(math.isfinite(float(v)) for v in result):
        raise AssertionError(f'nonfinite reference: {case["id"]}')
    return [float(v) for v in result]
