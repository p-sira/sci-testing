"""Generate Sci accuracy pages from measured outputs and independent references."""
import argparse
import cmath
from collections import defaultdict
from contextlib import redirect_stdout
import hashlib
import io
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import warnings

import filterpy
from filterpy.kalman import KalmanFilter
import mpmath as mp
import numpy as np
import scipy
from scipy import cluster, constants, differentiate, fft, integrate, interpolate
from scipy import linalg, optimize, signal, spatial, special, stats

ROOT = Path(__file__).resolve().parent
EPS = mp.mpf(2) ** -52
TINY = mp.mpf(2) ** -1022
MODULES = ('cluster complex constants differentiate fft filter integrate interpolate '
           'linalg optimize signal spatial special stats').split()
BOOST = 'https://www.boost.org/doc/libs/latest/libs/math/doc/html/math_toolkit/'
NOTES = {
    'cluster': ('Nearest-codebook distances use Euclidean geometry and first-minimum tie '
                'breaking. Code indices precede distances; whitening divides each column '
                'by its population standard deviation without centering. References use '
                'direct high-precision distances and moments. Hierarchy, k-means and other '
                'clustering functions remain unmeasured.'),
    'complex': ('Principal square roots and complex exponentials are compared component '
                'by component, real then imaginary. Inputs approach both sides of the '
                'negative-real branch cut without signed-zero ambiguity. Other complex '
                'operations and nonfinite inputs remain unmeasured.'),
    'constants': ('Vacuum wavelength is in metres and frequency in hertz; both conversions '
                  'use the exact SI value c = 299792458 m/s. NumPy evaluates c/x with '
                  'binary64 division. CODATA lookups, uncertainty tables, temperature '
                  'conversions, zero and nonfinite arguments remain unmeasured.'),
    'differentiate': ('The scalar callback is sin(x), with analytic reference cos(x). The '
                      'gradient callback is x² + 3xy + 2y², with reference [2x+3y, 3x+4y]. '
                      'Sci uses an absolute step of 1e-5. SciPy uses its adaptive order-8 '
                      'default algorithm with atol = rtol = 1e-10, so this compares achieved '
                      'accuracy at different work budgets. Truncation and subtraction '
                      'error can be much larger than epsilon. Hessians, vector Jacobians '
                      'and invalid steps remain unmeasured.'),
    'fft': ('Complex 1D transforms use the unnormalized negative-exponent forward DFT '
            'and positive-exponent inverse DFT divided by N. Results retain frequency '
            'order and flatten real then imaginary components for each bin. References '
            'are direct high-precision O(N²) sums. NumPy and SciPy may share PocketFFT '
            'ancestry. Real, multidimensional, cosine, sine, Hankel and periodic operators '
            'remain unmeasured.'),
    'filter': ('A correlated two-state Gaussian is predicted with F and additive Q, '
               'or updated from a scalar position observation with H and R. Output order '
               'is mean then row-major covariance. FilterPy uses its Joseph covariance '
               'update; Sci can differ at roundoff scale. References use high-precision '
               'matrix equations and Joseph covariance. EKF, UKF, particle filters, '
               'innovation diagnostics and invalid covariances remain unmeasured.'),
    'integrate': ('Trapezoid references evaluate the discrete rule on the supplied '
                  'binary64 samples, so sampling error against a continuous integral is '
                  'excluded. Fixed Gaussian quadrature uses order 16 on exp(x), compared '
                  'with the exact integral exp(b)-exp(a), including reversed and equal '
                  'bounds. SciPy receives a vectorized exponential callback. Other '
                  'integration rules and domain errors remain unmeasured.'),
    'interpolate': ('Distinct nodes sample x⁴-2x+3 exactly for the selected dyadic '
                    'coordinates. Queries include nodes, interior points and extrapolation. '
                    'The reference evaluates that polynomial in high precision. SciPy '
                    'barycentric weights use rng=0. Krogh uses no repeated-node derivative '
                    'data here. Splines, PCHIP, rational approximation and multidimensional '
                    'interpolation remain unmeasured.'),
    'linalg': ('General dense real systems use one right-hand side; solutions preserve '
               'coordinate order. Matrices are diagonally dominant with sizes 2, 3, 5 '
               'and 8. References use high-precision LU solve and determinant. Sci uses '
               'its default backend; NumPy and SciPy may share BLAS/LAPACK. Ill-conditioned '
               'or singular problems, decompositions and other solvers remain unmeasured.'),
    'optimize': ('The Rosenbrock objective is Σ[100(xᵢ₊₁-xᵢ²)²+(1-xᵢ)²]; its analytic '
                 'gradient preserves coordinate order. Cases include the minimizer and '
                 'nearby points. References directly evaluate these polynomials in high '
                 'precision. These results cover objective helpers only. Minimizers, '
                 'convergence status and finite-difference diagnostics remain unmeasured.'),
    'signal': ('Real convolution and correlation use full mode; correlation orders lags '
               'from -(len(h)-1) through len(x)-1. SciPy is forced to method="direct" '
               'and NumPy uses its direct API. References use high-precision discrete '
               'sums. Spectral analysis, filtering, windows, filter design, peaks and '
               'other signal functions remain unmeasured.'),
    'spatial': ('Unweighted Euclidean distance and cosine distance are measured on real '
                'vectors, including equal and nearly equal vectors and large finite '
                'coordinates. Cosine distance is 1-dot(u,v)/(norm(u)norm(v)), clipped to '
                '[0,2]. References use high-precision dot products and norms. Near-zero '
                'cosine distances amplify relative error; inspect absolute error too. '
                'The SciPy cosine evaluation on large coordinates emits an intermediate '
                'overflow warning and returns a finite inaccurate result; this is retained '
                'as an error measurement, with the warning in the raw data. Other '
                'distances, spatial structures and geometry remain unmeasured.'),
    'special': ('Real erf and gamma are compared with mpmath references. Gamma includes '
                'negative noninteger arguments but excludes poles. The stdlib column '
                'uses Python math and the host C math library. Sci and SciPy can share '
                'underlying special-function algorithms; agreement alone does not '
                'establish accuracy. Other special functions, poles and nonfinite '
                'arguments remain unmeasured.'),
    'stats': ('Normal CDF uses mean 0 and standard deviation 1; references use '
              'erfc(-x/√2)/2 to retain lower tails. Variance uses ddof=1, with constant, '
              'offset and very small samples, evaluated against high-precision central '
              'moments of the exact binary64 inputs. The stdlib uses statistics.variance '
              'for variance and erfc for CDF. Other distributions, inference, sampling '
              'and invalid inputs remain unmeasured.'),
}


def command(*args, cwd=ROOT):
    return subprocess.check_output(args, cwd=cwd, text=True).strip()


def high_precision(case):
    """References depend on mathematical definitions, never Sci or its backends."""
    api, a = case['api'], case['args']
    def num(x): return mp.mpf(x)  # Exact conversion of the parsed binary64 input.
    def vec(key): return list(map(num, a[key]))
    def mat(key): return mp.matrix([[num(x) for x in r] for r in a[key]])
    def flat(m): return list(m)
    def gauss(mean, cov): return flat(mean) + flat(cov)
    if api.startswith('complex.'):
        z = mp.mpc(num(a['x']), num(a['y']))
        r = mp.exp(z) if api.endswith('exp') else mp.sqrt(z)
        return [r.real, r.imag]
    if api.startswith('constants.'): return [mp.mpf(299792458)/num(a['x'])]
    if api == 'differentiate.derivative': return [mp.cos(num(a['x']))]
    if api == 'differentiate.gradient':
        x, y = vec('point'); return [2*x+3*y, 3*x+4*y]
    if api.startswith('fft.'):
        z = [mp.mpc(*map(num, r)) for r in a['z']]
        n = len(z); inverse = api == 'fft.ifft'; result = []
        for k in range(n):
            terms = []
            for j, value in enumerate(z):
                residue = (k*j) % n
                # Cardinal roots are exact, avoiding artificial sin(pi) residues.
                if (4*residue) % n == 0:
                    root = [mp.mpc(1),mp.mpc(0,1),mp.mpc(-1),mp.mpc(0,-1)][(4*residue//n)%4]
                else: root = mp.exp(2*mp.pi*mp.j*residue/n)
                terms.append(value * (root if inverse else root.conjugate()))
            value = mp.fsum(terms)/(n if inverse else 1)
            result.extend([value.real, value.imag])
        return result
    if api.startswith('filter.'):
        mean, cov = mp.matrix(vec('mean')), mat('covariance')
        if api.endswith('predict'):
            f, q = mat('f'), mat('q'); return gauss(f*mean, f*cov*f.T+q)
        h, r, z = mat('h'), mat('r'), mp.matrix(vec('z'))
        k = cov*h.T*(h*cov*h.T+r)**-1
        residual = mp.eye(len(mean))-k*h
        return gauss(mean+k*(z-h*mean),residual*cov*residual.T+k*r*k.T)
    if api == 'integrate.trapezoid':
        x, y = vec('coordinates'), vec('y')
        return [mp.fsum((x[i+1]-x[i])*(y[i]+y[i+1])/2 for i in range(len(y)-1))]
    if api == 'integrate.fixed_quad': return [mp.exp(num(a['b']))-mp.exp(num(a['a']))]
    if api.startswith('interpolate.'): return [x**4-2*x+3 for x in vec('query')]
    if api == 'linalg.solve': return flat(mp.lu_solve(mat('a'),mp.matrix(vec('b'))))
    if api == 'linalg.det': return [mp.det(mat('a'))]
    if api.startswith('optimize.'):
        x = vec('point')
        if api.endswith('rosen'): return [mp.fsum(100*(v-u*u)**2+(1-u)**2 for u,v in zip(x,x[1:]))]
        r = [mp.mpf(0)]*len(x)
        for i in range(len(x)-1):
            r[i] += -400*x[i]*(x[i+1]-x[i]**2)+2*(x[i]-1)
            r[i+1] += 200*(x[i+1]-x[i]**2)
        return r
    if api.startswith('signal.'):
        x,h = vec('x'),vec('h')
        if api.endswith('correlate'): h = h[::-1]
        return [mp.fsum(x[j]*h[k-j] for j in range(len(x)) if 0 <= k-j < len(h)) for k in range(len(x)+len(h)-1)]
    if api.startswith('spatial.'):
        u,v = vec('u'),vec('v')
        if api.endswith('euclidean'): return [mp.sqrt(mp.fsum((x-y)**2 for x,y in zip(u,v)))]
        dot = mp.fsum(x*y for x,y in zip(u,v))
        norm_product = mp.sqrt(mp.fsum(x*x for x in u)*mp.fsum(y*y for y in v))
        if dot > 0:
            # Lagrange's identity avoids cancellation for nearly parallel vectors.
            wedge = mp.fsum((u[i]*v[j]-u[j]*v[i])**2 for i in range(len(u)) for j in range(i+1,len(u)))
            distance = wedge/(norm_product*(norm_product+dot))
        else:
            distance = 1-dot/norm_product
        return [max(mp.mpf(0),min(mp.mpf(2),distance))]
    if api == 'special.erf': return [mp.erf(num(a['x']))]
    if api == 'special.gamma': return [mp.gamma(num(a['x']))]
    if api == 'stats.normal_cdf': return [mp.erfc(-num(a['x'])/mp.sqrt(2))/2]
    if api == 'stats.variance':
        v = vec('values'); mean = mp.fsum(v)/len(v)
        return [mp.fsum((x-mean)**2 for x in v)/(len(v)-1)]
    if api.startswith('cluster.'):
        obs = [[num(x) for x in r] for r in a['obs']]
        if api.endswith('vq'):
            codebook = [[num(x) for x in r] for r in a['codebook']]
            ds = [[mp.sqrt(mp.fsum((x-y)**2 for x,y in zip(row,c))) for c in codebook] for row in obs]
            codes = [min(range(len(row)),key=row.__getitem__) for row in ds]
            return list(map(mp.mpf,codes))+[row[c] for row,c in zip(ds,codes)]
        columns = list(zip(*obs)); std = []
        for col in columns:
            mean = mp.fsum(col)/len(col)
            std.append(mp.sqrt(mp.fsum((x-mean)**2 for x in col)/len(col)))
        return [x/std[j] for row in obs for j,x in enumerate(row)]
    raise ValueError(api)


def flatten(value):
    a = np.asarray(value)
    if np.iscomplexobj(a): return np.stack([a.real,a.imag],axis=-1).ravel().tolist()
    return a.astype(float).ravel().tolist()


def competitors(case):
    """Public library APIs, or explicitly labelled elementary NumPy expressions."""
    api, a = case['api'], case['args']
    v = lambda key: np.asarray(a[key],dtype=float)
    calls = {}
    if api.startswith('complex.'):
        z = complex(a['x'],a['y']); name = api.split('.')[1]
        calls = {'NumPy':lambda:getattr(np,name)(z),'stdlib':lambda:getattr(cmath,name)(z)}
    elif api.startswith('constants.'):
        name = api.split('.')[1]
        calls = {'SciPy':lambda:getattr(constants,name)(a['x']),'NumPy':lambda:np.divide(299792458.,a['x'])}
    elif api == 'differentiate.derivative':
        def derivative():
            r = differentiate.derivative(np.sin,a['x'],tolerances={'atol':1e-10,'rtol':1e-10})
            if not r.success: raise ValueError(f'SciPy status {r.status}')
            return r.df
        calls = {'SciPy':derivative}
    elif api == 'differentiate.gradient':
        def gradient():
            r = differentiate.jacobian(lambda p:np.asarray([p[0]**2+3*p[0]*p[1]+2*p[1]**2]),v('point'),tolerances={'atol':1e-10,'rtol':1e-10})
            if not np.all(r.success): raise ValueError(f'SciPy status {r.status}')
            return r.df
        calls = {'SciPy':gradient}
    elif api.startswith('fft.'):
        name = api.split('.')[1]; z = v('z')[:,0]+1j*v('z')[:,1]
        calls = {'SciPy':lambda:getattr(fft,name)(z),'NumPy':lambda:getattr(np.fft,name)(z)}
    elif api.startswith('filter.'):
        def kalman():
            k = KalmanFilter(dim_x=2,dim_z=1)
            k.x = v('mean').reshape(-1,1); k.P = v('covariance')
            if api.endswith('predict'):
                k.F=v('f'); k.Q=v('q'); k.predict()
            else:
                k.H=v('h'); k.R=v('r'); k.update(v('z'))
            return np.concatenate([k.x.ravel(),k.P.ravel()])
        calls = {'FilterPy':kalman}
    elif api == 'integrate.trapezoid':
        calls = {'SciPy':lambda:integrate.trapezoid(v('y'),x=v('coordinates')),'NumPy':lambda:np.trapezoid(v('y'),x=v('coordinates'))}
    elif api == 'integrate.fixed_quad': calls = {'SciPy':lambda:integrate.fixed_quad(np.exp,a['a'],a['b'],n=16)[0]}
    elif api.startswith('interpolate.'):
        name = api.split('.')[1]
        calls = {'SciPy':lambda:getattr(interpolate,name)(v('nodes'),v('values'),v('query'),**({'rng':0} if name.startswith('barycentric') else {}))}
    elif api.startswith('linalg.'):
        name = api.split('.')[1]; args = (v('a'),v('b')) if name == 'solve' else (v('a'),)
        calls = {'SciPy':lambda:getattr(linalg,name)(*args),'NumPy':lambda:getattr(np.linalg,name)(*args)}
    elif api.startswith('optimize.'): calls = {'SciPy':lambda:getattr(optimize,api.split('.')[1])(v('point'))}
    elif api.startswith('signal.'):
        name = api.split('.')[1]
        calls = {'SciPy':lambda:getattr(signal,name)(v('x'),v('h'),mode='full',method='direct'),'NumPy':lambda:getattr(np,name)(v('x'),v('h'),mode='full')}
    elif api.startswith('spatial.'): calls = {'SciPy':lambda:getattr(spatial.distance,api.split('.')[1])(v('u'),v('v'))}
    elif api.startswith('special.'):
        name = api.split('.')[1]
        calls = {'SciPy':lambda:getattr(special,name)(a['x']),'stdlib':lambda:getattr(math,name)(a['x'])}
    elif api == 'stats.normal_cdf': calls = {'SciPy':lambda:stats.norm.cdf(a['x']),'stdlib':lambda:math.erfc(-a['x']/math.sqrt(2))/2}
    elif api == 'stats.variance':
        import statistics
        calls = {'SciPy':lambda:stats.describe(v('values'),ddof=1).variance,'NumPy':lambda:np.var(v('values'),ddof=1),'stdlib':lambda:statistics.variance(a['values'])}
    elif api == 'cluster.vq': calls = {'SciPy':lambda:np.concatenate(cluster.vq.vq(v('obs'),v('codebook')))}
    elif api == 'cluster.whiten': calls = {'SciPy':lambda:cluster.vq.whiten(v('obs'))}
    else: raise ValueError(api)
    results = {}
    for library, call in calls.items():
        try:
            with warnings.catch_warnings(record=True) as captured:
                warnings.simplefilter('always')
                values = flatten(call())
            if not all(math.isfinite(x) for x in values): raise ValueError('nonfinite output')
            results[library] = {'values':values}
            if captured:
                results[library]['warnings'] = [str(w.message) for w in captured]
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as error:
            results[library] = {'error':f'{type(error).__name__}: {error}'}
    return results


def symmetric_error(actual, reference):
    """Boost's symmetric denominator with an explicit zero/subnormal adaptation."""
    actual, reference = mp.mpf(actual), mp.mpf(reference)
    if not mp.isfinite(actual) or not mp.isfinite(reference):
        raise ValueError('nonfinite accuracy sample')
    # All subnormals are treated as effectively zero, as in Boost's test reports.
    a = mp.mpf(0) if abs(actual) < TINY else actual
    b = mp.mpf(0) if abs(reference) < TINY else reference
    if a == b: return mp.mpf(0)
    if not a or not b: return 1/EPS
    return abs(a-b)/min(abs(a),abs(b))/EPS


def digest_sources():
    names = ['Cargo.toml','Cargo.lock','pyproject.toml','uv.lock','src/main.rs','report.py','test_report.py','cases.json','README.md']
    h = hashlib.sha256()
    for name in names:
        h.update(name.encode()+b'\0'+(ROOT/name).read_bytes()+b'\0')
    return h.hexdigest()


def environment(sci_path):
    configs = io.StringIO()
    with redirect_stdout(configs):
        np.show_config(); scipy.show_config()
    return dict(sci_commit=command('git','rev-parse','HEAD',cwd=sci_path),
                sci_source_dirty=bool(command('git','status','--porcelain','--','src','crates','Cargo.toml','Cargo.lock',cwd=sci_path)),
                testing_sources_sha256=digest_sources(),cases_sha256=hashlib.sha256((ROOT/'cases.json').read_bytes()).hexdigest(),
                rust=command('rustc','-Vv'),cargo=command('cargo','-V'),
                python=sys.version.split()[0],numpy=np.__version__,scipy=scipy.__version__,
                mpmath=mp.__version__,filterpy=filterpy.__version__,platform=platform.platform(),
                machine=platform.machine(),processor=next((l.split(':',1)[1].strip() for l in Path('/proc/cpuinfo').read_text().splitlines() if l.startswith('model name')),'unknown'),
                profile='dev; opt-level=0; no custom RUSTFLAGS',reference_decimal_digits=[100,160],
                library_config=configs.getvalue())


def number(x): return mp.nstr(x,5,min_fixed=-3,max_fixed=5)


def render(cases, records, meta, output):
    output.mkdir(parents=True,exist_ok=True)
    for module in MODULES:
        selected = [c for c in cases if c['api'].split('.')[0] == module]
        libraries = [lib for lib in ['Sci','SciPy','NumPy','stdlib','FilterPy'] if any(lib in records[c['id']]['libraries'] for c in selected)]
        groups = defaultdict(list)
        for case in selected: groups[(case['api'],case['group'])].append(case)
        rows = []
        for (api,group),batch in groups.items():
            cells = []
            for lib in libraries:
                errors=[]; absolute=[]; failures=0; unavailable=0
                for case in batch:
                    r = records[case['id']]; result = r['libraries'].get(lib)
                    if result is None: unavailable+=1; continue
                    if 'error' in result: failures+=1; continue
                    if len(result['values']) != len(r['reference']): raise ValueError(f'{case["id"]}: shape mismatch in {lib}')
                    skip = len(case['args']['obs']) if case['api']=='cluster.vq' else 0
                    for actual,ref in zip(result['values'][skip:],r['reference'][skip:]):
                        reference=mp.mpf(ref)
                        errors.append(symmetric_error(actual,reference)); absolute.append(abs(mp.mpf(actual)-reference))
                if unavailable == len(batch): cells.append('—'); continue
                if unavailable: raise ValueError('partially missing implementation')
                text = f'Max = {number(max(errors))}ε (Mean = {number(mp.fsum(errors)/len(errors))}ε); ∣Δ∣max = {number(max(absolute))}' if errors else 'No finite results'
                cells.append(text+f'; failures = {failures}/{len(batch)}')
            n = sum(len(records[c['id']]['reference'])-(len(c['args']['obs']) if c['api']=='cluster.vq' else 0) for c in batch)
            rows.append(f'| `{api.split(".",1)[1]}`: {group} | {len(batch)} / {n} | '+' | '.join(cells)+' |')
        text = f'''# {module}: accuracy comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-{module})

Measured binary64 (`f64`) results against independent mpmath references at
160 decimal digits, checked against a separate 100-digit evaluation. This page
covers {len(selected)} cases and {len(set(c['api'] for c in selected))} public APIs;
coverage is limited to the operations and inputs below.

## Accuracy

Each library cell follows Boost.Math's **Max = …ε (Mean = …ε)** format.
Lower values indicate smaller errors on this dataset. `∣Δ∣max` is the maximum
absolute error in output units. Failures count cases; statistics include only
successful finite scalar outputs. A dash means no comparator was run for that
operation. The mean is arithmetic across scalar components, including exact zeros.

| Function and test group | Cases / scalar outputs | {' | '.join(libraries)} |
| --- | ---: | {' | '.join('---' for _ in libraries)} |
'''+ '\n'.join(rows)+f'''

## Conventions and coverage gaps

{NOTES[module]}

Inputs are fixed in `sci-testing/cases.json`; raw per-case outputs, high-precision
references and the environment are stored in that repository's `results/`.
NaN, infinity, signed-zero contracts, invalid domains and work-limit enforcement
are outside this finite-input accuracy sample. Existing Sci API tests cover
additional behavior. No claims of full inventory coverage or worst-case bounds
follow from these tables. SciPy aliases, compatibility namespaces and Cython
variants are not counted as independent implementations.

## Reproduction

- Sci source: `{meta['sci_commit']}` (numerical sources clean: {str(not meta['sci_source_dirty']).lower()}).
- Testing source SHA-256: `{meta['testing_sources_sha256']}`.
- NumPy {meta['numpy']}; SciPy {meta['scipy']}; FilterPy {meta['filterpy']}; mpmath {meta['mpmath']}.
- {meta['rust'].splitlines()[0]}; Python {meta['python']}; Cargo dev profile, opt-level 0.
- {meta['platform']}; {meta['processor']}.

See the [shared metric, provenance and reproduction instructions](README.md).
'''
        (output/f'{module}.md').write_text(text)
    index=f'''# Module accuracy comparisons

Measured comparisons for all 14 public numerical modules. Each page uses
Boost.Math's [accuracy table layout]({BOOST}sf_erf/error_function.html):
functions and input regions form rows, libraries form columns, and cells show
maximum and mean error in units of machine epsilon. These are representative
samples, with explicit coverage gaps on each page.

| Module | Compared APIs | Cases |
| --- | ---: | ---: |
'''
    for module in MODULES:
        batch=[c for c in cases if c['api'].split('.')[0]==module]
        index+=f'| [{module}]({module}.md) | {len(set(c["api"] for c in batch))} | {len(batch)} |\n'
    index+='''
## Error definition

For finite nonzero output `v` and reference `r`, use the symmetric relative error

$$
E(v,r) = \\frac{|v-r|}{\\min(|v|,|r|)}, \\qquad
E_\\varepsilon = \\frac{E(v,r)}{\\varepsilon}, \\qquad
\\varepsilon = 2^{-52}.
$$

This is the symmetric denominator described by
[Boost.Math's relative error report]('''+BOOST+'''relative_error.html).
The arithmetic mean and maximum aggregate scalar errors within each row;
complex numbers contribute real and imaginary components separately, and
matrices flatten in row-major order. Vector/matrix lengths must match references.
Codebook indices are checked exactly and excluded from the scalar-output counts
and distance-error aggregates.

The actual binary64 result is promoted exactly to high precision and compared
with the unrounded reference. Fractional epsilon errors therefore remain visible;
zero error indicates exact agreement on the sampled output, not a general
rounding guarantee. Epsilon-scaled relative error is not a count of ULPs.

For the tables' zero convention, magnitudes below the smallest positive normal
binary64 value (`2^-1022`) are treated as zero. Two effective zeros have zero
relative error; exactly one effective zero has relative error 1 (`2^52` epsilon).
Opposite finite signs use the stated formula. These conventions adapt Boost's
[test-report treatment of zeros]('''+BOOST+'''relative_error.html);
they do not reproduce every sentinel rule of its `relative_difference` utility.
Absolute error always uses the original values, including subnormals. Near zeros
or cancellation, use the absolute error alongside the relative error.

Nonfinite outputs, exceptions and Sci errors count as failures and remain in raw
results. Comparator warnings are retained with each raw result, including warnings
that accompany finite inaccurate values. They cannot silently improve a table: every cell reports its failure
count. Unexpected output shapes or incorrect discrete code assignments stop report
generation. This suite currently samples finite valid inputs; special-value and
invalid-input contracts require separate API tests.

## References and shared algorithms

References come from mpmath 1.3.0 and independent mathematical definitions at
100 and 160 decimal digits. Parsed binary64 inputs are converted exactly, so the
comparison measures computation error rather than decimal input conversion.
The two evaluations must agree within `1e-90` times the largest reference
magnitude in that case (at least the smallest normal binary64 value). This
precision check is empirical evidence, not a rigorous error certificate.

References include analytic derivatives and integrals, direct DFT and convolution
sums, polynomial evaluation, high-precision LU solves, and Gaussian update
identities. Neither Sci nor its implementation backends generates references.
Python standard-library results use `math`, `cmath` or `statistics`, as described
on each page. NumPy constant conversions are explicitly labelled arithmetic
expressions. Versions and BLAS/LAPACK build details are recorded in `results/environment.json`.
Libraries can share algorithms or dependencies: matching Sci/SciPy results, or
NumPy/SciPy results, do not establish independence or a universal ranking.

## Reproduce in the separate repository

The harness is in the independent local Git repository `sci-testing`, next to
`sci`. It contains the Rust public-API runner, fixed `cases.json` inputs, Python
reference/comparison code, metric tests, Cargo and uv lockfiles, and raw results.
There is no comparison runner or additional dependency in the Sci workspace.
The testing repository has no published remote yet.

With sibling checkouts of `sci` and `sci-testing`, use:

```sh
cd ../sci-testing
uv sync --frozen
uv run --frozen python -m unittest -v
uv run --frozen python report.py --output ../sci/docs/accuracy
```

The runner requires Rust, a C++17 compiler and libclang, as does Sci's default
feature set. Its Cargo manifest uses `sci = { path = "../sci" }` with all default
modules. `Cargo.lock` pins the harness's resolved backend versions independently
of Sci's workspace lockfile. The current run uses the Cargo dev profile with no
custom compiler flags or target overrides; these settings are enforced by the
report command. The command builds the runner, retains raw outputs and references,
and writes all module pages and this index. To reuse an existing build directory,
set `CARGO_TARGET_DIR`; it does not change the selected Cargo profile.

For the published snapshot, check out the Sci commit recorded below and the
matching testing sources; `results/environment.json` records a SHA-256 over all
harness sources, inputs and lockfiles. This digest excludes generated results
and pages to avoid a circular dependency. Each module page repeats the source
revision, digest and platform. Regeneration refuses modified numerical sources.
Documentation edits are permitted. The schema and exact digest file list live
in `sci-testing/report.py`.

## Snapshot

'''
    index+=f'- Sci commit: `{meta["sci_commit"]}`.\n- Harness source SHA-256: `{meta["testing_sources_sha256"]}`.\n- Input SHA-256: `{meta["cases_sha256"]}`.\n- {len(cases)} cases; {len(set(c["api"] for c in cases))} measured public APIs.\n'
    (output/'README.md').write_text(index)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT.parent/'sci/docs/accuracy')
    args=parser.parse_args()
    if (np.__version__,scipy.__version__,mp.__version__,filterpy.__version__) != ('2.5.4','1.18.1','1.3.0','1.4.5'):
        raise RuntimeError('Run uv sync --frozen with the pinned library versions')
    if any(k.startswith(('RUSTFLAGS','CARGO_ENCODED_RUSTFLAGS','CARGO_PROFILE_')) or k == 'CARGO_BUILD_TARGET' for k in os.environ):
        raise RuntimeError('Unset custom compiler/profile/target overrides for this snapshot')
    cases=json.loads((ROOT/'cases.json').read_text())
    ids=[c['id'] for c in cases]
    if len(ids) != len(set(ids)) or {c['api'].split('.')[0] for c in cases} != set(MODULES):
        raise ValueError('case identifiers must be unique and cover every module')
    meta=environment(ROOT.parent/'sci')
    if meta['sci_source_dirty']: raise RuntimeError('Refusing to measure modified Sci numerical sources')
    subprocess.run(['cargo','build','--locked'],cwd=ROOT,check=True)
    target=Path(os.environ.get('CARGO_TARGET_DIR',ROOT/'target')).resolve()
    run=subprocess.run([str(target/'debug/sci-testing')],input=''.join(json.dumps(c)+'\n' for c in cases),text=True,capture_output=True,check=True)
    rust=[json.loads(line) for line in run.stdout.splitlines()]
    if [r['id'] for r in rust] != ids: raise ValueError('runner lost, duplicated or reordered cases')
    records={}
    for case,result in zip(cases,rust):
        with mp.workdps(100): low=high_precision(case)
        with mp.workdps(160):
            high=high_precision(case)
            if len(low) != len(high): raise ValueError('unstable reference shape')
            scale=max([TINY]+[abs(v) for v in high])
            if any(abs(x-y)>mp.mpf('1e-90')*scale for x,y in zip(low,high)):
                raise ValueError(f'unstable reference: {case["id"]}')
            libs={'Sci':{k:v for k,v in result.items() if k!='id'},**competitors(case)}
            for lib,res in libs.items():
                if 'values' in res:
                    if len(res['values']) != len(high): raise ValueError(f'{case["id"]}: {lib} shape mismatch')
                    if case['api']=='cluster.vq' and any(mp.mpf(x)!=r for x,r in zip(res['values'][:len(case['args']['obs'])],high)):
                        raise ValueError(f'incorrect code assignment: {case["id"]}: {lib}')
            records[case['id']]={'id':case['id'],'reference':[mp.nstr(x,165) for x in high],'libraries':libs}
    results=ROOT/'results'; results.mkdir(exist_ok=True)
    (results/'environment.json').write_text(json.dumps(meta,indent=2)+'\n')
    (results/'measurements.json').write_text(json.dumps(list(records.values()),indent=2,allow_nan=False)+'\n')
    with mp.workdps(160): render(cases,records,meta,args.output.resolve())
    print(f'Wrote {len(MODULES)} module pages, {len(cases)} cases to {args.output}')


if __name__ == '__main__': main()
