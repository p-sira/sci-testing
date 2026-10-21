"""Generate library comparison pages using measured Sci outputs as the baseline."""
import argparse
import cmath
import gc
from collections import defaultdict
from contextlib import redirect_stdout
import hashlib
import io
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
import warnings

from .native_compare import NativeComparisons

from .generate_cases import build_cases
from .case_csv import (case_files, disk_files, files_digest, measurement_files,
                      read_measurements, write_files, MeasurementStore, directory_digest, iter_case_files, write_cases,
                      csv_paths, write_measurements)

import filterpy
from filterpy.kalman import KalmanFilter
import mpmath as mp
import numpy as np
import scipy
from scipy import cluster, constants, differentiate, fft, integrate, interpolate
from scipy import linalg, optimize, signal, spatial, special, stats

from . import ROOT
EPS = mp.mpf(2) ** -52
TINY = mp.mpf(2) ** -1022
MODULES = ('cluster complex constants differentiate fft filter integrate interpolate '
           'linalg optimize signal sparse spatial special stats').split()
BOOST = 'https://www.boost.org/doc/libs/latest/libs/math/doc/html/math_toolkit/'
NOTES = {
    'sparse':'Canonical CSR and COO arrays retain explicit zeros. Matrix outputs include dimensions; graph distances retain infinities. Eigenpairs are compared by ordered values with independent residual checks; matching cardinality, permutation validity, flow conservation and random shape/count use invariants where outputs are nonunique. Sparse LU, iterative solves and matrix functions use pinned public SciPy references. Each scalar uses (class,value): class 0 for finite data, ±1 for signed infinity and 2 for NaN, with zero nonfinite payloads. This preserves unreachable vertices; exact tags contribute to component means. RNG streams differ. Dense samples preserve graph structure, integer capacities and dependent dimensions; fixed discrete profiles may repeat. Rust signature adaptations and bounded dense fallbacks are documented in Sci sparse coverage.',
    'cluster':'Vector quantization uses Euclidean distances and first-minimum ties; whitening divides by population standard deviations. Hierarchy fixtures also check counts, leaf ordering and partitions. K-means receives an explicit initial codebook. Returned discrete indices and tree metadata are recorded alongside distances.',
    'complex':'Complex outputs flatten real then imaginary components. Functions use principal branches; powers preserve their explicit real, complex or integer exponent conventions. The original exp/sqrt recipes use high-precision mpmath references; additional complex comparisons use the Python standard library as a binary64 reference. Signed-zero and nonfinite contracts remain in the Sci API tests.',
    'constants':'Vacuum conversions use metres and hertz and the exact SI value c = 299792458 m/s. Temperature recipes convert Celsius to Fahrenheit. CODATA values, uncertainties, unit strings and searches use the pinned SciPy constants table. Text metadata is encoded as UTF-8 lengths and bytes; Sci dtype f64 corresponds to NumPy float64.',
    'differentiate':'Callbacks are supplied explicitly: original sin and quadratic recipes, a quadratic Hessian and a two-output Jacobian. Fixed-step Sci stencils and adaptive SciPy differences have different work budgets. Analytic identities supply the callback references; raw inputs retain the selected steps and points.',
    'fft':'Transforms use unnormalized forward and inverse-divided-by-length conventions. Real transforms, multidimensional shapes, transform types, FFTLog and periodic operators follow the pinned fixture conventions. FFTLog uses SciPy fft; distinct periodic operators use its legacy fftpack implementation. Real and imaginary components are interleaved, including exact imaginary zeros.',
    'filter':'Gaussian outputs place means before row-major covariance; update and particle diagnostics follow the pinned filter fixture ordering. References use independent matrix identities, supplied deterministic particle draws and nonlinear callbacks. Original linear Kalman comparisons additionally call FilterPy.',
    'integrate':'Discrete-rule references integrate the supplied binary64 samples, excluding continuous sampling error. The dense added rules use dx = 0.25 and zero initial cumulative integrals. Newton-Cotes records weights followed by its error coefficient. Fixed Gaussian quadrature uses order 16 on exp(x), with analytic references.',
    'interpolate':'Original barycentric/Krogh recipes sample an exactly represented quartic. Added fitting and spline recipes compare evaluation, knots, coefficients and degree in the established fixture order. Queries, boundary modes, repeated knots and work sizes remain explicit in the CSV inputs. Equivalent spline representations can differ at roundoff scale.',
    'linalg':'Nonunique decompositions are compared through ordered eigenvalues, projectors and reconstructions, retaining independent structural invariants from the Sci fixture dispatcher. Interpolative decomposition compares spectral residual and rank; its SVD compares reconstruction and singular values. Type metadata uses the Sci f64 spelling. Dense sample sizes and constructor orders stay bounded.',
    'optimize':'Objective helpers use their pinned Rosenbrock and finite-difference conventions. Minimizers use 3+(x-1)^2+(y+2)^2, a 512-iteration budget and tolerance 1e-8; only converged points and costs are compared. Iteration trajectories are not compared. Analytic minimizers supply the reference; SciPy runs Nelder-Mead or BFGS separately.',
    'signal':'The pinned signal fixture adapters cover waveform, filtering, design, spectral, window, peak, spline, convolution and system functions through the public Sci facade. Output tuples flatten in the documented fixture order. Parameters, modes, shapes and normalization are retained in each function CSV; original real convolution/correlation also call NumPy.',
    'spatial':'Metrics, pairwise arrays, geometry and validation helpers use explicit weights, variance/covariance data and shapes from the pinned spatial fixtures. Output ordering follows those fixtures. Near-zero distances and cancellation can amplify relative error; the CSV retains full outputs for absolute-error inspection.',
    'special':'All public special-function implementations have typed registrations or an explicit policy recipe. Original erf/gamma recipes use high-precision mpmath; added recipes use independent public SciPy calls, identities or quadrature. Polynomial constructors return descending coefficients; sequence and tuple ordering follows the pinned fixtures. Boost.Math/GSL are called where the bridge has matching scalar APIs, converting elliptic parameter m to modulus sqrt(m).',
    'stats':'Distribution parameters, tail conventions, alternatives, ddof, explicit resampling draws and output ordering follow the pinned statistics fixtures. SciPy supplies binary64 references for added APIs, including documented adaptations and independently seeded probabilistic references. Weighted tau compares its statistic; its always-NaN p-value is left to the ordinary API tests. Normality, inverse and joint-tail approximations have reference and algorithm accuracy limits.',
}



def command(*args, cwd=ROOT):
    return subprocess.check_output(args, cwd=cwd, text=True).strip()


def binary64_reference(case):
    a=case['args']
    return any(key in a for key in ['p0','p','arguments','raw']) or bool(a.get('explicit'))


def high_precision(case):
    """References depend on mathematical definitions, never Sci or its backends."""
    api, a = case['api'], case['args']
    if binary64_reference(case):
        from .registered_cases import reference
        return list(map(mp.mpf, reference(case)))
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
    if 'p0' in a or 'p' in a or 'arguments' in a or 'raw' in a or a.get('explicit'):
        from .registered_cases import reference, reference_result
        library='stdlib' if api.startswith('complex.') else 'NumPy' if a.get('explicit') and api.startswith('linalg.') else 'SciPy'
        if api in ['special.geterr','special.seterr']:return {}
        try:
            with warnings.catch_warnings(record=True) as captured:
                warnings.simplefilter('always')
                if a.get('explicit'):
                    from .explicit_cases import comparison
                    values=comparison(case)
                else:
                    saved=reference_result(case);values=saved['values']
            if not all(math.isfinite(x) for x in values):raise ValueError('nonfinite output')
            result={'values':list(map(float,values))}
            observed=[str(w.message) for w in captured]
            if not a.get('explicit'):observed+=saved['warnings']
            if observed:result['warnings']=observed
            return {library:result}
        except (ValueError,ArithmeticError,TypeError,np.linalg.LinAlgError) as error:
            return {library:{'error':f'{type(error).__name__}: {error}'}}
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
    elif api == 'cluster.vq.vq': calls = {'SciPy':lambda:np.concatenate(cluster.vq.vq(v('obs'),v('codebook')))}
    elif api == 'cluster.vq.whiten': calls = {'SciPy':lambda:cluster.vq.whiten(v('obs'))}
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


ZERO=mp.mpf(0)


def component_error(actual, reference, *, binary64):
    # Equality of two binary64 values proves exact agreement with a binary64
    # oracle. Every unequal value and every high-precision oracle still uses mp.
    if binary64 and math.isfinite(actual) and actual==float(reference):return ZERO
    return symmetric_error(actual,mp.mpf(reference))


def digest_sources():
    names={'Cargo.toml','Cargo.lock','pyproject.toml','uv.lock','sci-revision.txt','api_inventory.csv','README.md'}
    names.update(p.relative_to(ROOT).as_posix() for pattern in ['*.py','sci_testing/**/*.py','tests/**/*.py','*.rs','src/*.rs','native/*.cpp','.github/workflows/*.yml','cases/**/*.csv'] for p in ROOT.glob(pattern))
    digest=hashlib.sha256()
    for name in sorted(names):digest.update(name.encode()+b'\0'+(ROOT/name).read_bytes()+b'\0')
    return digest.hexdigest()


def environment(sci_path):
    configs = io.StringIO()
    with redirect_stdout(configs):
        np.show_config(); scipy.show_config()
    return dict(sci_commit=command('git','rev-parse','HEAD',cwd=sci_path),
                sci_source_dirty=bool(command('git','status','--porcelain','--','src','crates','tests','Cargo.toml','Cargo.lock',cwd=sci_path)),
                testing_sources_sha256=digest_sources(),cases_sha256=directory_digest(ROOT/'cases/generated'),
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
        if hasattr(cases,'function'):
            names=[api for api in cases.apis if api.split('.')[0]==module]
            api_count,case_count=cases.module_counts(module)
            groups=((api,cases.function(api)) for api in names)
            available=set(meta['module_libraries'][module])
        else:
            selected=[c for c in cases if c['api'].split('.')[0]==module]
            case_count=len(selected);api_count=len({c['api'] for c in selected})
            grouped=defaultdict(list)
            for case in selected:grouped[case['api']].append(case)
            groups=grouped.items()
            available={lib for c in selected for lib in records[c['id']]['libraries']}
        available.discard('Sci')
        preferred=['SciPy','Boost.Math','GSL','NumPy','stdlib','FilterPy']
        libraries=sorted(available,key=lambda lib:(preferred.index(lib) if lib in preferred else len(preferred),lib))
        rows=[];failed={}
        for api,batch in groups:
            cells=[]
            baseline_failures=0
            for case in batch:
                record=records[case['id']]
                baseline=record['libraries']['Sci']
                if 'error' in baseline:
                    baseline_failures+=1
                    entry=failed.setdefault((api,'Sci'),[0,baseline['error'],case['id']]);entry[0]+=1
                elif len(baseline['values'])!=len(record['reference']):
                    raise ValueError(f'{case["id"]}: shape mismatch in Sci')
            for lib in libraries:
                sums=[];maximum=None;scalar_count=0;failures=0;unavailable=0;baseline_missing=0
                for case in batch:
                    record=records[case['id']];result=record['libraries'].get(lib)
                    if result is None:unavailable+=1;continue
                    if 'error' in result:
                        failures+=1
                        entry=failed.setdefault((api,lib),[0,result['error'],case['id']]);entry[0]+=1
                        continue
                    if len(result['values'])!=len(record['reference']):raise ValueError(f'{case["id"]}: shape mismatch in {lib}')
                    baseline=record['libraries']['Sci']
                    if 'error' in baseline:baseline_missing+=1;continue
                    skip=len(case['args']['obs']) if api=='cluster.vq.vq' else 0
                    errors=[component_error(actual,ref,binary64=True) for actual,ref in zip(result['values'][skip:],baseline['values'][skip:])]
                    if errors:
                        current=max(errors);maximum=current if maximum is None else max(maximum,current)
                        sums.append(mp.fsum(errors));scalar_count+=len(errors)
                if unavailable==len(batch):cells.append('—');continue
                if unavailable:raise ValueError('partially missing implementation')
                text=f'{number(maximum)} / {number(mp.fsum(sums)/scalar_count)}' if scalar_count else '—'
                cells.append(text+(f'; failed {failures}/{len(batch)}' if failures else '')+
                             (f'; Sci unavailable {baseline_missing}/{len(batch)}' if baseline_missing else ''))
            count=str(len(batch))+(f'; Sci failed {baseline_failures}/{len(batch)}' if baseline_failures else '')
            rows.append(f'| `{api.split(".",1)[1]}` | {count} | '+' | '.join(cells)+' |')
        failures=[f'- `{api}` ({lib}): {count} failed; {reason}. Example: `{identifier}`.' for (api,lib),(count,reason,identifier) in failed.items()]
        failure_section='\n## Evaluation failures\n\n'+'\n'.join(failures)+'\n' if failures else ''
        text = f'''# {module}: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-{module})

Measured binary64 (`f64`) results using **Sci-{module} as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers {case_count} cases and {api_count} public API names;
coverage is limited to the operations and inputs below.

## Agreement with Sci

Each library cell is **maximum / mean difference from Sci in ε** over cases
where both evaluations succeed. `Cases` counts input cases; `Sci failed N/T`
identifies missing baseline evaluations. `failed N/T` counts that library's own
failures; `Sci unavailable N/T` counts its successful evaluations excluded
because Sci failed. A dash means unavailable or no paired finite results.
Means include all compared scalar components, including exact zeros.
Full outputs, references, warnings and errors remain in the per-function CSVs.

| Function | Cases | {' | '.join(libraries)} |
| --- | ---: | {' | '.join('---' for _ in libraries)} |
'''+ '\n'.join(rows)+'\n'+failure_section+f'''

## Conventions and limitations

{NOTES[module]}

Inputs are generated deterministically into `sci-testing/cases/generated/<module>/<function>.csv`; raw outputs, high-precision
references and the environment are stored in that repository's `results/`.
NaN, infinity, signed-zero contracts, invalid domains and work-limit enforcement
are outside this finite-input accuracy sample. Existing Sci API tests cover
additional behavior. The [coverage audit](../../results/coverage.csv) checks every public function
against the compiler-resolved inventory. These samples do not establish
worst-case bounds or full domain coverage. SciPy aliases, compatibility namespaces and Cython
variants are not counted as independent implementations.

## Reproduction

- Sci source: `{meta['sci_commit']}` (numerical sources clean: {str(not meta['sci_source_dirty']).lower()}).
- Testing source SHA-256: `{meta['testing_sources_sha256']}`.
- NumPy {meta['numpy']}; SciPy {meta['scipy']}; FilterPy {meta['filterpy']}; mpmath {meta['mpmath']}.
- Native comparison versions: {meta.get('native_libraries',{})}.
- {meta['rust'].splitlines()[0]}; Python {meta['python']}; Cargo dev profile, opt-level 0.
- {meta['platform']}; {meta['processor']}.

See the [shared metric, provenance and reproduction instructions](README.md).
'''
        if 'reference_environment' in meta:
            text+='\nSci-only refresh: independent references and other library outputs are retained from the saved snapshot. Their original environment is recorded in `results/environment.json` under `reference_environment`.\n'
        if 'sci_measurement_commits' in meta:
            text+='\nSci results were refreshed by namespace. Per-function Sci revisions are recorded in `results/environment.json` under `sci_measurement_commits`.\n'
        (output/f'{module}.md').write_text(text)
    from .coverage_audit import inventory
    exports=inventory();canonical=len({row['canonical'] for row in exports})
    index=f"""# Sci library comparisons

Comparisons using each Sci module as the reference for all {canonical} public function implementations in
all 14 modules. The compiler inventory contains {len(exports)} exported names;
aliases share implementations. The [coverage audit](../../results/coverage.csv)
records sample counts, unique inputs, successful comparisons and failures.

Tables list **Function**, **Cases**, then relevant libraries (SciPy, Boost.Math,
GSL and others). Each library cell is **maximum / mean symmetric relative
difference from Sci in machine epsilon**, over cases where both evaluations
succeed. `failed N/T` counts that library's failures; `Sci failed N/T` in Cases
counts baseline failures, and `Sci unavailable N/T` marks otherwise successful
evaluations excluded from comparison. Unavailable comparisons use a dash. The per-function
CSV preserves every reference, result, warning and exception.

| Module | Recipe names | Cases |
| --- | ---: | ---: |
"""
    for module in MODULES:
        if hasattr(cases,'module_counts'):api_count,case_count=cases.module_counts(module)
        else:
            batch=[c for c in cases if c['api'].split('.')[0]==module];api_count=len({c['api'] for c in batch});case_count=len(batch)
        index+=f'| [{module}]({module}.md) | {api_count} | {case_count} |\n'
    index+=r"""

## Inputs and references

Each recipe has 1,000 endpoint-inclusive NumPy linspace cases and 1,000 uniform
random cases using seed `20261020` with an independent stream per API, plus
preserved boundary profiles. Broad original recipes complement local sweeps
around valid profiles. Integer-only functions, metadata queries and constrained
shapes use bounded discrete profiles, so some inputs repeat; the coverage CSV
reports unique inputs explicitly. Identical deterministic inputs share one
measured evaluation; all repeated case records remain in the CSV and totals.
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

Keep sibling `sci` and `sci-testing` checkouts, with Sci at `sci-revision.txt`.
Install a C++17 compiler, libclang, Boost.Math headers and GSL development files.
The Cargo and uv lockfiles pin Rust/Python dependencies separately from Sci.

```sh
cd ../sci-testing
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
per-function Sci revisions under `sci_measurement_commits`. Update `sci-revision.txt`
and, if necessary, `Cargo.lock` before measuring a newer Sci commit.
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

"""
    index+=f'- Sci commit: `{meta["sci_commit"]}`.\n- Harness source SHA-256: `{meta["testing_sources_sha256"]}`.\n- Input SHA-256: `{meta["cases_sha256"]}`.\n- {len(cases)} cases; {canonical} distinct public functions.\n- Native versions: {meta.get("native_libraries",{})}.\n'
    if 'reference_environment' in meta:
        index+='\nSci-only refresh: independent references and other library outputs are retained from the saved snapshot. See `reference_environment` in `results/environment.json` for their original provenance.\n'
    if 'sci_measurement_commits' in meta:
        index+='\nPer-function Sci revisions for partial refreshes are recorded under `sci_measurement_commits` in `results/environment.json`.\n'
    (output/'README.md').write_text(index)


def load_snapshot(cases, *, root=None):
    """Require identical inputs and intact measurements before reusing a snapshot."""
    root=ROOT if root is None else root
    expected_names=set()
    for name,content in iter_case_files(cases):
        path=root/'cases/generated'/name
        expected_names.add(name)
        if not path.exists() or path.read_text()!=content:raise ValueError('Stale case CSV files; rerun python -m sci_testing.report --full')
    actual_names={p.relative_to(root/'cases/generated').as_posix() for p in (root/'cases/generated').glob('*/*.csv')}
    if actual_names!=expected_names:raise ValueError('Stale case CSV files; rerun python -m sci_testing.report --full')
    meta=json.loads((root/'results/environment.json').read_text())
    if meta['cases_sha256'] != directory_digest(root/'cases/generated'):
        raise ValueError('Input digest mismatch; rerun python -m sci_testing.report --full')
    if meta.get('measurements_sha256') != directory_digest(root/'results/measurements'):
        raise ValueError('Raw measurement digest mismatch; rerun python -m sci_testing.report --full')
    records=MeasurementStore(root/'results/measurements')
    from itertools import groupby
    for api,iterator in groupby(cases,key=lambda c:c['api']):
        batch=list(iterator)
        records[batch[0]['id']]
        if set(records.records)!={c['id'] for c in batch}:raise ValueError('Measurement IDs differ from inputs; rerun python -m sci_testing.report --full')
    actual_names={p.relative_to(root/'results/measurements').as_posix().removesuffix('.gz') for p in csv_paths(root/'results/measurements')}
    if actual_names!=expected_names:raise ValueError('Measurement file set differs from inputs')
    return meta, records


def check_generated(cases, output, *, root=None):
    """Check derived inputs and pages against the source-pinned raw snapshot."""
    root=ROOT if root is None else root
    meta,records=load_snapshot(cases,root=root)
    if meta['testing_sources_sha256'] != digest_sources():
        raise ValueError('Testing sources changed; rerun python -m sci_testing.report')
    if 'coverage_sha256' in meta:
        from .coverage_audit import audit
        content=audit(cases,records)
        if (root/'results/coverage.csv').read_text()!=content or hashlib.sha256(content.encode()).hexdigest()!=meta['coverage_sha256']:raise ValueError('Stale coverage audit; rerun python -m sci_testing.report')
    with tempfile.TemporaryDirectory() as directory, mp.workdps(160):
        expected=Path(directory)
        render(cases,records,meta,expected)
        names={p.name for p in expected.glob('*.md')}
        if names != {p.name for p in output.glob('*.md')}:
            raise ValueError('Generated page set differs; rerun python -m sci_testing.report')
        for name in sorted(names):
            if (output/name).read_bytes() != (expected/name).read_bytes():
                raise ValueError(f'Stale accuracy page {name}; rerun python -m sci_testing.report')
    print(f'Checked {len(cases)} inputs, raw measurements and {len(names)} generated pages')


def measure_sci(cases, executable):
    """Evaluate distinct deterministic inputs once, preserving every case record."""
    unique={};keys=[]
    for case in cases:
        key=json.dumps([case['api'],case['args']],sort_keys=True,separators=(',',':'))
        if case['api'] in ['special.geterr','special.seterr']:key=case['id']
        keys.append(key)
        unique.setdefault(key,case)
    run=subprocess.run([str(executable)],input=''.join(json.dumps(c)+'\n' for c in unique.values()),
                       text=True,capture_output=True,check=True)
    rows=[json.loads(line) for line in run.stdout.splitlines()]
    if [r['id'] for r in rows]!=[c['id'] for c in unique.values()]:
        raise ValueError('runner lost, duplicated or reordered cases')
    measured=dict(zip(unique,rows))
    return [{**measured[key],'id':case['id']} for key,case in zip(keys,cases)]


def validate_result(case, result, reference):
    """Apply the same output checks to fresh and reused-reference measurements."""
    if 'values' in result:
        if len(result['values']) != len(reference):
            result['error']=f'shape mismatch: expected {len(reference)} scalar outputs, received {len(result["values"])}'
        elif case['api']=='cluster.vq.vq' and any(mp.mpf(x)!=mp.mpf(r) for x,r in zip(result['values'][:len(case['args']['obs'])],reference)):
            result['error']='incorrect code assignment'


def refresh_sci(cases, records, executable, destination, apis):
    """Replace only Sci evaluations; retain saved oracles and competitor results."""
    from itertools import groupby
    for api,iterator in groupby(cases,key=lambda c:c['api']):
        if api not in apis:continue
        batch=list(iterator)
        print(f'Measuring Sci only: {api}: {len(batch)} cases',flush=True)
        rust=measure_sci(batch,executable)
        updated={}
        for case,result in zip(batch,rust):
            saved=records[case['id']]
            if 'Sci' not in saved['libraries']:
                raise ValueError(f'{case["id"]}: saved Sci measurement missing')
            result={k:v for k,v in result.items() if k!='id'}
            if 'values' in result:result['values']=list(map(float,result['values']))
            with mp.workdps(160):validate_result(case,result,saved['reference'])
            updated[case['id']]={**saved,'libraries':{**saved['libraries'],'Sci':result}}
        write_measurements(destination,measurement_files(batch,updated))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'docs/accuracy')
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--check',action='store_true',help='validate inputs, raw results and generated pages without writing files')
    mode.add_argument('--full',action='store_true',help='recompute all Sci, reference and comparison-library results')
    mode.add_argument('--sci-only',action='store_true',help='refresh only Sci outputs (the default)')
    parser.add_argument('--module',action='append',default=[],metavar='NAMESPACE',help='refresh Sci only for this module or nested namespace; repeat to select several')
    args=parser.parse_args()
    if args.module and (args.full or args.check):
        parser.error('--module is available only for Sci-only refreshes')
    if (np.__version__,scipy.__version__,mp.__version__,filterpy.__version__) != ('2.5.4','1.18.1','1.3.0','1.4.5'):
        raise RuntimeError('Run uv sync --frozen with the pinned library versions')
    if not args.check and any(k.startswith(('RUSTFLAGS','CARGO_ENCODED_RUSTFLAGS','CARGO_PROFILE_')) or k == 'CARGO_BUILD_TARGET' for k in os.environ):
        raise RuntimeError('Unset custom compiler/profile/target overrides for this snapshot')
    cases=build_cases()
    if args.check:
        check_generated(cases,args.output.resolve())
        return
    selected=set(cases.apis)
    if args.module:
        selected=set()
        for module in args.module:
            matching={api for api in cases.apis if api.startswith(module+'.')}
            if not matching:parser.error(f'unknown or empty module: {module}')
            selected.update(matching)
    if not args.full:
        previous,records=load_snapshot(cases)
    else:
        write_cases(ROOT/'cases/generated', cases)
    if {api.split('.')[0] for api in cases.apis} != set(MODULES):
        raise ValueError('case identifiers must be unique and cover every module')
    pinned=(ROOT/'sci-revision.txt').read_text().strip()
    if command('git','rev-parse','HEAD',cwd=ROOT.parent/'sci') != pinned:
        raise RuntimeError('Sci checkout differs from sci-revision.txt; update the pin intentionally before measuring')
    meta=environment(ROOT.parent/'sci')
    if not args.full:
        meta['reference_environment']=previous.get('reference_environment',previous)
        if args.module:
            commits=previous.get('sci_measurement_commits',{api:previous['sci_commit'] for api in cases.apis}).copy()
            commits.update(dict.fromkeys(selected,meta['sci_commit']))
            meta['sci_measurement_commits']=commits
        for key in ['native_libraries','cxx','native_build_flags','module_libraries']:
            meta[key]=previous[key]
    else:
        native=NativeComparisons()
        meta['native_libraries']=native.versions
        meta['cxx']=command('c++','--version')
        meta['native_build_flags']='-std=c++17 -O2 -shared -fPIC -lgsl -lgslcblas'
    if meta['sci_source_dirty']: raise RuntimeError('Refusing to measure modified Sci numerical sources')
    subprocess.run(['cargo','build','--locked'],cwd=ROOT,check=True)
    target=Path(os.environ.get('CARGO_TARGET_DIR',ROOT/'target')).resolve()
    from itertools import groupby
    results=ROOT/'results'; results.mkdir(exist_ok=True)
    if not args.full:
        from .coverage_audit import audit
        with tempfile.TemporaryDirectory(dir=results) as directory:
            staged=Path(directory)
            shutil.copytree(results/'measurements',staged,dirs_exist_ok=True)
            refresh_sci(cases,records,target/'debug/sci-testing',staged,selected)
            audit(cases,MeasurementStore(staged))
            for path in csv_paths(staged):
                destination=results/'measurements'/path.relative_to(staged)
                destination.parent.mkdir(parents=True,exist_ok=True)
                path.replace(destination)
                alternate=destination.with_suffix('') if destination.suffix=='.gz' else destination.with_suffix('.csv.gz')
                alternate.unlink(missing_ok=True)
    else:
        files={};module_libraries={module:set() for module in MODULES}
        for api, iterator in groupby(cases, key=lambda c:c['api']):
            batch=list(iterator)
            identifiers=[c['id'] for c in batch]
            if len(identifiers)!=len(set(identifiers)) or any(not identifier.startswith(api+'/') for identifier in identifiers):raise ValueError('invalid or duplicate case ID')
            print(f'Measuring {api}: {len(batch)} cases', flush=True)
            rust=measure_sci(batch,target/'debug/sci-testing')
            for row in rust:
                if 'values' in row:row['values']=list(map(float,row['values']))
            if [r['id'] for r in rust] != [c['id'] for c in batch]: raise ValueError('runner lost, duplicated or reordered cases')
            records={}
            for case,result in zip(batch,rust):
                with mp.workdps(100): low=high_precision(case)
                with mp.workdps(160):
                    high=high_precision(case)
                    if len(low) != len(high): raise ValueError('unstable reference shape')
                    scale=max([TINY]+[abs(v) for v in high])
                    if any(abs(x-y)>mp.mpf('1e-90')*scale for x,y in zip(low,high)):
                        raise ValueError(f'unstable reference: {case["id"]}')
                    libs={'Sci':{k:v for k,v in result.items() if k!='id'},**competitors(case),**native.evaluate(case)}
                    for res in libs.values():validate_result(case,res,high)
                    records[case['id']]={'id':case['id'],'reference':[mp.nstr(x,165) for x in high],'libraries':libs}
                module_libraries[api.split('.')[0]].update(libs)
            output=measurement_files(batch, records)
            files.update(dict.fromkeys(write_measurements(results/'measurements',output)))
            from .registered_cases import cached_reference
            cached_reference.cache_clear()
            del batch, rust, records, output
            gc.collect()
        for path in csv_paths(results/'measurements'):
            if path.relative_to(results/'measurements').as_posix() not in files: path.unlink()
        meta['module_libraries']={module:sorted(names) for module,names in module_libraries.items()}
    meta['measurements_sha256']=directory_digest(results/'measurements')
    records=MeasurementStore(results/'measurements')
    from .coverage_audit import audit
    coverage=audit(cases,records)
    (results/'coverage.csv').write_text(coverage)
    meta['coverage_sha256']=hashlib.sha256(coverage.encode()).hexdigest()
    (results/'environment.json').write_text(json.dumps(meta,indent=2)+'\n')
    with mp.workdps(160): render(cases,records,meta,args.output.resolve())
    print(f'Wrote {len(MODULES)} module pages, {len(cases)} cases to {args.output}')


if __name__ == '__main__': main()
