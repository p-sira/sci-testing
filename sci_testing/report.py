"""Generate library comparison pages using measured Sci outputs as the baseline."""
import argparse
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
from .accuracy_cache import AccuracyCache, FILENAME as SUMMARY_FILENAME

from .generate_cases import build_cases
from .case_csv import (case_files, case_paths, case_directory_digest, csv_text,
                      disk_files, files_digest, measurement_files,
                      read_measurements, read_measurement_content, write_files, MeasurementStore, directory_digest, iter_case_files, write_cases,
                      csv_paths, write_measurements)

import filterpy
from filterpy.kalman import KalmanFilter
import mpmath as mp
import numpy as np
import scipy
from scipy import cluster, constants, differentiate, fft, integrate, interpolate
from scipy import linalg, optimize, signal, spatial, special, stats

from . import ROOT, SCI
EPS = mp.mpf(2) ** -52
TINY = mp.mpf(2) ** -1022
MODULES = ('cluster complex constants core differentiate fft filter integrate interpolate io '
           'linalg optimize signal sparse spatial special stats').split()
BOOST = 'https://www.boost.org/doc/libs/latest/libs/math/doc/html/math_toolkit/'
NOTES = {
    'sparse':'Canonical CSR and COO arrays retain explicit zeros. Matrix outputs include dimensions; graph distances retain infinities. Eigenpairs are compared by ordered values with independent residual checks; matching cardinality, permutation validity, flow conservation and random shape/count use invariants where outputs are nonunique. Sparse LU, iterative solves and matrix functions use pinned public SciPy references. Each scalar uses (class,value): class 0 for finite data, ±1 for signed infinity and 2 for NaN, with zero nonfinite payloads. This preserves unreachable vertices; exact tags contribute to component means. RNG streams differ. Dense samples preserve graph structure, integer capacities and dependent dimensions; fixed discrete profiles may repeat. Rust signature adaptations and bounded dense fallbacks are documented in Sci sparse coverage.',
    'cluster':'Vector quantization uses Euclidean distances and first-minimum ties; whitening divides by population standard deviations. Hierarchy fixtures also check counts, leaf ordering and partitions. K-means receives an explicit initial codebook. Returned discrete indices and tree metadata are recorded alongside distances.',
    'complex':'Complex outputs flatten real then imaginary components. Functions use principal branches; powers preserve their explicit real, complex or integer exponent conventions. All 34 complex APIs use high-precision mpmath mathematical references, with explicit signed-zero branch banks. NumPy provides the public binary64 complex operations; SciPy, Boost.Math and GSL do not expose equivalent general-purpose complex arithmetic. Signed-zero and nonfinite contracts remain in the Sci API tests.',
    'constants':'Vacuum conversions use metres and hertz and the exact SI value c = 299792458 m/s. Temperature recipes convert Celsius to Fahrenheit. CODATA values, uncertainties, unit strings and searches use the pinned SciPy constants table. Text metadata is encoded as UTF-8 lengths and bytes; Sci dtype f64 corresponds to NumPy float64.',
    'core':'CSR constructors, validators and conversion routines are exercised against independent SciPy sparse and explicit storage invariants. Canonical pointer arrays, index ordering, zeros and rejection behavior are retained in the per-case measurements.',
    'differentiate':'Callbacks are supplied explicitly: original sin and quadratic recipes, a quadratic Hessian and a two-output Jacobian. Fixed-step Sci stencils and adaptive SciPy differences have different work budgets. Analytic identities supply the callback references; raw inputs retain the selected steps and points.',
    'fft':'Transforms use unnormalized forward and inverse-divided-by-length conventions. Real transforms, multidimensional shapes, transform types, FFTLog and periodic operators follow the pinned fixture conventions. FFTLog uses SciPy fft; distinct periodic operators use its legacy fftpack implementation. Real and imaginary components are interleaved, including exact imaginary zeros.',
    'filter':'Gaussian outputs place means before row-major covariance; update and particle diagnostics follow the pinned filter fixture ordering. References use independent matrix identities, supplied deterministic particle draws and nonlinear callbacks. Original linear Kalman comparisons additionally call FilterPy.',
    'integrate':'Discrete-rule references integrate the supplied binary64 samples, excluding continuous sampling error. The dense added rules use dx = 0.25 and zero initial cumulative integrals. Newton-Cotes records weights followed by its error coefficient. Fixed Gaussian quadrature uses order 16 on exp(x), with analytic references.',
    'interpolate':'Original barycentric/Krogh recipes sample an exactly represented quartic. Added fitting and spline recipes compare evaluation, knots, coefficients and degree in the established fixture order. Queries, boundary modes, repeated knots and work sizes remain explicit in the CSV inputs. Equivalent spline representations can differ at roundoff scale.',
    'io':'Matrix Market, WAV and MATLAB comparisons use deterministic format fixtures and independent library references. Structured outputs follow the documented fixture ordering; raw input cases retain the shapes and format parameters.',
    'linalg':'Nonunique decompositions are compared through ordered eigenvalues, projectors and reconstructions, retaining independent structural invariants from the Sci fixture dispatcher. Interpolative decomposition compares spectral residual and rank; its SVD compares reconstruction and singular values. Type metadata uses the Sci f64 spelling. Dense sample sizes and constructor orders stay bounded.',
    'optimize':'Objective helpers use their pinned Rosenbrock and finite-difference conventions. Minimizers use 3+(x-1)^2+(y+2)^2, a 512-iteration budget and tolerance 1e-8; only converged points and costs are compared. Iteration trajectories are not compared. Analytic minimizers supply the reference; SciPy runs Nelder-Mead or BFGS separately.',
    'signal':'The pinned signal fixture adapters cover waveform, filtering, design, spectral, window, peak, spline, convolution and system functions through the public Sci facade. Output tuples flatten in the documented fixture order. Parameters, modes, shapes and normalization are retained in each function CSV; original real convolution/correlation also call NumPy.',
    'spatial':'Metrics, pairwise arrays, geometry and validation helpers use explicit weights, variance/covariance data and shapes from the pinned spatial fixtures. Output ordering follows those fixtures. Near-zero distances and cancellation can amplify relative error; the CSV retains full outputs for absolute-error inspection.',
    'special':'All public special-function implementations have typed registrations or an explicit policy recipe. Original erf/gamma recipes use high-precision mpmath; added recipes use independent public SciPy calls, identities or quadrature. Polynomial constructors return descending coefficients; sequence and tuple ordering follows the pinned fixtures. Boost.Math/GSL are called where the bridge has matching scalar APIs, converting elliptic parameter m to modulus sqrt(m). Public GSL provides complex digamma and complex Spence; its Spence on negative real inputs with negative imaginary zero selects the upper bank, which can disagree with the independent lower-bank reference.',
    'stats':'Distribution parameters, tail conventions, alternatives, ddof, explicit resampling draws and output ordering follow the pinned statistics fixtures. SciPy supplies binary64 references for added APIs, including documented adaptations and independently seeded probabilistic references. Weighted tau compares its statistic; its always-NaN p-value is left to the ordinary API tests. Normality, inverse and joint-tail approximations have reference and algorithm accuracy limits.',
}



def command(*args, cwd=ROOT):
    return subprocess.check_output(args, cwd=cwd, text=True).strip()


def binary64_reference(case):
    a=case['args']
    return any(key in a for key in ['p0','p','arguments','raw']) or bool(
        a.get('explicit') or a.get('io_case') or a.get('structural03_case'))


def high_precision(case):
    """References depend on mathematical definitions, never Sci or its backends."""
    api, a = case['api'], case['args']
    if api.startswith('complex.'):
        from .complex_accuracy import reference
        return reference(case)
    if a.get('special03_case'):
        from .special03_cases import reference
        return reference(case)
    if a.get('special_complex03_case'):
        from .special_complex03_cases import reference
        return reference(case)
    if a.get('structural03_case'):
        from .structural03_cases import reference
        return list(map(mp.mpf,reference(case)))
    if binary64_reference(case):
        from .registered_cases import reference
        return list(map(mp.mpf, reference(case)))
    def num(x): return mp.mpf(x)  # Exact conversion of the parsed binary64 input.
    def vec(key): return list(map(num, a[key]))
    def mat(key): return mp.matrix([[num(x) for x in r] for r in a[key]])
    def flat(m): return list(m)
    def gauss(mean, cov): return flat(mean) + flat(cov)
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


def reference_precisions(case):
    """Escalate precision around subnormal complex branch points as needed."""
    if not case['api'].startswith('complex.'):
        return ((100, 160),)
    a=case['args']
    components=[float(v) for value in a.values()
                for v in (value if isinstance(value,list) else [value])]
    # Low precision can agree on a false zero when a tiny term is cancelled
    # against a huge one (for example Im atan(1e308 i) ~ 1e-308).
    if any(abs(v)>1e150 or 0<abs(v)<1e-150 for v in components):
        return ((400, 600), (650, 1000))
    return ((200, 320), (400, 600), (650, 1000))


def stable_reference(case):
    for low_dps,high_dps in reference_precisions(case):
        with mp.workdps(low_dps):low=high_precision(case)
        with mp.workdps(high_dps):
            high=high_precision(case)
            scale=max([TINY]+[abs(value) for value in high])
            if len(low)==len(high) and all(abs(x-y)<=mp.mpf('1e-90')*scale
                                           for x,y in zip(low,high)):
                return high
    raise ValueError(f'Unstable independent reference: {case["id"]}')


def flatten(value):
    a = np.asarray(value)
    if np.iscomplexobj(a): return np.stack([a.real,a.imag],axis=-1).ravel().tolist()
    return a.astype(float).ravel().tolist()


def competitors(case):
    """Public library APIs, or explicitly labelled elementary NumPy expressions."""
    api, a = case['api'], case['args']
    if api.startswith('complex.'):
        from .complex_accuracy import numpy_comparison
        return numpy_comparison(case)
    if a.get('special_complex03_case'):
        from .special_complex03_cases import scipy_comparison
        return scipy_comparison(case)
    if a.get('special03_case'):
        # High-precision mathematical oracles are stored separately, without
        # presenting the same formula as an independent competitor library.
        return {}
    if a.get('structural03_case'):
        from .structural03_cases import reference
        return {'SciPy':{'values':list(map(float,reference(case)))}}
    if 'p0' in a or 'p' in a or 'arguments' in a or 'raw' in a or a.get('explicit') or a.get('io_case'):
        from .registered_cases import reference, reference_result
        library='NumPy' if a.get('explicit') and api.startswith('linalg.') else 'SciPy'
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
    if api.startswith('constants.'):
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
        calls = {'SciPy':lambda:getattr(special,name)(a['x'])}
    elif api == 'stats.normal_cdf': calls = {'SciPy':lambda:stats.norm.cdf(a['x'])}
    elif api == 'stats.variance':
        calls = {'SciPy':lambda:stats.describe(v('values'),ddof=1).variance,'NumPy':lambda:np.var(v('values'),ddof=1)}
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


def native_competitors(native, case):
    """Use only matching public Boost.Math/GSL interfaces and their true ABI."""
    if native is None:
        return {}
    if case['api'].startswith('special.complex.'):
        return native.evaluate_complex(case)
    return native.evaluate(case)


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
    names={'.gitmodules','Cargo.toml','Cargo.lock','pyproject.toml','uv.lock','sci-revision.txt','api_inventory.csv','README.md'}
    # Generated inputs already have their own logical-content digest. Keep the
    # source fingerprint focused on code and authored fixtures, including gzip.
    patterns=['*.py','sci_testing/**/*.py','tests/**/*.py','*.rs','src/*.rs',
              'native/*.cpp','.github/workflows/*.yml',
              'cases/base/**/*.csv','cases/base/**/*.csv.gz',
              'cases/recipes/**/*.csv','cases/recipes/**/*.csv.gz']
    names.update(p.relative_to(ROOT).as_posix() for pattern in patterns for p in ROOT.glob(pattern))
    digest=hashlib.sha256()
    for name in sorted(names):digest.update(name.encode()+b'\0'+(ROOT/name).read_bytes()+b'\0')
    return digest.hexdigest()


def environment(sci_path):
    configs = io.StringIO()
    with redirect_stdout(configs):
        np.show_config(); scipy.show_config()
    return dict(sci_commit=command('git','rev-parse','HEAD',cwd=sci_path),
                sci_source_dirty=bool(command('git','status','--porcelain','--','src','crates','tests','Cargo.toml','Cargo.lock',cwd=sci_path)),
                testing_sources_sha256=digest_sources(),cases_sha256=case_directory_digest(ROOT/'cases/generated'),
                rust=command('rustc','-Vv'),cargo=command('cargo','-V'),
                python=sys.version.split()[0],numpy=np.__version__,scipy=scipy.__version__,
                mpmath=mp.__version__,filterpy=filterpy.__version__,platform=platform.platform(),
                machine=platform.machine(),processor=next((l.split(':',1)[1].strip() for l in Path('/proc/cpuinfo').read_text().splitlines() if l.startswith('model name')),'unknown'),
                profile='dev; opt-level=0; debug=line-tables-only; no custom RUSTFLAGS',reference_decimal_digits=[100,160],
                library_config=configs.getvalue())


def number(x): return mp.nstr(x,5,min_fixed=-3,max_fixed=5)


def comparison_submodule(api, canonical):
    """Group historical aliases under their Sci 0.3 public namespace."""
    path=(canonical or api).split('.')
    return path[1] if len(path)>2 else None


def render(cases, records, meta, output, *, modules=None, cache=None):
    output.mkdir(parents=True,exist_ok=True)
    from .coverage_audit import inventory
    canonical={row['api']:row['canonical'] for row in inventory()}
    for module in (MODULES if modules is None else modules):
        if hasattr(cases,'function'):
            names=[api for api in cases.apis if api.split('.')[0]==module]
            api_count,case_count=cases.module_counts(module)
            groups=((api,cases.function(api)) for api in names)
            available=set(meta.get('module_libraries',{}).get(module,[]))
        else:
            selected=[c for c in cases if c['api'].split('.')[0]==module]
            case_count=len(selected);api_count=len({c['api'] for c in selected})
            grouped=defaultdict(list)
            for case in selected:grouped[case['api']].append(case)
            groups=grouped.items()
            available={lib for c in selected for lib in records[c['id']]['libraries']}
        available.discard('Sci')
        preferred=['SciPy','Boost.Math','GSL','NumPy','FilterPy']
        libraries=sorted(available,key=lambda lib:(preferred.index(lib) if lib in preferred else len(preferred),lib))
        sections=defaultdict(list);section_cases=defaultdict(int);failed={}
        for api,batch in groups:
            saved=cache.get(api,libraries) if cache is not None else None
            if saved is not None:
                for library,count,reason,identifier in saved['failures']:
                    failed[(api,library)]=[count,reason,identifier]
                section=comparison_submodule(api,canonical.get(api))
                section_cases[section]+=len(batch)
                sections[section].append(saved['row'])
                continue
            cells=[]
            baseline_failures=0
            reference_maximum=None;reference_sum=ZERO;reference_count=0
            for case in batch:
                record=records[case['id']]
                baseline=record['libraries']['Sci']
                if 'error' in baseline:
                    baseline_failures+=1
                    entry=failed.setdefault((api,'Sci'),[0,baseline['error'],case['id']]);entry[0]+=1
                elif len(baseline['values'])!=len(record['reference']):
                    raise ValueError(f'{case["id"]}: shape mismatch in Sci')
                else:
                    errors=[component_error(actual,reference,binary64=False)
                            for actual,reference in zip(baseline['values'],record['reference'])]
                    if errors:
                        current=max(errors)
                        reference_maximum=current if reference_maximum is None else max(reference_maximum,current)
                        reference_sum+=mp.fsum(errors)
                        reference_count+=len(errors)
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
                text=f'{number(maximum)} / {number(mp.fsum(sums)/scalar_count)}' if scalar_count else '—'
                cells.append(text+(f'; failed {failures}/{len(batch)}' if failures else '')+
                             (f'; unavailable {unavailable}/{len(batch)}' if unavailable else '')+
                             (f'; Sci unavailable {baseline_missing}/{len(batch)}' if baseline_missing else ''))
            count=str(len(batch))+(f'; Sci failed {baseline_failures}/{len(batch)}' if baseline_failures else '')
            reference_cell=(f'{number(reference_maximum)} / {number(reference_sum/reference_count)}'
                            if reference_count else '—')
            section=comparison_submodule(api,canonical.get(api))
            section_cases[section]+=len(batch)
            row=f'| `{api.split(".",1)[1]}` | {count} | {reference_cell} | '+' | '.join(cells)+' |'
            sections[section].append(row)
            if cache is not None:
                cache.put(api,libraries,{'row':row,'failures':[
                    [lib,*failure] for (name,lib),failure in failed.items() if name==api]})
        failures=[f'- `{api}` ({lib}): {count} failed; {reason}. Example: `{identifier}`.' for (api,lib),(count,reason,identifier) in failed.items()]
        failure_section='\n## Evaluation failures\n\n'+'\n'.join(failures)+'\n' if failures else ''
        headers=['Function','Cases','Sci vs reference',*libraries]
        table_header='| '+' | '.join(headers)+' |\n| '+' | '.join(['---','---:',*('---' for _ in headers[2:])])+' |\n'
        tables=[]
        for section in sorted(sections,key=lambda name:(name is None,name or '')):
            namespace=f'sci::{module}'+(f'::{section}' if section else '')
            title=f'### `{namespace}` ({len(sections[section])} functions, {section_cases[section]} cases)'
            tables.append(title+'\n\n'+table_header+'\n'.join(sections[section]))
        tables_text='\n\n'.join(tables)
        text = f'''# {module}: library comparison

[All modules and methodology](README.md) · [Public API](https://docs.rs/sci-{module})

Measured binary64 (`f64`) results using **Sci-{module} as the reference**.
Each library is compared with Sci on the same inputs. These differences measure
agreement with Sci; independent accuracy references remain in the raw CSVs. This page
covers {case_count} cases and {api_count} public API names;
coverage is limited to the operations and inputs below.

## Agreement with Sci

`Sci vs reference` is **maximum / mean difference between Sci and the saved
independent mathematical reference in ε**. Other library cells are **maximum /
mean difference from Sci in ε** over cases
where both evaluations succeed. `Cases` counts input cases; `Sci failed N/T`
identifies missing baseline evaluations. `failed N/T` counts that library's own
failures; `Sci unavailable N/T` counts its successful evaluations excluded
because Sci failed. `unavailable N/T` means the library has no compatible
public operation for those inputs or distribution families; they are excluded
from the reported mean and are not counted as numerical failures. A dash means
no saved comparison or no paired finite results.
Means include all compared scalar components, including exact zeros.
Full outputs, references, warnings and errors remain in the per-function CSVs.
Tables are grouped by Sci 0.3 public submodule using the compiler-resolved
canonical API inventory. Historical flat recipe aliases appear under their
current public namespace, with their measured recipe names preserved.

{tables_text}
'''+failure_section+f'''

## Conventions and limitations

{NOTES[module]}

Inputs are generated deterministically into `sci-testing/cases/generated/<module>/<function>.csv` (optionally `.gz`); raw outputs, high-precision
references and the environment are stored in that repository's `results/`.
NaN, infinity, signed-zero contracts, invalid domains and work-limit enforcement
are outside this finite-input accuracy sample. Existing Sci API tests cover
additional behavior. The [coverage audit](../../results/coverage.csv) checks every public function
against the compiler-resolved inventory. These samples do not establish
worst-case bounds or full domain coverage. SciPy aliases, compatibility namespaces and Cython
variants are not counted as independent implementations.

## Reproduction

The [shared methodology and environment](README.md) and
`results/environment.json` record source revisions, hashes, toolchain versions,
compiler flags and per-function measurement provenance.
'''
        (output/f'{module}.md').write_text(text)
    exports=inventory();canonical=len({row['canonical'] for row in exports if row['canonical']})
    sampled_apis=set(cases.apis) if hasattr(cases,'apis') else {case['api'] for case in cases}
    new_apis={row['api'] for row in exports if row.get('status')=='new_api_untested'}
    try:
        from .special_complex03_cases import APIS as complex_apis
        new_apis.update(complex_apis)
    except ImportError:
        pass  # The optional new special overload provider is not installed.
    pending_new=len(new_apis-sampled_apis)
    pending_other=sum(row.get('status')=='untested' and row['api'] not in sampled_apis
                      and row['api'] not in new_apis for row in exports)
    unresolved=sum(row.get('status')=='unresolved_legacy' for row in exports)
    index=f"""# Sci library comparisons

Historical and newly measured comparisons use each Sci module as the reference on
the same saved inputs. The compiler inventory lists {canonical} mapped canonical
names across {len(MODULES)} modules and {len(exports)} exported names; these are not
all independently verified against Sci 0.3. There are {pending_new} currently
unmeasured new API names, {pending_other} other untested exports and
{unresolved} unresolved legacy mappings.
The [coverage audit](../../results/coverage.csv) explicitly marks pending names
alongside measured sample counts, unique inputs, successes and failures.

Tables are grouped by Sci 0.3 public submodule (including historical aliases)
and list **Function**, **Cases**, **Sci vs reference**, then relevant libraries
(SciPy, Boost.Math, GSL and others). The reference column compares Sci with
the saved independent mathematical oracle; library columns compare each
competitor with Sci. Each cell is **maximum / mean symmetric relative
difference in machine epsilon**, over cases where both evaluations
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
module recipes use mpmath at 100 and 160 decimal digits, while complex facade
operations use 200/320 digits (with higher precision around subnormal branch
points); additional recipes use pinned public SciPy APIs, independent identities and quadrature. These
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
uv run --frozen python -m sci_testing.report --check --exhaustive
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
per-case identities, the authenticated coverage audit and cached per-function
numerical summaries, and every generated page without compiling or writing files.
`--check --exhaustive` independently recomputes coverage and every numerical
summary metric from saved measurements; it is also scheduled nightly in CI.
Content-addressed accuracy summaries live in `results/accuracy-summaries.json`.
They are indexed by input and measurement SHA-256, comparison columns, and
metric implementation fingerprint. A changed or missing summary is never silently
trusted; `--render` regenerates it from verified saved measurements without
rerunning the numerical libraries. Both check modes validate the recorded platform
snapshot; neither requires identical
floating-point measurements on another platform. Set `CARGO_TARGET_DIR` to reuse
an existing build directory. The runner uses Cargo's dev profile without custom
compiler flags. Numerical sources and fixture providers must be clean; Sci
documentation edits are permitted.

The GitHub workflow checks the committed snapshot and uploads the rendered reports,
coverage audit and provenance metadata (not duplicate copies of raw measurements).
An identical run on the same platform reproduces the snapshot without timestamps.
All reports and data belong to this repository; Sci links to them externally.

## Snapshot

"""
    index+=f'- Sci commit: `{meta["sci_commit"]}`.\n- Harness source SHA-256: `{meta["testing_sources_sha256"]}`.\n- Input SHA-256: `{meta["cases_sha256"]}`.\n- {len(cases)} saved cases covering {len(sampled_apis)} recipe API names; {canonical} mapped canonical public functions in the inventory.\n- Native versions: {meta.get("native_libraries",{})}.\n'
    if 'reference_environment' in meta:
        index+='\nSci-only refresh: independent references and other library outputs are retained from the saved snapshot. See `reference_environment` in `results/environment.json` for their original provenance.\n'
    if 'sci_measurement_commits' in meta:
        index+='\nPer-function Sci revisions for partial refreshes are recorded under `sci_measurement_commits` in `results/environment.json`.\n'
    if 'reference_extensions' in meta:
        index+='\nNew API references were evaluated separately from the historical snapshot. Their provenance is recorded in `reference_extensions` in `results/environment.json`.\n'
    (output/'README.md').write_text(index)


def render_cached(cases,records,meta,output,*,modules=None,root=None,inputs=None,
                  exhaustive=False,require_hits=False):
    """Stage the rendered pages and their authenticated summary index together."""
    root=ROOT if root is None else Path(root)
    prior=root/'results'/SUMMARY_FILENAME
    checksum=meta.get('accuracy_summaries_sha256')
    cache=AccuracyCache(root,records.directory,Path(__file__),inputs=inputs,
                        previous=prior if checksum else None,
                        previous_sha256=checksum,
                        reuse=not exhaustive,require_hits=require_hits)
    # Selective rendering can reuse unaffected API summaries only when the
    # stored metric implementation fingerprint still matches. If that
    # fingerprint changes, all previous summaries are untrusted and must be
    # rebuilt together before the full index can be authenticated.
    if modules is not None and set(cache.entries) != set(cases.apis):
        modules = None
    render(cases,records,meta,output,modules=modules,cache=cache)
    content=cache.serialize(cases.apis)
    output.mkdir(parents=True,exist_ok=True)
    (output/SUMMARY_FILENAME).write_bytes(content)
    meta['accuracy_summaries_sha256']=hashlib.sha256(content).hexdigest()


def publish_summary(pages,root=None):
    root=ROOT if root is None else Path(root)
    (pages/SUMMARY_FILENAME).replace(root/'results'/SUMMARY_FILENAME)


def load_snapshot(cases, *, root=None, cached=False, verify_apis=None):
    """Verify the complete saved digests and selected per-case identities.

    Comparator refreshes can check only their selected API IDs because neither
    the other inputs nor other Sci results change. Whole-directory digests and
    exact file inventories are checked regardless of this selection.
    """
    root=ROOT if root is None else root
    inputs=root/'cases/generated'
    paths=case_paths(inputs)
    expected_names={api.rsplit('.',1)[0]+'/'+api.rsplit('.',1)[1]+'.csv' for api in cases.apis} if cached else set()
    if not cached:
        for name,content in iter_case_files(cases):
            expected_names.add(name)
            if name not in paths or csv_text(paths[name])!=content:
                raise ValueError('Stale case CSV files; rerun python -m sci_testing.report --full')
    actual_names=set(paths)
    if actual_names!=expected_names:raise ValueError('Stale case CSV files; rerun python -m sci_testing.report --full')
    meta=json.loads((root/'results/environment.json').read_text())
    if meta['cases_sha256'] != case_directory_digest(inputs):
        raise ValueError('Input digest mismatch; rerun python -m sci_testing.report --full')
    if meta.get('measurements_sha256') != directory_digest(root/'results/measurements'):
        raise ValueError('Raw measurement digest mismatch; rerun python -m sci_testing.report --full')
    records=MeasurementStore(root/'results/measurements')
    if hasattr(cases,'function'):
        batches=((api,cases.function(api)) for api in
                 (cases.apis if verify_apis is None else sorted(verify_apis)))
    else:
        from itertools import groupby
        batches=((api,list(group)) for api,group in groupby(cases,key=lambda c:c['api'])
                 if verify_apis is None or api in verify_apis)
    for api,batch in batches:
        records[batch[0]['id']]
        if set(records.records)!={c['id'] for c in batch}:raise ValueError('Measurement IDs differ from inputs; rerun python -m sci_testing.report --full')
    actual_names={p.relative_to(root/'results/measurements').as_posix().removesuffix('.gz') for p in csv_paths(root/'results/measurements')}
    if actual_names!=expected_names:raise ValueError('Measurement file set differs from inputs')
    return meta, records


def load_stored_snapshot(*, root=None):
    """Verify the historical input and measurement snapshot without regeneration."""
    root=ROOT if root is None else root
    meta=json.loads((root/'results/environment.json').read_text())
    inputs=case_paths(root/'cases/generated')
    if meta['cases_sha256']!=case_directory_digest(root/'cases/generated'):
        raise ValueError('Input digest mismatch; cannot extend the saved references')
    if meta.get('measurements_sha256')!=directory_digest(root/'results/measurements'):
        raise ValueError('Raw measurement digest mismatch; cannot extend the saved references')
    measured={p.relative_to(root/'results/measurements').as_posix().removesuffix('.gz')
              for p in csv_paths(root/'results/measurements')}
    if measured!=set(inputs):
        raise ValueError('Measurement files differ from saved inputs')
    return meta,MeasurementStore(root/'results/measurements'),inputs


def preserved_coverage(previous, *, cases, records, root=None):
    """Reuse a signed-off Sci coverage audit when Sci results and inputs are unchanged.

    Call only after load_snapshot or load_stored_snapshot has verified the saved
    input and measurement digests. Comparator-only/reference-only transformations
    retain all Sci results and original inputs, so their coverage is invariant.
    """
    root=ROOT if root is None else root
    path=root/'results/coverage.csv'
    if not path.exists() or 'coverage_sha256' not in previous:
        from .coverage_audit import audit
        return audit(cases,records)
    content=path.read_text()
    if hashlib.sha256(content.encode()).hexdigest()!=previous['coverage_sha256']:
        raise ValueError('Coverage audit digest mismatch; rerun a full validation')
    return content


def selective_render_modules(output, dirty_modules):
    """Create missing pages on the first run; otherwise render only changed modules.

    Module pages contain only module-local data; environment metadata lives in
    the overview, which render() always regenerates.
    """
    existing={path.name for path in output.glob('*.md')}
    expected={'README.md',*(f'{module}.md' for module in MODULES)}
    return None if existing!=expected else sorted(dirty_modules)


def extend_references(cases, added, executable, destination):
    """Evaluate only new cases; preserve all older reference strings and results."""
    native=NativeComparisons() if any(api.startswith(('special.','stats.'))
                                      for api in added) else None
    libraries={}
    for api,batch in added.items():
        print(f'Measuring new inputs: {api}: {len(batch)} cases',flush=True)
        rust=measure_sci(batch,executable)
        if [result.get('id') for result in rust]!=[case['id'] for case in batch]:
            raise ValueError(f'{api}: runner lost, duplicated or reordered cases')
        module,function=api.rsplit('.',1)
        updated={}
        for case,result in zip(batch,rust):
            if 'values' not in result and 'error' not in result:
                raise ValueError(f'{case["id"]}: runner returned neither values nor an error')
            result={k:v for k,v in result.items() if k!='id'}
            if 'values' in result:result['values']=list(map(float,result['values']))
            high=stable_reference(case)
            with mp.workdps(160):
                libs={'Sci':result,**competitors(case),**native_competitors(native,case)}
                for value in libs.values():validate_result(case,value,high)
                updated[case['id']]={'id':case['id'],'reference':[mp.nstr(x,165) for x in high],
                                     'libraries':libs}
                libraries.setdefault(api.split('.')[0],set()).update(libs)
        complete=cases.function(api)
        filename=f'{module}/{function}.csv'
        for path in [destination/filename,destination/(filename+'.gz')]:path.unlink(missing_ok=True)
        write_measurements(destination,measurement_files(complete,updated))
    return libraries


def check_generated(cases, output, *, root=None, cached=False, exhaustive=False):
    """Check derived inputs and pages against the source-pinned raw snapshot."""
    root=ROOT if root is None else root
    meta,records=load_snapshot(cases,root=root,cached=cached)
    if meta['testing_sources_sha256'] != digest_sources():
        raise ValueError('Testing sources changed; rerun python -m sci_testing.report')
    if 'coverage_sha256' in meta:
        saved=(root/'results/coverage.csv').read_text()
        if hashlib.sha256(saved.encode()).hexdigest()!=meta['coverage_sha256']:
            raise ValueError('Coverage audit digest mismatch')
        if exhaustive or 'accuracy_summaries_sha256' not in meta:
            from .coverage_audit import audit
            if audit(cases,records)!=saved:
                raise ValueError('Stale coverage audit; rerun python -m sci_testing.report')
    with tempfile.TemporaryDirectory() as directory, mp.workdps(160):
        expected=Path(directory)
        original_checksum=meta.get('accuracy_summaries_sha256')
        if original_checksum:
            render_cached(cases,records,meta,expected,root=root,
                          exhaustive=exhaustive,require_hits=not exhaustive)
            if meta['accuracy_summaries_sha256']!=original_checksum or (
                    expected/SUMMARY_FILENAME).read_bytes()!=(root/'results'/SUMMARY_FILENAME).read_bytes():
                raise ValueError('Stale accuracy summaries; rerun python -m sci_testing.report --render')
        else:
            render(cases,records,meta,expected)
        names={p.name for p in expected.glob('*.md')}
        if names != {p.name for p in output.glob('*.md')}:
            raise ValueError('Generated page set differs; rerun python -m sci_testing.report')
        for name in sorted(names):
            if (output/name).read_bytes() != (expected/name).read_bytes():
                raise ValueError(f'Stale accuracy page {name}; rerun python -m sci_testing.report')
    print(f'Checked {len(cases)} inputs, raw measurements and {len(names)} generated pages'+
          (' exhaustively' if exhaustive else ''))


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
    changed=set()
    for api,iterator in groupby(cases,key=lambda c:c['api']):
        if api not in apis:continue
        batch=list(iterator)
        print(f'Measuring Sci only: {api}: {len(batch)} cases',flush=True)
        rust=measure_sci(batch,executable)
        if [row.get('id') for row in rust] != [case['id'] for case in batch]:
            raise ValueError(f'{api}: runner lost, duplicated or reordered cases')
        updated={}
        for case,result in zip(batch,rust):
            saved=records[case['id']]
            if 'Sci' not in saved['libraries']:
                raise ValueError(f'{case["id"]}: saved Sci measurement missing')
            result={k:v for k,v in result.items() if k!='id'}
            if 'error' not in result and 'values' not in result:
                raise ValueError(f'{case["id"]}: runner returned neither values nor an error')
            if 'values' in result:result['values']=list(map(float,result['values']))
            with mp.workdps(160):validate_result(case,result,saved['reference'])
            updated[case['id']]={**saved,'libraries':{**saved['libraries'],'Sci':result}}
        # The staging tree may share hard links with the saved snapshot. Break
        # those links before write_measurements opens a file for truncation.
        files=measurement_files(batch,updated)
        for name in files:
            (destination/name).unlink(missing_ok=True)
            (destination/(name+'.gz')).unlink(missing_ok=True)
        changed.update(write_measurements(destination,files))
    return changed


def refresh_comparators(cases, records, destination, apis, *, native=None):
    """Add measured competitors without evaluating Sci or rewriting stored oracles."""
    from itertools import groupby
    changed=set();libraries=defaultdict(set)
    for api,iterator in groupby(cases,key=lambda c:c['api']):
        if api not in apis:continue
        batch=list(iterator)
        print(f'Measuring comparison libraries only: {api}: {len(batch)} cases',flush=True)
        updated={}
        for case in batch:
            saved=records[case['id']]
            new={**competitors(case),**native_competitors(native,case)}
            for library,result in new.items():
                if library=='Sci':raise ValueError('Competitor cannot replace Sci measurement')
                validate_result(case,result,saved['reference'])
                libraries[api.split('.')[0]].add(library)
            updated[case['id']]={**saved,'libraries':{**{lib:value for lib,value in saved['libraries'].items()
                                                      if lib != 'stdlib'},**new}}
        files=measurement_files(batch,updated)
        # A staged tree can share hard links with the source snapshot.
        for name in files:
            (destination/name).unlink(missing_ok=True)
            (destination/(name+'.gz')).unlink(missing_ok=True)
        changed.update(write_measurements(destination,files))
    return changed,libraries


def upgrade_scientific_comparisons(cases, records, destination):
    """Replace cmath-derived complex oracles and remove legacy stdlib comparisons.

    Existing Sci outputs, cases, and non-stdlib library results remain pinned.
    The caller stages files with hard links, so always unlink before writing.
    """
    changed=set()
    for api in cases.apis:
        if not api.startswith('complex.'):
            continue
        batch=cases.function(api)
        print(f'Upgrading independent complex references: {api}: {len(batch)} cases',flush=True)
        updated={}
        for case in batch:
            saved=records[case['id']]
            high=stable_reference(case)
            with mp.workdps(160):
                libraries={lib:value for lib,value in saved['libraries'].items()
                           if lib not in ('stdlib','NumPy')}
                libraries.update(competitors(case))
                for result in libraries.values(): validate_result(case,result,high)
                updated[case['id']]={'id':case['id'],
                                    'reference':[mp.nstr(x,165) for x in high],
                                    'libraries':libraries}
        files=measurement_files(batch,updated)
        for name in files:
            (destination/name).unlink(missing_ok=True)
            (destination/(name+'.gz')).unlink(missing_ok=True)
        changed.update(write_measurements(destination,files))

    # Other modules already use SciPy and Boost.Math/GSL where available.
    # Only obsolete stdlib rows need to be stripped from their saved reports.
    for path in csv_paths(ROOT/'results/measurements'):
        if path.parent.name not in ('special','stats'):
            continue
        content=csv_text(path)
        if ',stdlib,' not in content:
            continue
        batch=read_measurement_content(content)
        for record in batch.values():record['libraries'].pop('stdlib',None)
        api=next(iter(batch)).split('/')[0]
        files=measurement_files([{'id':id,'api':api} for id in batch],batch)
        for name in files:
            (destination/name).unlink(missing_ok=True)
            (destination/(name+'.gz')).unlink(missing_ok=True)
        changed.update(write_measurements(destination,files))
        print(f'Removed Python stdlib comparison: {api}',flush=True)
    return changed


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'docs/accuracy')
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--check',action='store_true',help='validate inputs, raw results and generated pages without writing files')
    parser.add_argument('--exhaustive',action='store_true',help='with --check, independently recompute coverage and all numerical summary metrics')
    mode.add_argument('--render',action='store_true',help='regenerate coverage and accuracy pages from verified saved measurements, without running Sci or reference libraries')
    mode.add_argument('--extend',action='store_true',help='measure only new API names and retain every historical independent reference')
    mode.add_argument('--full',action='store_true',help='recompute all Sci, reference and comparison-library results')
    mode.add_argument('--sci-only',action='store_true',help='refresh only Sci outputs (the default)')
    mode.add_argument('--comparators-only',action='store_true',help='backfill comparison libraries from verified saved inputs, preserving Sci and independent reference results')
    mode.add_argument('--upgrade-scientific-comparisons',action='store_true',
                      help='replace cmath complex references with mpmath/NumPy and remove all stdlib comparisons, preserving Sci outputs')
    parser.add_argument('--module',action='append',default=[],metavar='NAMESPACE',help='refresh Sci only for this module or nested namespace; repeat to select several')
    args=parser.parse_args()
    if args.exhaustive and not args.check:
        parser.error('--exhaustive requires --check')
    if args.module and (args.full or args.check or args.extend or args.upgrade_scientific_comparisons):
        parser.error('--module is available only for Sci-only, comparator and render refreshes')
    if args.comparators_only and not args.module:
        parser.error('--comparators-only requires explicit --module selection')
    if (np.__version__,scipy.__version__,mp.__version__,filterpy.__version__) != ('2.5.4','1.18.1','1.3.0','1.4.5'):
        raise RuntimeError('Run uv sync --frozen with the pinned library versions')
    if not args.check and any(k.startswith(('RUSTFLAGS','CARGO_ENCODED_RUSTFLAGS','CARGO_PROFILE_')) or k == 'CARGO_BUILD_TARGET' for k in os.environ):
        raise RuntimeError('Unset custom compiler/profile/target overrides for this snapshot')
    cached=False
    if args.full:
        cases=build_cases()
    elif args.extend:
        from .case_set import historical_cases
        from .generate_cases import extend_historical_cases
        previous,records,_=load_stored_snapshot()
        historic=historical_cases(ROOT/'cases/generated',
                                  expected_cases_sha256=previous['cases_sha256'],
                                  measurements_directory=ROOT/'results/measurements',
                                  expected_measurements_sha256=previous['measurements_sha256'])
        cases=extend_historical_cases(historic)
        added=cases.added
        if not added:raise ValueError('No new public API cases to extend')
    else:
        from .case_set import cached_cases, historical_cases
        from .generate_cases import case_signature
        cases=cached_cases(ROOT/'cases/generated',case_signature())
        if cases is None and (ROOT/'results/environment.json').is_file():
            previous=json.loads((ROOT/'results/environment.json').read_text())
            cases=historical_cases(ROOT/'cases/generated',
                                   expected_cases_sha256=previous['cases_sha256'],
                                   measurements_directory=ROOT/'results/measurements',
                                   expected_measurements_sha256=previous['measurements_sha256'])
        cached=cases is not None
        if cases is None:cases=build_cases()
    if args.check:
        check_generated(cases,args.output.resolve(),cached=cached,exhaustive=args.exhaustive)
        return
    if args.render:
        render_apis=set(cases.apis)
        if args.module:
            render_apis=set()
            for module in args.module:
                matching={api for api in cases.apis if api==module or api.startswith(module+'.')}
                if not matching:parser.error(f'unknown or empty module: {module}')
                render_apis.update(matching)
        previous,records=load_snapshot(cases,cached=cached,
                                       verify_apis=render_apis if args.module else None)
        meta=previous.copy()
        if meta['testing_sources_sha256'] != digest_sources():
            meta['reference_environment']=meta.get('reference_environment',previous)
            meta['testing_sources_sha256']=digest_sources()
        coverage=preserved_coverage(previous,cases=cases,records=records)
        meta['coverage_sha256']=hashlib.sha256(coverage.encode()).hexdigest()
        with tempfile.TemporaryDirectory(dir=ROOT/'results') as directory:
            pages=Path(directory)
            modules=(selective_render_modules(args.output.resolve(),
                                              {api.split('.')[0] for api in render_apis})
                     if args.module else None)
            with mp.workdps(160):render_cached(cases,records,meta,pages,modules=modules)
            (ROOT/'results/coverage.csv').write_text(coverage)
            publish_summary(pages)
            (ROOT/'results/environment.json').write_text(json.dumps(meta,indent=2)+'\n')
            args.output.resolve().mkdir(parents=True,exist_ok=True)
            for page in pages.glob('*.md'):
                page.replace(args.output.resolve()/page.name)
        print(f'Rendered {len(MODULES)} module pages, {len(cases)} saved cases to {args.output}')
        return
    if args.upgrade_scientific_comparisons:
        previous,records=load_snapshot(cases,cached=cached)
        meta=previous.copy()
        meta['reference_environment']=meta.get('reference_environment',previous)
        meta['testing_sources_sha256']=digest_sources()
        meta['module_libraries']={module:sorted((set(libraries)-{'stdlib'})|
                                 ({'NumPy'} if module=='complex' else set()))
                                  for module,libraries in previous['module_libraries'].items()}
        meta.setdefault('reference_refreshes',[]).append({
            'apis':[api for api in cases.apis if api.startswith('complex.')],
            'reference':'mpmath, 200/320 decimal digits; escalated to 400/600 or 650/1000 near subnormal branch points; signed-zero banks',
            'replaced':'binary64 Python cmath and stdlib comparisons',
            'testing_sources_sha256':meta['testing_sources_sha256']})
        with tempfile.TemporaryDirectory(dir=ROOT/'results') as directory:
            workspace=Path(directory)
            staged=workspace/'measurements'
            shutil.copytree(ROOT/'results/measurements',staged,copy_function=os.link)
            changed=upgrade_scientific_comparisons(cases,records,staged)
            staged_records=MeasurementStore(staged)
            coverage=preserved_coverage(previous,cases=cases,records=staged_records)
            meta['measurements_sha256']=directory_digest(staged)
            meta['coverage_sha256']=hashlib.sha256(coverage.encode()).hexdigest()
            pages=workspace/'pages'
            with mp.workdps(160):
                render_cached(cases,staged_records,meta,pages,
                       modules=selective_render_modules(args.output.resolve(),
                                                        {'complex','special','stats'}))
            for name in sorted(changed):
                src=staged/name;dest=ROOT/'results/measurements'/name
                src.replace(dest)
                alternate=dest.with_suffix('') if dest.suffix=='.gz' else dest.with_suffix('.csv.gz')
                alternate.unlink(missing_ok=True)
            (ROOT/'results/coverage.csv').write_text(coverage)
            publish_summary(pages)
            (ROOT/'results/environment.json').write_text(json.dumps(meta,indent=2)+'\n')
            args.output.resolve().mkdir(parents=True,exist_ok=True)
            for page in pages.glob('*.md'):page.replace(args.output.resolve()/page.name)
        print(f'Upgraded {len(changed)} saved function reports without rerunning Sci')
        return
    selected=set(cases.apis)
    if args.module:
        selected=set()
        for module in args.module:
            # A full public API name selects just that function; namespace
            # prefixes retain their existing multi-function selection behavior.
            matching={api for api in cases.apis if api == module or api.startswith(module+'.')}
            if not matching:parser.error(f'unknown or empty module: {module}')
            selected.update(matching)
    if not (args.full or args.extend):
        previous,records=load_snapshot(cases,cached=cached,
                                       verify_apis=selected if args.comparators_only else None)
    elif args.full:
        from .generate_cases import case_signature
        write_cases(ROOT/'cases/generated',cases,signature=case_signature())
    if {api.split('.')[0] for api in cases.apis} != set(MODULES):
        raise ValueError('case identifiers must be unique and cover every module')
    if args.comparators_only:
        native=NativeComparisons() if any(api.startswith(('special.','stats.'))
                                          for api in selected) else None
        meta=previous.copy()
        meta['module_libraries']={module:list(names) for module,names in previous['module_libraries'].items()}
        if meta['testing_sources_sha256']!=digest_sources():
            meta['reference_environment']=meta.get('reference_environment',previous)
            meta['testing_sources_sha256']=digest_sources()
        if native:
            meta['native_libraries']=native.versions
            meta['cxx']=command(os.environ.get('CXX','c++'),'--version')
            meta['native_build_flags']='-std=c++17 -O2 -shared -fPIC -lgsl -lgslcblas'
        with tempfile.TemporaryDirectory(dir=ROOT/'results') as directory:
            workspace=Path(directory)
            staged=workspace/'measurements'
            shutil.copytree(ROOT/'results/measurements',staged,copy_function=os.link)
            changed,libraries=refresh_comparators(cases,records,staged,selected,native=native)
            staged_records=MeasurementStore(staged)
            coverage=preserved_coverage(previous,cases=cases,records=staged_records)
            meta['measurements_sha256']=directory_digest(staged)
            meta['coverage_sha256']=hashlib.sha256(coverage.encode()).hexdigest()
            for module,names in libraries.items():
                meta['module_libraries'][module]=sorted(set(meta['module_libraries'][module])|names)
            meta.setdefault('comparison_refreshes',[]).append({
                'apis':sorted(selected),'testing_sources_sha256':meta['testing_sources_sha256'],
                'scipy':meta['scipy'],'native_libraries':meta.get('native_libraries',{})})
            pages=workspace/'pages'
            with mp.workdps(160):
                render_cached(cases,staged_records,meta,pages,
                       modules=selective_render_modules(args.output.resolve(),
                                                        {api.split('.')[0] for api in selected}))
            for name in sorted(changed):
                src=staged/name;dest=ROOT/'results/measurements'/name
                dest.parent.mkdir(parents=True,exist_ok=True)
                src.replace(dest)
                alternate=dest.with_suffix('') if dest.suffix=='.gz' else dest.with_suffix('.csv.gz')
                alternate.unlink(missing_ok=True)
            (ROOT/'results/coverage.csv').write_text(coverage)
            publish_summary(pages)
            (ROOT/'results/environment.json').write_text(json.dumps(meta,indent=2)+'\n')
            args.output.resolve().mkdir(parents=True,exist_ok=True)
            for page in pages.glob('*.md'):page.replace(args.output.resolve()/page.name)
        print(f'Backfilled {len(selected)} APIs from saved inputs and references')
        return
    pinned=(ROOT/'sci-revision.txt').read_text().strip()
    if command('git','rev-parse','HEAD',cwd=SCI) != pinned:
        raise RuntimeError('Sci checkout differs from sci-revision.txt; update the pin intentionally before measuring')
    meta=environment(SCI)
    if not args.full:
        meta['reference_environment']=previous.get('reference_environment',previous)
        # Append-only reference extensions document independently evaluated
        # new APIs. A Sci-only refresh must not erase that provenance.
        if 'reference_extensions' in previous:
            meta['reference_extensions']=previous['reference_extensions']
        if 'comparison_refreshes' in previous:
            meta['comparison_refreshes']=previous['comparison_refreshes']
        if 'reference_refreshes' in previous:
            meta['reference_refreshes']=previous['reference_refreshes']
        if args.module or args.extend:
            commits=previous.get('sci_measurement_commits',{api:previous['sci_commit'] for api in cases.apis}).copy()
            if args.extend:
                commits={api:commits.get(api,previous['sci_commit']) for api in cases.apis}
                commits.update(dict.fromkeys(added,meta['sci_commit']))
            else:commits.update(dict.fromkeys(selected,meta['sci_commit']))
            meta['sci_measurement_commits']=commits
        for key in ['native_libraries','cxx','native_build_flags','module_libraries']:
            meta[key]=previous[key]
    else:
        native=NativeComparisons()
        meta['native_libraries']=native.versions
        meta['cxx']=command(os.environ.get('CXX','c++'),'--version')
        meta['native_build_flags']='-std=c++17 -O2 -shared -fPIC -lgsl -lgslcblas'
    if meta['sci_source_dirty']: raise RuntimeError('Refusing to measure modified Sci numerical sources')
    subprocess.run(['cargo','build','--locked'],cwd=ROOT,check=True)
    # Cargo may inherit a workspace target-dir from .cargo/config.toml even
    # when CARGO_TARGET_DIR is unset. Resolve the actual artifact directory.
    cargo_metadata=json.loads(subprocess.check_output(
        ['cargo','metadata','--format-version','1','--no-deps'],cwd=ROOT,text=True))
    target=Path(cargo_metadata['target_directory'])
    from itertools import groupby
    results=ROOT/'results'; results.mkdir(exist_ok=True)
    if args.extend:
        from .coverage_audit import audit
        from .generate_cases import case_signature
        from .case_csv import stage_case_extension, finalize_extended_cache
        new_cases=[case for api in sorted(added) for case in added[api]]
        with tempfile.TemporaryDirectory(dir=results) as directory:
            workspace=Path(directory)
            staged_inputs=workspace/'cases'
            stage_case_extension(historic,new_cases,staged_inputs)
            staged=workspace/'measurements'
            shutil.copytree(results/'measurements',staged,copy_function=os.link)
            new_libraries=extend_references(cases,added,target/'debug/sci-testing',staged)
            meta['module_libraries']={module:list(names) for module,names in previous['module_libraries'].items()}
            for module,names in new_libraries.items():
                meta['module_libraries'].setdefault(module,[])
                meta['module_libraries'][module]=sorted(set(meta['module_libraries'][module])|names)
            meta['cases_sha256']=finalize_extended_cache(
                staged_inputs,historic,new_cases,signature=case_signature(),
                historical_cases_sha256=previous['cases_sha256'],
                historical_measurements_directory=results/'measurements',
                historical_measurements_sha256=previous['measurements_sha256'],
                measurements_directory=staged)
            meta['measurements_sha256']=directory_digest(staged)
            extensions=list(previous.get('reference_extensions',[]))
            extensions.append(dict(apis=sorted(added),sci_commit=meta['sci_commit'],
                                   testing_sources_sha256=meta['testing_sources_sha256'],
                                   cases_sha256=meta['cases_sha256'],
                                   numpy=meta['numpy'],scipy=meta['scipy'],filterpy=meta['filterpy'],
                                   mpmath=meta['mpmath'],python=meta['python'],rust=meta['rust']))
            meta['reference_extensions']=extensions
            staged_records=MeasurementStore(staged)
            coverage=audit(cases,staged_records,required_apis=added)
            meta['coverage_sha256']=hashlib.sha256(coverage.encode()).hexdigest()
            pages=workspace/'pages'
            with mp.workdps(160):render_cached(cases,staged_records,meta,pages,inputs=staged_inputs)
            # All existing inputs were checked before preparing these files.
            old_inputs=ROOT/'cases/generated'
            backup=workspace/'previous-cases'
            old_inputs.replace(backup)
            try:staged_inputs.replace(old_inputs)
            except Exception:
                backup.replace(old_inputs)
                raise
            for api in added:
                module,function=api.rsplit('.',1)
                for suffix in ('.csv','.csv.gz'):
                    path=staged/module/(function+suffix)
                    if path.exists():
                        destination=results/'measurements'/module/(function+suffix)
                        destination.parent.mkdir(parents=True,exist_ok=True)
                        path.replace(destination)
                        other=destination.with_suffix('') if destination.suffix=='.gz' else destination.with_suffix('.csv.gz')
                        other.unlink(missing_ok=True)
            (results/'coverage.csv').write_text(coverage)
            publish_summary(pages)
            (results/'environment.json').write_text(json.dumps(meta,indent=2)+'\n')
            args.output.resolve().mkdir(parents=True,exist_ok=True)
            for page in pages.glob('*.md'):page.replace(args.output.resolve()/page.name)
        print(f'Extended {len(added)} functions, {sum(map(len,added.values()))} new cases')
        return
    if not args.full:
        from .coverage_audit import audit
        with tempfile.TemporaryDirectory(dir=results) as directory:
            staged=Path(directory)/'measurements'
            # Hard links are safe because refresh_sci removes selected links
            # before writing. Reusing unchanged files avoids copying the whole
            # (potentially multi-gigabyte) reference snapshot on each refresh.
            shutil.copytree(results/'measurements',staged,copy_function=os.link)
            changed=refresh_sci(cases,records,target/'debug/sci-testing',staged,selected)
            staged_records=MeasurementStore(staged)
            coverage=audit(cases,staged_records)
            meta['measurements_sha256']=directory_digest(staged)
            meta['coverage_sha256']=hashlib.sha256(coverage.encode()).hexdigest()
            pages=Path(directory)/'pages'
            with mp.workdps(160):render_cached(cases,staged_records,meta,pages)
            # Finish all audits and page generation before publishing changes.
            for name in sorted(changed):
                path=staged/name
                destination=results/'measurements'/name
                destination.parent.mkdir(parents=True,exist_ok=True)
                path.replace(destination)
                alternate=destination.with_suffix('') if destination.suffix=='.gz' else destination.with_suffix('.csv.gz')
                alternate.unlink(missing_ok=True)
            (results/'coverage.csv').write_text(coverage)
            publish_summary(pages)
            (results/'environment.json').write_text(json.dumps(meta,indent=2)+'\n')
            args.output.resolve().mkdir(parents=True,exist_ok=True)
            for page in pages.glob('*.md'):
                page.replace(args.output.resolve()/page.name)
        print(f'Wrote {len(MODULES)} module pages, {len(cases)} cases to {args.output}')
        return
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
                high=stable_reference(case)
                with mp.workdps(160):
                    native_results=native_competitors(native,case)
                    libs={'Sci':{k:v for k,v in result.items() if k!='id'},**competitors(case),**native_results}
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
    with mp.workdps(160): render_cached(cases,records,meta,args.output.resolve())
    publish_summary(args.output.resolve())
    (results/'environment.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(f'Wrote {len(MODULES)} module pages, {len(cases)} cases to {args.output}')


if __name__ == '__main__': main()
