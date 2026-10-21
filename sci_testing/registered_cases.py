"""Decorator-based reference and recipe registration for additional API families."""
from functools import lru_cache
import hashlib
import math
import json
import ast
import types
import warnings
import random

import numpy as np
from scipy import special

from .api_registry import SCI, registered_functions
from .case_csv import read_cases

from . import ROOT
REFERENCES = {}


def reference_family(namespace):
    def register(function):
        REFERENCES[namespace] = function
        return function
    return register


@lru_cache(maxsize=1)
def scipy_special_reference():
    return fixture_definitions('special')


@reference_family('special')
def special_reference(case):
    name = case['api'].split('.')[1]
    args = [case['args'][f'p{i}'] for i in range(len(case['args']))]
    function = specs()[case['api']]
    args = [int(value) if kind in ['usize','isize','u32','i32','u64','i64'] else
            complex(*value) if kind == 'Complex' else value
            for value, (_, kind) in zip(args, function.arguments)]
    if name == 'harmonic': value = float(special.digamma(args[0]+1)+np.euler_gamma) if args[0] else 0.
    elif name == 'gen_harmonic': value = math.fsum(j**(-args[1]) for j in range(1,args[0]+1))
    elif name == 'ln_factorial': value = special.gammaln(args[0]+1)
    elif name == 'ln_binomial': value = special.gammaln(args[0]+1)-special.gammaln(args[1]+1)-special.gammaln(args[0]-args[1]+1)
    elif name == 'seterr': value = args[0]
    elif name in ['multinomial','ln_multinomial']:
        n,counts=args
        value=special.gammaln(n+1)-math.fsum(special.gammaln(c+1) for c in counts)
        if name=='multinomial':value=math.exp(value)
    elif name in ['softmax','log_softmax','logsumexp']:
        from scipy import special as sp
        value=getattr(sp,name)(args[0])
    else: value = scipy_special_reference().reference(name, args)
    if hasattr(value, 'coeffs'): value = value.coeffs
    value = scipy_special_reference().flatten(value)
    if 'Complex' in function.returns:
        value = np.column_stack([value.real, value.imag]).ravel()
    elif np.iscomplexobj(value):
        if np.any(value.imag != 0): raise ValueError(f'{name}: nonreal reference')
        value = value.real
    return np.asarray(value, dtype=float).tolist()


@reference_family('complex')
def complex_reference(case):
    import cmath
    name=case['api'].split('.')[1]
    args=[case['args'][f'p{i}'] for i in range(len(case['args']))]
    if name=='from_polar': value=cmath.rect(*args)
    else:
        z=complex(*args[0])
        w=args[1] if len(args)>1 else None
        if name=='add':value=z+complex(*w)
        elif name=='subtract':value=z-complex(*w)
        elif name=='multiply':value=z*complex(*w)
        elif name=='divide':value=z/complex(*w)
        elif name=='negate':value=-z
        elif name=='conjugate':value=z.conjugate()
        elif name=='reciprocal':value=1/z
        elif name=='real':return [z.real]
        elif name=='imag':return [z.imag]
        elif name=='abs':return [abs(z)]
        elif name=='norm_sqr':return [z.real*z.real+z.imag*z.imag]
        elif name=='arg':return [cmath.phase(z)]
        elif name=='to_polar':return list(cmath.polar(z))
        elif name in ['powf','powi','powu']:value=z**w
        elif name=='powc':value=z**complex(*w)
        elif name=='cbrt':value=cmath.rect(abs(z)**(1/3),cmath.phase(z)/3)
        else:value=getattr(cmath, {'ln':'log'}.get(name,name))(z)
    return [value.real,value.imag]


@lru_cache(maxsize=1)
def specs():
    return {f.api:f for f in registered_functions()}


@lru_cache(maxsize=8192)
def cached_reference(api, arguments):
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter('always')
        values=REFERENCES[api.split('.')[0]]({'api':api,'args':json.loads(arguments)})
    return {'values':values,'warnings':[str(w.message) for w in captured]}


def reference(case):
    if case['args'].get('explicit'):
        from .explicit_cases import reference as explicit_reference
        return explicit_reference(case)
    return cached_reference(case['api'], json.dumps(case['args'],sort_keys=True))['values']


def scaled(value, factor):
    return [scaled(x,factor) for x in value] if isinstance(value,list) else value*factor


def has_nonzero(value):
    if isinstance(value,dict):return any(has_nonzero(x) for x in value.values())
    if isinstance(value,list):return any(has_nonzero(x) for x in value)
    return isinstance(value,(int,float)) and not isinstance(value,bool) and value!=0


@lru_cache(maxsize=512)
def density_distribution(name, parameters):
    return fixture_definitions('statistics').distribution(name,parameters)


def expanded_registered(count, seed, apis=None):
    """Local sweeps around valid fixture profiles; integer work stays bounded.

    Continuous arguments vary within ±0.1% of one valid profile. This narrow
    family complements the broader module recipes; it is not a domain-wide
    accuracy claim. Discrete-only functions repeat bounded profiles explicitly.
    """
    from .case_set import profiles as profile_map
    profiles=profile_map()
    definitions = specs()
    grid = np.linspace(-1., 1., count)
    for api, seeds in sorted(profiles.items()):
        if apis is not None and api not in apis:continue
        yield from seeds
        preferred=seeds[0]
        for candidate in seeds:
            arguments=candidate['args'];raw=arguments.get('raw',{})
            values=raw.get('args',raw.get('x',arguments.get('arguments',arguments.get('x',[]))))
            if has_nonzero(values):preferred=candidate;break
        if api=='stats.truncate':
            preferred=next(c for c in seeds if c['args']['p'][0]<c['args']['x'][0]<c['args']['p'][1])
        if api.startswith('complex.'):
            for candidate in seeds:
                values=candidate['args']
                z=values.get('p0')
                if isinstance(z,list) and any(z) and values.get('p1',1)!=0:
                    preferred=candidate;break
        if api=='special.gen_harmonic':
            preferred=next((c for c in seeds if c['args']['p0']>1 and .5<c['args']['p1']<4),preferred)
        rng = random.Random(int.from_bytes(hashlib.sha256(f'{seed}:{api}'.encode()).digest(), 'big'))
        kinds = dict((f'p{i}',kind) for i,(_,kind) in enumerate(definitions[api].arguments)) if api in definitions else {}
        for method in ['linspace','random']:
            for i in range(count):
                # A fixed valid seed keeps parameter dependencies and dimensions intact.
                selected=seeds[i%len(seeds)] if api in ['stats.pdf','stats.logpdf','stats.cdf','stats.sf','stats.ppf'] else preferred
                args = json.loads(json.dumps(selected['args']))
                varied=False
                t = float(grid[i]) if method == 'linspace' else rng.uniform(-1,1)
                if 'raw' in args:
                    raw=args['raw']
                    if api.startswith('sparse.'):
                        from .sparse_cases import perturb
                        args=perturb(api,args['suite_name'],raw,t,i)
                        raw=args['raw']
                    elif api.startswith('filter.'):
                        def shift(value):return [shift(x) for x in value] if isinstance(value,list) else value+.001*t
                        raw['args'][0]=shift(raw['args'][0])
                    elif api.startswith('linalg.'):
                        name=raw['api']
                        if name in ['hilbert','invhilbert','pascal','invpascal','helmert','dft']:
                            raw['args'][0]=1+i%16 if method=='linspace' else rng.randint(1,16)
                        elif name=='hadamard':raw['args'][0]=1 << (i%5 if method=='linspace' else rng.randrange(5))
                        elif name=='cossin':
                            # A Givens rotation preserves the required orthogonal input.
                            matrix=raw['args'][0]
                            first,second=matrix[:2]
                            c,s=math.cos(.001*t),math.sin(.001*t)
                            matrix[0]=[c*x-s*y for x,y in zip(first,second)]
                            matrix[1]=[s*x+c*y for x,y in zip(first,second)]
                        elif name not in ['qr_update','qr_insert','qr_delete']:
                            for j,arg in enumerate(raw['args']):
                                if isinstance(arg,list) and arg and isinstance(arg[0],list) and all(isinstance(x,(int,float)) for row in arg for x in row):
                                    raw['args'][j]=[[x*(1+.001*t) for x in row] for row in arg]
                                    break
                            else:
                                for j,arg in enumerate(raw['args']):
                                    if isinstance(arg,list) and arg and all(isinstance(x,(int,float)) for x in arg):
                                        raw['args'][j]=scaled(arg,1+.001*t);break
                        else:
                            raw['args'][1]=scaled(raw['args'][1],1+.001*t)
                    elif api.startswith('fft.') and raw.get('re'):
                        raw['re']=[x*(1+.001*t) for x in raw['re']]
                    elif api.startswith('interpolate.'):
                        # Keep knots, degrees and fit constraints from the valid profile.
                        # Query-only variation is encoded in the imported profile cases.
                        query_indices={'make_interp_spline':3,'make_lsq_spline':5,'make_splrep':5,'splrep':5,'make_smoothing_spline':4,'splev':0,'splder':2,'splantider':2,'insert':3,'spalde':0}
                        name=raw['api']
                        j=query_indices.get(name)
                        if j is not None and isinstance(raw['args'][j],list):
                            raw['args'][j]=[x+.001*t for x in raw['args'][j]]
                            if name=='spalde':
                                s=raw['args'][1];degree=s['degree'];left=s['knots'][degree];right=s['knots'][-degree-1]
                                raw['args'][j]=[min(right,max(left,x)) for x in raw['args'][j]]
                        else:
                            index={'bisplrep':2,'generate_knots':1,'griddata':1,'interpn':1,'make_splprep':0,'splprep':0}.get(name)
                            if index is not None:raw['args'][index]=scaled(raw['args'][index],1+.001*t)
                            elif name in ['bisplev','splint','sproot']:
                                spline=raw['args'][0 if name=='sproot' else 2]
                                spline['coefficients']=scaled(spline['coefficients'],1+.001*t)
                            elif name=='approximate_taylor_polynomial':raw['args'][1]+=.001*t
                    elif api.startswith('optimize.') and raw.get('x'):
                        raw['x']=[x*(1+.001*t) for x in raw['x']]
                    elif api.startswith('cluster.'):
                        if raw.get('y'):raw['y']=scaled(raw['y'],1+.001*t)
                        elif raw.get('x'):raw['x']=scaled(raw['x'],1+.001*t)
                        if raw.get('z') and raw['api']!='is_valid_im':
                            for row in raw['z']:row[2]*=1+.001*t
                    elif api.startswith('spatial.'):
                        boolean={'dice','hamming','jaccard','rogerstanimoto','russellrao','sokalsneath','yule'}
                        name=raw['api']
                        key='w' if name=='geometric_slerp' else 'v' if name=='tsearch' else 'u'
                        if name not in boolean and raw.get('param') not in boolean and raw.get(key):
                            raw[key]=':'.join(';'.join(str(float(x)*(1+.001*t)) for x in row.split(';')) for row in raw[key].split(':'))
                    args=materialize_case(args['suite_name'],raw)
                elif api.startswith('stats.') and args.get('x'):
                    if api in ['stats.fisher_exact','stats.barnard_exact','stats.boschloo_exact']:
                        maximum=min(16,64-int(sum(args['x'][1:])))
                        args['x'][0]=1+min(maximum-1,int((t+1)*maximum/2))
                    elif api in ['stats.multinomial_pmf','stats.multinomial_logpmf']:
                        args['y'][0]*=1+.001*t
                        total=math.fsum(args['y']);args['y']=[x/total for x in args['y']]
                    elif api=='stats.pointbiserialr':
                        args['y']=[x*(1+.001*t) for x in args['y']]
                    else:
                        args['x'] = [x*(1+.001*t) if isinstance(x,(int,float)) else x for x in args['x']]
                    if api in ['stats.pdf','stats.logpdf']:
                        # Quantiles preserve support, including integer PMF arguments.
                        dist=density_distribution(args['d'],tuple(args['p']))
                        # The p=1 geometric law is concentrated at one; SciPy's
                        # legacy inverse formula divides by log(0) in this case.
                        quantile=1. if args['d']=='Geometric' and args['p']==[1] else float(dist.ppf(.5+.49*t))
                        args['x']=[quantile for _ in args['x']]
                    if api=='stats.ppf':args['x']=[min(.999,max(.001,x)) for x in args['x']]
                    if api=='stats.truncate':
                        left,right=args['p'][:2];margin=.001*(right-left)
                        args['x']=[min(right-margin,max(left+margin,x+.001*t*(right-left))) for x in selected['args']['x']]
                    if api in ['stats.qmc.update_discrepancy','stats.false_discovery_control']:
                        args['x']=[min(1.,max(0.,x)) if isinstance(x,(int,float)) else x for x in args['x']]
                    if api in ['stats.chisquare','stats.power_divergence'] and args.get('y'):
                        args['y']=[x*(1+.001*t) for x in args['y']]
                    if api in ['stats.dirichlet_pdf','stats.dirichlet_logpdf']:
                        args['x'][0]*=1+.001*t
                        total=math.fsum(args['x']);args['x']=[x/total for x in args['x']]
                    if api.startswith('stats.contingency.'):
                        args['x']=list(seeds[0]['args']['x'])
                elif api=='stats.dunnett':
                    args['g'][0]=[x+.001*t for x in args['g'][0]]
                elif api=='stats.multiscale_graphcorr':
                    n=int(args['p'][0])
                    for j in range(n):args['g'][j]=[x+.001*t for x in args['g'][j]]
                elif api.startswith('stats.'):
                    if args.get('g') and not api.startswith('stats.qmc.'):
                        args['g'][0]=scaled(args['g'][0],1+.001*t)
                    else:
                        for j,value in enumerate(args.get('p',[])):
                            if isinstance(value,float) and value!=0 and not value.is_integer():
                                args['p'][j]*=1+.001*t;break
                elif api.startswith('signal.'):
                    values=args['arguments']
                    value=values[0] if values else None
                    if api.startswith('signal.windows.') and values and isinstance(value,int):
                        minimum=2*math.ceil(values[1])+1 if api=='signal.windows.dpss' else 2
                        length=minimum+min(30,int((t+1)*15))
                        values[0]=2*length if api=='signal.windows.kaiser_bessel_derived' else length
                    elif api in ['signal.freqz_sos','signal.sosfreqz','signal.sos2tf','signal.sos2zpk','signal.sosfilt','signal.sosfiltfilt','signal.sosfilt_zi']:
                        for row in value:row[:3]=scaled(row[:3],1+.001*t)
                    elif isinstance(value,list) and len(value)==3 and isinstance(value[-1],(int,float)) and all(isinstance(x,list) for x in value[:2]):
                        value[-1]*=1+.001*t
                    elif api!='signal.firwin_2d' and isinstance(value,list) and value and all(isinstance(x,(int,float)) and not isinstance(x,bool) for x in value):
                        values[0]=[x*(1+.001*t) for x in value]
                    elif isinstance(value,list) and value and all(isinstance(row,list) and all(isinstance(x,(int,float)) for x in row) for row in value):
                        values[0]=scaled(value,1+.001*t)
                    else:
                        for j,value in enumerate(values[1:],1):
                            if isinstance(value,float) and value!=0:values[j]*=1+.001*t;break
                            if isinstance(value,list) and value and all(isinstance(x,float) for x in value):values[j]=scaled(value,1+.001*t);break
                for key, kind in (kinds.items() if api in ['special.ellip_harm','special.ellip_harm_2','special.ellip_normal'] else reversed(kinds.items())):
                    if kind == 'f64':
                        value = args[key]
                        if value != 0 and math.isfinite(value):
                            args[key] = value*(1+.001*t)
                            varied=True
                            break
                    elif kind == '&[f64]':
                        args[key] = [x+.001*t*(j+1) for j,x in enumerate(args[key])]
                        varied=True
                        break
                    elif kind == 'Complex':
                        re, im = args[key]
                        args[key] = [re+.001*t, im+.0005*t]
                        varied=True
                        break
                if not varied and kinds:
                    integer_keys=[k for k,kind in kinds.items() if kind in ['usize','isize','u32','i32','u64','i64']]
                    if len(integer_keys)==1 and all(kind in ['usize','isize','u32','i32','u64','i64','bool'] for kind in kinds.values()):
                        args[integer_keys[0]]=1+i%16 if method=='linspace' else rng.randint(1,16)
                    if api=='special.ln_binomial':args['p0']=int(args['p1'])+1+(i if method=='linspace' else rng.randrange(count))
                    if api in ['special.multinomial','special.ln_multinomial']:
                        first=1+(i if method=='linspace' else rng.randrange(count));args['p1']=[first,2];args['p0']=first+2
                    if api in ['special.factorialk','special.jnyn_zeros','special.stirling2']:
                        n=1+min(15,int((t+1)*8))
                        args['p0']=n+(args['p1'] if api=='special.stirling2' else 0)
                if api=='special.xlogy':args={'p0':1+.001*t,'p1':1+.002*t}
                if args==selected['args']:
                    # Sweep valid discrete/structural profiles when no numerical
                    # perturbation applies, including empty and metadata inputs.
                    index=min(len(seeds)-1,int((t+1)*len(seeds)/2))
                    args=json.loads(json.dumps(seeds[index]['args']))
                yield dict(id=f'{api}/{method}-{i:04d}', api=api,
                           group=f'{method}, local fixture sweep; binary64 SciPy reference', args=args)


@lru_cache(maxsize=4)
def fixture_definitions(module):
    """Load only definitions, excluding the upstream generator's file-writing body."""
    path=SCI/f'tests/scipy/generate_{module}.py'
    tree=ast.parse(path.read_text(),filename=str(path))
    definitions=[]
    for node in tree.body:
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='cases' for t in node.targets):break
        definitions.append(node)
    namespace={'__file__':str(path),'__name__':f'fixture_reference_{module}'}
    exec(compile(ast.Module(body=definitions,type_ignores=[]),str(path),'exec'),namespace)
    return types.SimpleNamespace(**namespace)


@reference_family('stats')
def statistics_reference(case):
    module=fixture_definitions('statistics')
    values=np.asarray(module.reference({'api':case['args'].get('dispatch_api',case['api'].split('.',1)[1]),**case['args']}),dtype=float).ravel().tolist()
    return values[:1] if case['api']=='stats.weightedtau' else values


@reference_family('signal')
def signal_reference(case):
    module=fixture_definitions('signal')
    name=case['api'].split('.',1)[1]
    rows=module.output_rows(module.reference(name,case['args']['arguments']))
    return [value for row in rows for value in row]


@reference_family('interpolate')
@reference_family('fft')
@reference_family('cluster')
@reference_family('spatial')
@reference_family('differentiate')
@reference_family('optimize')
def extracted_reference(case):
    from .fixture_program import evaluate
    return evaluate(case['args']['suite_name'],case['args']['raw'])[1]


def materialize_case(suite, raw):
    from .fixture_program import input_fields
    if suite in ('sparse','sparse_parity'):
        from .sparse_cases import materialize
        return materialize(suite,raw)
    if suite=='linalg':return linalg_case(raw)
    if suite=='filter':return filter_case(raw)
    fields=input_fields(suite,raw)
    adapter={'fft_periodic':'periodic','geometry':'spatial'}.get(suite,suite)
    return dict(fields=fields, raw=raw, suite=adapter, suite_name=suite)


@reference_family('linalg')
def linalg_reference(case):
    raw=case['args']['raw']
    module=fixture_definitions('linalg')
    rows=module.evaluate(raw['api'],raw['args'])
    if raw['api']=='find_best_blas_type':rows[1]='f64'
    if raw['api']=='interpolative.svd': rows=rows[:2]
    if raw['api']=='interpolative.interp_decomp':
        from scipy import linalg
        rows=[linalg.norm(np.asarray(raw['args'][0])-rows[0],2),rows[1]]
    values=[]
    for row in rows:
        if isinstance(row,str):values.extend([len(row.encode()),*row.encode()])
        else:values.extend(np.asarray(row,dtype=float).ravel().tolist())
    return values


def linalg_case(raw):
    module=fixture_definitions('linalg')
    fields=['0',raw['api'],*[module.encode(a) for a in raw['args']],*['']*max(0,5-len(raw['args']))]
    return dict(fields=fields,raw=raw,suite='linalg',suite_name='linalg')


@reference_family('filter')
def filter_reference(case):
    module=fixture_definitions('filter')
    raw=case['args']['raw']
    result=module.reference(raw['function'],raw['args'])
    return [float(x) for part in result for x in np.asarray(part).ravel()]


def filter_case(raw):
    module=fixture_definitions('filter')
    fields=[raw['id'],raw['function'],*[module.encode(v) for v in raw['args']]]
    return dict(fields=fields,raw=raw,suite='filter',suite_name='filter')


def reference_result(case):
    return cached_reference(case['api'],json.dumps(case['args'],sort_keys=True))


@reference_family('sparse')
def sparse_reference(case):
    from .sparse_cases import reference
    return reference(case)
