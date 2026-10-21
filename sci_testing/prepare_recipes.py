"""Import valid input profiles from Sci's pinned, independently checked fixtures.

Only inputs are imported; expectations are evaluated independently at run time.
"""
import csv
import json
import math
import warnings
from .rust_source import formatted
from .api_registry import SCI, registered_functions, rust_registrations
from .case_csv import case_files, write_files

from . import ROOT


def prepare():
    source=[]
    for name in ['special','dependency_extras']:
        source.extend(json.loads((SCI/f'tests/scipy/{name}_cases.json').read_text()))
    manual={'harmonic':[5],'gen_harmonic':[5,2.],'ln_factorial':[5],'ln_binomial':[5,2],
            'seterr':[.5,1.], 'ln_multinomial':[5,[2,3]], 'multinomial':[5,[2,3]],
            'cbrt':[.5], 'logsumexp':[[1.,2.,3.]], 'softmax':[[1.,2.,3.]], 'log_softmax':[[1.,2.,3.]]}
    cases=[]
    for function in registered_functions():
        module,name=function.api.split('.')
        if name in ['erf','gamma'] and module=='special' or name in ['exp','sqrt'] and module=='complex': continue
        if module=='special':
            seeds=[case['args'] for case in source if case['name']==name and len(case['args'])==len(function.arguments)]
            if not seeds: seeds=[manual[name]]
        else:
            seeds=[]
            with (SCI/'tests/fixtures/complex_inputs.csv').open() as stream:
                for row in csv.DictReader(stream):
                    if row['api']!=name:continue
                    values=list(row.values())
                    z=[float(values[2]),float(values[3])];w=[float(values[4]),float(values[5])]
                    if name=='from_polar': args=z
                    elif name in ['add','subtract','multiply','divide','powc']:args=[z,w]
                    elif name=='powf':args=[z,w[0]]
                    elif name in ['powi','powu']:args=[z,int(w[0])]
                    else:args=[z]
                    seeds.append(args)
            if not seeds: seeds=[[[.5,.25],3]] if name=='powu' else [[[.5,.25]]]
        def finite(value):
            return all(map(finite,value)) if isinstance(value,list) else math.isfinite(value) if isinstance(value,(int,float)) else True
        seeds=[args for args in seeds if finite(args)]
        for i,args in enumerate(seeds):
            if len(args)!=len(function.arguments):raise ValueError((function.api,args,function.arguments))
            typed=[bool(v) if t=='bool' else float(v) if t=='f64' else v for v,(_,t) in zip(args,function.arguments)]
            cases.append(dict(id=f'{function.api}/recipe-{i:04d}',api=function.api,
                group='reference recipe; binary64 reference',args={f'p{j}':v for j,v in enumerate(typed)}))
    for module,filename in [('stats','statistics'),('signal','signal')]:
        for index,case in enumerate(json.loads((SCI/f'tests/scipy/{filename}_cases.json').read_text())):
            name=case['api'] if module=='stats' else case['function']
            if (module,name) in [('stats','normal_cdf'),('stats','variance'),('signal','convolve'),('signal','correlate')]:continue
            args={key:value for key,value in case.items() if key in ['p','x','y','g','d']} if module=='stats' else {'arguments':case['arguments']}
            try:json.dumps(args,allow_nan=False)
            except ValueError:continue
            cases.append(dict(id=f'{module}.{name}/recipe-{index:04d}',api=f'{module}.{name}',group='reference recipe; binary64 SciPy reference',args=args))
    from .registered_cases import materialize_case
    existing={c['api'] for c in cases}|{'fft.fft','fft.ifft','optimize.rosen','optimize.rosen_der','spatial.euclidean','spatial.cosine'}
    for suite,filename,namespace in [('interpolation','interpolation','interpolate'),('fft','fft','fft'),('fft_periodic','fft_periodic','fft'),('cluster','cluster','cluster.hierarchy'),('optimize','optimize','optimize')]:
        for index,raw in enumerate(json.loads((SCI/f'tests/scipy/{filename}_cases.json').read_text())):
            api=f'{namespace}.{raw["api"]}'
            if api in existing:continue
            args=materialize_case(suite,raw)
            cases.append(dict(id=f'{api}/recipe-{index:04d}',api=api,group='reference recipe; binary64 SciPy reference',args=args))
    with (SCI/'tests/fixtures/geometry_inputs.csv').open() as stream:
        for index,raw in enumerate(csv.DictReader(stream)):
            api=f'spatial.{raw["api"]}'
            if api in existing:continue
            cases.append(dict(id=f'{api}/recipe-{index:04d}',api=api,group='reference recipe; binary64 SciPy reference',args=materialize_case('geometry',raw)))
    from .registered_cases import linalg_case
    for index,raw in enumerate(json.loads((SCI/'tests/scipy/linalg_cases.json').read_text())):
        api=f'linalg.{raw["api"]}'
        if api in ['linalg.solve','linalg.det']:continue
        cases.append(dict(id=f'{api}/recipe-{index:04d}',api=api,group='reference recipe; binary64 SciPy reference',args=linalg_case(raw)))
    from .registered_cases import filter_case
    for index,raw in enumerate(json.loads((SCI/'tests/scipy/filter_cases.json').read_text())):
        api=f'filter.{raw["function"]}'
        if api in ['filter.kalman_predict','filter.kalman_update']:continue
        cases.append(dict(id=f'{api}/recipe-{index:04d}',api=api,group='reference recipe; binary64 SciPy reference',args=filter_case(raw)))
    from .sparse_cases import recipes as sparse_recipes
    cases.extend(sparse_recipes())
    # Resolve fixture dispatch labels to actual public Sci namespaces.
    aliases={
        'cluster.vq':'cluster.vq.vq','cluster.whiten':'cluster.vq.whiten',
        **{f'stats.{name}':f'stats.contingency.{name}' for name in ['association','crosstab','expected_freq','margins','odds_ratio','relative_risk']},
        **{f'stats.{name}':f'stats.qmc.{name}' for name in ['discrepancy','geometric_discrepancy','scale','update_discrepancy']},
    }
    for case in cases:
        original=case['api']
        if original in aliases:
            case['api']=aliases[original]
            case['id']=case['id'].replace(original+'/',case['api']+'/',1)
            if 'p' in case['args']:case['args']['dispatch_api']=original.split('.',1)[1]
    from .registered_cases import reference
    valid=[]
    for case in cases:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                values=reference(case)
            if all(map(math.isfinite,values)): valid.append(case)
        except (ValueError,ArithmeticError,TypeError): pass
    missing={c['api'] for c in cases}-{c['api'] for c in valid}
    if missing: raise ValueError(f'No finite reference profile: {sorted(missing)}')
    cases=valid
    write_files(ROOT/'cases/recipes',case_files(cases))
    (ROOT/'src/generated.rs').write_text(formatted(rust_registrations()))
    return cases


if __name__=='__main__':print(f'Prepared {len(prepare())} input profiles')
