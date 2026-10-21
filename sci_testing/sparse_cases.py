"""Sparse input recipes and independent pinned SciPy reference adapters.

The fixture generators supply public SciPy calls without their file writers.
Outputs follow sorted observable keys. Matrix results include their dimensions;
nonunique permutations compare size plus Rust-side permutation invariants.
Random samples compare shape and stored count; streams differ across libraries.
"""
import ast
import copy
from functools import lru_cache
import io
import json
import numpy as np
from .api_registry import SCI

GRAPH = {
    'shortest_path','dijkstra','bellman_ford','floyd_warshall','johnson',
    'breadth_first_order','depth_first_order','breadth_first_tree','depth_first_tree',
    'connected_components','minimum_spanning_tree','reverse_cuthill_mckee',
    'laplacian','maximum_bipartite_matching','structural_rank',
    'min_weight_full_bipartite_matching','maximum_flow','yen',
    'csgraph_from_dense','csgraph_masked_from_dense','csgraph_from_masked',
    'csgraph_to_dense','csgraph_to_masked','reconstruct_path','construct_dist_matrix',
}
LINALG = {
    'spsolve','splu','factorized','spsolve_triangular','norm','inv','matrix_power',
    'spbandwidth','is_sptriangular','cg','bicg','bicgstab','cgs','qmr','tfqmr',
    'gmres','lgmres','gcrotmk','minres','lsqr','lsmr','eigs','eigsh','lobpcg',
    'svds','expm','expm_multiply','funm_multiply_krylov','onenormest','spilu',
    'aslinearoperator','norm_axis',
}
MATRICES = {
    **{name:{'d'} for name in ('shortest_path','dijkstra','bellman_ford','floyd_warshall','johnson')},
    'breadth_first_tree':{'a'},'depth_first_tree':{'a'},'laplacian':{'a'},
    'csgraph_from_dense':{'a'},'csgraph_masked_from_dense':{'a'},
    'csgraph_from_masked':{'a'},'csgraph_to_dense':{'a'},'csgraph_to_masked':{'a'},
    'reconstruct_path':{'a'},'construct_dist_matrix':{'d'},
    'expm':{'a'},'expm_multiply':{'a'},'spdiags':{'a'},
    'load_npz':{'a'},'save_npz':{'a'},
}


def api_name(name):
    if name=='spectral_norm':name='norm'
    return 'sparse.'+('csgraph.' if name in GRAPH else 'linalg.' if name in LINALG else '')+name


def enc(values):
    return ';'.join(str(v) for v in np.asarray(values).ravel())


def rows(values):
    return ':'.join(enc(row) for row in values)


def materialize(suite, raw):
    # Match the fixture input encodings without computing any expected values.
    name=raw['api']
    if suite=='sparse':
        fields=[raw['id'],name,str(raw['rows']),str(raw['cols']),rows(raw.get('a',[])),
            str(raw.get('brows',0)),str(raw.get('bcols',0)),rows(raw.get('b',[])),
            rows(raw.get('x',[])) if name in ('diags','diags_array') else enc(raw.get('x',[])),enc(raw.get('k',0))]
    else:
        matrix=np.asarray(raw['a'],dtype=float)
        triples=[[r,c,matrix[r,c]] for r,c in zip(*np.nonzero(matrix))]
        if raw.get('zero_edge'):triples.append([0,1,0.])
        fields=[raw['id'],name,*map(str,matrix.shape),rows(triples),enc(raw['p']),
            rows(raw['x']) if name=='spdiags' else enc(raw['x']),str(int(raw.get('zero_edge',False)))]
    if 'maxiter' in raw:
        fields.extend(map(str,(raw['maxiter'],raw['atol'],raw['rtol'])))
    if name=='load_npz':
        # Independent NumPy archive input, rather than an archive made by Sci.
        ordered=sorted(triples,key=lambda v:(v[0],v[1]))
        pointers=np.zeros(matrix.shape[0]+1,dtype=np.int64)
        for r,_,_ in ordered:pointers[int(r)+1]+=1
        archive=io.BytesIO()
        writer=np.savez_compressed if raw['p'][0] else np.savez
        writer(archive,format=np.array(b'csr'),shape=np.array(matrix.shape,dtype=np.int64),
               data=np.array([v for _,_,v in ordered],dtype=float),
               indices=np.array([c for _,c,_ in ordered],dtype=np.int64),
               indptr=np.cumsum(pointers))
        fields.append(archive.getvalue().hex())
    return dict(fields=fields,raw=raw,suite=suite,suite_name=suite)


def recipes():
    cases=[]
    for suite in ('sparse','sparse_parity'):
        for i,raw in enumerate(json.loads((SCI/f'tests/scipy/{suite}_cases.json').read_text())):
            if raw['api']=='lu_solve':continue  # This record method is covered by splu.
            api=api_name(raw['api'])
            cases.append(dict(id=f'{api}/recipe-{suite}-{i:04d}',api=api,
                group='reference recipe; binary64 SciPy reference',args=materialize(suite,raw)))
    for scale in (1.,1e-20):
        raw=dict(id=f'gmres_scale_{scale}',api='gmres',
            a=[[4*scale,scale],[scale,3*scale]],x=[scale,2*scale],p=[],
            maxiter=2,atol=0.,rtol=1e-10)
        api=api_name('gmres')
        cases.append(dict(id=f'{api}/review-scale-{scale}',api=api,
            group='uniform operator/RHS scaling; two inner iterations',
            args=materialize('sparse_parity',raw)))
    return cases


@lru_cache(maxsize=2)
def program(suite):
    path=SCI/f'tests/scipy/generate_{suite}.py'
    tree=ast.parse(path.read_text(),filename=str(path))
    prefix=[]
    for node in tree.body:
        if isinstance(node,ast.With):
            loop=next(n for n in node.body if isinstance(n,ast.For))
            break
        if isinstance(node,(ast.Import,ast.ImportFrom,ast.FunctionDef,ast.Assert,ast.Assign)):
            prefix.append(node)
    namespace={'__file__':str(path),'__name__':'sparse_reference','io':io}
    exec(compile(ast.Module(body=prefix,type_ignores=[]),str(path),'exec'),namespace)
    body=[]
    for node in copy.deepcopy(loop.body):
        # Input/output writers are never executed. NPZ interchange stays in memory.
        if isinstance(node,ast.Expr) and isinstance(node.value,ast.Call) and isinstance(node.value.func,ast.Attribute) and isinstance(node.value.func.value,ast.Name) and node.value.func.value.id in ('iw','ew'):
            if node.value.func.value.id=='ew':
                node=ast.parse('return [float(result.shape[0]), float(result.shape[1]), *result.ravel().tolist()]').body[0]
            else:continue
        body.append(node)
    class MemoryArchive(ast.NodeTransformer):
        def visit_Assign(self,node):
            if any(isinstance(t,ast.Name) and t.id=='path' for t in node.targets):
                return ast.copy_location(ast.parse('path=io.BytesIO()').body[0],node)
            for target in node.targets:
                if isinstance(target,ast.Subscript) and isinstance(target.value,ast.Name) and target.value.id=='opts' and isinstance(target.slice,ast.Constant) and target.slice.value=='atol':
                    node.value=ast.parse("c.get('atol',1e-13)",mode='eval').body
            return self.generic_visit(node)
        def visit_Call(self,node):
            if isinstance(node.func,ast.Name) and node.func.id=='dict':
                for keyword in node.keywords:
                    if keyword.arg in ('maxiter','rtol'):
                        keyword.value=ast.Call(func=ast.Attribute(value=ast.Name(id='c',ctx=ast.Load()),attr='get',ctx=ast.Load()),args=[ast.Constant(keyword.arg),keyword.value],keywords=[])
            if isinstance(node.func,ast.Attribute) and node.func.attr=='load_npz':
                node.args=[ast.parse('(path.seek(0), path)[1]',mode='eval').body]
            return self.generic_visit(node)
    if suite=='sparse_parity':
        body=[MemoryArchive().visit(node) for node in body]
        body.append(ast.parse('return collected').body[0])
    function=ast.FunctionDef(name='evaluate',args=ast.arguments(posonlyargs=[],args=[ast.arg(arg='c')],kwonlyargs=[],kw_defaults=[],defaults=[]),body=body,decorator_list=[])
    exec(compile(ast.fix_missing_locations(ast.Module(body=[function],type_ignores=[])),str(path),'exec'),namespace)
    return namespace


def raw_reference(case):
    args=case['args'];raw=args['raw'];name=raw['api'];namespace=program(args['suite'])
    if args['suite']=='sparse':return namespace['evaluate'](raw)
    collected={}
    def emit(c,key,value):
        if namespace['s'].issparse(value):value=value.toarray()
        array=np.atleast_2d(np.asarray(value,dtype=float))
        values=array.ravel().tolist()
        if key in MATRICES.get(name,set()):values=[*map(float,array.shape),*values]
        collected[key]=values
    namespace['emit']=emit
    namespace['collected']=collected
    namespace['evaluate'](raw)
    collected.pop('pred',None) if name in ('reconstruct_path','construct_dist_matrix') else None
    if name=='reverse_cuthill_mckee':collected['permutation']=[float(len(raw['a']))]
    if name=='maximum_flow':collected.pop('a')  # Net flow is nonunique; conservation is checked in Rust.
    if name=='eigs':
        pairs=sorted(zip(collected.pop('wre'),collected.pop('wim')))
        collected['w']=[v for pair in pairs for v in pair]
    return [v for key in sorted(collected) for v in collected[key]]


def tagged(values):
    """Encode each scalar as (class,value): 0 finite, ±1 infinity, 2 NaN.

    This keeps unreachable graph vertices observable without treating them as
    failed accuracy samples. Finite payloads retain signed zero and precision.
    """
    import math
    return [v for x in values for v in ((0.,x) if math.isfinite(x) else
        (2.,0.) if math.isnan(x) else (math.copysign(1.,x),0.))]


def reference(case):
    return tagged(raw_reference(case))


def perturb(api, suite, raw, t, index):
    """Preserve shape/graph constraints while varying numerical inputs locally."""
    name=raw['api']
    factor=1+.001*t
    if name in ('random','random_array'):
        raw['p'][0]=.5+.49*t
    elif name=='get_index_dtype':
        raw['x']=[int(2**31+index-500)]
    elif name in ('eye','eye_array','identity'):
        raw['rows']=max(1,1-int(raw.get('k',0)))+index%16
        if name=='identity':raw['cols']=raw['rows']
        else:raw['cols']=max(raw['cols'],1+int(raw.get('k',0)))
    elif name in ('diags','diags_array','spdiags'):
        raw['x']=[[v*factor for v in row] for row in raw['x']]
    elif name in ('cg','bicg','bicgstab','cgs','qmr','tfqmr','gmres','lgmres',
                  'gcrotmk','minres','lsqr','lsmr','spilu','splu','factorized',
                  'spsolve','spsolve_triangular','matvec','aslinearoperator',
                  'funm_multiply_krylov','expm_multiply'):
        raw['x']=[v*factor for v in raw['x']]
    elif name=='maximum_flow':
        raw['a']=[[v*(1+index%8) for v in row] for row in raw['a']]
    elif name=='scale':raw['k']*=factor
    elif suite=='sparse':
        for key in ('a','b'):
            for entry in raw.get(key,[]):entry[2]*=factor
    else:
        raw['a']=[[v*factor for v in row] for row in raw['a']]
    return materialize(suite,raw)
