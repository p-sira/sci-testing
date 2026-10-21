"""Reuse independent SciPy fixture recipes without running their file writers.

The selected loop body executes with in-memory collectors for the input/output
CSV rows. Its numerical calls remain the upstream public SciPy calls.
"""
import ast
from functools import lru_cache
import types
from .api_registry import SCI


class CollectRows(ast.NodeTransformer):
    def visit_Expr(self,node):
        call=node.value
        if isinstance(call,ast.Call) and isinstance(call.func,ast.Attribute) and call.func.attr=='writerow':
            receiver=call.func.value
            if isinstance(receiver,ast.Name):
                name='fields' if receiver.id=='iw' else 'expected'
                return ast.copy_location(ast.Assign(targets=[ast.Name(id=name,ctx=ast.Store())],value=call.args[0]),node)
        return self.generic_visit(node)


@lru_cache(maxsize=8)
def program(module):
    path=SCI/f'tests/scipy/generate_{module}.py'
    tree=ast.parse(path.read_text(),filename=str(path))
    namespace={'__file__':str(path),'__name__':'fixture_recipe'}
    prefix=[];loop=None
    for node in tree.body:
        if isinstance(node,ast.With):
            loops=[v for v in ast.walk(node) if isinstance(v,ast.For) and ((isinstance(v.target,ast.Name) and v.target.id in ['c','case']) or isinstance(v.target,ast.Tuple))]
            if loops:loop=loops[0];break
        if isinstance(node,(ast.Import,ast.ImportFrom,ast.FunctionDef,ast.ClassDef,ast.Assert)):
            prefix.append(node)
        elif isinstance(node,ast.Assign) and not any(isinstance(t,ast.Name) and t.id=='cases' for t in node.targets):prefix.append(node)
    if loop is None:raise ValueError(f'{module}: no independent case loop')
    exec(compile(ast.Module(body=prefix,type_ignores=[]),str(path),'exec'),namespace)
    body=[CollectRows().visit(node) for node in loop.body]
    aliases=ast.parse('case=c\nindex=0').body
    fields_init=ast.Assign(targets=[ast.Name(id='fields',ctx=ast.Store())],value=ast.Constant(None))
    result=ast.Return(value=ast.Tuple(elts=[ast.Name(id='fields',ctx=ast.Load()),ast.Name(id='expected',ctx=ast.Load())],ctx=ast.Load()))
    function=ast.FunctionDef(name='evaluate',args=ast.arguments(posonlyargs=[],args=[ast.arg(arg='c')],kwonlyargs=[],kw_defaults=[],defaults=[]),body=[*aliases,fields_init,*body,result],decorator_list=[])
    generated=ast.fix_missing_locations(ast.Module(body=[function],type_ignores=[]))
    exec(compile(generated,str(path),'exec'),namespace)
    input_body=[]
    for node in body:
        input_body.append(node)
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='fields' for t in node.targets):break
    else:input_body=[ast.Assign(targets=[ast.Name(id='fields',ctx=ast.Store())],value=ast.Call(func=ast.Name(id='list',ctx=ast.Load()),args=[ast.Call(func=ast.Attribute(value=ast.Name(id='c',ctx=ast.Load()),attr='values',ctx=ast.Load()),args=[],keywords=[])],keywords=[]))]
    field_function=ast.FunctionDef(name='input_fields',args=function.args,body=[*aliases,*input_body,ast.Return(value=ast.Name(id='fields',ctx=ast.Load()))],decorator_list=[])
    exec(compile(ast.fix_missing_locations(ast.Module(body=[field_function],type_ignores=[])),str(path),'exec'),namespace)
    return namespace['evaluate'],namespace['input_fields']


def evaluate(module,raw):
    fields,expected=program(module)[0](raw)
    if fields is not None:fields=['' if value is None else str(value) for value in fields]
    else:fields=list(raw.values())
    def vector(text):return [float(x) for x in str(text).replace(':',';').split(';')] if str(text) else []
    if module in ['fft','fft_periodic']:
        re,im=vector(expected[2]),vector(expected[3])
        values=[x for pair in zip(re,im) for x in pair] if module=='fft' else re
    elif module=='interpolation':values=[x for row in expected[1:] for x in vector(row)]
    else:values=vector(expected[-1])
    return fields,values


def input_fields(module,raw):
    return ['' if value is None else str(value) for value in program(module)[1](raw)]
