"""Build and call public Boost.Math/GSL APIs through a small checked C ABI."""
import ctypes as ct
import subprocess

from . import ROOT


class NativeComparisons:
    def __init__(self):
        target = ROOT/'target/native-comparisons.so'
        target.parent.mkdir(exist_ok=True)
        subprocess.run(['c++', '-std=c++17', '-O2', '-shared', '-fPIC',
                        str(ROOT/'native/comparisons.cpp'), '-lgsl', '-lgslcblas',
                        '-o', str(target)], check=True)
        self.library = ct.CDLL(str(target))
        signatures = {
            'comparison_version': ([ct.c_int], ct.c_char_p),
            'comparison_error': ([], ct.c_char_p),
            'comparison_supported': ([ct.c_char_p], ct.c_int),
            'comparison_evaluate': ([ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
                                     ct.POINTER(ct.c_double)], ct.c_int),
        }
        for name, (arguments, result) in signatures.items():
            function = getattr(self.library, name)
            function.argtypes, function.restype = arguments, result
        self.versions = {name: self.library.comparison_version(i).decode()
                         for i, name in enumerate(['Boost.Math', 'GSL'])}

    def evaluate(self, case):
        module, name = case['api'].split('.',1)
        if module not in ['special', 'stats']: return {}
        count = self.library.comparison_supported(name.encode())
        if not count: return {}
        a = case['args']
        values = [a[f'p{i}'] for i in range(len(a))] if 'p0' in a else [a['x']] if 'x' in a else []
        if len(values) != count: raise ValueError(f'{name}: native arity mismatch')
        arguments = (ct.c_double*count)(*values)
        results = {}
        for i, library in enumerate(self.versions):
            output = ct.c_double()
            status = self.library.comparison_evaluate(i, name.encode(), arguments, count, ct.byref(output))
            result={'values':[output.value]} if status >= 0 else {}
            if status:result['error']=self.library.comparison_error().decode()
            results[library]=result
        return results
