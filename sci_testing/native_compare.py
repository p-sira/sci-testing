"""Build and call public Boost.Math/GSL APIs through a small checked C ABI."""
import ctypes as ct
import os
import subprocess
from functools import lru_cache

from . import ROOT


@lru_cache(maxsize=1)
def _native_library():
    """Compile and load the native bridge once per Python process."""
    target = ROOT/'target/native-comparisons.so'
    target.parent.mkdir(exist_ok=True)
    subprocess.run([os.environ.get('CXX', 'c++'), '-std=c++17', '-O2', '-shared', '-fPIC',
                        str(ROOT/'native/comparisons.cpp'),
                        str(ROOT/'native/stats_distributions.cpp'),
                        str(ROOT/'native/stats_summary.cpp'),
                        str(ROOT/'native/probability_extra.cpp'),
                        str(ROOT/'native/special_extra.cpp'),
                        str(ROOT/'native/stats_extra.cpp'),
                        str(ROOT/'native/complex_extra.cpp'),
                        str(ROOT/'native/vector_special.cpp'),
                        str(ROOT/'native/combinatorics_extra.cpp'),
                        str(ROOT/'native/spherical_extra.cpp'), '-lgsl', '-lgslcblas',
                        '-o', str(target)], check=True)
    return ct.CDLL(str(target))


class NativeComparisons:
    def __init__(self):
        self.library = _native_library()
        signatures = {
            'comparison_version': ([ct.c_int], ct.c_char_p),
            'comparison_error': ([], ct.c_char_p),
            'comparison_supported': ([ct.c_char_p], ct.c_int),
            'comparison_library_supported': ([ct.c_int, ct.c_char_p], ct.c_int),
            'comparison_output_count': ([ct.c_char_p], ct.c_int),
            'comparison_complex_supported': ([ct.c_char_p], ct.c_int),
            'comparison_complex_evaluate': ([ct.c_int, ct.c_char_p,
                                             ct.POINTER(ct.c_double), ct.c_int,
                                             ct.POINTER(ct.c_double), ct.c_int], ct.c_int),
            'comparison_evaluate': ([ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
                                     ct.POINTER(ct.c_double)], ct.c_int),
            'comparison_evaluate_many': ([ct.c_int, ct.c_char_p, ct.POINTER(ct.c_double), ct.c_int,
                                          ct.POINTER(ct.c_double), ct.c_int], ct.c_int),
        }
        for name, (arguments, result) in signatures.items():
            function = getattr(self.library, name)
            function.argtypes, function.restype = arguments, result
        self.versions = {name: self.library.comparison_version(i).decode()
                         for i, name in enumerate(['Boost.Math', 'GSL'])}
        from .stats_distributions import StatsDistributionComparisons
        self.distributions = StatsDistributionComparisons(self.library)
        from .stats_summary import StatsSummaryComparisons
        self.summaries = StatsSummaryComparisons(self.library)
        from .probability_extra import ProbabilityExtraComparisons
        self.probability_extra = ProbabilityExtraComparisons(self.library)
        from .special_extra import SpecialExtraComparisons
        self.special_extra = SpecialExtraComparisons(self.library)
        from .stats_extra import StatsExtraComparisons
        self.stats_extra = StatsExtraComparisons(self.library)
        from .complex_extra import ComplexExtraComparisons
        self.complex_extra = ComplexExtraComparisons(self.library)
        from .vector_special import VectorSpecialComparisons
        self.vector_special = VectorSpecialComparisons(self.library)
        from .combinatorics_extra import CombinatoricsExtraComparisons
        self.combinatorics_extra = CombinatoricsExtraComparisons(self.library)
        from .spherical_extra import SphericalExtraComparisons
        self.spherical_extra = SphericalExtraComparisons(self.library)

    def evaluate_complex(self, case):
        """Compare supported Sci complex-special inputs with public GSL functions."""
        api = case['api']
        if not api.startswith('special.complex.'):
            return {}
        extra = self.complex_extra.evaluate(case)
        if extra:
            return extra
        name = api.rsplit('.', 1)[-1]
        encoded = name.encode()
        if not self.library.comparison_complex_supported(encoded):
            return {}
        z = case['args'].get('z')
        if not isinstance(z, (list, tuple)) or len(z) != 2:
            raise ValueError(f'{name}: expected a complex [real, imag] pair')
        arguments = (ct.c_double * 2)(*z)
        output = (ct.c_double * 2)()
        status = self.library.comparison_complex_evaluate(1, encoded, arguments,
                                                         2, output, 2)
        if status:
            return {'GSL': {'error': self.library.comparison_error().decode()}}
        return {'GSL': {'values': list(output)}}

    def evaluate(self, case):
        module, _, name = case['api'].partition('.')
        if module not in ['special', 'stats']: return {}
        if case['api'] in ('stats.pdf', 'stats.logpdf', 'stats.cdf', 'stats.sf', 'stats.ppf'):
            return self.distributions.evaluate(case)
        if module == 'stats':
            extra = self.stats_extra.evaluate(case)
            if extra:
                return extra
            summarized = self.summaries.evaluate(case)
            if summarized:
                return summarized
        if module == 'special' and not name.startswith('complex.'):
            for adapter in (self.vector_special, self.combinatorics_extra, self.spherical_extra):
                extra = adapter.evaluate(case)
                if extra:
                    return extra
            extra = self.probability_extra.evaluate(case)
            if extra:
                return extra
            extra = self.special_extra.evaluate(case)
            if extra:
                return extra
        # The native ABI accepts real doubles, not complex-valued submodules.
        if name.startswith('complex.'): return {}
        # 0.3 uses mathematical-family modules. Native Boost/GSL identifiers
        # remain the final function name, irrespective of Sci's namespace.
        name = name.rsplit('.',1)[-1]
        encoded = name.encode()
        count = self.library.comparison_supported(encoded)
        if not count: return {}
        output_count = self.library.comparison_output_count(encoded)
        a = case['args']
        if 'p0' in a:
            values = [a[f'p{i}'] for i in range(count)]
        elif name == 'gamma_ratio' and a.get('special03_case'):
            values = [a['a'], a['b']]
        elif count == 1 and 'x' in a:
            values = [a['x']]
        else:
            return {}
        if len(values) != count: raise ValueError(f'{name}: native arity mismatch')
        arguments = (ct.c_double*count)(*values)
        results = {}
        for i, library in enumerate(self.versions):
            if not self.library.comparison_library_supported(i, encoded):
                continue
            output = (ct.c_double*output_count)()
            if output_count == 1:
                status = self.library.comparison_evaluate(i, encoded, arguments, count, output)
            else:
                status = self.library.comparison_evaluate_many(i, encoded, arguments, count,
                                                                 output, output_count)
            result={'values':list(output)} if status >= 0 else {}
            if status:result['error']=self.library.comparison_error().decode()
            results[library]=result
        return results
