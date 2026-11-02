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
                        str(ROOT/'native/spherical_extra.cpp'),
                        str(ROOT/'native/bessel_ratio.cpp'),
                        str(ROOT/'native/bessel_riccati.cpp'),
                        str(ROOT/'native/bessel_zeros.cpp'),
                        str(ROOT/'native/bessel_hankel.cpp'),
                        str(ROOT/'native/bessel_integrals.cpp'),
                        str(ROOT/'native/bessel_hypergeom.cpp'),
                        str(ROOT/'native/bessel_kelvin_k.cpp'),
                        str(ROOT/'native/bessel_struve_integrals.cpp'),
                        str(ROOT/'native/gsl_bessel_zeros.cpp'),
                        str(ROOT/'native/boost_stats_more.cpp'),
                        str(ROOT/'native/boost_special_transforms.cpp'),
                        str(ROOT/'native/gsl_special_derived.cpp'),
                        str(ROOT/'native/gsl_stats_more.cpp'),
                        str(ROOT/'native/gsl_quadrature_roots.cpp'),
                        str(ROOT/'native/boost_special_arrays.cpp'),
                        str(ROOT/'native/gsl_kelvin_ber.cpp'),
                        str(ROOT/'native/boost_stats_inference.cpp'),
                        str(ROOT/'native/gsl_stats_aggregates.cpp'),
                        '-lgsl', '-lgslcblas',
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
        from .bessel_ratio import BesselRatioComparisons
        self.bessel_ratio = BesselRatioComparisons(self.library)
        from .bessel_riccati import BesselRiccatiComparisons
        self.bessel_riccati = BesselRiccatiComparisons(self.library)
        from .bessel_zeros import BesselZerosComparisons
        self.bessel_zeros = BesselZerosComparisons(self.library)
        from .bessel_hankel import BesselHankelComparisons
        self.bessel_hankel = BesselHankelComparisons(self.library)
        from .bessel_integrals import BesselIntegralComparisons
        self.bessel_integrals = BesselIntegralComparisons(self.library)
        from .bessel_hypergeom import BesselHypergeomComparisons
        self.bessel_hypergeom = BesselHypergeomComparisons(self.library)
        from .bessel_kelvin_k import BesselKelvinKComparisons
        self.bessel_kelvin_k = BesselKelvinKComparisons(self.library)
        from .bessel_struve_integrals import BesselStruveIntegralsComparisons
        self.bessel_struve_integrals = BesselStruveIntegralsComparisons(self.library)
        from .gsl_bessel_zeros import GSLBesselZerosComparisons
        self.gsl_bessel_zeros = GSLBesselZerosComparisons(self.library)
        from .boost_stats_more import BoostStatsMoreComparisons
        self.boost_stats_more = BoostStatsMoreComparisons(self.library)
        from .boost_special_transforms import BoostSpecialTransformsComparisons
        self.boost_special_transforms = BoostSpecialTransformsComparisons(self.library)
        from .gsl_special_derived import GslSpecialDerivedComparisons
        self.gsl_special_derived = GslSpecialDerivedComparisons(self.library)
        from .gsl_stats_more import StatsMoreGSLComparisons
        self.gsl_stats_more = StatsMoreGSLComparisons(self.library)
        from .gsl_quadrature_roots import GSLQuadratureRootsComparisons
        self.gsl_quadrature_roots = GSLQuadratureRootsComparisons(self.library)
        from .boost_special_arrays import BoostSpecialArraysComparisons
        self.boost_special_arrays = BoostSpecialArraysComparisons(self.library)
        from .gsl_kelvin_ber import GSLKelvinBerComparisons
        self.gsl_kelvin_ber = GSLKelvinBerComparisons(self.library)
        from .boost_stats_inference import BoostStatsInferenceComparisons
        self.boost_stats_inference = BoostStatsInferenceComparisons(self.library)
        from .gsl_stats_aggregates import StatsAggregatesGSLComparisons
        self.gsl_stats_aggregates = StatsAggregatesGSLComparisons(self.library)

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
        # Independent providers may implement the same Sci API through distinct
        # native primitives. Preserve existing dispatch precedence within a
        # provider, while retaining supplementary library comparisons.
        supplemental = {}
        for adapter in (self.gsl_bessel_zeros, self.boost_stats_more,
                        self.boost_special_transforms, self.gsl_special_derived,
                        self.gsl_stats_more, self.gsl_quadrature_roots,
                        self.boost_special_arrays, self.gsl_kelvin_ber,
                        self.boost_stats_inference, self.gsl_stats_aggregates):
            supplemental.update(adapter.evaluate(case))
        return {**supplemental, **self._evaluate_existing(case)}

    def _evaluate_existing(self, case):
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
            for adapter in (self.bessel_ratio, self.bessel_riccati, self.bessel_zeros,
                            self.bessel_hankel, self.bessel_integrals,
                            self.bessel_hypergeom, self.bessel_kelvin_k,
                            self.bessel_struve_integrals,
                            self.vector_special, self.combinatorics_extra, self.spherical_extra):
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
