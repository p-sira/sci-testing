"""Build and call public Boost.Math/GSL APIs through a small checked C ABI."""
import ctypes as ct
import fcntl
import hashlib
import json
import os
import platform
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path

from . import ROOT


# Preserve this explicit list: not every C++ file in native/ is a separate TU,
# and some are included by other source files.
_NATIVE_SOURCES = (
    'comparisons.cpp', 'stats_distributions.cpp', 'stats_summary.cpp',
    'probability_extra.cpp', 'special_extra.cpp', 'stats_extra.cpp',
    'complex_extra.cpp', 'vector_special.cpp', 'combinatorics_extra.cpp',
    'spherical_extra.cpp', 'bessel_ratio.cpp', 'bessel_riccati.cpp',
    'bessel_zeros.cpp', 'bessel_hankel.cpp', 'bessel_integrals.cpp',
    'bessel_hypergeom.cpp', 'bessel_kelvin_k.cpp',
    'bessel_struve_integrals.cpp', 'gsl_bessel_zeros.cpp',
    'boost_stats_more.cpp', 'boost_special_transforms.cpp',
    'gsl_special_derived.cpp', 'gsl_stats_more.cpp',
    'gsl_quadrature_roots.cpp', 'boost_special_arrays.cpp',
    'gsl_kelvin_ber.cpp', 'boost_stats_inference.cpp',
    'gsl_stats_aggregates.cpp', 'gsl_numerical_core.cpp',
    'boost_complex_multiprecision.cpp', 'boost_graph.cpp',
    'boost_numerical_core.cpp', 'boost_ublas.cpp', 'gsl_linalg.cpp',
    'boost_spatial.cpp', 'gsl_spatial.cpp', 'boost_interpolation.cpp',
    'boost_units.cpp', 'gsl_sparse.cpp', 'gsl_signal_stats.cpp',
    'gsl_interpolation.cpp', 'gsl_optimize.cpp', 'gsl_fft_simple.cpp',
    'extra_special_scalar.cpp', 'extra_stats_scalar.cpp',
)
_NATIVE_FLAGS = ('-std=c++17', '-O2', '-shared', '-fPIC')
_NATIVE_LIBS = ('-lgsl', '-lgslcblas')
_BUILD_ENV = ('CPATH', 'CPLUS_INCLUDE_PATH', 'C_INCLUDE_PATH', 'LIBRARY_PATH',
              'LD_LIBRARY_PATH', 'GCC_EXEC_PREFIX', 'COMPILER_PATH')


def _compiler_command():
    # CXX may be a compiler wrapper, e.g. "ccache g++".
    command = shlex.split(os.environ.get('CXX', 'c++'))
    if not command:
        raise ValueError('CXX must specify a C++ compiler')
    compiler = shutil.which(command[0])
    if compiler is None:
        raise FileNotFoundError(f'C++ compiler not found: {command[0]}')
    return (str(Path(compiler).resolve()), *command[1:])


def _compiler_and_libraries(compiler):
    """Discover the actual Boost/GSL headers selected by this C++ compiler."""
    version = subprocess.run([*compiler, '--version'], check=True,
                             capture_output=True, text=True).stdout
    headers = subprocess.run(
        [*compiler, '-E', '-P', '-x', 'c++', '-'],
        input=('#include <boost/version.hpp>\n'
               '#include <gsl/gsl_version.h>\n'
               'BOOST_LIB_VERSION\nGSL_VERSION\n'),
        check=True, capture_output=True, text=True).stdout
    libraries = []
    # GCC and Clang both implement -print-file-name. Record the real linked
    # libraries' inode metadata as well as header versions, when discoverable.
    for name in ('libgsl.so', 'libgslcblas.so'):
        result = subprocess.run([*compiler, '-print-file-name='+name],
                                check=True, capture_output=True, text=True)
        path = Path(result.stdout.strip())
        if path.is_file():
            resolved = path.resolve()
            s = resolved.stat()
            libraries.append((name, str(resolved), s.st_size, s.st_mtime_ns))
        else:
            libraries.append((name, result.stdout.strip()))
    return {'compiler_version': version, 'library_headers': headers,
            'linked_libraries': libraries}


def _native_fingerprint(compiler, sources):
    metadata = {
        'schema': 1,
        'compiler': compiler,
        'toolchain': _compiler_and_libraries(compiler),
        'compile_flags': _NATIVE_FLAGS,
        'link_flags': _NATIVE_LIBS,
        'environment': {key: os.environ.get(key) for key in _BUILD_ENV},
        'platform': (platform.system(), platform.machine()),
        'sources': [],
    }
    for source in sources:
        data = source.read_bytes()
        metadata['sources'].append(
            (str(source.resolve()), hashlib.sha256(data).hexdigest()))
    return hashlib.sha256(json.dumps(
        metadata, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def _file_sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _cache_valid(binary, manifest, fingerprint):
    try:
        metadata = json.loads(manifest.read_text(encoding='utf-8'))
        return (metadata.get('fingerprint') == fingerprint and
                binary.is_file() and binary.stat().st_size > 4 and
                metadata.get('sha256') == _file_sha256(binary))
    except (OSError, ValueError, TypeError):
        return False


def _publish_native_binary(compiler, sources, binary, manifest, fingerprint):
    """Compile privately and atomically publish only a verified ELF binary."""
    temporary = None
    metadata_temp = None
    published = False
    try:
        fd, filename = tempfile.mkstemp(dir=binary.parent,
                                        prefix='.native-comparisons-',
                                        suffix='.so')
        os.close(fd)
        temporary = Path(filename)
        subprocess.run([*compiler, *_NATIVE_FLAGS, *map(str, sources),
                        *_NATIVE_LIBS, '-o', str(temporary)], check=True)
        with temporary.open('rb') as stream:
            elf_magic = stream.read(4)
        if temporary.stat().st_size < 4 or elf_magic != b'\x7fELF':
            raise RuntimeError('native compiler did not produce a shared ELF library')
        # Protect against source/header/toolchain changes during the long build.
        if _native_fingerprint(compiler, sources) != fingerprint:
            raise RuntimeError('native build inputs changed during compilation; retry')
        digest = _file_sha256(temporary)
        fd, filename = tempfile.mkstemp(dir=binary.parent,
                                        prefix='.native-comparisons-',
                                        suffix='.json')
        os.close(fd)
        metadata_temp = Path(filename)
        metadata_temp.write_text(json.dumps(
            {'fingerprint': fingerprint, 'sha256': digest},
            sort_keys=True) + '\n', encoding='utf-8')
        os.replace(temporary, binary)
        published = True
        os.replace(metadata_temp, manifest)
    except BaseException:
        # A binary lacking its matching manifest must never be considered
        # cache-valid. Clean up even on interruption after the binary rename.
        if published:
            binary.unlink(missing_ok=True)
            manifest.unlink(missing_ok=True)
        raise
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
        if metadata_temp is not None:
            metadata_temp.unlink(missing_ok=True)


def _native_library():
    """Reuse a checked, content-addressed .so across processes and runs.

    The global flock covers cache inspection, compilation and publication.
    Fingerprint-specific paths avoid dlopen reusing an old inode if the source
    changes within the current Python process. No new Python dependencies.
    """
    cache = ROOT/'target'
    cache.mkdir(parents=True, exist_ok=True)
    sources = tuple(ROOT/'native'/filename for filename in _NATIVE_SOURCES)
    compiler = _compiler_command()
    with (cache/'.native-comparisons.lock').open('a+b') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            fingerprint = _native_fingerprint(compiler, sources)
            binary = cache/f'native-comparisons-{fingerprint}.so'
            manifest = cache/f'native-comparisons-{fingerprint}.json'
            if not _cache_valid(binary, manifest, fingerprint):
                _publish_native_binary(compiler, sources, binary,
                                       manifest, fingerprint)
            return ct.CDLL(str(binary))
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)


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
        # All Boost sublibraries use the same Boost installation, but are
        # separately identified in the saved comparator columns.
        boost_version = self.library.comparison_version(0).decode()
        self.versions = {'Boost.Math': boost_version,
                         'Boost.Multiprecision': boost_version,
                         'Boost.Numeric.uBLAS': boost_version,
                         'Boost.Graph': boost_version,
                         'Boost.Units': boost_version,
                         'GSL': self.library.comparison_version(1).decode()}
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
        from .gsl_numerical_core import GSLNumericalCoreComparisons
        self.gsl_numerical_core = GSLNumericalCoreComparisons(self.library)
        from .boost_complex_multiprecision import BoostComplexMultiprecisionComparisons
        self.boost_complex_multiprecision = BoostComplexMultiprecisionComparisons(self.library)
        from .boost_graph import BoostGraphComparisons
        self.boost_graph = BoostGraphComparisons(self.library)
        from .boost_numerical_core import BoostNumericalCoreComparisons
        self.boost_numerical_core = BoostNumericalCoreComparisons(self.library)
        from .boost_ublas import BoostUblasComparisons
        self.boost_ublas = BoostUblasComparisons(self.library)
        from .gsl_linalg import GslLinalgComparisons
        self.gsl_linalg = GslLinalgComparisons(self.library)
        from .boost_spatial import BoostSpatialComparisons
        self.boost_spatial = BoostSpatialComparisons(self.library)
        from .gsl_spatial import GSLSpatialComparisons
        self.gsl_spatial = GSLSpatialComparisons(self.library)
        from .boost_interpolation import BoostInterpolationComparisons
        self.boost_interpolation = BoostInterpolationComparisons(self.library)
        from .boost_units import BoostUnitsComparisons
        self.boost_units = BoostUnitsComparisons(self.library)
        from .gsl_sparse import GslSparseComparisons
        self.gsl_sparse = GslSparseComparisons(self.library)
        from .gsl_signal_stats import GSLSignalStatsComparisons
        self.gsl_signal_stats = GSLSignalStatsComparisons(self.library)
        from .gsl_interpolation import GslInterpolationComparisons
        self.gsl_interpolation = GslInterpolationComparisons(self.library)
        from .gsl_optimize import GSLOptimizeComparisons
        self.gsl_optimize = GSLOptimizeComparisons(self.library)
        from .gsl_fft_simple import GSLFFTSimpleComparisons
        self.gsl_fft_simple = GSLFFTSimpleComparisons(self.library)
        from .extra_native_comparisons import ExtraNativeComparisons
        self.extra_native = ExtraNativeComparisons(self.library)
        from .extra_native_stats import ExtraNativeStatsComparisons
        self.extra_native_stats = ExtraNativeStatsComparisons(self.library)

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
        for adapter in (self.extra_native, self.extra_native_stats,
                        self.gsl_bessel_zeros, self.boost_stats_more,
                        self.boost_special_transforms, self.gsl_special_derived,
                        self.gsl_stats_more, self.gsl_quadrature_roots,
                        self.boost_special_arrays, self.gsl_kelvin_ber,
                        self.boost_stats_inference, self.gsl_stats_aggregates,
                        self.gsl_numerical_core,
                        self.boost_complex_multiprecision, self.boost_graph,
                        self.boost_numerical_core, self.boost_ublas,
                        self.gsl_linalg, self.boost_spatial, self.gsl_spatial,
                        self.boost_interpolation, self.boost_units,
                        self.gsl_sparse, self.gsl_signal_stats,
                        self.gsl_interpolation, self.gsl_optimize,
                        self.gsl_fft_simple):
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
        for i, library in enumerate(('Boost.Math', 'GSL')):
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
