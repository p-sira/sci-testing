"""Optional high-precision WolframEngine mathematical references.

This is an intentionally curated subset of Sci signatures.  ``reference`` and
``references`` return decimal strings, not rounded binary64 measurements.
WolframEngine is launched through ``wolframscript -code`` without a shell; no
caller-provided expression text is ever evaluated.  Every Python binary64 float
is transmitted as its *exact* integer ratio, rather than a decimal approximation.

Wolfram's standard principal branches apply.  These particular APIs have real
outputs in their supported domains; a nonreal, unevaluated, nonfinite or
insufficient-precision result is an explicit failure.  More signatures can be
added only after checking their argument order, conventions and output shape.
"""

from __future__ import annotations

import json
import math
import os
import re
import shutil
import subprocess
from collections.abc import Mapping, Sequence
from decimal import Decimal, InvalidOperation, localcontext
from fractions import Fraction
from functools import lru_cache
from typing import cast


class WolframReferenceError(RuntimeError):
    """An available Wolfram reference could not be evaluated faithfully."""


class WolframUnavailableError(WolframReferenceError):
    """WolframEngine is missing, unlicensed or did not finish."""


class WolframEvaluationError(WolframReferenceError):
    """A requested expression yielded no valid finite real reference."""


class UnsupportedReferenceError(ValueError):
    """The Sci API, signature or input domain has not been curated."""


# Formatted, fixed Wolfram Language function calls; the only interpolations are
# exact rational literals constructed from validated Python numeric values.
# The index order is the Sci input order, not necessarily the Wolfram order.
_SPECS: dict[str, tuple[str, int]] = {
    'agm': ('ArithmeticGeometricMean[{0},{1}]', 2),
    'beta': ('Beta[{0},{1}]', 2),
    'binom': ('Binomial[{0},{1}]', 2),
    'digamma': ('PolyGamma[{0}]', 1),
    'psi': ('PolyGamma[{0}]', 1),
    'trigamma': ('PolyGamma[1,{0}]', 1),
    'polygamma': ('PolyGamma[{0},{1}]', 2),
    'gamma': ('Gamma[{0}]', 1),
    'rgamma': ('1/Gamma[{0}]', 1),
    'erf': ('Erf[{0}]', 1),
    'erfc': ('Erfc[{0}]', 1),
    'erfi': ('Erfi[{0}]', 1),
    'erfinv': ('InverseErf[{0}]', 1),
    'exp1': ('ExpIntegralE[1,{0}]', 1),
    'expi': ('ExpIntegralEi[{0}]', 1),
    'expn': ('ExpIntegralE[{0},{1}]', 2),
    'log1p': ('Log[1+{0}]', 1),
    'expm1': ('Exp[{0}]-1', 1),
    'i0': ('BesselI[0,{0}]', 1),
    'i1': ('BesselI[1,{0}]', 1),
    'j0': ('BesselJ[0,{0}]', 1),
    'j1': ('BesselJ[1,{0}]', 1),
    'k0': ('BesselK[0,{0}]', 1),
    'k1': ('BesselK[1,{0}]', 1),
    'y0': ('BesselY[0,{0}]', 1),
    'y1': ('BesselY[1,{0}]', 1),
    'iv': ('BesselI[{0},{1}]', 2),
    'jv': ('BesselJ[{0},{1}]', 2),
    'kv': ('BesselK[{0},{1}]', 2),
    'yv': ('BesselY[{0},{1}]', 2),
    'jn': ('BesselJ[{0},{1}]', 2),
    'kn': ('BesselK[{0},{1}]', 2),
    'yn': ('BesselY[{0},{1}]', 2),
    'hyp0f1': ('Hypergeometric0F1[{0},{1}]', 2),
    'hyp1f1': ('Hypergeometric1F1[{0},{1},{2}]', 3),
    'hyp2f1': ('Hypergeometric2F1[{0},{1},{2},{3}]', 4),
    'hyperu': ('HypergeometricU[{0},{1},{2}]', 3),
    'elliprc': ('CarlsonRC[{0},{1}]', 2),
    'elliprf': ('CarlsonRF[{0},{1},{2}]', 3),
    'elliprd': ('CarlsonRD[{0},{1},{2}]', 3),
    'elliprj': ('CarlsonRJ[{0},{1},{2},{3}]', 4),
    'elliprg': ('CarlsonRG[{0},{1},{2}]', 3),
    'ellipk': ('EllipticK[{0}]', 1),
    'ellipe': ('EllipticE[{0}]', 1),
    'ellipf': ('EllipticF[{0},{1}]', 2),
    'ellipeinc': ('EllipticE[{0},{1}]', 2),
    'ellipkinc': ('EllipticF[{0},{1}]', 2),
    'ellippi': ('EllipticPi[{0},{1}]', 2),
    'ellippiinc': ('EllipticPi[{1},{0},{2}]', 3),
    'ellippiinc_bulirsch': ('EllipticPi[{1},{0},{2}]', 3),
    'ellipkm1': ('EllipticK[1-{0}]', 1),
    'ellipd': ('CarlsonRD[0,1-{0},1]/3', 1),
    'eval_legendre': ('LegendreP[{0},{1}]', 2),
    'eval_chebyt': ('ChebyshevT[{0},{1}]', 2),
    'eval_chebyu': ('ChebyshevU[{0},{1}]', 2),
    'eval_hermite': ('HermiteH[{0},{1}]', 2),
    'eval_laguerre': ('LaguerreL[{0},{1}]', 2),
    'eval_genlaguerre': ('LaguerreL[{0},{1},{2}]', 3),
    'eval_jacobi': ('JacobiP[{0},{1},{2},{3}]', 4),
    'eval_gegenbauer': ('GegenbauerC[{0},{1},{2}]', 3),
    'zeta': ('Zeta[{0}]', 1),
    'zetac': ('Zeta[{0}]-1', 1),
    'sinc': ('Sinc[Pi*{0}]', 1),
    'spence': ('PolyLog[2,1-{0}]', 1),
    'harmonic': ('HarmonicNumber[{0}]', 1),
    'lambertw': ('ProductLog[{0}]', 1),
    'lambertwm1': ('ProductLog[-1,{0}]', 1),
    'ndtr': ('Erfc[-{0}/Sqrt[2]]/2', 1),
    'ndtri': ('Sqrt[2]*InverseErf[2*{0}-1]', 1),
    'erfcx': ('Exp[{0}^2]*Erfc[{0}]', 1),
    'dawsn': ('Sqrt[Pi]/2*Exp[-{0}^2]*Erfi[{0}]', 1),
    'expit': ('1/(1+Exp[-{0}])', 1),
    'logit': ('Log[{0}/(1-{0})]', 1),
    'cbrt': ('Surd[{0},3]', 1),
    'exp2': ('2^{0}', 1),
    'exp10': ('10^{0}', 1),
}

_MARKER = '__SCI_WOLFRAM_REFERENCE_V1__:'
_NUMBER = re.compile(r'^([+-]?(?:\d+(?:\.\d*)?|\.\d+))(?:`[0-9.]*)?(?:\*\^([+-]?\d+))?$')
# Wolfram's 0``a denotes zero with *absolute accuracy* a, not a
# significant-digit precision.  It can arise from Im[real-valued result].
_IMAGINARY_ACCURACY_ZERO = re.compile(r'^0(?:\.0*)?``(\d+(?:\.\d*)?|\.\d+)$')
_DECIMAL = re.compile(r'^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$')
_MAX_ARGV = 100_000  # Linux individual argv strings are limited to ~128 KiB.


def _rational(value: object) -> str:
    if isinstance(value, bool):
        raise UnsupportedReferenceError('Boolean values are not numeric inputs')
    if isinstance(value, int):
        fraction = Fraction(value)
    elif isinstance(value, float):
        if not math.isfinite(value):
            raise UnsupportedReferenceError('Nonfinite binary64 input')
        # WL converts -0. to 0.; silently discarding the sign could change
        # branch-side semantics for functions with branch cuts at the origin.
        if value == 0.0 and math.copysign(1.0, value) < 0:
            raise UnsupportedReferenceError('Wolfram does not preserve signed negative zero')
        fraction = Fraction(*value.as_integer_ratio())
    elif isinstance(value, (str, Decimal)):
        raw = str(value)
        if len(raw) > 1_000 or not _DECIMAL.fullmatch(raw):
            raise UnsupportedReferenceError('Only finite decimal numeric strings are supported')
        try:
            decimal_value = Decimal(raw)
            if decimal_value.is_zero() and decimal_value.is_signed():
                raise UnsupportedReferenceError('Wolfram does not preserve signed negative zero')
            fraction = Fraction(decimal_value)
        except (InvalidOperation, ValueError, OverflowError) as exc:
            raise UnsupportedReferenceError('Invalid decimal numeric input') from exc
    else:
        raise UnsupportedReferenceError(f'Unsupported numeric input type {type(value).__name__}')
    numerator, denominator = fraction.numerator, fraction.denominator
    if len(str(numerator)) > 2_000 or len(str(denominator)) > 2_000:
        raise UnsupportedReferenceError('Input rational is too large')
    return f'Rational[{numerator},{denominator}]'


def _expression(case: Mapping[str, object]) -> str:
    api = case.get('api')
    args_value = case.get('args')
    if not isinstance(api, str) or not isinstance(args_value, dict):
        raise UnsupportedReferenceError('Case needs an API string and argument mapping')
    args = cast(dict[str, object], args_value)
    if api in ('special.erf', 'special.gamma') and set(args) == {'x'}:
        spec, _ = _SPECS[api.removeprefix('special.')]
        return spec.format(_rational(args['x']))
    if api == 'stats.normal_cdf' and set(args) == {'x'}:
        return f'Erfc[-{_rational(args["x"])}/Sqrt[2]]/2'
    if api == 'integrate.fixed_quad' and set(args) == {'a', 'b'}:
        a, b = _rational(args['a']), _rational(args['b'])
        return f'Exp[{b}]-Exp[{a}]'
    if api == 'differentiate.derivative' and set(args) == {'x'}:
        return f'Cos[{_rational(args["x"])}]'
    if api == 'stats.variance' and set(args) == {'values'}:
        raw_values = args['values']
        if not isinstance(raw_values, (list, tuple)):
            raise UnsupportedReferenceError('Sample variance requires two or more values')
        values = cast(Sequence[object], raw_values)
        if len(values) < 2:
            raise UnsupportedReferenceError('Sample variance requires two or more values')
        return 'Variance[{' + ','.join(_rational(v) for v in values) + '}]'
    if not api.startswith('special.') or api.count('.') != 1:
        raise UnsupportedReferenceError(f'Unsupported API: {api}')
    name = api.removeprefix('special.')
    if name == 'zeta' and set(args) == {'p0', 'p1'}:
        return f'Zeta[{_rational(args["p0"])},{_rational(args["p1"])}]'
    spec = _SPECS.get(name)
    if spec is None:
        raise UnsupportedReferenceError(f'Unsupported API: {api}')
    expression, count = spec
    if set(args) != {f'p{i}' for i in range(count)}:
        raise UnsupportedReferenceError(f'Unsupported argument signature for {api}: {sorted(args)}')
    return expression.format(*(_rational(args[f'p{i}']) for i in range(count)))


def supports(case: Mapping[str, object]) -> bool:
    """Whether this exact API signature and numeric input representation is curated."""
    try:
        _ = _expression(case)
    except (UnsupportedReferenceError, OverflowError):
        return False
    return True


def definition(case: Mapping[str, object]) -> str:
    """Exact, input-free Wolfram Language template for a supported case.

    Placeholders are numbered by the Sci argument order (``p0``, ``p1``, ...)
    or the documented named input order.  This records the chosen mathematical
    convention without placing case-specific, user-controlled source in a
    provenance record.  Unsupported signatures raise explicitly.
    """
    _ = _expression(case)
    api = case.get('api')
    args_value = case.get('args')
    if not isinstance(api, str) or not isinstance(args_value, dict):
        raise UnsupportedReferenceError('Case needs an API string and argument mapping')
    args = cast(dict[str, object], args_value)
    if api == 'stats.normal_cdf':
        return 'Erfc[-{0}/Sqrt[2]]/2'
    if api == 'stats.variance':
        return 'Variance[{values...}]'
    if api == 'integrate.fixed_quad':
        return 'Exp[{1}]-Exp[{0}]  (arguments: a, b)'
    if api == 'differentiate.derivative':
        return 'Cos[{0}]'
    name = api.removeprefix('special.')
    if name == 'zeta' and set(args) == {'p0', 'p1'}:
        return 'Zeta[{0},{1}]'
    return _SPECS[name][0]


def _executable() -> str:
    requested = os.environ.get('WOLFRAM_REFERENCE_EXECUTABLE', 'wolframscript')
    executable = shutil.which(requested)
    if executable is None:
        raise WolframUnavailableError(f'Wolfram executable unavailable: {requested}')
    return executable


def _execute(code: str, *, timeout: float, executable: str | None = None) -> str:
    if len(code.encode('utf-8')) >= _MAX_ARGV:
        raise WolframReferenceError('Wolfram batch code exceeds subprocess argument limit')
    try:
        process = subprocess.run(
            [executable or _executable(), '-code', code], check=False,
            capture_output=True,
            text=True, encoding='utf-8', errors='replace', timeout=timeout,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise WolframUnavailableError(f'WolframEngine did not complete: {exc}') from exc
    if process.returncode:
        raise WolframUnavailableError(
            f'WolframEngine exited {process.returncode}: {process.stderr[-2000:]}')
    if process.stdout.count(_MARKER) != 1:
        raise WolframReferenceError(
            f'Unexpected WolframEngine output: {process.stdout[-2000:]} {process.stderr[-1000:]}')
    return process.stdout.split(_MARKER, 1)[1].lstrip()


def _number(text: str) -> str:
    match = _NUMBER.fullmatch(text.strip())
    if match is None:
        raise WolframEvaluationError(f'Invalid Wolfram numeric representation: {text[:120]}')
    significand, exponent = match.groups()
    return significand + (('e' + exponent) if exponent else '')


def _normalize(row: object, *, digits: int, label: str) -> list[str]:
    if not isinstance(row, list):
        raise WolframReferenceError(f'Invalid result shape for {label}: {row!r}')
    values = cast(list[object], row)
    if len(values) != 4 or not all(isinstance(x, str) for x in values):
        raise WolframReferenceError(f'Invalid result shape for {label}: {row!r}')
    status, real, imaginary, precision = cast(list[str], values)
    if status != 'OK':
        raise WolframEvaluationError(f'{label}: WolframEngine returned {status}: {real[:160]}')
    a = _number(real)
    real_value = Decimal(a)
    zero_accuracy = _IMAGINARY_ACCURACY_ZERO.fullmatch(imaginary.strip())
    if zero_accuracy is not None:
        # An imaginary zero known to absolute accuracy a only certifies
        # a + floor(log10(abs(real))) relative decimal digits.  Do not let
        # this notation bypass the existing real-result precision guard.
        accuracy = Decimal(zero_accuracy.group(1))
        if (not real_value.is_finite() or not real_value
                or accuracy + real_value.adjusted() < digits - 3):
            raise WolframEvaluationError(
                f'{label}: insufficient imaginary-zero accuracy: {imaginary}')
        b = '0'
    else:
        b = _number(imaginary)
    if not real_value.is_finite() or not Decimal(b).is_finite():
        raise WolframEvaluationError(f'{label}: nonfinite numerical result')
    if Decimal(b) != 0:
        raise WolframEvaluationError(f'{label}: nonreal principal-branch value ({b}i)')
    if precision != 'Infinity':
        try:
            if Decimal(_number(precision)) < digits - 3:
                raise WolframEvaluationError(f'{label}: insufficient numerical precision: {precision}')
        except InvalidOperation as exc:
            raise WolframReferenceError(f'Invalid precision for {label}: {precision}') from exc
    # InputForm can expose guard digits beyond the backtick-tracked precision.
    # Publish only the requested significant decimal digits; otherwise those
    # internal mantissa digits would be mistaken for a validated reference.
    try:
        with localcontext() as context:
            context.prec = digits
            rounded = +Decimal(a)
    except InvalidOperation as exc:
        raise WolframEvaluationError(f'{label}: cannot round decimal reference') from exc
    return [str(rounded)]


def _batch_code(expressions: Sequence[str], *, digits: int) -> str:
    """Produce one bounded, static-expression, machine-readable WL program."""
    entries: list[str] = []
    for expression in expressions:
        entries.append(
            ''.join((
                f'Module[{{r=Quiet[Check[N[{expression},{digits + 20}],$Failed]]}},',
                'If[NumberQ[r]&&FreeQ[r,ComplexInfinity|Indeterminate|_DirectedInfinity],',
                f'If[Precision[r]<{digits - 3},{{"PRECISION",',
                'ToString[Precision[r],InputForm],"0","0"},',
                f'r=SetPrecision[r,{digits}];',
                '{"OK",ToString[Re[r],InputForm],ToString[Im[r],InputForm],',
                'ToString[Precision[r],InputForm]}],',
                '{"ERROR",ToString[r,InputForm],"0","0"}]]',
            )))
    return 'Print["' + _MARKER + '",ExportString[{' + ','.join(entries) + '},"RawJSON"]]'


def references(cases: Sequence[Mapping[str, object]], *, digits: int = 160,
               batch_size: int = 256, timeout: float = 90.0) -> list[list[str]]:
    """Evaluate curated scalar references in batches, preserving order.

    Unsupported signatures and evaluation failures raise explicit exceptions.
    Calling with an empty sequence never starts WolframEngine.
    """
    if isinstance(digits, bool) or not 20 <= digits <= 1000:
        raise ValueError('Reference precision must be 20..1000 decimal digits')
    if isinstance(batch_size, bool) or batch_size < 1:
        raise ValueError('batch_size must be positive')
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError('timeout must be positive and finite')
    expressions = [_expression(case) for case in cases]
    if not expressions:
        return []
    results: list[list[str]] = []
    offset = 0
    while offset < len(expressions):
        count = min(batch_size, len(expressions)-offset)
        # N works on exact rational inputs throughout: the extra digits are
        # guard digits, and the final N specifies the published decimal dps.
        # Quiet/Check isolates individual domain errors inside a batch.
        while True:
            batch = expressions[offset:offset+count]
            code = _batch_code(batch, digits=digits)
            if len(code.encode('utf-8')) < _MAX_ARGV:
                break
            if count == 1:
                raise WolframReferenceError('A single Wolfram reference exceeds subprocess argument limit')
            count = max(1, count // 2)
        raw = _execute(code, timeout=timeout)
        try:
            decoded, consumed = cast(tuple[object, int], json.JSONDecoder().raw_decode(raw))
        except (json.JSONDecodeError, ValueError) as exc:
            raise WolframReferenceError(f'Invalid WolframEngine JSON: {raw[:500]}') from exc
        if raw[consumed:].strip() not in ('', 'Null') or not isinstance(decoded, list):
            raise WolframReferenceError(f'Unexpected WolframEngine batch response: {raw[:500]}')
        rows = cast(list[object], decoded)
        if len(rows) != len(batch):
            raise WolframReferenceError(f'Unexpected WolframEngine batch length: {len(rows)}')
        for index, row in enumerate(rows):
            case = cases[offset + index]
            results.append(_normalize(row, digits=digits,
                                      label=str(case.get('id', case.get('api', offset+index)))))
        offset += count
    return results


def reference(case: Mapping[str, object], *, digits: int = 160) -> list[str]:
    """Evaluate one high-precision mathematical reference (decimal strings)."""
    return references([case], digits=digits)[0]


@lru_cache(maxsize=8)
def _version(executable: str) -> str:
    raw = _execute('Print["' + _MARKER + '",$Version]', timeout=30.0,
                   executable=executable)
    value = raw.splitlines()[0].strip()
    if not value or not re.match(r'^\d+\.\d+\.', value):
        raise WolframReferenceError(f'Unexpected Wolfram version: {value!r}')
    return value


def version() -> str:
    """Actual installed WolframEngine version; raises when unavailable."""
    return _version(_executable())


def provenance(*, digits: int = 160) -> dict[str, object]:
    """Reproducible evaluator metadata for the generated reference snapshot."""
    return {
        'provider': 'WolframEngine',
        'version': version(),
        'digits': digits,
        'float_input': 'exact IEEE 754 binary64 integer ratio',
        'integer_and_decimal_input': 'exact rational',
        'branches': 'Wolfram Language principal branches; real results required',
        'output': 'one full-precision decimal string per scalar result',
        'apis': sorted('special.' + key for key in _SPECS),
        'extra_apis': ['stats.normal_cdf', 'stats.variance',
                       'integrate.fixed_quad', 'differentiate.derivative'],
    }
