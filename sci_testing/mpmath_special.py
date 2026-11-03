"""Public mpmath comparisons for saved Sci special-function cases.

The historical ``p0``, ``p1``, ... inputs use the published mpmath functions
and analytic identities.  The Sci 0.3 typed and complex case generators already
define their independent mpmath mathematical functions and signed-zero branch
conventions, so those can be evaluated with the same precision here.

Only supported, shape-compatible, finite binary64 results are reported.  In
particular, polynomial *constructors*, distribution inverse-parameter solvers,
and Sci-specific multi-output zero searches are intentionally not claimed.
"""

import math
from functools import lru_cache

import mpmath as mp

PROVIDER = 'mpmath'
PRECISION = 80

# The direct public functions that match the saved real-input signatures.
_DIRECT = {
    'agm': mp.agm, 'beta': mp.beta, 'binom': mp.binomial,
    'comb': mp.binomial, 'digamma': mp.digamma, 'psi': mp.digamma,
    'trigamma': lambda x: mp.polygamma(1, x),
    'polygamma': mp.polygamma, 'gamma': mp.gamma, 'rgamma': mp.rgamma,
    'erfc': mp.erfc, 'erfi': mp.erfi, 'erfinv': mp.erfinv,
    'exp1': mp.e1, 'expi': mp.ei, 'expn': mp.expint,
    'log1p': mp.log1p, 'expm1': mp.expm1,
    'i0': lambda x: mp.besseli(0, x), 'i1': lambda x: mp.besseli(1, x),
    'j0': lambda x: mp.besselj(0, x), 'j1': lambda x: mp.besselj(1, x),
    'k0': lambda x: mp.besselk(0, x), 'k1': lambda x: mp.besselk(1, x),
    'y0': lambda x: mp.bessely(0, x), 'y1': lambda x: mp.bessely(1, x),
    'iv': mp.besseli, 'jv': mp.besselj, 'kv': mp.besselk,
    'yv': mp.bessely, 'jn': mp.besselj, 'kn': mp.besselk,
    'yn': mp.bessely, 'hyp0f1': mp.hyp0f1, 'hyp1f1': mp.hyp1f1,
    'hyp2f1': mp.hyp2f1, 'hyperu': mp.hyperu,
    'elliprc': mp.elliprc, 'elliprf': mp.elliprf,
    'elliprd': mp.elliprd, 'elliprj': mp.elliprj, 'elliprg': mp.elliprg,
    'ellipk': mp.ellipk, 'ellipe': mp.ellipe,
    'ellipf': mp.ellipf, 'ellipeinc': mp.ellipe,
    'ellipkinc': mp.ellipf, 'ellippi': mp.ellippi,
    'ellipkm1': lambda x: mp.ellipk(1 - x),
    'ellippiinc': lambda phi, n, m: mp.ellippi(n, phi, m),
    'ellippiinc_bulirsch': lambda phi, n, m: mp.ellippi(n, phi, m),
    'eval_gegenbauer': mp.gegenbauer,
    'eval_genlaguerre': mp.laguerre,
    'eval_jacobi': mp.jacobi,
    'eval_legendre': mp.legendre, 'legendre_p': mp.legendre,
    'eval_chebyt': mp.chebyt, 'eval_chebyu': mp.chebyu,
    'eval_hermite': mp.hermite,
    'eval_laguerre': lambda n, x: mp.laguerre(n, 0, x),
    'eval_sh_chebyt': lambda n, x: mp.chebyt(n, 2*x - 1),
    'eval_sh_chebyu': lambda n, x: mp.chebyu(n, 2*x - 1),
    'eval_sh_legendre': lambda n, x: mp.legendre(n, 2*x - 1),
    'hankel1': mp.hankel1, 'hankel2': mp.hankel2,
    'struve': mp.struveh, 'modstruve': mp.struvel,
    'poch': mp.rf, 'zeta': mp.zeta,
    'sindg': lambda x: mp.sin(x*mp.pi/180),
    'cosdg': lambda x: mp.cos(x*mp.pi/180),
    'tandg': lambda x: mp.tan(x*mp.pi/180),
    'cotdg': lambda x: mp.cot(x*mp.pi/180),
    'sinc': lambda x: mp.sinc(mp.pi*x),
    'spence': lambda x: mp.polylog(2, 1-x),
    'zetac': lambda x: mp.zeta(x) - 1,
    'harmonic': mp.harmonic,
    'lambertw': mp.lambertw,
    'lambertwm1': lambda x: mp.lambertw(x, -1),
    'wrightomega': lambda x: mp.lambertw(mp.exp(x)),
    'cbrt': lambda x: mp.sign(x) * mp.root(abs(x), 3),
    'exp10': lambda x: mp.power(10, x),
    'exp2': lambda x: mp.power(2, x),
    'softplus': lambda x: mp.log1p(mp.exp(x)),
    'exprel': lambda x: mp.expm1(x)/x if x else mp.mpf(1),
    'expit': lambda x: 1/(1 + mp.exp(-x)),
    'logit': lambda x: mp.log(x/(1-x)),
    'log_expit': lambda x: -mp.log1p(mp.exp(-x)),
    'log_ndtr': lambda x: mp.log(mp.erfc(-x/mp.sqrt(2))/2),
    'ndtr': lambda x: mp.erfc(-x/mp.sqrt(2))/2,
    'ndtri': lambda x: -mp.sqrt(2)*mp.erfinv(1-2*x),
    'powm1': lambda x, y: mp.expm1(y*mp.log(x)),
    'cosm1': lambda x: -2*mp.sin(x/2)**2,
    'boxcox': lambda x, a: mp.expm1(a*mp.log(x))/a if a else mp.log(x),
    'boxcox1p': lambda x, a: mp.expm1(a*mp.log1p(x))/a if a else mp.log1p(x),
    'inv_boxcox': lambda y, a: mp.exp(mp.log1p(a*y)/a) if a else mp.exp(y),
    'inv_boxcox1p': lambda y, a: mp.expm1(mp.log1p(a*y)/a) if a else mp.expm1(y),
    'xlogy': lambda x, y: x*mp.log(y) if x else mp.mpf(0),
    'xlog1py': lambda x, y: x*mp.log1p(y) if x else mp.mpf(0),
    'entr': lambda x: -x*mp.log(x) if x else mp.mpf(0),
    'rel_entr': lambda x, y: x*mp.log(x/y) if x else mp.mpf(0),
    'kl_div': lambda x, y: x*mp.log(x/y)-x+y if x else y,
    'huber': lambda delta, r: (r*r/2 if abs(r) <= delta
                               else delta*(abs(r)-delta/2)),
    'pseudo_huber': lambda delta, r: delta*delta*(mp.sqrt(1+(r/delta)**2)-1),
    'radian': lambda d, m, s: (d + m/60 + s/3600)*mp.pi/180,
    'gammaln': lambda x: mp.log(abs(mp.gamma(x))),
    'loggamma': lambda x: mp.loggamma(x),
    'gammasgn': lambda x: mp.sign(mp.gamma(x)),
    'betaln': lambda a, b: mp.log(abs(mp.beta(a, b))),
    'multigammaln': lambda a, d: d*(d-1)/4*mp.log(mp.pi)
                     + mp.fsum(mp.loggamma(a+(1-j)/2) for j in range(1, int(d)+1)),
    'factorial': lambda x: mp.gamma(x+1),
    'factorial2': lambda n: mp.fprod(range(int(n), 0, -2)),
    'factorialk': lambda n, k: mp.fprod(range(int(n), 0, -int(k))),
    'ln_factorial': lambda n: mp.loggamma(n+1),
    'ln_binomial': lambda n, k: mp.loggamma(n+1)-mp.loggamma(k+1)-mp.loggamma(n-k+1),
    'ndtri_exp': lambda x: -mp.sqrt(2)*mp.erfinv(1-2*mp.exp(x)),
    'wofz': lambda z: mp.exp(-z*z)*mp.erfc(-mp.j*z),
    'voigt_profile': lambda x, sigma, gamma: mp.re(
        mp.exp(-((x+mp.j*gamma)/(sigma*mp.sqrt(2)))**2)
        * mp.erfc(-mp.j*(x+mp.j*gamma)/(sigma*mp.sqrt(2))))
        / (sigma*mp.sqrt(2*mp.pi)),
    'erfcx': lambda x: mp.exp(x*x)*mp.erfc(x),
    'dawsn': lambda x: mp.sqrt(mp.pi)/2*mp.exp(-x*x)*mp.erfi(x),
    'assoc_laguerre': lambda x, n, alpha: mp.laguerre(n, alpha, x),
    'assoc_legendre_p': lambda n, m, x: mp.legenp(n, m, x),
    'lpmv': lambda m, n, x: mp.legenp(n, m, x),
    'eval_hermitenorm': lambda n, x: mp.hermite(n, x/mp.sqrt(2))/mp.power(2, n/2),
    'eval_chebyc': lambda n, x: 2*mp.chebyt(n, x/2),
    'eval_chebys': lambda n, x: mp.chebyu(n, x/2),
    'ber': lambda x: mp.ber(0,x), 'bei': lambda x: mp.bei(0,x),
    'ker': lambda x: mp.ker(0,x), 'kei': lambda x: mp.kei(0,x),
    'stirling2': lambda n, k: mp.fsum((-1)**(k-j)*mp.binomial(k,j)*j**n
                                          for j in range(int(k)+1))/mp.factorial(k),
}


def _ellipd(m):
    return (mp.ellipk(m)-mp.ellipe(m))/m if m else mp.pi/4


def _ellipdinc(phi, m):
    return (mp.ellipf(phi, m)-mp.ellipe(phi, m))/m if m else (phi-mp.sin(phi)*mp.cos(phi))/2


def _bessel_derivative(name, order, x, n):
    """Derivative from the public Bessel order-shift identities.

    A numeric mp.diff evaluates a large number of neighboring arguments,
    especially for complex Hankel functions.  At 80 dps the finite
    difference-in-order identity is exact for each supported family.
    """
    funcs = {'ivp': mp.besseli, 'jvp': mp.besselj, 'kvp': mp.besselk,
             'yvp': mp.bessely, 'h1vp': mp.hankel1, 'h2vp': mp.hankel2}
    n = int(n)
    if n > 24:
        return mp.diff(lambda t: funcs[name](order, t), x, n)
    sign = -1 if name in ('kvp',) and n % 2 else 1
    alternating = name in ('jvp', 'yvp', 'h1vp', 'h2vp')
    with mp.extradps(15 + n):
        result = mp.fsum(((-1 if alternating and k % 2 else 1)
                          * mp.binomial(n, k)
                          * funcs[name](order - n + 2*k, x))
                         for k in range(n+1))
        return +(sign*result / 2**n)


def _spherical(name, n, x):
    # All orders and arguments in the saved suite are real and nonsingular.
    funcs = {'spherical_jn': mp.besselj, 'spherical_yn': mp.bessely,
             'spherical_in': mp.besseli, 'spherical_kn': mp.besselk}
    return mp.sqrt(mp.pi/(2*x))*funcs[name](n+mp.mpf('.5'), x)


def _legendre_q(n, m, x):
    if abs(x)<1:
        return mp.legenq(n,m,x,type=2)
    # mpmath's type=2 includes an imaginary branch-cut contribution outside
    # [-1,1].  Sci's real Q sequences use the real type=3 continuation,
    # with the associated-order parity on the negative outer interval.
    return (-1 if x<0 and m%2 else 1)*mp.re(mp.legenq(n,m,x,type=3))


def _legendre_degrees(n, x, *, second_kind=False):
    """Base P_l or real Q_l values from the three-term degree recurrence."""
    if second_kind:
        q0 = mp.log((1+x)/(1-x))/2 if abs(x)<1 else mp.log((x+1)/(x-1))/2
        result = [q0]
        if n:
            result.append(x*q0-1)
    else:
        result = [mp.mpf(1)]
        if n:
            result.append(x)
    for l in range(1, n):
        result.append(((2*l+1)*x*result[l] - l*result[l-1])/(l+1))
    return result


def _associated_legendre_p_all(n, m, x):
    """Ferrers P through degree n and orders 0..m, -m..-1.

    Every positive order is generated from the previous two orders in
    O(n*m) arithmetic, avoiding the hypergeometric call in every legenp.
    """
    n, m = int(n), int(m)
    if not -1 < x < 1:
        # Near/at the branch boundary preserve public mpmath conventions.
        orders = [*range(m+1), *range(-m, 0)]
        return [[mp.legenp(k,j,x) if abs(j)<=k else mp.mpf(0)
                 for j in orders] for k in range(n+1)]
    with mp.extradps(20):
        s = mp.sqrt((1-x)*(1+x))
        positives = [_legendre_degrees(n, x)]
        if m:
            positives.append([mp.mpf(0)] + [
                -l*(positives[0][l-1]-x*positives[0][l])/s
                for l in range(1,n+1)])
        for j in range(1,m):
            positives.append([
                -2*j*x*positives[j][l]/s
                -(l+j)*(l-j+1)*positives[j-1][l]
                if l>j else mp.mpf(0)
                for l in range(n+1)])
        result = []
        for l in range(n+1):
            row = [+(positives[j][l]) for j in range(m+1)]
            negative = []
            ratio = mp.mpf(1)
            for j in range(1,min(l,m)+1):
                ratio = -ratio/((l+j)*(l-j+1))
                negative.append(+(ratio*positives[j][l]))
            row.extend([mp.mpf(0)]*(m-l) + negative[::-1])
            result.append(row)
        return result


def _associated_legendre_q_all(m, n, x):
    """Real associated Legendre Q and x-derivatives, order-major.

    The two-order recurrence computes all Q_l^j in O(n*m) at high
    precision.  The raising identity supplies derivatives algebraically:
    Q'_l^j = (j*x*Q_l^j/s + Q_l^(j+1))/s times -1 (inside)
    or +1 (outside), with s = sqrt(abs(1-x*x)).
    """
    m, n = int(m), int(n)
    if abs(x)==1:
        # Singular endpoints are outside the saved finite measurement domain.
        raise ValueError('Legendre Q singular at x=+/-1')
    extra = 20 + max(m,n)//2
    if abs(x)>1:
        # The degree recurrence has a growing and a decaying solution.
        # When computing Q_n, roundoff in the growing solution is amplified
        # by (|x|+sqrt(x*x-1))**(2*n).  Supply those extra decimal digits
        # before the recurrence, then round back to the caller's 80 dps.
        growth = abs(x)+mp.sqrt(x*x-1)
        extra += math.ceil((2*n+2)*float(mp.log10(growth))) + 2*m
    with mp.extradps(extra):
        inside = abs(x) < 1
        s = mp.sqrt((1-x)*(1+x) if inside else (x-1)*(x+1))
        orders = [_legendre_degrees(n,x,second_kind=True)]
        orders.append([-1/s] + [
            -l*(orders[0][l-1]-x*orders[0][l])/s
            for l in range(1,n+1)])
        for j in range(1,m+1):
            coefficient = -1 if inside else 1
            orders.append([
                -2*j*x*orders[j][l]/s
                +coefficient*(l+j)*(l-j+1)*orders[j-1][l]
                for l in range(n+1)])
        derivatives = [[
            (-1 if inside else 1) *
            (orders[j+1][l]+j*x*orders[j][l]/s)/s
            for l in range(n+1)]
            for j in range(m+1)]
        return ([[+v for v in row] for row in orders[:m+1]],
                [[+v for v in row] for row in derivatives])


def _quadrature(name, values):
    n=int(values[0])
    # scipy.special roots_* return the concatenation of nodes and weights.
    rules={
        'roots_legendre':'legendre','roots_sh_legendre':'legendre',
        'roots_chebyt':'chebyshev1','roots_chebyc':'chebyshev1',
        'roots_sh_chebyt':'chebyshev1','roots_chebyu':'chebyshev2',
        'roots_chebys':'chebyshev2','roots_sh_chebyu':'chebyshev2',
        'roots_hermite':'hermite','roots_hermitenorm':'hermite',
        'roots_laguerre':'laguerre','roots_genlaguerre':'glaguerre',
        'roots_jacobi':'jacobi','roots_gegenbauer':'jacobi',
    }
    alpha=mp.mpf(0)
    beta=mp.mpf(0)
    if name=='roots_genlaguerre':alpha=values[1]
    if name=='roots_jacobi':alpha,beta=values[1:]
    if name=='roots_gegenbauer':alpha=beta=values[1]-mp.mpf('.5')
    nodes,weights=mp.gauss_quadrature(n,rules[name],alpha=alpha,beta=beta)
    nodes,weights=list(nodes),list(weights)
    if name in ('roots_chebyc','roots_chebys'):
        nodes=[x*2 for x in nodes]
        weights=[w*2 for w in weights]
    if name in ('roots_sh_chebyt','roots_sh_chebyu','roots_sh_legendre'):
        nodes=[(x+1)/2 for x in nodes]
        if name=='roots_sh_chebyu': weights=[w/4 for w in weights]
        if name=='roots_sh_legendre': weights=[w/2 for w in weights]
    if name=='roots_hermitenorm':
        nodes=[x*mp.sqrt(2) for x in nodes]
        weights=[w*mp.sqrt(2) for w in weights]
    # Symmetric weights have exactly zero center nodes.  Asymmetric Jacobi
    # weights and Laguerre weights do not, even when their degree is odd.
    symmetric=name not in ('roots_laguerre','roots_genlaguerre') and (
        name!='roots_jacobi' or alpha==beta)
    if n%2 and symmetric:
        nodes[n//2]=mp.mpf('.5') if name.startswith('roots_sh_') else mp.mpf(0)
    return [*nodes,*weights]


def _historic(name, v):
    if name in _DIRECT:
        return _DIRECT[name](*v)
    if name.startswith('roots_'):
        return _quadrature(name,v)
    if name in ('jn_zeros','jnp_zeros','yn_zeros','ynp_zeros'):
        order,n=v
        fn=mp.besseljzero if name.startswith('j') else mp.besselyzero
        derivative=int('p' in name)
        # mpmath counts J_0'(0)=0 as its first derivative root; SciPy's
        # jnp_zeros returns only strictly positive zeros.
        first=2 if name=='jnp_zeros' and order==0 else 1
        return [fn(order,k,derivative=derivative)
                for k in range(first,first+int(n))]
    if name=='jnyn_zeros':
        order,n=v
        return [[fn(order,k,derivative=d) for k in range(1,int(n)+1)]
                for fn,d in ((mp.besseljzero,0),(mp.besseljzero,1),
                             (mp.besselyzero,0),(mp.besselyzero,1))]
    if name in ('ai_zeros','bi_zeros'):
        n=int(v[0]);fn=mp.airyaizero if name=='ai_zeros' else mp.airybizero
        airy=mp.airyai if name=='ai_zeros' else mp.airybi
        zeros=[fn(k) for k in range(1,n+1)]
        derivatives=[fn(k,derivative=1) for k in range(1,n+1)]
        return [zeros,derivatives,[airy(x) for x in derivatives],
                [airy(x,1) for x in zeros]]
    if name=='euler':return [mp.eulernum(k) for k in range(int(v[0])+1)]
    if name in ('berp','beip','kerp','keip'):
        fn={'berp':mp.ber,'beip':mp.bei,'kerp':mp.ker,'keip':mp.kei}[name]
        return mp.diff(lambda x:fn(0,x),v[0])
    if name=='kelvin':
        x=v[0]
        return [mp.ber(0,x)+mp.j*mp.bei(0,x),
                mp.ker(0,x)+mp.j*mp.kei(0,x),
                mp.diff(lambda t:mp.ber(0,t)+mp.j*mp.bei(0,t),x),
                mp.diff(lambda t:mp.ker(0,t)+mp.j*mp.kei(0,t),x)]
    if name in ('ivp','jvp','kvp','yvp','h1vp','h2vp'):
        return _bessel_derivative(name, *v)
    if name in ('spherical_jn','spherical_yn','spherical_in','spherical_kn'):
        n,x,*optional = v
        fn = lambda t: _spherical(name, n, t)
        return mp.diff(fn,x) if optional and optional[0] else fn(x)
    if name in ('i0e','i1e','k0e','k1e','jve','yve','ive','kve',
                'hankel1e','hankel2e'):
        if name in ('i0e','i1e','k0e','k1e'):
            n,x = (0 if name[1]=='0' else 1), v[0]
            return (mp.besseli(n,x)*mp.exp(-abs(x)) if name[0]=='i'
                    else mp.besselk(n,x)*mp.exp(x))
        nu,x = v
        if name == 'ive': return mp.besseli(nu,x)*mp.exp(-abs(x))
        if name == 'kve': return mp.besselk(nu,x)*mp.exp(x)
        if name == 'jve': return mp.besselj(nu,x)*mp.exp(-abs(mp.im(x)))
        if name == 'yve': return mp.bessely(nu,x)*mp.exp(-abs(mp.im(x)))
        if name == 'hankel1e': return mp.hankel1(nu,x)*mp.exp(-mp.j*x)
        return mp.hankel2(nu,x)*mp.exp(mp.j*x)
    if name == 'ellipke': return [mp.ellipk(v[0]),mp.ellipe(v[0])]
    if name == 'ellipd': return _ellipd(*v)
    if name == 'ellipdinc': return _ellipdinc(*v)
    if name == 'ellipj':
        u,m=v
        sn=mp.ellipfun('sn',u,m)
        # Amplitude is continuous across the quarter-periods; asin(sn)
        # alone loses the accumulated phase beyond [-K(m), K(m)].
        periods=int(mp.floor((u+mp.ellipk(m))/(2*mp.ellipk(m))))
        amplitude=periods*mp.pi+(-1)**periods*mp.asin(sn)
        return [sn, mp.ellipfun('cn',u,m), mp.ellipfun('dn',u,m),amplitude]
    if name == 'legendre_p_all':
        n,x=v
        return _legendre_degrees(int(n),x)
    if name == 'assoc_legendre_p_all':
        n,m,x=v
        return _associated_legendre_p_all(n,m,x)
    if name == 'sph_harm_y_all':
        n,m,theta,phi=v
        orders=[*range(int(m)+1),*range(-int(m),0)]
        return [[mp.spherharm(k,j,theta,phi) if abs(j)<=k else mp.mpc(0)
                 for j in orders] for k in range(int(n)+1)]
    if name == 'lqn':
        n,x=v
        q,dq=_associated_legendre_q_all(0,n,x)
        return [q[0],dq[0]]
    if name == 'lqmn':
        m,n,x=v
        return _associated_legendre_q_all(m,n,x)
    if name == 'betainc': return mp.betainc(v[0],v[1],0,v[2],regularized=True)
    if name == 'betaincc': return mp.betainc(v[0],v[1],v[2],1,regularized=True)
    if name == 'gammainc': return mp.gammainc(v[0],0,v[1],regularized=True)
    if name == 'gammaincc': return mp.gammainc(v[0],v[1],mp.inf,regularized=True)
    if name == 'fresnel': return [mp.fresnels(v[0]),mp.fresnelc(v[0])]
    if name == 'sici': return [mp.si(v[0]),mp.ci(v[0])]
    if name == 'shichi': return [mp.shi(v[0]),mp.chi(v[0])]
    if name in ('airy','airye'):
        x=v[0]
        result=[mp.airyai(x),mp.airyai(x,1),mp.airybi(x),mp.airybi(x,1)]
        if name == 'airye':
            phase=mp.mpf(2)/3*x*mp.sqrt(x)
            return [result[0]*mp.exp(phase),result[1]*mp.exp(phase),
                    result[2]*mp.exp(-abs(mp.re(phase))),
                    result[3]*mp.exp(-abs(mp.re(phase)))]
        return result
    if name == 'gen_harmonic':
        n,a=v
        if n>10000: raise ValueError('large harmonic sums are not sampled')
        return mp.fsum(mp.power(k,-a) for k in range(1,int(n)+1))
    if name == 'multinomial':
        n,counts=v
        return mp.exp(mp.loggamma(n+1)-mp.fsum(mp.loggamma(c+1) for c in counts))
    if name == 'ln_multinomial':
        n,counts=v
        return mp.loggamma(n+1)-mp.fsum(mp.loggamma(c+1) for c in counts)
    if name == 'logsumexp':
        xs=v[0]
        a=max(xs)
        return a+mp.log(mp.fsum(mp.exp(x-a) for x in xs))
    if name == 'softmax':
        xs=v[0]; a=max(xs); weights=[mp.exp(x-a) for x in xs]
        s=mp.fsum(weights)
        return [x/s for x in weights]
    if name == 'log_softmax':
        xs=v[0]; a=max(xs)
        scale=a+mp.log(mp.fsum(mp.exp(x-a) for x in xs))
        return [x-scale for x in xs]
    if name == 'bernoulli':
        return [mp.bernoulli(i) for i in range(int(v[0])+1)]
    if name in ('sph_harm_y',):
        n,m,theta,phi=v
        return mp.spherharm(n,m,theta,phi)
    if name in ('bdtr','bdtrc'):
        k,n,p=v
        result=mp.betainc(n-k,k+1,0,1-p,regularized=True)
        return 1-result if name=='bdtrc' else result
    if name in ('nbdtr','nbdtrc'):
        k,n,p=v
        result=mp.betainc(n,k+1,0,p,regularized=True)
        return 1-result if name=='nbdtrc' else result
    if name in ('chdtr','chdtrc'):
        df,x=v
        return mp.gammainc(df/2,0 if name=='chdtr' else x/2,
                           x/2 if name=='chdtr' else mp.inf,regularized=True)
    if name in ('gdtr','gdtrc'):
        a,b,x=v
        return mp.gammainc(b,0 if name=='gdtr' else a*x,
                           a*x if name=='gdtr' else mp.inf,regularized=True)
    if name in ('fdtr','fdtrc'):
        a,b,x=v
        z=a*x/(b+a*x)
        return mp.betainc(a/2,b/2,0 if name=='fdtr' else z,
                          z if name=='fdtr' else 1,regularized=True)
    if name == 'pdtr': return mp.gammainc(v[0]+1,v[1],mp.inf,regularized=True)
    if name == 'pdtrc': return mp.gammainc(v[0]+1,0,v[1],regularized=True)
    if name == 'stdtr':
        df,t=v
        u=df/(df+t*t)
        h=mp.betainc(df/2,mp.mpf('.5'),0,u,regularized=True)/2
        return h if t<0 else 1-h
    if name == 'zeta' and len(v)==2: return mp.zeta(v[0],v[1])
    raise LookupError(name)


@lru_cache(maxsize=512)
def _discrete_sequence(name, values):
    """Reuse expensive high-precision roots for repeated saved integer inputs.

    All calls originate inside evaluate()'s fixed 80-digit precision context.
    A bounded cache keeps the results stable and avoids repeated root searches
    when 2,000 saved cases contain only a few distinct counts and orders.
    """
    return _historic(name, list(values))


def _historic_evaluation(name, values):
    if (name in {'ai_zeros', 'bi_zeros', 'jn_zeros', 'jnp_zeros',
                 'yn_zeros', 'ynp_zeros', 'jnyn_zeros', 'euler'} or
            name.startswith('roots_')):
        return _discrete_sequence(name, tuple(values))
    return _historic(name, values)


# These are additional historical signatures with verified formulas above.
_COMPOSITE = frozenset({
    'airy','airye','ellipke','ellipd','ellipdinc','ellipj',
    'euler','ai_zeros','bi_zeros','jn_zeros','jnp_zeros','yn_zeros',
    'ynp_zeros','jnyn_zeros','berp','beip','kerp','keip','kelvin',
    'roots_legendre','roots_sh_legendre','roots_chebyt','roots_chebyc',
    'roots_sh_chebyt','roots_chebyu','roots_chebys','roots_sh_chebyu',
    'roots_hermite','roots_hermitenorm','roots_laguerre',
    'roots_genlaguerre','roots_jacobi','roots_gegenbauer',
    'legendre_p_all','assoc_legendre_p_all','sph_harm_y_all','lqn','lqmn',
    'betainc','betaincc','gammainc','gammaincc','fresnel','sici','shichi',
    'gen_harmonic','multinomial','ln_multinomial','logsumexp','softmax',
    'log_softmax','bernoulli','sph_harm_y',
    'bdtr','bdtrc','nbdtr','nbdtrc','chdtr','chdtrc','gdtr','gdtrc',
    'fdtr','fdtrc','pdtr','pdtrc','stdtr',
    'spherical_jn','spherical_yn','spherical_in','spherical_kn',
    'ivp','jvp','kvp','yvp','h1vp','h2vp',
    'i0e','i1e','k0e','k1e','ive','jve','kve','yve','hankel1e','hankel2e',
})

SUPPORTED_NAMES = frozenset(_DIRECT) | _COMPOSITE


def _number(value):
    if isinstance(value, bool): return value
    if isinstance(value, (list, tuple)):
        if len(value)==2 and all(isinstance(x,(int,float)) for x in value):
            return mp.mpc(*[mp.mpf(x) for x in value])
        return [_number(v) for v in value]
    return mp.mpf(value)


_HISTORICAL_COMPLEX = frozenset({
    'hankel1','hankel2','hankel1e','hankel2e',
    'h1vp','h2vp',
    'sph_harm_y','sph_harm_y_all','wofz','kelvin',
})


def _flatten(value, *, complex_output=False):
    if isinstance(value,(list,tuple)):
        return [part for x in value for part in _flatten(x,complex_output=complex_output)]
    if isinstance(value, mp.mpc):
        if complex_output:
            return [float(mp.re(value)),float(mp.im(value))]
        if mp.im(value):
            raise ValueError('complex result for real Sci function')
        return [float(mp.re(value))]
    if mp.im(value):
        # A real-returning Sci signature cannot accept a complex continuation.
        raise ValueError('nonreal mpmath result for real Sci function')
    return [float(mp.re(value))]


def evaluate(case):
    """Return ``{'mpmath': {'values': [...]}}``, or empty when unavailable.

    Use mp.workdps so high-precision intermediate results never depend on the
    report's current process-wide precision. Conversion to floats happens only
    once, after mathematical evaluation and output-shape normalization.
    """
    api = case.get('api','')
    if not api.startswith('special.'):
        return {}
    args = case.get('args',{})
    with mp.workdps(PRECISION):
        try:
            if args.get('special03_case') is True:
                from .special03_cases import APIS, reference
                if api not in APIS: return {}
                result = reference(case)
            elif args.get('special_complex03_case') is True:
                from .special_complex03_cases import APIS, reference
                if api not in APIS: return {}
                result = reference(case)
            elif api in ('special.erf', 'special.gamma') and set(args)=={'x'}:
                # The original dense Sci 0.1 recipes have named input x.
                # Later historical registrations use positional p0.
                result = (mp.erf if api=='special.erf' else mp.gamma)(
                    mp.mpf(args['x']))
            else:
                name = api.rsplit('.',1)[-1]
                if name not in SUPPORTED_NAMES or set(args) != {
                    f'p{i}' for i in range(len(args))
                }:
                    return {}
                v=[_number(args[f'p{i}']) for i in range(len(args))]
                # Arrays such as softmax/logsumexp are lists of scalar numbers,
                # not mpmath complex pairs.
                if name in ('softmax','log_softmax','logsumexp',
                            'multinomial','ln_multinomial'):
                    v=[[mp.mpf(t) for t in x] if isinstance(x,list) else
                       mp.mpf(x) for x in (args[f'p{i}'] for i in range(len(args)))]
                # The two-argument Hurwitz zeta needs the same name as the
                # one-argument Riemann counterpart.
                if name=='zeta' and len(v)==2:
                    result=mp.zeta(*v)
                else:
                    result=_historic_evaluation(name,v)
            values=_flatten(result,complex_output=(not args.get('special03_case')
                                                  and not args.get('special_complex03_case')
                                                  and api.rsplit('.',1)[-1] in _HISTORICAL_COMPLEX))
            if not values or not all(math.isfinite(x) for x in values):
                return {}
            return {PROVIDER:{'values':values}}
        except (ValueError,TypeError,ZeroDivisionError,OverflowError,
                ArithmeticError,mp.libmp.libhyper.NoConvergence):
            return {}
