"""Deterministic finite-input accuracy cases; no numerical reference dependencies.

Retain the original hand-picked cases, then add domain grids, adjacent floating
point values, scaled problems, and samples from a fixed seed. This file generates
inputs only. Every generated problem has bounded dimensions and work.
"""
from collections import Counter
import hashlib
import numpy as np
import math
import random

from .case_csv import case_files, read_cases, write_files, write_cases

from . import ROOT
CASES_PER_METHOD = 1000
SEED = 20261020


def build_cases():
    cases = read_cases(ROOT / 'cases/base')
    counters = Counter()
    rng = random.Random(SEED)

    def add(api, group, **args):
        i = counters[api]
        counters[api] += 1
        cases.append(dict(id=f'{api}/expanded-{i:04d}', api=api, group=group, args=args))

    def sample(n, scale=1.0):
        return [rng.uniform(-2, 2) * scale for _ in range(n)]

    # Elementary complex functions: all quadrants, large magnitudes, branch sides.
    for api in ['complex.exp', 'complex.sqrt']:
        for _ in range(48):
            add(api, 'seeded quadrants', x=rng.uniform(-20, 20), y=rng.uniform(-20, 20))
        for exponent in [-100, -20, 0, 20, 100]:
            scale = 10.0**exponent
            if api.endswith('sqrt'):
                for sign in [-1, 1]:
                    add(api, 'scaled quadrants', x=-scale, y=sign * scale / 8)
        if api.endswith('exp'):
            for x in [-700., -100., -1., 1., 100., 700.]:
                for y in [-100., .25, 100.]:
                    add(api, 'finite exponential range', x=x, y=y)
        else:
            for x in [-10., -1., -.125]:
                for y in [-1e-14, 1e-14]:
                    add(api, 'branch-cut approaches', x=x, y=y)

    for api in ['constants.lambda2nu', 'constants.nu2lambda']:
        for exponent in range(-280, 281, 10):
            for sign in [-1, 1]:
                add(api, 'signed logarithmic range', x=sign * 1.25 * 10.0**exponent)
        for x in [1., 299792458., 1e14]:
            for direction in [0., math.inf]:
                add(api, 'adjacent representable inputs', x=math.nextafter(x, direction))

    for step in [1e-2, 1e-5, 1e-7]:
        for i in range(-32, 33):
            add('differentiate.derivative', f'sin grid, h = {step:g}', x=i / 4, step=step)
        for _ in range(24):
            add('differentiate.gradient', f'quadratic samples, h = {step:g}', point=sample(2), step=step)

    # Powers of two, composites, primes, and exact-zero spectra. Use small dyadic
    # data so fixture generation is independent of host transcendental functions.
    for api in ['fft.fft', 'fft.ifft']:
        for n in [1, 2, 3, 5, 6, 8, 9, 11, 12, 17, 32, 48]:
            for kind in ['zero', 'impulse', 'seeded dyadic']:
                if kind == 'zero': z = [[0., 0.] for _ in range(n)]
                elif kind == 'impulse': z = [[1., -.5]] + [[0., 0.] for _ in range(n - 1)]
                else: z = [[rng.randrange(-32, 33) / 16, rng.randrange(-32, 33) / 16] for _ in range(n)]
                add(api, f'{kind}, N = {n}', z=z)

    # SPD two-state priors constructed analytically; observations stay scalar.
    for _ in range(64):
        scale = 10.0**rng.choice([-8, -3, 0, 3, 8])
        rho = rng.uniform(-.8, .8)
        prior = dict(mean=sample(2), covariance=[[2*scale, rho*scale], [rho*scale, scale]])
        add('filter.kalman_predict', 'seeded SPD priors', **prior,
            f=[[1., rng.uniform(-1, 1)], [0., 1.]], q=[[scale/10, 0.], [0., scale/5]])
        add('filter.kalman_update', 'seeded SPD observations', **prior,
            z=sample(1), h=[[1., rng.uniform(-1, 1)]], r=[[scale*rng.uniform(.1, 2)]])

    for _ in range(64):
        n = rng.choice([2, 3, 5, 9, 17, 33])
        coordinates = sorted(sample(n))
        if rng.randrange(2): coordinates.reverse()
        if rng.randrange(4) == 0: coordinates[1] = coordinates[0]
        add('integrate.trapezoid', 'uneven, reversed or repeated coordinates',
            coordinates=coordinates, y=sample(n, 10.0**rng.choice([-100, 0, 100])))
        a = rng.uniform(-10, 10)
        b = a + rng.choice([-1., 0., 1.]) * 10.0**rng.choice([-6, -2, 0])
        add('integrate.fixed_quad', 'short, reversed or zero-width intervals', a=a, b=b)

    # Dyadic integer nodes make polynomial samples exact in binary64. Distinct
    # sample values, rather than a rounded transcendental source, define the data.
    for api in ['interpolate.barycentric_interpolate', 'interpolate.krogh_interpolate']:
        for _ in range(48):
            nodes = sorted(rng.sample(range(-32, 33), rng.choice([5, 6, 8, 10])))
            nodes = [x/8 for x in nodes]
            query = [nodes[0], nodes[-1], *sample(6), nodes[0]-.25, nodes[-1]+.25]
            add(api, 'dyadic quartic grids and extrapolation', nodes=nodes,
                values=[x**4 - 2*x + 3 for x in nodes], query=query)

    for n in [1, 2, 3, 4, 6, 8]:
        for _ in range(8):
            a = [sample(n) for _ in range(n)]
            for i in range(n): a[i][i] += 4*n
            add('linalg.solve', f'seeded dominant matrices, N = {n}', a=a, b=sample(n))
            add('linalg.det', f'seeded dominant matrices, N = {n}', a=a)
        # Hilbert matrices are a separate conditioning group, without singular
        # inputs or tolerances that hide the expected sensitivity.
        a = [[1/(i+j+1) for j in range(n)] for i in range(n)]
        add('linalg.solve', f'Hilbert matrix, N = {n}', a=a, b=[1.]*n)
        add('linalg.det', f'Hilbert matrix, N = {n}', a=a)

    for api in ['optimize.rosen', 'optimize.rosen_der']:
        for n in [2, 3, 5, 10, 32]:
            for _ in range(8): add(api, f'seeded objectives, N = {n}', point=sample(n))
            for displacement in [0., 1e-8, -1e-8]:
                add(api, f'near minimum, N = {n}', point=[1.+displacement*(i+1) for i in range(n)])

    for api in ['signal.convolve', 'signal.correlate']:
        for n in [1, 2, 7, 16, 33, 64]:
            for _ in range(8):
                scale = 10.0**rng.choice([-100, 0, 100])
                add(api, f'scaled sequences, N = {n}', x=sample(n, scale), h=sample(rng.choice([1, 2, 5, 9])))
            add(api, f'zero sequence, N = {n}', x=[0.]*n, h=[1., -1.])

    for n in [2, 3, 8, 32]:
        for _ in range(12):
            scale = 10.0**rng.choice([-50, 0, 50])
            u, v = sample(n, scale), sample(n, scale)
            for api in ['spatial.euclidean', 'spatial.cosine']:
                add(api, f'seeded vectors, N = {n}', u=u, v=v)
        for delta in [0., 1e-4, 1e-8]:
            u = [float(i+1) for i in range(n)]
            v = u.copy(); v[0] += delta
            add('spatial.cosine', f'nearly parallel, N = {n}', u=u, v=v)
            add('spatial.euclidean', f'nearly equal, N = {n}', u=u, v=v)

    for i in range(-64, 65): add('special.erf', 'dense grid [-8, 8]', x=i/8)
    for x in [-1e-300, -1e-100, 1e-100, 1e-300]: add('special.erf', 'tiny magnitudes', x=x)
    for _ in range(64): add('special.gamma', 'positive logarithmic samples', x=10.0**rng.uniform(-6, 2.2))
    for pole in [-1., -2., -10., -30.]:
        for direction in [-math.inf, math.inf]:
            add('special.gamma', 'adjacent to negative integer poles', x=math.nextafter(pole, direction))
    for x in [1e-250, 1e-100, 170., 171., 171.5]: add('special.gamma', 'finite range extremes', x=x)

    for i in range(-76, 25): add('stats.normal_cdf', 'central and tail grid [-38, 12]', x=i/2)
    for n in [2, 3, 8, 32, 64]:
        for _ in range(12):
            scale = 10.0**rng.choice([-150, -50, 0, 50, 150])
            add('stats.variance', f'scaled samples, N = {n}', values=sample(n, scale))
        for offset in [0., 1e8, 1e12]:
            add('stats.variance', f'offset samples, N = {n}', values=[offset+rng.uniform(-1, 1) for _ in range(n)])

    for n in [3, 5, 8, 16]:
        for _ in range(12):
            scale = 10.0**rng.choice([-6, 0, 6])
            obs = [sample(2, scale) for _ in range(n)]
            add('cluster.vq.vq', f'seeded codebooks, N = {n}', obs=obs, codebook=[sample(2, scale) for _ in range(3)])
            add('cluster.vq.whiten', f'seeded whitening, N = {n}', obs=obs)
    # Exact first-minimum tie behavior and exact codebook matches are also tested.
    add('cluster.vq.vq', 'ties and exact matches', obs=[[0., 0.], [-1., 0.], [1., 0.]], codebook=[[-1., 0.], [1., 0.]])
    cases.extend(dense_cases(sorted({case['api'] for case in cases})))
    from .case_set import CaseSet
    return CaseSet(cases,CASES_PER_METHOD,SEED)


def dense_arguments(api, i, value):
    """Map normalized coordinates into bounded, valid numerical problems.

    Grid coordinates are permutations of the same endpoint-inclusive linspace;
    random coordinates are independent uniform samples. Sizes cycle separately
    so both methods exercise each supported shape with distinct numerical data.
    """
    def vector(n, offset=0): return [2*value(offset+j) for j in range(n)]
    n = (2, 3, 5, 8)[i % 4]
    t = value(0)
    if api.startswith('complex.'): return dict(x=20*t, y=20*value(1))
    if api.startswith('constants.'): return dict(x=1e-6+(t+1)/2*1e6)
    if api == 'differentiate.derivative': return dict(x=8*t, step=(1e-2, 1e-5, 1e-7)[i % 3])
    if api == 'differentiate.gradient': return dict(point=vector(2), step=(1e-2, 1e-5, 1e-7)[i % 3])
    if api.startswith('fft.'): return dict(z=[vector(2, 2*j) for j in range(n)])
    if api.startswith('filter.'):
        scale = 10**(3*t)
        rho = .8*value(1)
        args = dict(mean=vector(2, 2), covariance=[[2*scale, rho*scale], [rho*scale, scale]])
        if api.endswith('predict'):
            args.update(f=[[1., value(4)], [0., 1.]], q=[[scale/10, 0.], [0., scale/5]])
        else:
            args.update(z=vector(1, 4), h=[[1., value(5)]], r=[[scale*(1.1+value(6))]])
        return args
    if api == 'integrate.trapezoid':
        return dict(coordinates=sorted(vector(n)), y=vector(n, n))
    if api == 'integrate.fixed_quad': return dict(a=5*t, b=5*value(1))
    if api.startswith('interpolate.'):
        # Dyadic nodes and values exactly represent the quartic used by the reference.
        nodes = [-2., -1., 0., 1., 2.]
        return dict(nodes=nodes, values=[x**4-2*x+3 for x in nodes], query=vector(n))
    if api.startswith('linalg.'):
        a = [vector(n, n*j) for j in range(n)]
        for j in range(n): a[j][j] += 4*n
        return dict(a=a, b=vector(n, n*n)) if api.endswith('solve') else dict(a=a)
    if api.startswith('optimize.'): return dict(point=vector(n))
    if api.startswith('signal.'): return dict(x=vector(n), h=vector(3, n))
    if api.startswith('spatial.'): return dict(u=vector(n), v=vector(n, n))
    if api == 'special.erf': return dict(x=8*t)
    if api == 'special.gamma': return dict(x=.05+(t+1)/2*169.95)
    if api == 'stats.normal_cdf': return dict(x=-38+(t+1)*25)
    if api == 'stats.variance': return dict(values=vector(n))
    if api.startswith('cluster.'):
        args = dict(obs=[vector(2, 2*j) for j in range(n)])
        if api.endswith('vq'): args['codebook'] = [vector(2, 2*n+2*j) for j in range(3)]
        return args
    raise ValueError(f'Missing dense recipe: {api}')


def dense_cases(apis,count=CASES_PER_METHOD):
    grid = np.linspace(-1., 1., count).tolist()
    for api in apis:
        # Independent streams prevent additions or reordering from changing other APIs.
        seed = int.from_bytes(hashlib.sha256(f'{SEED}:{api}'.encode()).digest(), 'big')
        rng = random.Random(seed)
        for method in ['linspace', 'random']:
            for i in range(count):
                coordinates = {}
                def value(j):
                    if j not in coordinates:
                        coordinates[j] = grid[(i+37*j) % count] if method == 'linspace' else rng.uniform(-1, 1)
                    return coordinates[j]
                yield dict(id=f'{api}/{method}-{i:04d}', api=api,
                           group=f'{method}, {count} cases', args=dense_arguments(api, i, value))


if __name__ == '__main__':
    cases = build_cases()
    write_cases(ROOT / 'cases/generated', cases)
    print(f'Generated {len(cases)} cases with seed {SEED}')
