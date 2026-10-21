"""Bounded callback, metadata, and small-rule recipes with independent oracles."""
import hashlib
import random
import numpy as np
from scipy import constants, integrate, interpolate, stats, cluster

APIS='''linalg.dot linalg.matvec linalg.matmul linalg.transpose special.geterr stats.normal_pdf stats.normal_ppf constants.value constants.unit
constants.precision constants.find constants.convert_temperature differentiate.jacobian
differentiate.hessian integrate.simpson integrate.cumulative_trapezoid
integrate.cumulative_simpson integrate.romb integrate.newton_cotes
interpolate.pchip_interpolate interpolate.lagrange interpolate.pade
cluster.vq.kmeans cluster.vq.kmeans2 optimize.minimize optimize.try_minimize
optimize.minimize_with_gradient optimize.try_minimize_with_gradient'''.split()


def cases(count,seed,apis=None):
    grid=np.linspace(-1.,1.,count)
    keys=sorted(constants.physical_constants)
    for api in sorted(set(APIS) if apis is None else set(apis)&set(APIS)):
        rng=random.Random(int.from_bytes(hashlib.sha256(f'{seed}:{api}'.encode()).digest(),'big'))
        for method in ['linspace','random']:
            for i in range(count):
                t=float(grid[i]) if method=='linspace' else rng.uniform(-1,1)
                args={'explicit':True}
                if api.startswith('linalg.'):
                    args.update(a=[[1.+t,.5],[.25,2.]],b=[[2.,1.],[0.,1.+t/2]],u=[t,1.],v=[2.,1.+t])
                elif api.startswith('constants.') and api.rsplit('.',1)[1] in ['value','unit','precision','find']:
                    key=keys[i*len(keys)//count] if method=='linspace' else rng.choice(keys)
                    args['key']=key[:min(len(key),8)] if api.endswith('find') else key
                elif api=='stats.normal_ppf':args['x']=.001+.998*(t+1)/2
                elif api in ['stats.normal_pdf','constants.convert_temperature']:args['x']=8*t
                elif api.startswith('differentiate.'):args['point']=[t,1+t/2]
                elif api.startswith('optimize.'):args['start']=[3*t,1+2*t]
                elif api.startswith('cluster.'):
                    args.update(obs=[[-2.+t/10,1.],[-1.+t/10,2.],[3.+t/10,4.],[4.+t/10,5.]],initial=[[-2.,1.],[4.,5.]])
                elif api.startswith('interpolate.'):
                    args.update(nodes=[-2.,-1.,0.,1.,2.],values=[1.,2.+t/10,3.,2.,1.],query=[t,t/2])
                elif api=='integrate.newton_cotes':args['values']=[0.,1.+t/100,2.,3.,4.]
                elif api.startswith('integrate.'):args['values']=[float(j*j)/16+t*(j+1)/10 for j in range(9)]
                yield dict(id=f'{api}/{method}-{i:04d}',api=api,group=f'{method}, bounded explicit recipe',args=args)


def reference(case):
    api,a=case['api'],case['args']
    def text(value):return [len(value.encode()),*value.encode()]
    if api.startswith('linalg.'):
        name=api.rsplit('.',1)[1]
        value=np.dot(a['u'],a['v']) if name=='dot' else np.asarray(a['a'])@a['u'] if name=='matvec' else np.asarray(a['a'])@a['b'] if name=='matmul' else np.asarray(a['a']).T
        return np.asarray(value).ravel().tolist()
    if api=='special.geterr':return [0.]
    if api=='stats.normal_pdf':return [stats.norm.pdf(a['x'])]
    if api=='stats.normal_ppf':return [stats.norm.ppf(a['x'])]
    if api=='constants.value':return [constants.value(a['key'])]
    if api=='constants.precision':return [constants.precision(a['key'])]
    if api=='constants.unit':return text(constants.unit(a['key']))
    if api=='constants.find':return [x for value in constants.find(a['key']) for x in text(value)]
    if api=='constants.convert_temperature':return [constants.convert_temperature(a['x'],'Celsius','Fahrenheit')]
    if api=='differentiate.hessian':return [2.,3.,3.,4.]
    if api=='differentiate.jacobian':
        x,y=a['point'];return [y,x,1.,1.]
    if api.startswith('optimize.'):return [1.,-2.,3.]
    if api.startswith('cluster.'):
        obs=np.asarray(a['obs']);initial=np.asarray(a['initial'])
        if api.endswith('kmeans2'):
            centers,labels=cluster.vq.kmeans2(obs,initial,iter=8,minit='matrix')
            return [*centers.ravel(),*labels.astype(float)]
        centers,distortion=cluster.vq.kmeans(obs,initial,thresh=1e-8)
        return [*centers.ravel(),float(distortion)]
    if api.startswith('interpolate.'):
        name=api.rsplit('.',1)[1]
        if name=='pade':
            p,q=interpolate.pade(a['values'],2,2);return [*p.c,*q.c]
        if name=='lagrange':return interpolate.lagrange(a['nodes'],a['values']).c.tolist()
        return interpolate.pchip_interpolate(a['nodes'],a['values'],a['query']).tolist()
    name=api.rsplit('.',1)[1]
    if name=='newton_cotes':
        w,b=integrate.newton_cotes(np.asarray(a['values']));return [*w,float(b)]
    kwargs={'dx':.25}
    if name.startswith('cumulative_'):kwargs['initial']=0.
    return np.asarray(getattr(integrate,name)(a['values'],**kwargs)).reshape(-1).tolist()


def comparison(case):
    api,a=case['api'],case['args']
    if api.startswith('optimize.'):
        from scipy import optimize
        f=lambda p:3.+(p[0]-1.)**2+(p[1]+2.)**2
        g=lambda p:np.array([2*(p[0]-1),2*(p[1]+2)])
        gradient='gradient' in api
        result=optimize.minimize(f,a['start'],method='BFGS' if gradient else 'Nelder-Mead',jac=g if gradient else None,tol=1e-8,options={'maxiter':512})
        if not result.success:raise ValueError(result.message)
        return [*result.x,float(result.fun)]
    if api.startswith('differentiate.'):
        from scipy import differentiate
        point=np.asarray(a['point'])
        if api.endswith('hessian'):
            r=differentiate.hessian(lambda p:p[0]**2+3*p[0]*p[1]+2*p[1]**2,point,tolerances={'atol':1e-10,'rtol':1e-10})
            values=r.ddf
        else:
            r=differentiate.jacobian(lambda p:np.array([p[0]*p[1],p[0]+p[1]]),point,tolerances={'atol':1e-10,'rtol':1e-10})
            values=r.df
        if not np.all(r.success):raise ValueError('SciPy adaptive differences did not converge')
        return np.asarray(values).ravel().tolist()
    return reference(case)
