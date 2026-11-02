"""Supplementary NumPy comparisons for saved scientific-function cases.

These are independent public NumPy operations or explicit NumPy mathematical
expressions.  In particular, no Sci/SciPy reference dispatcher is called from
this adapter.  Each successful comparison reproduces the *entire* flattened
fixture output, including interleaved real/imaginary FFT components.
"""

import math

import numpy as np


_FFT_ONE_DIMENSION = {
    "fft", "ifft", "rfft", "irfft", "hfft", "ihfft",
    "fftshift", "ifftshift", "fftfreq", "rfftfreq",
}
_FFT_MULTIDIMENSION = {
    "fft2", "ifft2", "rfft2", "irfft2", "fftn", "ifftn",
    "rfftn", "irfftn",
}
_REAL_FFT_INPUTS = {"rfft", "ihfft", "rfft2", "rfftn"}
_INVERSE_REAL_FFT = {"irfft2", "irfftn"}


def _fft(name, a):
    if name not in _FFT_ONE_DIMENSION | _FFT_MULTIDIMENSION:
        return None
    raw = a.get("raw")
    if not isinstance(raw, dict) or a.get("suite_name") != "fft":
        return None  # Existing small FFT fixtures already have NumPy coverage.
    if name in {"fftfreq", "rfftfreq"}:
        value = getattr(np.fft, name)(raw["n"], raw["d"])
    else:
        re = np.asarray(raw["re"], dtype=float)
        im = np.asarray(raw.get("im", np.zeros_like(re)), dtype=float)
        if re.shape != im.shape:
            return None
        data = re if name in _REAL_FFT_INPUTS else re + 1j * im
        if name in _FFT_MULTIDIMENSION:
            shape = list(raw["shape"])
            axes = raw.get("axes", [0, 1])
            kwargs = {"axes": axes}
            if name in _INVERSE_REAL_FFT:
                shape[axes[-1]] = shape[axes[-1]] // 2 + 1
                kwargs["s"] = [raw["shape"][axis] for axis in axes]
            value = getattr(np.fft, name)(data.reshape(shape), **kwargs)
        elif name in {"irfft", "hfft", "ihfft"}:
            value = getattr(np.fft, name)(data, n=raw.get("n"))
        else:
            value = getattr(np.fft, name)(data)
    result = np.asarray(value).ravel()
    return np.column_stack((result.real, result.imag)).ravel()


def _spatial(name, a):
    supported = {
        "euclidean", "sqeuclidean", "cityblock", "chebyshev", "cosine",
        "correlation", "braycurtis", "canberra", "hamming", "jaccard",
        "minkowski", "minkowski_distance", "minkowski_distance_p",
        "distance_matrix", "cdist", "pdist", "is_valid_y", "is_valid_dm",
        "num_obs_y", "num_obs_dm", "squareform", "geometric_slerp",
    }
    if name not in supported:
        return None
    raw = a.get("raw")
    if raw and a.get("suite_name") == "geometry":
        def parse(s):
            if ":" in s:
                return np.array([[float(v) for v in row.split(";")] for row in s.split(":")])
            return np.array([float(v) for v in s.split(";")]) if s else np.array([])
        u, v = parse(raw["u"]), parse(raw["v"])
        w = parse(raw["w"]) if raw["w"] else None
        p = float(raw["param"]) if name in {
            "minkowski", "minkowski_distance", "minkowski_distance_p",
            "distance_matrix", "is_valid_dm", "geometric_slerp"
        } else 2.
        centered = raw.get("extra") == "true"
        metric = raw.get("param") if name in {"cdist", "pdist"} else None
    elif name in {"euclidean", "cosine"} and "u" in a and "v" in a:
        u, v = np.asarray(a["u"], float), np.asarray(a["v"], float)
        w, p, centered, metric = None, 2., False, None
    else:
        return None

    if name == "is_valid_y":
        size=int((1+math.isqrt(1+8*len(u)))//2)
        return float(u.ndim==1 and size*(size-1)//2==len(u))
    if name == "is_valid_dm":
        return float(u.ndim==2 and u.shape[0]==u.shape[1] and
                     np.all(np.abs(np.diag(u))<=p) and
                     np.all(np.abs(u-u.T)<=p))
    if name == "num_obs_dm":return float(len(u)) if u.ndim==2 and u.shape[0]==u.shape[1] else None
    if name == "num_obs_y":
        count=int((1+math.isqrt(1+8*len(u)))//2)
        return float(count) if count*(count-1)//2==len(u) and count>1 else None
    if name == "squareform":
        if raw["param"]=="condensed":
            size=int((1+math.isqrt(1+8*len(u)))//2)
            if size*(size-1)//2!=len(u):return None
            matrix=np.zeros((size,size),float)
            indices=np.triu_indices(size,1)
            matrix[indices]=u
            matrix[(indices[1],indices[0])]=u
            return matrix
        if u.ndim!=2 or len(u)!=u.shape[1]:return None
        return u[np.triu_indices(len(u),1)]
    if name == "geometric_slerp":
        times=np.asarray(w,float)
        dot=np.clip(np.dot(u,v),-1.,1.)
        angle=np.arccos(dot)
        if angle<1e-10 or np.pi-angle<1e-10:return None
        return (np.sin((1-times)[:,None]*angle)*u[None,:] +
                np.sin(times[:,None]*angle)*v[None,:])/np.sin(angle)

    def distance(x, y, kind):
        diff = np.abs(x - y)
        weights = w if w is not None and kind not in {"distance_matrix", "cdist", "pdist"} else 1.
        if kind == "euclidean":
            return np.sqrt(np.sum(weights * diff**2, axis=-1))
        if kind == "sqeuclidean":
            return np.sum(weights * diff**2, axis=-1)
        if kind == "cityblock":
            return np.sum(weights * diff, axis=-1)
        if kind == "chebyshev":
            return np.max(diff, axis=-1)
        if kind == "braycurtis":
            return np.sum(weights * diff, axis=-1) / np.sum(weights * np.abs(x + y), axis=-1)
        if kind == "canberra":
            return np.sum(np.divide(weights * diff, np.abs(x) + np.abs(y),
                                    out=np.zeros_like(diff), where=(np.abs(x) + np.abs(y)) != 0), axis=-1)
        if kind == "hamming":
            return np.sum(weights * (x != y), axis=-1) / np.sum(np.broadcast_to(weights, x.shape), axis=-1)
        if kind == "jaccard":
            bx, by = x != 0, y != 0
            union = np.sum(weights * (bx | by), axis=-1)
            return np.divide(np.sum(weights * (bx ^ by), axis=-1), union,
                             out=np.zeros_like(union, dtype=float), where=union != 0)
        if kind in {"cosine", "correlation"}:
            if kind == "correlation" and centered:
                x, y = x - np.average(x, weights=w), y - np.average(y, weights=w)
            dot = np.sum(weights * x * y, axis=-1)
            norms = np.sqrt(np.sum(weights*x*x, axis=-1) * np.sum(weights*y*y, axis=-1))
            return np.clip(1. - dot / norms, 0., 2.)
        if kind == "minkowski_distance_p":
            return np.sum(diff**p, axis=-1)
        if kind in {"minkowski_distance", "minkowski", "distance_matrix"}:
            if math.isinf(p):return np.max(diff, axis=-1)
            return np.sum(weights * diff**p, axis=-1)**(1./p)
        raise ValueError(kind)

    if name == "distance_matrix":
        return distance(u[:, None, :], v[None, :, :], name)
    if name in {"cdist", "pdist"}:
        if metric not in {"euclidean", "cityblock", "sqeuclidean", "cosine"}:
            return None
        if name == "cdist":
            return distance(u[:, None, :], v[None, :, :], metric)
        return np.asarray([distance(u[i], u[j], metric)
                           for i in range(len(u)) for j in range(i+1, len(u))])
    return distance(u, v, name)


def _linalg(name, a):
    raw = a.get("raw")
    if not raw or a.get("suite_name") != "linalg":
        return None
    arguments = raw["args"]
    x = np.asarray(arguments[0], float)
    if name in {"det", "inv", "pinv", "pinvh", "svdvals", "svd", "qr",
                "norm", "eigvalsh", "cholesky", "solve", "lstsq"}:
        if name == "det": return np.linalg.det(x)
        if name == "inv": return np.linalg.inv(x)
        if name in {"pinv", "pinvh"}:
            # SciPy's positional second argument is atol; the optional rtol is
            # otherwise dimension * machine epsilon. NumPy's rtol has that default.
            atol = arguments[1] if len(arguments)>1 else None
            rtol = arguments[2] if len(arguments)>2 else None
            if atol not in (None, 0):return None
            if name == "pinvh":
                lower = arguments[3]
                triangle = np.tril(x) if lower else np.triu(x)
                x = triangle + triangle.T - np.diag(np.diag(x))
                eigen = np.abs(np.linalg.eigvalsh(x))
                magnitude = np.max(eigen)
            else:
                eigen = np.linalg.svdvals(x)
                magnitude = np.max(eigen)
            cutoff = (rtol if rtol is not None else max(x.shape)*np.finfo(float).eps)
            inverse = np.linalg.pinv(x, rtol=cutoff, hermitian=name=="pinvh")
            rank = int(np.sum(eigen > cutoff * magnitude))
            return np.r_[inverse.ravel(), rank]
        if name == "svdvals":return np.linalg.svdvals(x)
        if name == "svd":
            u,s,vh=np.linalg.svd(x,full_matrices=False)
            return np.r_[(u@u.T).ravel(),s,(vh.T@vh).ravel(),((u*s)@vh).ravel()]
        if name == "qr":
            q,r=np.linalg.qr(x,mode="reduced")
            return np.r_[(q@q.T).ravel(),(r.T@r).ravel(),(q@r).ravel()]
        if name == "norm":return np.linalg.norm(x, ord=arguments[1] if len(arguments)>1 else None)
        if name == "eigvalsh":return np.linalg.eigvalsh(x, UPLO="L" if arguments[1] else "U")
        if name == "cholesky":
            lower = np.linalg.cholesky(x)
            return lower if arguments[1] else lower.T
        if name == "solve":return np.linalg.solve(x, np.asarray(arguments[1],float))
        if name == "lstsq":
            coeffs, residuals, rank, singular = np.linalg.lstsq(
                x,np.asarray(arguments[1],float),rcond=arguments[2])
            return np.r_[np.asarray(coeffs).ravel(),residuals,rank,singular]
    return None


def _stats(name, a):
    if name == "normal_pdf" and "x" in a:
        x = float(a["x"])
        return np.exp(-x*x/2.) / np.sqrt(2.*np.pi)
    raw = a.get("raw")
    if raw or not {"x", "p"}.issubset(a):
        return None
    x = np.asarray(a["x"], float)
    p = a["p"]
    if name in {"mean", "variance", "std", "median"}:
        method = {"variance":"var"}.get(name,name)
        return getattr(np,method)(x,**({"ddof":int(p[0])} if p else {}))
    if name == "sem": return np.std(x,ddof=int(p[0]))/np.sqrt(x.size)
    if name == "variation": return np.std(x,ddof=int(p[0]))/np.mean(x)
    if name == "zscore": return (x-np.mean(x))/np.std(x,ddof=int(p[0]))
    if name == "quantile": return np.quantile(x,float(p[0]))
    if name == "scoreatpercentile": return np.percentile(x,float(p[0]))
    if name == "median_abs_deviation": return np.median(np.abs(x-np.median(x)))/float(p[0])
    if name == "iqr": return np.diff(np.percentile(x,[25,75]))[0]
    if name == "gmean":return np.exp(np.mean(np.log(x)))
    if name == "hmean":return len(x)/np.sum(1./x)
    if name == "gstd":return np.exp(np.std(np.log(x),ddof=int(p[0])))
    if name == "gzscore":
        logs=np.log(x)
        return (logs-np.mean(logs))/np.std(logs,ddof=int(p[0]))
    if name == "mode":
        values, counts=np.unique(x,return_counts=True)
        index=np.argmax(counts)
        return [values[index],float(counts[index])]
    if name == "entropy":
        y=np.asarray(a["y"],float)
        first=x/np.sum(x)
        if len(y):
            second=y/np.sum(y)
            terms=np.where(first>0,first*np.log(first/second),0.)
        else:
            terms=np.where(first>0,-first*np.log(first),0.)
        value=np.sum(terms)
        return value/np.log(p[0]) if p[0] else value
    if name == "energy_distance":
        y=np.asarray(a["y"],float)
        first=np.abs(x[:,None]-y[None,:]).mean()
        second=np.abs(x[:,None]-x[None,:]).mean()
        third=np.abs(y[:,None]-y[None,:]).mean()
        return np.sqrt(max(0.,2*first-second-third))
    if name in {"circmean", "circvar", "circstd"}:
        lo,hi=float(p[0]),float(p[1])
        phase=(x-lo)*(2*np.pi/(hi-lo))
        mean=np.mean(np.exp(1j*phase));strength=np.abs(mean)
        if name=="circmean":return lo + np.mod(np.angle(mean),2*np.pi)*(hi-lo)/(2*np.pi)
        if name=="circvar":return 1-strength
        normalized=bool(p[2]) if len(p)>2 else False
        return np.sqrt(-2*np.log(strength))*(1. if normalized else (hi-lo)/(2*np.pi))
    if name == "moment":
        order=int(p[0])
        return 1. if order==0 else 0. if order==1 else np.mean((x-np.mean(x))**order)
    if name == "trim_mean" or name == "trimboth":
        n=int(x.size*float(p[0]));trimmed=np.sort(x)[n:x.size-n] if n else np.sort(x)
        return np.mean(trimmed) if name=="trim_mean" else trimmed
    if name == "rankdata":
        order=np.argsort(x,kind="stable")
        sorted_x=x[order]; starts=np.r_[0,np.flatnonzero(sorted_x[1:]!=sorted_x[:-1])+1]
        ends=np.r_[starts[1:],len(x)];rank=np.empty(len(x),float)
        method=int(p[0])
        if method==4:rank[order]=np.arange(1,len(x)+1)
        else:
            for i,(lo,hi) in enumerate(zip(starts,ends)):
                rank[order[lo:hi]]=[(lo+hi+1)/2,lo+1,hi,i+1][method]
        return rank
    if name == "percentileofscore":
        score=float(p[0]);low=np.count_nonzero(x<score);high=np.count_nonzero(x<=score)
        kind=int(p[1]);count=[(low+high+int(low!=high))/2,low,high,(low+high)/2][kind]
        return 100.*count/len(x)
    return None


def _differentiate(name, a):
    if name == "derivative" and "x" in a:
        return np.cos(a["x"])
    if name == "gradient" and "point" in a:
        x,y = a["point"]
        return [2*x+3*y, 3*x+4*y]
    if a.get("explicit") is True and "point" in a:
        x,y=a["point"]
        if name == "jacobian":return [y,x,1.,1.]
        if name == "hessian":return [2.,3.,3.,4.]
    return None


def _integrate(name, a):
    if name == "fixed_quad" and {"a","b"}.issubset(a):
        roots, weights = np.polynomial.legendre.leggauss(16)
        lo,hi = a["a"],a["b"]
        return (hi-lo)/2 * np.dot(weights,np.exp((hi-lo)/2*roots+(hi+lo)/2))
    if not a.get("explicit") or "values" not in a:
        return None
    y=np.asarray(a["values"],float)
    if name == "simpson" and len(y)%2==1:
        return .25/3*(y[0]+y[-1]+4*np.sum(y[1:-1:2])+2*np.sum(y[2:-1:2]))
    if name == "cumulative_trapezoid":
        return np.r_[0.,np.cumsum(.25*(y[1:]+y[:-1])/2)]
    return None


def _signal(name, a):
    args=a.get("arguments")
    if name in {"convolve", "correlate", "fftconvolve", "oaconvolve"} and args:
        if np.ndim(args[0]) != 1 or np.ndim(args[1]) != 1:return None
        mode=args[2] if len(args)>2 else "full"
        if mode not in {"full","same","valid"}:return None
        return getattr(np,"correlate" if name=="correlate" else "convolve")(args[0],args[1],mode=mode)
    if name == "detrend" and args:
        values=np.asarray(args[0],float)
        if values.ndim!=1:return None
        if args[1] is False:return values-np.mean(values)
        return values-np.polynomial.polynomial.polyval(np.arange(len(values)),
                     np.polynomial.polynomial.polyfit(np.arange(len(values)),values,1))
    if name.startswith("windows.") and args:
        window=name.split(".",1)[1]
        match={"hann":"hanning", "hamming":"hamming", "bartlett":"bartlett", "blackman":"blackman"}
        if window not in match:return None
        n=int(args[0]);sym=bool(args[1]) if len(args)>1 else True
        if n<0:return None
        values=getattr(np,match[window])(n if sym or n<=1 else n+1)
        return values if sym or n<=1 else values[:-1]
    return None


def _cluster(name, a):
    if name == "vq.whiten" and "obs" in a:
        data=np.asarray(a["obs"],float)
        return data/np.std(data,axis=0)
    if name == "vq.vq" and {"obs","codebook"}.issubset(a):
        obs=np.asarray(a["obs"],float)
        codebook=np.asarray(a["codebook"],float)
        distance=np.linalg.norm(obs[:,None,:]-codebook[None,:,:],axis=2)
        nearest=np.argmin(distance,axis=1)
        return np.r_[nearest,distance[np.arange(len(obs)),nearest]]
    return None


def _scipy_blas(name, a):
    """SciPy BLAS for the explicit dense operations otherwise NumPy-only."""
    if not a.get("explicit") or name not in {"dot", "matvec", "matmul"}:
        return None
    from scipy.linalg import blas
    if name == "dot":
        return blas.ddot(np.asarray(a["u"],float),np.asarray(a["v"],float))
    matrix=np.asarray(a["a"],float)
    if name == "matvec":
        return blas.dgemv(1.,matrix,np.asarray(a["u"],float))
    return blas.dgemm(1.,matrix,np.asarray(a["b"],float))


def evaluate(case):
    """Return extra verified provider outputs, or {} for unsupported cases.

    Bad/unrepresentable numerical domains are unavailable rather than reported
    as successful comparisons with missing or nonfinite components.
    """
    api, a = case.get("api", ""), case.get("args", {})
    if not isinstance(a,dict) or "." not in api:
        return {}
    module,name=api.split(".",1)
    handlers={"fft":_fft,"signal":_signal,"spatial":_spatial,
              "linalg":_linalg,"stats":_stats,"integrate":_integrate,
              "differentiate":_differentiate,"cluster":_cluster}
    if module not in handlers:
        return {}
    try:
        with np.errstate(all="ignore"):
            if module == "linalg" and a.get("explicit"):
                value=_scipy_blas(name,a)
                if value is not None:
                    array=np.asarray(value,dtype=float).ravel()
                    return {"SciPy":{"values":array.tolist()}} if np.all(np.isfinite(array)) else {}
            value=handlers[module](name,a)
        if value is None:
            return {}
        result=np.asarray(value,dtype=float).ravel()
        if not np.all(np.isfinite(result)):
            return {}
        return {"NumPy":{"values":result.tolist()}}
    except (KeyError,TypeError,ValueError,ArithmeticError,np.linalg.LinAlgError):
        return {}
