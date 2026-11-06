import refine_entropy as r
import numpy as np,time,concurrent.futures,json
from flint import fmpz_poly,ctx
ctx.prec=128
old=np.load(r.ROOT/'entropy_boundary_log_3000.npz');ids=np.flatnonzero(old['fill_mask']|(abs(old['h64']-old['h32'])>.001));u,v=old['u'],old['v'];d=np.load(r.ROOT/'entropy_refinement_512.npz');values=d['values']
def job(k):
 t=time.time();r.init(512);i,j=divmod(int(ids[k]),3000);a=(r.ctypes.c_longlong*513)();r.lib.refine_coeff(u[i],v[j],512,a)
 p=fmpz_poly(list(a));factors=p.factor_squarefree()[1]
 if not np.isfinite(values[k]):
  estimates=[]
  for f,e in factors:
   aa=(r.ctypes.c_longlong*(f.degree()+1))(*[int(c) for c in f.coeffs()]);estimates.append(r.lib.refine_entropy(aa,f.degree()))
  return {'k':int(k),'fallback_degrees':[f.degree() for f,e in factors],'fallback_values':estimates}
 roots=p.complex_roots();real=[z.real for z,m in roots if z.imag.is_zero() and z.real>0 and z.real<1]
 h=-r.math.log(float(min(real))) if real else 0.
 return {'k':int(k),'reference':h,'computed':float(values[k]),'error':abs(h-values[k]),'seconds':time.time()-t}
if __name__=='__main__':
 rng=np.random.default_rng(20261105);sample=list(rng.choice(np.flatnonzero(np.isfinite(values)),8,replace=False));sample+=list(np.flatnonzero(~np.isfinite(values))[:3])
 results=[]
 with concurrent.futures.ProcessPoolExecutor(max_workers=3) as pool:
  for out in pool.map(job,sample):results.append(out);print(out,flush=True)
 (r.ROOT/'entropy-refined-validation.json').write_text(json.dumps(results,indent=2))
