import refine_entropy as r
import numpy as np,ctypes,json,time,concurrent.futures
from flint import fmpz_poly,ctx
ctx.prec=256
d=np.load(r.ROOT/'entropy_boundary_refined_3000.npz');ids=d['ids'];u,v=d['u'],d['v']
def job(pair):
 k,n=pair;i,j=divmod(int(ids[k]),3000);a=(ctypes.c_longlong*(n+1))();r.lib.refine_coeff(u[i],v[j],n,a)
 roots=fmpz_poly(list(a)).complex_roots();real=[z.real for z,m in roots if z.imag.is_zero() and z.real>0 and z.real<1];h=-r.math.log(float(min(real))) if real else 0
 return k,n,h
if __name__=='__main__':
 pairs=[(k,n) for k in [91241,91617,92000,92386,93154,93540,95104,97016,97539] for n in [1024,1025]];out=[]
 with concurrent.futures.ProcessPoolExecutor(max_workers=6) as pool:
  for x in pool.map(job,pairs):out.append(x);print(x,flush=True)
 (r.ROOT/'entropy_ambiguous_orders.json').write_text(json.dumps(out))
