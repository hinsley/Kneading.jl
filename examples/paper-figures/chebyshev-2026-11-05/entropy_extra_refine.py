import refine_entropy as r
import numpy as np,ctypes,json,concurrent.futures
from flint import fmpz_poly,ctx
ctx.prec=256
d=np.load(r.ROOT/'entropy_boundary_log_3000.npz');u,v=d['u'],d['v'];points=np.load(r.ROOT/'entropy_extra_points.npy')
def job(row):
 i,j,n=map(int,row);ctx.prec=256;a=(ctypes.c_longlong*(n+1))();r.lib.refine_coeff(u[i],v[j],n,a);rr=fmpz_poly(list(a)).complex_roots();real=[z.real for z,m in rr if z.imag.is_zero() and z.real>0 and z.real<1]
 return i,j,n,-r.math.log(float(min(real))) if real else 0.
if __name__=='__main__':
 out=[]
 with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
  for rec in pool.map(job,[(i,j,n) for i,j in points for n in [512,513,1024,1025]]):out.append(rec);print(rec,flush=True)
 (r.ROOT/'entropy-extra-arb.json').write_text(json.dumps(out))
