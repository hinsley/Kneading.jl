import refine_entropy as r
import numpy as np,ctypes,json,concurrent.futures
from flint import fmpz_poly,ctx
ctx.prec=128
lib=ctypes.CDLL(str(r.ROOT/'symmetric_entropy.so'));lib.symmetric_factors.argtypes=[ctypes.c_double,ctypes.c_int,ctypes.POINTER(ctypes.c_longlong),ctypes.POINTER(ctypes.c_longlong)]
d=np.load(r.ROOT/'entropy_boundary_log_3000.npz');u=d['u']
def solve(row):
 i,n=map(int,row);B=(ctypes.c_longlong*(n+1))();Q=(ctypes.c_longlong*(n+1))();lib.symmetric_factors(u[i],n,B,Q);vals=[]
 for a in [B,Q]:
  p=fmpz_poly(list(a));rr=p.complex_roots();real=[z.real for z,m in rr if z.imag.is_zero() and z.real>0 and z.real<1];vals.append(-r.math.log(float(min(real))) if real else 0.)
 return i,n,max(vals)
if __name__=='__main__':
 refined=np.load(r.ROOT/'entropy_boundary_refined_3000.npz');ids=refined['ids'];indices=set((ids[(ids//3000+ids%3000)==2999]//3000).tolist());indices.update(np.load(r.ROOT/'entropy_extra_points.npy')[:,0].tolist());out=[]
 with concurrent.futures.ProcessPoolExecutor(max_workers=6) as pool:
  for rec in pool.map(solve,[(i,n) for i in sorted(indices) for n in [256,512]],chunksize=1):out.append(rec);print(rec,flush=True)
 (r.ROOT/'entropy-symmetric-factors.json').write_text(json.dumps(out))
