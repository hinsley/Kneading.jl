import validate_refined_entropy as q
from flint import fmpz_poly,ctx
import concurrent.futures,math,json,time
ctx.prec=128
def solve(k):
 i,j=divmod(int(q.ids[k]),3000);a=(q.r.ctypes.c_longlong*513)();q.r.lib.refine_coeff(q.u[i],q.v[j],512,a)
 roots=fmpz_poly(list(a)).complex_roots();real=[z.real for z,m in roots if z.imag.is_zero() and z.real>0 and z.real<1]
 return int(k),-math.log(float(min(real))) if real else 0.
if __name__=='__main__':
 ids=q.np.flatnonzero(~q.np.isfinite(q.values));out=[]
 with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
  for item in pool.map(solve,ids):out.append(item);print(item,flush=True)
 (q.r.ROOT/'entropy-nan512-arb.json').write_text(json.dumps(out))
