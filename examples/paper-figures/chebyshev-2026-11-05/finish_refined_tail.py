import refine_entropy as r
import numpy as np,ctypes,json,pickle,concurrent.futures,time
from flint import fmpz_poly,ctx
ctx.prec=256
d=np.load(r.ROOT/'entropy_boundary_refined_3000.npz');ext=np.load(r.ROOT/'entropy_extended_8192.npz');ids=d['ids'];u,v=d['u'],d['v'];best=ext['best'].copy();orders=ext['orders'].copy();conv=ext['converged'].copy()
records=pickle.load(open(r.ROOT/'entropy_periodic_tails.pkl','rb'));roots=pickle.load(open(r.ROOT/'entropy_periodic_roots.pkl','rb'));periodic=np.zeros(len(ids),bool)
for k,n,p,meta in records:
 if p is not None:best[k]=roots[p];conv[k]=True;orders[k]=32768;periodic[k]=True

def arbjob(pair):
 k,n=pair;ctx.prec=256;i,j=divmod(int(ids[k]),3000);a=(ctypes.c_longlong*(n+1))();r.lib.refine_coeff(u[i],v[j],n,a)
 rr=fmpz_poly(list(a)).complex_roots();real=[z.real for z,m in rr if z.imag.is_zero() and z.real>0 and z.real<1];h=-r.math.log(float(min(real))) if real else 0
 return k,n,h
if __name__=='__main__':
 # Symmetry creates nearly double roots; compare adjacent truncations with Arb.
 diag=np.flatnonzero(((ids//3000+ids%3000)==2999)&~periodic)
 cache=r.ROOT/'entropy-diagonal-arb.json';out=json.loads(cache.read_text()) if cache.exists() else []
 done={(k,n) for k,n,h in out};pairs=[(int(k),n) for k in diag for n in [512,513,1024,1025] if (int(k),n) not in done]
 t=time.time()
 with concurrent.futures.ProcessPoolExecutor(max_workers=6) as pool:
  for count,res in enumerate(pool.map(arbjob,pairs,chunksize=1),1):
   out.append(res);cache.write_text(json.dumps(out))
   if count%10==0:print(f'Arb {count}/{len(pairs)}, {time.time()-t:.1f}s',flush=True)
 vals={(k,n):h for k,n,h in out}
 for k in diag:
  a=max(vals[k,512],vals[k,513]);b=max(vals[k,1024],vals[k,1025]);best[k]=b;orders[k]=1025;conv[k]=abs(a-b)<=.001
 print('After tails and high-precision symmetric points, remaining',sum(~conv),flush=True)
 # Remaining nonsymmetric points: double the already long calculation once more.
 active=np.flatnonzero(~conv);print('remaining ids',active.tolist(),flush=True)
 np.savez(r.ROOT/'entropy_refinement_combined.npz',ids=ids,best=best,orders=orders,converged=conv,periodic_tail=periodic)
