import refine_entropy as r
import numpy as np,ctypes,json,time,concurrent.futures,pickle
from flint import fmpz_poly,ctx
lib=ctypes.CDLL(str(r.ROOT/'periodic_tail.so'));lib.periodic_numerator.argtypes=[ctypes.c_double,ctypes.c_double,ctypes.c_int,ctypes.c_int,ctypes.POINTER(ctypes.c_longlong),ctypes.POINTER(ctypes.c_int)]
ctx.prec=192
def job(row):
 k,u,v=row;a=(ctypes.c_longlong*8192)();meta=(ctypes.c_int*4)();n=lib.periodic_numerator(u,v,32768,8192,a,meta)
 if n<0:return k,n,None,None
 p=fmpz_poly([a[i] for i in range(n)])
 # Cancel all denominator factors to reduce the entropy numerator exactly.
 p0,p1=meta[1],meta[3];q0=fmpz_poly([1]+[0]*(p0-1)+[-1]);q1=fmpz_poly([1]+[0]*(p1-1)+[-1]);den=q0*q1*fmpz_poly([1,-1]);p=p//p.gcd(den)
 return k,n,tuple(int(c) for c in p.coeffs()),tuple(meta)
def rootjob(p):
 ctx.prec=192
 roots=fmpz_poly(list(p)).complex_roots();real=[z.real for z,m in roots if z.imag.is_zero() and z.real>0 and z.real<1]
 return p,-r.math.log(float(min(real))) if real else 0.
if __name__=='__main__':
 d=np.load(r.ROOT/'entropy_boundary_refined_3000.npz');ids=d['ids'];u,v=d['u'],d['v'];cache=r.ROOT/'entropy_periodic_tails.pkl'
 if cache.exists():records=pickle.loads(cache.read_bytes())
 else:
  params=[(k,float(u[ix//3000]),float(v[ix%3000])) for k,ix in enumerate(ids)];records=[];t=time.time()
  with concurrent.futures.ProcessPoolExecutor(max_workers=10) as pool:
   for count,rec in enumerate(pool.map(job,params,chunksize=64),1):
    records.append(rec)
    if count%10000==0:print(f'Long tail check {count}/{len(ids)} in {time.time()-t:.1f}s',flush=True)
  cache.write_bytes(pickle.dumps(records))
 groups={};status={}
 for k,n,p,meta in records:
  status[str(n if n<0 else 1)]=status.get(str(n if n<0 else 1),0)+1
  if p is not None:groups.setdefault(p,[]).append(k)
 print('status',status,'unique',len(groups),'degrees',np.quantile([len(p)-1 for p in groups],[0,.5,.9,1]),flush=True)
 results={};rootcache=r.ROOT/'entropy_periodic_roots.pkl'
 if rootcache.exists():results=pickle.loads(rootcache.read_bytes())
 todo=[p for p in groups if p not in results]
 with concurrent.futures.ProcessPoolExecutor(max_workers=5) as pool:
  for count,(p,h) in enumerate(pool.map(rootjob,todo,chunksize=1),1):
   results[p]=h
   if count%25==0:rootcache.write_bytes(pickle.dumps(results));print(f'Arb rational roots {count}/{len(todo)}',flush=True)
 rootcache.write_bytes(pickle.dumps(results))
 r.ROOT.joinpath('entropy_periodic_summary.json').write_text(json.dumps({'status':status,'unique':len(groups),'verified_iterates':32768},indent=2))
