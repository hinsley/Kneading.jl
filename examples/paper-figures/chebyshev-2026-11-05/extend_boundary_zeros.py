import refine_entropy as r
import numpy as np,concurrent.futures,time
if __name__=='__main__':
 d=np.load(r.ROOT/'entropy_boundary_zeros_final.npz');base=np.load(r.ROOT/'entropy_boundary_refined_final_3000.npz');ids=d['ids'];best=d['best'].copy();conv=d['converged'].copy();orders=d['orders'].copy();u,v=base['u'],base['v']
 for n in [16384,32768]:
  active=np.flatnonzero(~conv);params=[(int(k),float(u[ids[k]//3000]),float(v[ids[k]%3000])) for k in active];t=time.time()
  with concurrent.futures.ProcessPoolExecutor(max_workers=1,initializer=r.init,initargs=(n,)) as pool:
   for k,h in pool.map(r.job,params,chunksize=1):
    conv[k]=np.isfinite(h)&(abs(h-best[k])<=.001);best[k]=h;orders[k]=n;print(n,k,h,bool(conv[k]),time.time()-t,flush=True)
  np.savez(r.ROOT/f'entropy_boundary_zeros_extended_{n}.npz',ids=ids,best=best,orders=orders,converged=conv,periodic=d['periodic'])
  if conv.all():break
 np.savez(r.ROOT/'entropy_boundary_zeros_final.npz',ids=ids,best=best,orders=orders,converged=conv,periodic=d['periodic'])
 print('Remaining',sum(~conv),flush=True)
