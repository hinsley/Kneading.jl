import refine_entropy as r
import numpy as np,concurrent.futures,time
if __name__=='__main__':
 d=np.load(r.ROOT/'entropy_boundary_refined_3000.npz');x=np.load(r.ROOT/'entropy_extended_8192.npz');ids=d['ids'];u,v=d['u'],d['v'];active=np.load(r.ROOT/'entropy_remaining_nonsym.npy');best=x['best'][active].copy();conv=np.zeros(len(active),bool)
 for n in [16384,32768]:
  t=time.time();todo=np.flatnonzero(~conv);params=[]
  for ix in todo:
   k=active[ix];i,j=divmod(int(ids[k]),3000);params.append((int(ix),float(u[i]),float(v[j])))
  with concurrent.futures.ProcessPoolExecutor(max_workers=1,initializer=r.init,initargs=(n,)) as pool:
   for ix,h in pool.map(r.job,params,chunksize=1):
    conv[ix]=np.isfinite(h)&(abs(best[ix]-h)<=.001);best[ix]=h;print(n,int(active[ix]),h,conv[ix],time.time()-t,flush=True)
  np.savez(r.ROOT/f'entropy_last_edge_{n}.npz',active=active,best=best,converged=conv)
  if conv.all():break
