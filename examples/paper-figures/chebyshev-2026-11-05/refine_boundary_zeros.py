import refine_entropy as r
import periodic_tail_refine as p
import numpy as np,ctypes,pickle,json,concurrent.futures,time
if __name__=='__main__':
 d=np.load(r.ROOT/'entropy_boundary_refined_final_3000.npz');u,v=d['u'],d['v'];ids=np.unique(np.r_[np.load(r.ROOT/'entropy_boundary_zeros.npy'),d['ids'][(d['h_refined'].ravel()[d['ids']]==0)&(d['method']==1)]])
 params=[(k,float(u[ix//3000]),float(v[ix%3000])) for k,ix in enumerate(ids)];periodic=np.zeros(len(ids),bool);best=np.zeros(len(ids));conv=np.zeros(len(ids),bool);orders=np.zeros(len(ids),dtype='int32')
 cache=r.ROOT/'entropy_boundary_zeros_tails.pkl'
 if cache.exists():records=pickle.loads(cache.read_bytes())
 else:
  with concurrent.futures.ProcessPoolExecutor(max_workers=10) as pool:records=list(pool.map(p.job,params,chunksize=64))
  cache.write_bytes(pickle.dumps(records))
 roots=pickle.load(open(r.ROOT/'entropy_periodic_roots.pkl','rb'));todo={rec[2] for rec in records if rec[2] is not None and rec[2] not in roots}
 print('New boundary cells',len(ids),'periodic',sum(rec[2] is not None for rec in records),'new root calculations',len(todo),flush=True)
 with concurrent.futures.ProcessPoolExecutor(max_workers=6) as pool:
  for poly,h in pool.map(p.rootjob,todo,chunksize=1):roots[poly]=h
 for k,n,poly,meta in records:
  if poly is not None:periodic[k]=conv[k]=True;best[k]=roots[poly];orders[k]=32768
 previous=np.full(len(ids),np.nan)
 for n in [512,1024,2048,4096,8192]:
  active=np.flatnonzero(~conv);cache=r.ROOT/f'entropy_boundary_zeros_{n}.npz';t=time.time()
  if cache.exists():
   q=np.load(cache);assert np.array_equal(active,q['active']);vals=q['values']
  else:
   vals=[]
   with concurrent.futures.ProcessPoolExecutor(max_workers=3 if n==8192 else 10,initializer=r.init,initargs=(n,)) as pool:
    for k,h in pool.map(r.job,[params[k] for k in active],chunksize=16):vals.append(h)
   vals=np.array(vals);np.savez(cache,active=active,values=vals)
  good=np.isfinite(vals);delta=abs(vals-previous[active]);
  # Do not accept short zero/zero agreement at the boundary: require 4096 terms.
  if n>=1024:conv[active]=good&np.isfinite(previous[active])&(delta<=.001)&((vals>0)|(n>=4096))
  best[active[good]]=vals[good];orders[active[good]]=n;previous[active]=vals
  print('N',n,'remaining',sum(~conv),'elapsed',time.time()-t,flush=True)
  if conv.all():break
 np.savez(r.ROOT/'entropy_boundary_zeros_final.npz',ids=ids,best=best,orders=orders,converged=conv,periodic=periodic)
 print('Finished boundary refinement: unresolved',sum(~conv),flush=True)
