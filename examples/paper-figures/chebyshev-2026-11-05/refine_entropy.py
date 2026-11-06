from pathlib import Path
import numpy as np,ctypes,concurrent.futures,math,json,time
ROOT=Path(__file__).resolve().parent
lib=ctypes.CDLL(str(ROOT/'entropy_refine.so'))
lib.refine_init.argtypes=[ctypes.c_int]
lib.refine_coeff.argtypes=[ctypes.c_double,ctypes.c_double,ctypes.c_int,ctypes.POINTER(ctypes.c_longlong)]
lib.refine_entropy.argtypes=[ctypes.POINTER(ctypes.c_longlong),ctypes.c_int];lib.refine_entropy.restype=ctypes.c_double
ORDER=0
def init(n):
 global ORDER;ORDER=n;lib.refine_init(n)
def job(row):
 ix,u,v=row;a=(ctypes.c_longlong*(ORDER+1))();lib.refine_coeff(u,v,ORDER,a)
 h=lib.refine_entropy(a,ORDER)
 return ix,h
if __name__=='__main__':
 old=np.load(ROOT/'entropy_boundary_log_3000.npz');u,v=old['u'],old['v']
 # Tighten the old degree comparison for the edge as well as its original filled cells.
 target=old['fill_mask']|(abs(old['h64']-old['h32'])>.001)
 ids=np.flatnonzero(target);print('Refinement cells',len(ids),flush=True)
 previous=old['h64'].ravel()[ids].copy();best=previous.copy();orders=np.full(len(ids),64,dtype='int32');converged=np.zeros(len(ids),dtype=bool)
 # At least 512 terms, checked against 256; continue doubling where needed.
 for n in [256,512,1024,2048,4096]:
  active=np.flatnonzero(~converged);cache=ROOT/f'entropy_refinement_{n}.npz'
  if cache.exists():
   d=np.load(cache);assert np.array_equal(active,d['active']);vals=d['values']
  else:
   params=[(int(k),float(u[ids[k]//3000]),float(v[ids[k]%3000])) for k in active]
   t=time.time();vals=np.empty(len(active))
   with concurrent.futures.ProcessPoolExecutor(max_workers=10,initializer=init,initargs=(n,)) as pool:
    for count,(_,h) in enumerate(pool.map(job,params,chunksize=64),1):
     vals[count-1]=h
     if count%10000==0:print(f'N={n}: {count}/{len(active)} in {time.time()-t:.1f}s',flush=True)
   np.savez(cache,active=active,values=vals)
  good=np.isfinite(vals);delta=np.abs(vals-previous[active]);
  if n>=512:converged[active]=good&np.isfinite(previous[active])&(delta<=.001)
  best[active[good]]=vals[good];orders[active[good]]=n;previous[active]=vals
  print(f'N={n}: finite {good.sum()}/{len(vals)}, converged {converged.sum()}/{len(ids)}, max change {np.nanmax(delta):.6g}',flush=True)
  if converged.all():break
 raw=old['h64'].copy();raw.ravel()[ids]=best
 remaining=np.zeros_like(target);remaining.ravel()[ids[~converged]]=True
 display=raw.copy();display[remaining]=.01
 np.savez(ROOT/'entropy_boundary_refined_3000.npz',u=u,v=v,h_refined=raw,h_display=display,remaining_mask=remaining,target_mask=target,ids=ids,orders=orders,converged=converged)
 display.tofile(ROOT/'entropy_display_3000.bin')
 report={'target_cells':len(ids),'converged':int(converged.sum()),'remaining':int(remaining.sum()),'tolerance':.001,'minimum_order':512,'max_order':int(orders.max()),'order_counts':{str(n):int((orders==n).sum()) for n in np.unique(orders)}}
 (ROOT/'entropy-refinement.json').write_text(json.dumps(report,indent=2));print(report,flush=True)
