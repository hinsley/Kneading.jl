import refine_entropy as r
import numpy as np,ctypes,time,json,concurrent.futures
r.lib.refine_signature.argtypes=[ctypes.c_double,ctypes.c_double,ctypes.c_int,ctypes.POINTER(ctypes.c_ubyte)]
if __name__=='__main__':
 d=np.load(r.ROOT/'entropy_boundary_refined_3000.npz');u,v=d['u'],d['v'];ids=d['ids'];best=d['h_refined'].ravel()[ids].copy();orders=d['orders'].copy();conv=d['converged'].copy()
 for n in [8192,16384,32768]:
  active=np.flatnonzero(~conv);groups={};t=time.time()
  for k in active:
   i,j=divmod(int(ids[k]),3000);buf=(ctypes.c_ubyte*(2*n))();r.lib.refine_signature(u[i],v[j],n,buf);groups.setdefault(bytes(buf),[]).append(int(k))
  print(f'N={n}: {len(active)} cells, {len(groups)} distinct matrices',flush=True)
  keys=list(groups);params=[]
  for key in keys:
   k=groups[key][0];i,j=divmod(int(ids[k]),3000);params.append((k,float(u[i]),float(v[j])))
  vals=[]
  with concurrent.futures.ProcessPoolExecutor(max_workers=1 if n>=16384 else 3,initializer=r.init,initargs=(n,)) as pool:
   for count,(k,h) in enumerate(pool.map(r.job,params,chunksize=1),1):
    kk=groups[keys[count-1]];delta=abs(h-best[kk]);good=np.isfinite(h);conv[kk]=good&(delta<=.001)
    if good:best[kk]=h;orders[kk]=n
    vals.append((k,h));print(f'N={n} matrix {count}/{len(params)} h={h:.7g} elapsed {time.time()-t:.1f}s',flush=True)
  np.savez(r.ROOT/f'entropy_extended_{n}.npz',ids=ids,best=best,orders=orders,converged=conv)
  print(f'N={n}: remaining {(~conv).sum()}',flush=True)
  if conv.all():break
 raw=d['h_refined'].copy();raw.ravel()[ids]=best;remaining=np.zeros_like(raw,dtype=bool);remaining.ravel()[ids[~conv]]=True;display=raw.copy();display[remaining]=.01
 np.savez(r.ROOT/'entropy_boundary_refined_3000.npz',u=u,v=v,h_refined=raw,h_display=display,remaining_mask=remaining,target_mask=d['target_mask'],ids=ids,orders=orders,converged=conv)
 display.tofile(r.ROOT/'entropy_display_3000.bin')
 report={'target_cells':len(ids),'converged':int(conv.sum()),'remaining':int(remaining.sum()),'tolerance':.001,'minimum_order':512,'max_order':int(orders.max()),'order_counts':{str(n):int((orders==n).sum()) for n in np.unique(orders)}}
 (r.ROOT/'entropy-refinement.json').write_text(json.dumps(report,indent=2));print(report,flush=True)
