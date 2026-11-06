from pathlib import Path
import numpy as np,pickle,json
root=Path(__file__).resolve().parent
d=np.load(root/'entropy_boundary_refined_3000.npz');x=np.load(root/'entropy_extended_8192.npz');ids=d['ids'];best=x['best'].copy();orders=x['orders'].copy();conv=x['converged'].copy();method=np.ones(len(ids),dtype='uint8')
records=pickle.load(open(root/'entropy_periodic_tails.pkl','rb'));roots=pickle.load(open(root/'entropy_periodic_roots.pkl','rb'))
for k,n,p,m in records:
 if p is not None:best[k]=roots[p];orders[k]=32768;conv[k]=True;method[k]=2
last=np.load(root/'entropy_last_edge_32768.npz'); previous_last=np.load(root/'entropy_last_edge_16384.npz')
for j,k in enumerate(last['active']):best[k]=last['best'][j];conv[k]=last['converged'][j];orders[k]=16384 if previous_last['converged'][j] else 32768
sym=json.loads((root/'entropy-symmetric-factors.json').read_text());sv={(i,n):h for i,n,h in sym};imap={int(ix):k for k,ix in enumerate(ids)};extra=[]
for i in sorted(set(i for i,n,h in sym)):
 ix=i*3000+2999-i;a,b=sv[i,256],sv[i,512]
 if ix in imap:
  k=imap[ix];best[k]=b;conv[k]=abs(a-b)<=.001;orders[k]=512;method[k]=3
 else:extra.append((ix,b,abs(a-b)<=.001))
raw=d['h_refined'].copy();raw.ravel()[ids]=best;remaining=np.zeros(raw.shape,dtype=bool);remaining.ravel()[ids[~conv]]=True
for ix,h,ok in extra:raw.ravel()[ix]=h;remaining.ravel()[ix]=not ok
# Only the unresolved cells retain the requested small positive fallback.
display=raw.copy();display[remaining]=.01
np.savez(root/'entropy_boundary_refined_final_3000.npz',u=d['u'],v=d['v'],h_refined=raw,h_display=display,remaining_mask=remaining,ids=ids,orders=orders,converged=conv,method=method,extra_points=np.array(extra))
display.tofile(root/'entropy_display_3000.bin')
report={'target_cells':len(ids)+len(extra),'passed_check':int(conv.sum())+sum(ok for ix,h,ok in extra),'remaining':int(remaining.sum()),'remaining_grid_indices':np.argwhere(remaining).tolist(),'tolerance':.001,'maximum_iterates':32768,'periodic_tail_cells':int((method==2).sum()),'symmetry_factor_cells':int((method==3).sum())+len(extra),'fill_value':.01}
(root/'entropy-refinement-final.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
