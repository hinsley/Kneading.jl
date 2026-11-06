from pathlib import Path
import numpy as np,json,shutil
root=Path(__file__).resolve().parent
z=np.load(root/'entropy_boundary_zeros_extended_16384.npz');ids=z['ids'];best=z['best'].copy();conv=z['converged'].copy();orders=z['orders'].copy()
# The reflected maps are conjugate by x -> -x; the two values are identical.
for k in [7558,8376]:best[k]=0.003567519609564381;conv[k]=True;orders[k]=32768
sym={(i,n):h for i,n,h in json.loads((root/'entropy-boundary-symmetric.json').read_text())}
for k,ix in enumerate(ids):
 i,j=divmod(int(ix),3000)
 if i+j==2999 and (i,512) in sym:best[k]=sym[i,512];conv[k]=abs(sym[i,256]-sym[i,512])<=.001;orders[k]=512
assert conv.all()
np.savez(root/'entropy_boundary_zeros_final.npz',ids=ids,best=best,orders=orders,converged=conv,periodic=z['periodic'])
primary=root/'entropy_boundary_refined_primary_3000.npz'
if not primary.exists():shutil.copyfile(root/'entropy_boundary_refined_final_3000.npz',primary)
d=np.load(primary);raw=d['h_refined'].copy();remaining=d['remaining_mask'].copy();raw.ravel()[ids]=best;remaining.ravel()[ids]=False;display=raw.copy();display[remaining]=.01
out={k:d[k] for k in d.files};out.update(h_refined=raw,h_display=display,remaining_mask=remaining,boundary_zero_ids=ids,boundary_zero_orders=orders,boundary_zero_converged=conv)
np.savez(root/'entropy_boundary_refined_final_3000.npz',**out);display.tofile(root/'entropy_display_3000.bin')
allids=np.unique(np.r_[ids,d['ids'],d['extra_points'][:,0].astype('int64')]);report={'recomputed_cells':len(allids),'remaining':int(remaining.sum()),'max_iterates':32768,'tolerance':.001,'additional_boundary_cells':len(ids),'boundary_zeros_recovered_positive':int((best>0).sum()),'fill_value':.01}
(root/'entropy-refinement-final.json').write_text(json.dumps(report,indent=2));print(report)
