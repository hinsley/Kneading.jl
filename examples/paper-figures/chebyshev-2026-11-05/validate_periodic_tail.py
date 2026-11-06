import periodic_tail_refine as q
import numpy as np,ctypes,pickle,json
records=pickle.load(open(q.r.ROOT/'entropy_periodic_tails.pkl','rb'));d=np.load(q.r.ROOT/'entropy_boundary_refined_3000.npz');ids=d['ids'];rng=np.random.default_rng(20261105);good=[x[0] for x in records if x[2] is not None]
for k in rng.choice(good,32,replace=False):
 i,j=divmod(int(ids[k]),3000);u,v=d['u'][i],d['v'][j];a=(ctypes.c_longlong*8192)();m=(ctypes.c_int*4)();n=q.lib.periodic_numerator(u,v,32768,8192,a,m);s=[int(a[t]) if t<n else 0 for t in range(513)]
 for step in [1,m[1],m[3]]:
  for t in range(step,513):s[t]+=s[t-step]
 b=(ctypes.c_longlong*513)();q.r.lib.refine_coeff(u,v,512,b);assert s==list(b),(k,m[:])
print('32 periodic-tail rational determinants reproduce all 513 finite-series coefficients exactly.')
