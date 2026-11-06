from pathlib import Path
import numpy as np, ctypes, math,json,concurrent.futures
from flint import fmpz_poly,ctx
ROOT=Path(__file__).resolve().parent
lib=ctypes.CDLL(str(ROOT/'entropy_boundary.so'));lib.entropy_init()
lib.entropy_coefficients.argtypes=[ctypes.c_double,ctypes.c_double,ctypes.POINTER(ctypes.c_int)]
lib.entropy_value.argtypes=[ctypes.POINTER(ctypes.c_int),ctypes.c_int];lib.entropy_value.restype=ctypes.c_double
ctx.prec=128
N=3000
axes=np.fromfile(ROOT/f'entropy_axes_{N}.bin',dtype='float64');u,v=axes[:N],axes[N:]
h=np.fromfile(ROOT/f'entropy_raw_{N}.bin',dtype='float64').reshape(2,N,N) # [order,u,v]
def coeff(i,j):
 a=(ctypes.c_int*65)();lib.entropy_coefficients(u[i],v[j],a);return tuple(a)
def solve(p):
 ctx.prec=128
 poly=fmpz_poly(list(p))
 while poly.degree()>0 and poly(1)==0:poly=poly//fmpz_poly([-1,1])
 roots=poly.complex_roots()
 real=[r.real for r,m in roots if r.imag.is_zero() and r.real>0 and r.real<1]
 return -math.log(float(min(real))) if real else 0.
def job(p):
 try:return p,solve(p)
 except Exception as e:return p,float('nan')
if __name__=='__main__':
 bad=np.argwhere(~np.isfinite(h));tasks={}
 for order,i,j in bad:
  p=coeff(i,j)[:65 if order==0 else 33];tasks.setdefault(p,[]).append((order,i,j))
 print(f'{len(bad)} ambiguous estimates; {len(tasks)} distinct polynomials',flush=True)
 with concurrent.futures.ProcessPoolExecutor(max_workers=6) as pool:
  for count,(p,val) in enumerate(pool.map(job,tasks,chunksize=4),1):
   for ix in tasks[p]:h[ix]=val
   if count%100==0:print(f'Arb {count}/{len(tasks)}',flush=True)
 rng=np.random.default_rng(20261105);errors=[]
 # Independent Arb isolation validates 300 production polynomials at both orders.
 samples=[]
 for i,j in rng.integers(0,N,size=(300,2)):
  p=coeff(i,j)
  for k,n in [(0,65),(1,33)]:samples.append((p[:n],h[k,i,j]))
 with concurrent.futures.ProcessPoolExecutor(max_workers=6) as pool:
  for (p,val),(_,expected) in zip(pool.map(job,[p for p,_ in samples],chunksize=4),samples):
   if not math.isfinite(val):raise RuntimeError('Arb validation failed')
   errors.append(abs(val-expected))
 assert max(errors)<1e-9,max(errors)
 failure=~np.isfinite(h[0]);sensitive=(~np.isfinite(h[1]))|(np.abs(h[0]-h[1])>.01)
 mask=failure|sensitive
 filled=h[0].copy();filled[mask]=.01
 np.savez(ROOT/'entropy_boundary_log_3000.npz',u=u,v=v,h64=h[0],h32=h[1],failure_mask=failure,sensitivity_mask=sensitive,fill_mask=mask,h_display=filled,fill_value=.01)
 # For Julia: byte order [v,u], same as its cached Lyapunov matrices.
 filled.tofile(ROOT/'entropy_display_3000.bin')
 report={'grid':[N,N],'orders':[64,32],'ambiguous_estimates':len(bad),'distinct_arb_fallbacks':len(tasks),'unresolved':int(failure.sum()),'sensitive_cells':int(sensitive.sum()),'filled_cells':int(mask.sum()),'fill_value':.01,'validation_polynomials':600,'max_arb_difference':max(errors),'h64_range':[float(np.nanmin(h[0])),float(np.nanmax(h[0]))]}
 (ROOT/'entropy-validation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)
