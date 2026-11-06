import entropy_symmetric_driver as s
import concurrent.futures,json
if __name__=='__main__':
 out=[]
 with concurrent.futures.ProcessPoolExecutor(max_workers=4) as pool:
  for x in pool.map(s.solve,[(i,n) for i in [950,961,962,965] for n in [256,512]]):out.append(x);print(x,flush=True)
 (s.r.ROOT/'entropy-boundary-symmetric.json').write_text(json.dumps(out))
