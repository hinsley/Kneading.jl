using Kneading.OneDimensionalMaps, Random
rng=MersenneTwister(20261105)
params=[(0.,0.),(.2,-.2),(.5,-.5),(1.,-1.),(.5,-1.),(1.,-.5),(.8,-.8),(0.,-.9),(.99,0.),(.99,-.99)]
append!(params,[(rand(rng),-rand(rng)) for _ in 1:100])
open("entropy-reference-polynomials.txt","w") do io
 for (u,v) in params
  u==v&&continue
  f(x)=((u-v)/2)*(4*x*x*x-3*x)+(u+v)/2
  m=PartitionedIntervalMap(f,LapPartition((-Inf,Inf),[-.5,.5]);orientations=[Increasing,Decreasing,Increasing])
  kd=kneading_determinant(kneading_matrix(kneading_data(m,64));deleted_lap=1)
  p=Int.(kd.numerator.coefficients);resize!(p,65)
  for n in 2:65;p[n]+=p[n-1];end
  lines=split(read(`./entropy_boundary point $u $v`,String),'\n')
  q=parse.(Int,split(lines[1]));@assert p==q (u,v,p,q)
  println(io,u," ",v," ",join(p,' '))
 end
end
println("All 109 nonconstant degree-64 polynomials match Kneading.jl exactly, including exact critical hits.")
