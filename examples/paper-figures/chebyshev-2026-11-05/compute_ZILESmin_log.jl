using CUDA, JLD2, Random, Statistics, Printf

const N = parse(Int, get(ENV, "GRID_SIZE", "3000"))
const BURN = parse(Int, get(ENV, "LYAP_BURN", "8192"))
const STEPS = parse(Int, get(ENV, "LYAP_STEPS", "16384"))

# Same map as Kneading.jl's Chebyshev tutorial, without fast-math.
@inline cubic(u, v, x) = ((u-v)/2) * (4*x*x*x - 3*x) + (u+v)/2
@inline deriv(u, v, x) = (u-v) * (6*x*x - 1.5)

function lyapunov_kernel!(outleft, outright, us, vs, burn, steps)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    n = length(us)
    while i <= length(outleft)
        j = (i-1) % n + 1
        k = (i-1) ÷ n + 1
        u, v = us[k], vs[j]
        # First critical images are exactly u and v; never include f'(±1/2)=0.
        xl, xr = u, v
        for t in 1:burn
            xl = cubic(u,v,xl)
            xr = cubic(u,v,xr)
        end
        sl, sr = 0.0, 0.0
        for t in 1:steps
            sl += log(abs(deriv(u,v,xl)))
            sr += log(abs(deriv(u,v,xr)))
            xl = cubic(u,v,xl)
            xr = cubic(u,v,xr)
        end
        outleft[i] = sl / steps
        outright[i] = sr / steps
        i += stride
    end
    return
end

function cpu_pair(u,v; burn=BURN, steps=STEPS)
    xl,xr=u,v
    for _ in 1:burn
        xl=cubic(u,v,xl); xr=cubic(u,v,xr)
    end
    sl=sr=0.0
    for _ in 1:steps
        sl+=log(abs(deriv(u,v,xl))); sr+=log(abs(deriv(u,v,xr)))
        xl=cubic(u,v,xl); xr=cubic(u,v,xr)
    end
    (sl/steps,sr/steps)
end

function main()
    CUDA.device!(0)
    @assert occursin("4070 Ti", CUDA.name(CUDA.device()))
    println("GPU: ", CUDA.name(CUDA.device())); flush(stdout)
    u=exp.(range(log(0.01),0.0;length=N)); v=-reverse(u)
    left=CUDA.zeros(Float64,N,N); right=similar(left)
    du,dv=CuArray(u),CuArray(v)
    # Compile and validate on a tiny grid before starting the production scan.
    tinyu=CuArray([0.0,0.2,0.8,1.0]); tinyv=CuArray([-1.0,-0.8,-0.2,0.0])
    tl=CUDA.zeros(Float64,4,4); tr=similar(tl)
    CUDA.@sync @cuda threads=64 lyapunov_kernel!(tl,tr,tinyu,tinyv,BURN,STEPS)
    al,ar=Array(tl),Array(tr)
    @assert max(al[1,4],ar[1,4]) ≈ log(9.0) atol=1e-10
    @assert al[4,1] == -Inf && ar[4,1] == -Inf
    # Stable fixed point for (u,v)=(0.2,-0.2): x=0, f'(0)=-0.6.
    @assert al[3,2] ≈ log(0.6) atol=1e-10
    @assert ar[3,2] ≈ log(0.6) atol=1e-10
    println("Analytic checks passed: zero map, attracting fixed point, Chebyshev endpoint orbits.")
    t=@elapsed CUDA.@sync @cuda threads=256 blocks=min(cld(N*N,256),4096) lyapunov_kernel!(left,right,du,dv,BURN,STEPS)
    println("GPU scan seconds: ",t); flush(stdout)
    lambda_left=Array(left); lambda_right=Array(right)
    lambda_max=max.(lambda_left,lambda_right)
    @assert !any(isnan,lambda_max)
    @assert !any(==(Inf),lambda_max)
    jldsave("lyapunov_log_$(N).jld2";u,v,lambda_left,lambda_right,lambda_max,burn=BURN,steps=STEPS,gpu=CUDA.name(CUDA.device()),scan_seconds=t)
    println("Saved scan; finite extrema: ", extrema(filter(isfinite,lambda_max)))
    # Independent CPU comparisons and longer averages on a reproducible subset.
    rng=MersenneTwister(20261105)
    ids=rand(rng,1:N*N,256)
    errors=Float64[]; convergence=Float64[]
    for i in ids
        j=(i-1)%N+1; k=(i-1)÷N+1
        a,b=cpu_pair(u[k],v[j])
        x=min(a,b); y=min(lambda_left[i],lambda_right[i])
        if isfinite(x)&&isfinite(y); push!(errors,abs(x-y)); end
        c,d=cpu_pair(u[k],v[j];burn=2BURN,steps=2STEPS)
        z=min(c,d)
        if isfinite(x)&&isfinite(z); push!(convergence,abs(tanh(x/0.5)-tanh(z/0.5))); end
    end
    open("validation_log.txt","w") do io
        println(io,"Grid: $N x $N; burn-in $BURN; averaging $STEPS; Float64, natural log")
        println(io,"GPU: ", CUDA.name(CUDA.device()))
        println(io,"Analytic checks: passed (zero map, stable fixed point, Chebyshev endpoint orbits)")
        println(io,"CPU/GPU |delta lambda| on 256 sample points: median=$(median(errors)), 95th percentile=$(quantile(errors,.95)), max=$(maximum(errors))")
        println(io,"Doubling burn-in and averaging: |delta tanh(lambda/0.5)| median=$(median(convergence)), 95th percentile=$(quantile(convergence,.95)), max=$(maximum(convergence))")
        println(io,"Finite lambda_max extrema: ",extrema(filter(isfinite,lambda_max)))
        println(io,"Negative-infinite values: ",count(==(-Inf),lambda_max))
        println(io,"NaN/+Inf values: 0")
    end
    print(read("validation_log.txt",String))
end
main()
