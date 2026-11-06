using JLD2

# Solve the defining condition, rather than iterating a rounded parameter pair.
function compute_homoclinic()
    setprecision(384) do
        v = -BigFloat(1)/2
        f(x,u) = (u-v)/2*(4x^3-3x)+(u+v)/2
        condition(u) = f(f(u,u),u)-f(u,u)
        lo, hi = big"0.83", big"0.84"
        @assert condition(lo)*condition(hi)<0
        for _ in 1:390
            mid=(lo+hi)/2
            if sign(condition(mid))==sign(condition(lo)); lo=mid; else; hi=mid; end
        end
        u=(lo+hi)/2
        star=f(u,u)
        @assert abs(f(star,u)-star)<big"1e-110"
        @assert -BigFloat(1)/2 < star < BigFloat(1)/2
        slope=(u-v)*(6star^2-BigFloat(3)/2)
        @assert slope < -1
        # Invert only the middle (decreasing) lap. The first preimage is exactly
        # +1/2 since v=-1/2. Repeated inverse images then converge to the repeller.
        back=BigFloat[-BigFloat(1)/2,BigFloat(1)/2]
        for _ in 2:80
            target=back[end]
            a,b=-BigFloat(1)/2,BigFloat(1)/2
            for _ in 1:390
                mid=(a+b)/2
                if f(mid,u)>target; a=mid; else; b=mid; end
            end
            push!(back,(a+b)/2)
        end
        # A compact trapping interval for the inverse branch proves the omitted
        # backward tail converges: g maps [-1/5,3/10] into itself and |g'|<1.
        a,b=-BigFloat(1)/5,BigFloat(3)/10
        min_slope=(u-v)*(BigFloat(3)/2-6max(abs(a),abs(b))^2)
        @assert f(a,u)>b && f(b,u)<a && min_slope>1
        @assert any(x->a<x<b,back)
        hom=vcat(reverse(back),[u,star,star])
        residual=maximum(abs(f(hom[k],u)-hom[k+1]) for k in 1:length(hom)-1)
        @assert residual<big"1e-110"
        @assert abs(back[end]-star)<big"1e-22"
        # Same map, nearby initial condition; no artificial post-hit drift in blue.
        seed=-BigFloat(1)/2+big"1e-6"
        n=5000
        red=Vector{BigFloat}(undef,n+1);red[1]=seed
        for k in 1:n; red[k+1]=f(red[k],u); end
        @assert all(x->v<=x<=u,red)
        lyap=sum(log(abs((u-v)*(6x^2-BigFloat(3)/2))) for x in red[1001:end-1])/(n-1000)
        @assert lyap>0
        bins=length(unique(clamp(floor(Int,200*(x-v)/(u-v)),0,199) for x in red[1001:end]))
        @assert bins>=190
        output=joinpath(@__DIR__,"cubic_homoclinic_data.jld2")
        jldsave(output; u=Float64(u),v=Float64(v),xstar=Float64(star),slope=Float64(slope),
            homoclinic=Float64.(hom),chaotic=Float64.(red),u_decimal=string(u),
            xstar_decimal=string(star),seed_decimal=string(seed),precision_bits=384,
            residual=string(residual),backward_distance=string(abs(back[end]-star)),
            lyapunov=Float64(lyap),occupied_bins=bins,n)
        open(joinpath(@__DIR__,"cubic_homoclinic_validation.txt"),"w") do io
            println(io,"u = ",u,"\nv = ",v,"\nxstar = ",star,"\nmultiplier = ",slope)
            println(io,"maximum orbit residual = ",residual,"\nbackward-tail distance = ",abs(back[end]-star))
            println(io,"red seed = ",seed,"\nred finite-time Lyapunov exponent = ",lyap,"\noccupied bins = ",bins," / 200")
        end
        println(read(joinpath(@__DIR__,"cubic_homoclinic_validation.txt"),String))
    end
end
compute_homoclinic()
