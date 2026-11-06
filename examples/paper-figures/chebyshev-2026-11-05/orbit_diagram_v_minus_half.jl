using CairoMakie, LaTeXStrings, JLD2

# Fixed-initial-state orbit diagram, as documented by ChaosTools.orbitdiagram:
# https://juliadynamics.github.io/DynamicalSystemsDocs.jl/chaostools/stable/orbitdiagram/
# Direct scalar iteration avoids introducing additional package dependencies.
# At v=-1/2, f(+1/2)=-1/2, so both critical orbits have the same tail.
const V = -0.5
const NPARAM = 6001
const TRANSIENT = 32768
const NSAVE = 1024
@inline cubic(x, u, v) = ((u-v)/2)*(4*x*x*x-3*x)+(u+v)/2

function main()
    us = collect(range(0.0, 1.0; length=NPARAM))
    states = Matrix{Float64}(undef, NSAVE, NPARAM)
    Threads.@threads for j in eachindex(us)
        x = -0.5
        u = us[j]
        for _ in 1:TRANSIENT
            x = cubic(x,u,V)
        end
        for k in 1:NSAVE
            x = cubic(x,u,V)
            states[k,j] = x
        end
    end
    @assert all(isfinite, states)
    @assert minimum(states) >= -0.5-1e-12 && maximum(states) <= 1+1e-12
    @assert all(u -> isapprox(cubic(0.5,u,V),-0.5;atol=2eps(Float64)), us)
    @assert all(u -> isapprox(cubic(-0.5,u,V),u;atol=2eps(Float64)), us)
    @assert all(x -> abs(x) == 0.5, states[:,(NPARAM+1)÷2]) # superstable 2-cycle at u=1/2
    @assert all(==(1.0), states[:,end]) # u=1: -1/2 -> 1 -> 1.
    jldsave(joinpath(@__DIR__,"orbit_diagram_v_minus_half.jld2");
        u=us, v=V, states, initial_condition=-0.5, transient=TRANSIENT, saved=NSAVE)
    println("Computed $(length(states)) points; range $(extrema(states)). Rendering..."); flush(stdout)
    pts = Vector{Point2f}(undef, length(states))
    for j in eachindex(us), k in 1:NSAVE
        pts[(j-1)*NSAVE+k] = Point2f(us[j],states[k,j])
    end
    fig=Figure(size=(1500,900),backgroundcolor=:white,figure_padding=(16,24,16,16))
    ax=Axis(fig[1,1];xlabel=L"u",ylabel=L"x",
        xlabelsize=38,ylabelsize=38,xticklabelsize=28,yticklabelsize=28,
        xticks=([0,.2,.4,.6,.8,1],[L"0",L"0.2",L"0.4",L"0.6",L"0.8",L"1"]),
        yticks=([-.5,0,.5,1],[L"-0.5",L"0",L"0.5",L"1"]),
        xgridvisible=false,ygridvisible=false,
        xticklabelpad=12,yticklabelpad=12,xlabelpadding=8,ylabelpadding=8,
        xtickalign=1,ytickalign=1)
    scatter!(ax,pts;markersize=0.8,color=(:black,0.06),rasterize=true)
    xlims!(ax,0,1); ylims!(ax,-.52,1.02)
    save(joinpath(@__DIR__,"orbit_diagram_v_minus_half.png"),fig;px_per_unit=2)
    save(joinpath(@__DIR__,"orbit_diagram_v_minus_half.pdf"),fig)
    println("Saved PNG, PDF, data, and reproducible source.")
end
main()
