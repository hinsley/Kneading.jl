using CairoMakie, LaTeXStrings, JLD2

# Fixed-initial-state orbit diagram, as documented by ChaosTools.orbitdiagram:
# https://juliadynamics.github.io/DynamicalSystemsDocs.jl/chaostools/stable/orbitdiagram/
# Direct scalar iteration avoids introducing additional package dependencies.
# Log coordinates and cutoff match ZILESmin: -log(1-u), 0 <= u <= .99.
# Original orbit diagram retained 1024 iterates; this version retains 10240.
# At v=-1/2, f(+1/2)=-1/2, so both critical orbits have the same tail.
const V = -0.5
const NPARAM = 6001
const TRANSIENT = 32768
const NSAVE = 10240
@inline cubic(x, u, v) = ((u-v)/2)*(4*x*x*x-3*x)+(u+v)/2

function main()
    d=load(joinpath(@__DIR__,"orbit_diagram_v_minus_half_log_10x.jld2"))
    us,logus,states=d["u"],d["log_u"],d["states"]
    @assert size(states)==(NSAVE,NPARAM)
    println("Loaded all $(length(states)) points; rendering compressed upper panel at alpha=0.025.");flush(stdout)
    pts = Vector{Point2f}(undef, length(states))
    for j in eachindex(us), k in 1:NSAVE
        pts[(j-1)*NSAVE+k] = Point2f(logus[j],states[k,j])
    end
    fig=Figure(size=(1000,270),backgroundcolor=:white,figure_padding=(16,24,16,16))
    ax=Axis(fig;bbox=BBox(440.5/3,2448.5/3,70,255),xlabel=L"u",ylabel=L"x",
        xlabelvisible=false,xticksvisible=false,xticklabelsvisible=false,
        yaxisposition=:left,yticksvisible=true,yticklabelsvisible=true,
        xlabelsize=42,ylabelsize=42,xticklabelsize=30,yticklabelsize=30,
        xticks=(-log1p.(-[0.0,0.5,0.8,0.9,0.95,0.98,0.99]),[L"0",L"0.5",L"0.8",L"0.9",L"0.95",L"0.98",L"0.99"]),
        yticks=([-.5,0,.5,1],[L"-0.5",L"0",L"0.5",L"1"]),
        xgridvisible=false,ygridvisible=false,
        xticklabelpad=12,yticklabelpad=12,xlabelpadding=8,ylabelpadding=8,
        xtickalign=1,ytickalign=0)
    scatter!(ax,pts;markersize=0.8,color=(:black,0.025),rasterize=true)
    xlims!(ax,0,-log(0.01)); ylims!(ax,-.8,1.3)
    save(joinpath(@__DIR__,"Figure6_orbit_upper_panel.png"),fig;px_per_unit=3)
    println("Axis pixel viewport: ",ax.scene.viewport[])
    println("Saved compressed upper panel.")
end
main()
