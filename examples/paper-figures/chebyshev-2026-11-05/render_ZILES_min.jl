using CairoMakie, LaTeXStrings, JLD2
using Kneading.Diagrams
using Printf

const N = parse(Int, get(ENV, "GRID_SIZE", "3000"))
const ENTRIES = 20
const OUTPUT = get(ENV, "FIGURE_OUTPUT", "ZILESmin.png")
const DRAW_SIZE = 1000

function compute_contours(u,v)
    plane=ParameterPlane(u,v;xname="u",yname="v")
    diagram=KneadingDiagram(plane;metadata=(family=:chebyshev_cubic,critical_points=(-0.5,0.5),iterates=ENTRIES))
    values=Array{Float64}(undef,length(v),length(u),2)
    values[:,:,1].=-0.5; values[:,:,2].=0.5
    for iterate in 2:ENTRIES
        scan_plane!(values,plane) do vals,u,v
            vals[1]=((u-v)/2)*(4*vals[1]^3-3*vals[1])+(u+v)/2
            vals[2]=((u-v)/2)*(4*vals[2]^3-3*vals[2])+(u+v)/2
        end
        for (s,source) in enumerate((:left_critical,:right_critical)), level in (-0.5,0.5)
            add_level_contours!(diagram,view(values,:,:,s);source,iterate,level)
        end
        @printf("Orbit entry %d/%d; %d segments so far\n",iterate,ENTRIES,sum(length(l.segments) for l in diagram.layers)); flush(stdout)
    end
    jldsave("contours_$(N).jld2";diagram)
    diagram
end

function main()
    d=load("lyapunov_$(N).jld2")
    u,v=d["u"],d["v"]
    lambda=min.(d["lambda_left"],d["lambda_right"])
    @assert all(lambda .<= d["lambda_max"])
    @assert !any(isnan,lambda)
    @assert size(lambda)==(N,N)
    diagram=isfile("contours_$(N).jld2") ? load("contours_$(N).jld2","diagram") : compute_contours(u,v)
    fig=Figure(size=(DRAW_SIZE,DRAW_SIZE),backgroundcolor=:white,figure_padding=(12,12,10,24))
    ax=Axis(fig[1,1];aspect=DataAspect(),xlabel=L"u",ylabel=L"v",
        xlabelsize=42,ylabelsize=42,xticklabelsize=30,yticklabelsize=30,
        xticks=([0,.5,1],[L"0",L"0.5",L"1"]),
        yticks=([-1,-.5,0],[L"-1",L"-0.5",L"0"]),
        xgridvisible=false,ygridvisible=false,
        xticksize=5,yticksize=5,xtickalign=1,ytickalign=1,
        xticklabelpad=10,yticklabelpad=10,
        spinewidth=1,xtickwidth=1,ytickwidth=1,
        xlabelpadding=6,ylabelpadding=6)
    # The page uses tanh(lambda / 0.5): retain it, with true-lambda labels.
    cmap=cgrad([RGBf(0x24/255,0xc7/255,0x7a/255),RGBf(0,0,0),RGBf(0xf2/255,0xc9/255,0x4c/255)], [0,0.5,1])
    heatmap!(ax,u,v,permutedims(tanh.(lambda ./ 0.5));colormap=cmap,colorrange=(-1,1),interpolate=false,rasterize=true)
    colors=Dict(:left_critical=>RGBf(.85,.15,.12),:right_critical=>RGBf(.10,.30,.90))
    # Exactly the package extension's ordering, color, opacity, and width.
    for layer in sort(diagram.layers;by=l->l.iterate,rev=true)
        isempty(layer.segments)&&continue
        points=Vector{Point2f}(undef,2length(layer.segments))
        for (i,s) in enumerate(layer.segments)
            points[2i-1]=Point2f(s.x1,s.y1);points[2i]=Point2f(s.x2,s.y2)
        end
        linesegments!(ax,points;color=(colors[layer.source],Float64(layer.iterate-1)^(-1.2)),linewidth=2.0)
    end
    xlims!(ax,0,1);ylims!(ax,-1,0)
    Colorbar(fig[1,2];colormap=cmap,limits=(-1,1),width=20,
        label=L"\lambda_{\min}",labelsize=40,ticklabelsize=28,labelpadding=8,
        ticks=(tanh.([-0.5,-0.25,0,0.25,0.5] ./ 0.5),[L"-0.5",L"-0.25",L"0",L"0.25",L"0.5"]),
        ticksize=5,tickalign=1,spinewidth=1,tickwidth=1,ticklabelpad=10,
        height=Relative(1.0))
    colgap!(fig.layout,18)
    # Explicit layout keeps the colorbar aligned with the square plot panel.
    rowsize!(fig.layout,1,Aspect(1,1.0))
    resize_to_layout!(fig)
    println("Saving CairoMakie figure...");flush(stdout)
    save(OUTPUT,fig;px_per_unit=3)
    save("preview_ZILESmin.png",fig;px_per_unit=1)
    println("Saved ",OUTPUT,"; ",length(diagram.layers)," layers, ",sum(length(l.segments) for l in diagram.layers)," segments")
end
main()
