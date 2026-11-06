using CairoMakie, LaTeXStrings, JLD2

# Homoclinic and nearby chaotic cobwebs at the same parameter point.
# Run compute_cubic_homoclinic.jl once; subsequent style changes use cached data.
function main()
    data=load(joinpath(@__DIR__, "cubic_homoclinic_data.jld2"))
    u,v,xstar=data["u"],data["v"],data["xstar"]
    f(x) = (u-v)/2 * (4x^3-3x) + (u+v)/2
    function cobweb(xs)
        points=Point2d[Point2d(xs[1],xs[1])]
        for k in 1:length(xs)-1
            push!(points,Point2d(xs[k],xs[k+1]),Point2d(xs[k+1],xs[k+1]))
        end
        points
    end
    fig = Figure(size=(900,900), backgroundcolor=:white)
    ticks = ([-.5,0.,.5,1.], [L"-0.5",L"0",L"0.5",L"1"])
    ax = Axis(fig; bbox=BBox(140,875,110,845), aspect=DataAspect(),
        xlabel=L"x", ylabel=L"f_{u,v}(x)",
        xlabelsize=42, ylabelsize=42, xticklabelsize=36, yticklabelsize=36,
        xticks=ticks, yticks=ticks, xgridvisible=false, ygridvisible=false,
        xticksize=5, yticksize=5, xtickalign=1, ytickalign=1,
        xticklabelpad=10, yticklabelpad=10, xlabelpadding=6, ylabelpadding=6,
        spinewidth=1, xtickwidth=1, ytickwidth=1)
    lines!(ax, cobweb(data["chaotic"]); color=(:indianred,0.50), linewidth=0.16)
    lines!(ax, [-0.7,1.0], [-0.7,1.0]; color=:grey60, linestyle=:dash, linewidth=1.3)
    grid = range(-0.7,1.0; length=4001)
    lines!(ax, grid, f.(grid); color=:black, linewidth=2.2)
    lines!(ax, cobweb(data["homoclinic"]); color=(:royalblue3,1.0), linewidth=3.0)
    scatter!(ax, [xstar], [xstar]; color=:black, markersize=14)
    text!(ax, xstar-0.065, xstar+0.025; text=L"x_{u,v}^{\ast}",
        fontsize=36, color=:black, align=(:right,:center))
    limits!(ax, -0.7,1.0,-0.7,1.0)
    # Use the PNG in the paper: PDF viewers may enforce a minimum stroke width
    # and make the densely overlaid sub-pixel cobweb substantially too dark.
    save(joinpath(@__DIR__, "cubic_critical_cobweb.png"), fig; px_per_unit=3)
    save(joinpath(@__DIR__, "cubic_critical_cobweb.pdf"), fig)
    println("Rendered homoclinic and nearby chaotic cobwebs at u=$u, v=$v")
end

main()
