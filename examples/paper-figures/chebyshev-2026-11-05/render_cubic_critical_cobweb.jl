using CairoMakie, LaTeXStrings, JLD2

# Reproducible critical-orbit cobweb; no transient iterates are discarded.
function main()
    u, v = 0.83679757, -0.50003822
    n = 5000
    # High precision avoids an early numerical departure near the repelling point.
    orbit, xstar, slope = setprecision(256) do
        ub, vb = parse(BigFloat, "0.83679757"), parse(BigFloat, "-0.50003822")
        f(x) = (ub-vb)/2 * (4x^3-3x) + (ub+vb)/2
        lo, hi = BigFloat(-0.5), BigFloat(0.5)
        for _ in 1:260
            mid = (lo+hi)/2
            if f(mid) > mid
                lo = mid
            else
                hi = mid
            end
        end
        root = (lo+hi)/2
        @assert abs(f(root)-root) < big"1e-70"
        xs = Vector{BigFloat}(undef, n+1)
        xs[1] = BigFloat(-0.5)
        for k in 1:n
            xs[k+1] = f(xs[k])
        end
        @assert all(isfinite, xs)
        (Float64.(xs), Float64(root), Float64((ub-vb)*(6root^2-1.5)))
    end
    jldsave(joinpath(@__DIR__, "cubic_critical_cobweb_data.jld2"); u, v, n, orbit, xstar, slope, precision_bits=256)
    f(x) = (u-v)/2 * (4x^3-3x) + (u+v)/2
    cobweb = Point2d[Point2d(orbit[1], orbit[1])]
    for k in 1:n
        push!(cobweb, Point2d(orbit[k], orbit[k+1]), Point2d(orbit[k+1], orbit[k+1]))
    end

    fig = Figure(size=(900,900), backgroundcolor=:white)
    ticks = ([-1.,-.5,0.,.5,1.], [L"-1",L"-0.5",L"0",L"0.5",L"1"])
    ax = Axis(fig; bbox=BBox(140,875,110,845), aspect=DataAspect(),
        xlabel=L"x", ylabel=L"f_{u,v}(x)",
        xlabelsize=42, ylabelsize=42, xticklabelsize=30, yticklabelsize=30,
        xticks=ticks, yticks=ticks, xgridvisible=false, ygridvisible=false,
        xticksize=5, yticksize=5, xtickalign=1, ytickalign=1,
        xticklabelpad=10, yticklabelpad=10, xlabelpadding=6, ylabelpadding=6,
        spinewidth=1, xtickwidth=1, ytickwidth=1)
    lines!(ax, cobweb; color=(:indianred,0.50), linewidth=0.16)
    lines!(ax, [-1.25,1.25], [-1.25,1.25]; color=:grey60, linestyle=:dash, linewidth=1.3)
    grid = range(-1.25,1.25; length=4001)
    lines!(ax, grid, f.(grid); color=:black, linewidth=2.2)
    scatter!(ax, [xstar], [xstar]; color=:black, markersize=14)
    text!(ax, xstar-0.065, xstar+0.025; text=L"x_{u,v}^{\ast}",
        fontsize=36, color=:black, align=(:right,:center))
    limits!(ax, -1.25,1.25,-1.25,1.25)
    # Use the PNG in the paper: PDF viewers may enforce a minimum stroke width
    # and make the densely overlaid sub-pixel cobweb substantially too dark.
    save(joinpath(@__DIR__, "cubic_critical_cobweb.png"), fig; px_per_unit=3)
    save(joinpath(@__DIR__, "cubic_critical_cobweb.pdf"), fig)
    println("Rendered $n iterates; middle fixed point = $xstar; multiplier = $slope")
end

main()
