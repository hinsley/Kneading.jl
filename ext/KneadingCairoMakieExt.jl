"""
Kneading.jl's optional CairoMakie extension. It loads only when CairoMakie is
available and provides the plotting implementation for
`Kneading.Diagrams.plot_kneading_contours` without making CairoMakie a core
dependency.
"""
module KneadingCairoMakieExt

using CairoMakie
using Kneading.Diagrams

import Kneading.Diagrams: plot_kneading_contours, save_kneading_contours

function _segment_points(segments)
    points = Point2f[]
    sizehint!(points, 2 * length(segments))
    for segment in segments
        push!(
            points,
            Point2f(segment.x1, segment.y1),
            Point2f(segment.x2, segment.y2),
        )
    end
    return points
end

function plot_kneading_contours(
    diagram::KneadingDiagram;
    colors = Dict{Symbol,Any}(),
    opacity = _ -> 1.0,
    linewidth::Real = 2.0,
    size::Tuple{Int,Int} = (1000, 1000),
    xticks = nothing,
    yticks = nothing,
)
    figure = Figure(size = size, backgroundcolor = :white)
    axis = Axis(
        figure[1, 1];
        aspect = DataAspect(),
        xlabel = diagram.plane.xname,
        ylabel = diagram.plane.yname,
        xlabelsize = 30,
        ylabelsize = 30,
        xticklabelsize = 16,
        yticklabelsize = 16,
        xgridcolor = (:gray, 0.16),
        ygridcolor = (:gray, 0.16),
    )
    xlims!(axis, extrema(diagram.plane.x)...)
    ylims!(axis, extrema(diagram.plane.y)...)
    isnothing(xticks) || (axis.xticks = xticks)
    isnothing(yticks) || (axis.yticks = yticks)

    layers = sort(diagram.layers; by = layer -> layer.iterate, rev = true)
    for layer in layers
        isempty(layer.segments) && continue
        base_color = get(colors, layer.source, RGBf(0.1, 0.1, 0.1))
        linesegments!(
            axis,
            _segment_points(layer.segments);
            color = (base_color, clamp(opacity(layer.iterate), 0.0, 1.0)),
            linewidth,
        )
    end
    return figure
end

function save_kneading_contours(
    path::AbstractString,
    diagram::KneadingDiagram;
    px_per_unit::Real = 1,
    kwargs...,
)
    figure = plot_kneading_contours(diagram; kwargs...)
    save(path, figure; px_per_unit)
    return path
end

end
