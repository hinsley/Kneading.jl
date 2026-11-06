# Same CairoMakie projection, clipping and per-segment stroke, with temporary
# buffers bounded to 200,000 points rather than retained for every plot.
@eval CairoMakie begin
    function draw_atomic(::Scene, screen::Screen, plot::LineSegments)
        attr=plot.attributes
        Makie.compute_colors!(attr)
        Makie.add_computation!(attr,Val(:uniform_clip_planes),:clip)
        points=attr.positions_transformed_f32c[]
        transform=attr.projectionview[] * attr.model_f32c[]
        colors=attr.computed_color[]
        linewidth=attr.linewidth[]
        for first in 1:200000:length(points)
            last=min(first+199999,length(points))
            projected=map(view(points,first:last)) do point
                transform * to_ndim(Vec4d,to_ndim(Vec3d,point,0),1)
            end
            cc=colors isa AbstractArray ? view(colors,first:last) : colors
            sp,sc,sw=clip_line_points(projected,cc,linewidth,false,attr.uniform_clip_planes[],attr.resolution[])
            draw_lineplot(screen.context,(clipped_points=sp,clipped_colors=sc,clipped_linewidths=sw,
                linestyle=attr.linestyle[],linecap=attr.linecap[],joinstyle=nothing,
                miter_limit=nothing,is_lines_plot=false))
        end
        return
    end
end
