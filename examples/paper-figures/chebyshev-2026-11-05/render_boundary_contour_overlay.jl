using CairoMakie, JLD2, Kneading.Diagrams
include("cairo_bounded_segments.jl")
function main()
 diagram=load("contours_boundary_log_3000.jld2","diagram")
 # Same physical panel dimensions and output density as Figure 6's lower panel.
 fig=Figure(size=(669,669),backgroundcolor=:transparent,figure_padding=0)
 ax=Axis(fig[1,1];backgroundcolor=:transparent,width=669,height=669)
 hidedecorations!(ax);hidespines!(ax)
 colors=Dict(:left_critical=>RGBf(.85,.15,.12),:right_critical=>RGBf(.1,.3,.9))
 layers=sort(diagram.layers;by=l->l.iterate,rev=true);empty!(diagram.layers);diagram=nothing;GC.gc()
 while !isempty(layers)
  layer=popfirst!(layers);isempty(layer.segments)&&continue
  pts=Vector{Point2f}(undef,2length(layer.segments))
  for (i,s) in enumerate(layer.segments);pts[2i-1]=Point2f(s.x1,s.y1);pts[2i]=Point2f(s.x2,s.y2);end
  linesegments!(ax,pts;color=(colors[layer.source],Float64(layer.iterate-1)^(-.9)),linewidth=2*Float64(layer.iterate-1)^(-.5))
  layer=nothing;GC.gc();println("Layers remaining: ",length(layers));flush(stdout)
 end
 xlims!(ax,0,-log(.01));ylims!(ax,log(.01),0)
 println("Rendering transparent contour cache");flush(stdout)
 save("boundary_contours_opacity09.png",fig;px_per_unit=3)
end
main()
