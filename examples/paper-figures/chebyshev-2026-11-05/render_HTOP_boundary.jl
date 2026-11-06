using CairoMakie, LaTeXStrings, JLD2
function main()
 n=3000;d=load("lyapunov_boundary_log_3000.jld2");u=-log1p.(-d["u"]);v=log1p.(d["v"])
 h=reshape(reinterpret(Float64,read("entropy_display_3000.bin")),n,n)
 fig=Figure(size=(970,791),backgroundcolor=:white)
 ax=Axis(fig;bbox=BBox(147,816,98,767),xlabel=L"u",ylabel=L"v",
  xlabelsize=42,ylabelsize=42,xticklabelsize=30,yticklabelsize=30,
  xticks=(-log1p.(-[0.,.5,.8,.9,.95,.98,.99]),[L"0",L"0.5",L"0.8",L"0.9",L"0.95",L"0.98",L"0.99"]),
  yticks=(log1p.([-.99,-.98,-.95,-.9,-.8,-.5,0.]),[L"-0.99",L"-0.98",L"-0.95",L"-0.9",L"-0.8",L"-0.5",L"0"]),
  xgridvisible=false,ygridvisible=false,xticksize=5,yticksize=5,xtickalign=1,ytickalign=1,
  xticklabelpad=10,yticklabelpad=10,spinewidth=1,xtickwidth=1,ytickwidth=1,xlabelpadding=6,ylabelpadding=6)
 cmap=cgrad([RGBf(0,0,0),RGBf(0xf2/255,0xc9/255,0x4c/255)])
 heatmap!(ax,u,v,permutedims(h);colormap=cmap,colorrange=(0,log(3)),interpolate=false,rasterize=true)
 xlims!(ax,0,-log(.01));ylims!(ax,log(.01),0)
 Colorbar(fig;bbox=BBox(834,854,98,767),colormap=cmap,limits=(0,log(3)),
  label=L"h_{\mathrm{top}}",labelsize=40,ticklabelsize=28,labelpadding=4,
  ticks=([0.,.25,.5,.75,1.],[L"0",L"0.25",L"0.5",L"0.75",L"1"]),
  ticksize=5,tickalign=1,spinewidth=1,tickwidth=1,ticklabelpad=10)
 save("HTOP_boundary_base.png",fig;px_per_unit=3)
end
main()
