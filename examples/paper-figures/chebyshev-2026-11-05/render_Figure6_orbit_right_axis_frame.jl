using CairoMakie, LaTeXStrings
function main()
    fig=Figure(size=(1000,270),backgroundcolor=:white,figure_padding=(16,24,16,16))
    ax=Axis(fig;bbox=BBox(440.5/3,2448.5/3,70,255),xlabel=L"u",ylabel=L"x",
        xlabelvisible=false,xticksvisible=false,xticklabelsvisible=false,
        yaxisposition=:right,yticksvisible=true,yticklabelsvisible=true,
        xlabelsize=42,ylabelsize=42,xticklabelsize=30,yticklabelsize=30,
        xticks=(-log1p.(-[0.0,0.5,0.8,0.9,0.95,0.98,0.99]),[L"0",L"0.5",L"0.8",L"0.9",L"0.95",L"0.98",L"0.99"]),
        yticks=([0,.5,1],[L"0",L"0.5",L"1"]), # Avoid colliding labels at the shared panel boundary.
        xgridvisible=false,ygridvisible=false,
        xticklabelpad=12,yticklabelpad=12,xlabelpadding=8,ylabelpadding=8,
        xtickalign=1,ytickalign=0)
    xlims!(ax,0,-log(0.01)); ylims!(ax,-.52,1.02)
    save(joinpath(@__DIR__,"Figure6_orbit_right_axis_frame.png"),fig;px_per_unit=3)
end
main()
