using Test

using Kneading.Diagrams

include("chebyshev_cubic_scan.jl")

using .ChebyshevCubicScan

@testset "Chebyshev cubic scan" begin
    u = 0.37
    v = -0.81
    left, right = chebyshev_critical_points()
    @test chebyshev_cubic(u, v, left) ≈ u
    @test chebyshev_cubic(u, v, right) ≈ v

    diagram = compute_chebyshev_cubic_diagram(
        grid_size = 31,
        iterates = 6,
    )
    @test diagram.metadata.family == :chebyshev_cubic
    @test length(diagram.layers) == 2 * 2 * 5
    @test any(!isempty(layer.segments) for layer in diagram.layers)
    @test all(
        segment -> all(
            isfinite,
            (segment.x1, segment.y1, segment.x2, segment.y2),
        ),
        Iterators.flatten(layer.segments for layer in diagram.layers),
    )
end
