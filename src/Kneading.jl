"""
Symbolic-dynamics tools for maps and flow-derived return data.

Use [`Kneading.OneDimensionalMaps`](@ref) for interval-map kneading data and
entropy estimates. Use [`Kneading.Diagrams`](@ref) for parameter-plane scans
and kneading-diagram contours. Use [`Kneading.FlowKneading`](@ref) to initialize
critical orbits from real saddles or saddle-foci, capture extremum events,
and scan their orientation words across parameter planes.
"""
module Kneading

export Diagrams,
    FlowKneading,
    FlowNormalTangents,
    OneDimensionalMaps,
    RealSaddleInitialization

include("one_dimensional_maps.jl")
include("diagrams.jl")
include("flow_normal_tangents.jl")
include("real_saddle_initialization.jl")
include("flow_kneading.jl")

using .FlowKneading
using .FlowNormalTangents

export LocalMaximum, LocalMinimum, CoordinateComponent,
    RealSaddleInitializer, SaddleFocusInitializer, SaddleFocusSeed,
    SaddleFocusInitializationError, init_saddle_focus,
    FlowKneadingProblem, FlowKneadingEvent, FlowKneadingResult, flow_kneading,
    FlowKneadingDiagram, scan_flow_kneading, write_flow_scan,
    FlowNormalTolerances, init_flow_normal, solve_flow_normal!

end
