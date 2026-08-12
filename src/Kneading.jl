"""
Symbolic-dynamics tools for one-dimensional maps.

Use [`Kneading.OneDimensionalMaps`](@ref) for interval-map kneading data and
entropy estimates. Use [`Kneading.Diagrams`](@ref) for parameter-plane scans
and kneading-diagram contours.
"""
module Kneading

export Diagrams, OneDimensionalMaps

include("one_dimensional_maps.jl")
include("diagrams.jl")

end
