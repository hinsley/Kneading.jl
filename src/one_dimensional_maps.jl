module OneDimensionalMaps

export Decreasing,
    EntropyEstimate,
    Increasing,
    KneadingData,
    KneadingDeterminant,
    KneadingMatrix,
    LapOrientation,
    LapPartition,
    LeftSide,
    PartitionedIntervalMap,
    Polynomial,
    PowerSeriesJet,
    RootInterval,
    RightSide,
    TangentDirection,
    entropy_estimate,
    kneading_data,
    kneading_determinant,
    kneading_matrix,
    polynomial_approximation,
    real_root_intervals

include("interval_maps.jl")
include("kneading_algebra.jl")
include("entropy.jl")

end
