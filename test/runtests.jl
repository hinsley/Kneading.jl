import Kneading

using Kneading.Diagrams
using Kneading.OneDimensionalMaps
using Test

@testset "Public namespaces" begin
    @test :Diagrams in names(Kneading)
    @test :OneDimensionalMaps in names(Kneading)
    @test !(:ParameterPlane in names(Kneading))
    @test !(:PartitionedIntervalMap in names(Kneading))
end

function example_map(x)
    if x < 1 // 3
        return x / 2
    elseif x < 2 // 3
        return 3x - 1
    else
        return 3 - 3x
    end
end

partition = LapPartition((0 // 1, 1 // 1), [1 // 3, 2 // 3])
interval_map = PartitionedIntervalMap(example_map, partition)
data = kneading_data(interval_map, 10)
matrix = kneading_matrix(data)

@testset "Example workflow" begin
    @test interval_map.orientations == [Increasing, Increasing, Decreasing]
    @test size(data.itineraries) == (2, 2)
    @test data[1, LeftSide] == fill(1, 11)
    @test data[1, RightSide] == [2, fill(1, 10)...]
    @test data[2, LeftSide] == [2, 3, fill(1, 9)...]
    @test data[2, RightSide] == [3, 3, fill(1, 9)...]
    @test size(matrix) == (2, 3)
    @test all(length(entry) == 11 for entry in matrix)

    determinant = kneading_determinant(matrix)
    polynomial_determinant = kneading_determinant(matrix; mode = :polynomial)

    @test determinant.numerator.coefficients == BigInt[1, -2, zeros(BigInt, 9)...]
    @test polynomial_determinant.numerator.coefficients == BigInt[1, -2]
    @test denominator(determinant).coefficients == BigInt[1, -1]

    entropy = entropy_estimate(determinant)
    exact_entropy = log(BigFloat(2))
    @test entropy.entropy_interval[1] <= exact_entropy <= entropy.entropy_interval[2]
end

@testset "Interval-map regressions" begin
    tent_map(x) = x < 0.5 ? x : 1.0 - x
    tent_partition = LapPartition((0.0, 1.0), [0.5])
    tent_data = kneading_data(PartitionedIntervalMap(tent_map, tent_partition), 2)

    @test tent_data[1, LeftSide] == [1, 1, 1]
    @test tent_data[1, RightSide] == [2, 1, 1]

    function varied_map(x)
        x == 0.5 && return 2.0
        return x <= 0.5 ? sin(x) : x^2
    end

    varied_partition = LapPartition((0.0, 1.0), [0.5])
    varied_map_data = kneading_data(
        PartitionedIntervalMap(varied_map, varied_partition),
        1,
    )

    @test varied_map_data[1, LeftSide] == [1, 1]
    @test varied_map_data[1, RightSide] == [2, 1]
end

@testset "Parameter-plane scans" begin
    plane = ParameterPlane(
        collect(range(-1.0, 1.0; length = 7)),
        collect(range(-1.0, 1.0; length = 7));
        xname = "a",
        yname = "b",
    )
    fields = zeros(Float64, 7, 7, 2)
    scan_plane!(fields, plane) do destination, x, y
        destination[1] = x + y
        destination[2] = x - y
    end

    @test fields[:, :, 1] == [x + y for y in plane.y, x in plane.x]
    @test fields[:, :, 2] == [x - y for y in plane.y, x in plane.x]
    @test_throws ArgumentError ParameterPlane([0.0, 0.0], [0.0, 1.0])
    @test_throws DimensionMismatch level_contours(plane, zeros(2, 2))

    diagram = KneadingDiagram(plane; metadata = (source = :test,))
    add_level_contours!(
        diagram,
        view(fields, :, :, 1);
        source = :sum,
        iterate = 2,
        level = 0.0,
    )

    @test diagram.metadata.source == :test
    @test length(diagram.layers) == 1
    @test diagram.layers[1].source == :sum
    @test !isempty(diagram.layers[1].segments)
    for segment in diagram.layers[1].segments
        @test abs(segment.x1 + segment.y1) <= 32eps(Float64)
        @test abs(segment.x2 + segment.y2) <= 32eps(Float64)
    end
end

@testset "Continuation parameter-plane scans" begin
    plane = ParameterPlane([10.0, 20.0, 30.0], [1.0, 2.0])
    state_factory = (line_index, initial_x, initial_y) ->
        [Float64(line_index), 0.0, initial_x, initial_y]
    step! = function (output, state, x, y)
        state[2] += 1
        output[1] = state[1]
        output[2] = state[2]
        output[3] = state[3]
        output[4] = state[4]
        output[5] = x
        output[6] = y
    end

    column_forward = zeros(Float64, 2, 3, 6)
    scan_continuation!(step!, state_factory, column_forward, plane)
    @test column_forward[:, :, 1] == [
        1.0 2.0 3.0
        1.0 2.0 3.0
    ]
    @test column_forward[:, :, 2] == [
        1.0 1.0 1.0
        2.0 2.0 2.0
    ]
    @test column_forward[:, :, 3] == [
        10.0 20.0 30.0
        10.0 20.0 30.0
    ]
    @test column_forward[:, :, 4] == [
        1.0 1.0 1.0
        1.0 1.0 1.0
    ]
    @test column_forward[:, :, 5] == [
        10.0 20.0 30.0
        10.0 20.0 30.0
    ]
    @test column_forward[:, :, 6] == [
        1.0 1.0 1.0
        2.0 2.0 2.0
    ]

    column_backward = zeros(Float64, 2, 3, 6)
    scan_continuation!(
        step!,
        state_factory,
        column_backward,
        plane;
        along = :columns,
        direction = :backward,
    )
    @test column_backward[:, :, 1] == column_forward[:, :, 1]
    @test column_backward[:, :, 2] == [
        2.0 2.0 2.0
        1.0 1.0 1.0
    ]
    @test column_backward[:, :, 3] == column_forward[:, :, 3]
    @test column_backward[:, :, 4] == [
        2.0 2.0 2.0
        2.0 2.0 2.0
    ]
    @test column_backward[:, :, 5] == column_forward[:, :, 5]
    @test column_backward[:, :, 6] == column_forward[:, :, 6]

    row_forward = zeros(Float64, 2, 3, 6)
    scan_continuation!(
        step!,
        state_factory,
        row_forward,
        plane;
        along = :rows,
    )
    @test row_forward[:, :, 1] == [
        1.0 1.0 1.0
        2.0 2.0 2.0
    ]
    @test row_forward[:, :, 2] == [
        1.0 2.0 3.0
        1.0 2.0 3.0
    ]
    @test row_forward[:, :, 3] == fill(10.0, 2, 3)
    @test row_forward[:, :, 4] == [
        1.0 1.0 1.0
        2.0 2.0 2.0
    ]
    @test row_forward[:, :, 5] == column_forward[:, :, 5]
    @test row_forward[:, :, 6] == column_forward[:, :, 6]

    row_backward = zeros(Float64, 2, 3, 6)
    scan_continuation!(
        step!,
        state_factory,
        row_backward,
        plane;
        along = :rows,
        direction = :backward,
    )
    @test row_backward[:, :, 1] == row_forward[:, :, 1]
    @test row_backward[:, :, 2] == [
        3.0 2.0 1.0
        3.0 2.0 1.0
    ]
    @test row_backward[:, :, 3] == fill(30.0, 2, 3)
    @test row_backward[:, :, 4] == row_forward[:, :, 4]
    @test row_backward[:, :, 5] == column_forward[:, :, 5]
    @test row_backward[:, :, 6] == column_forward[:, :, 6]

    @test_throws ArgumentError scan_continuation!(
        step!,
        state_factory,
        zeros(Float64, 2, 3, 6),
        plane;
        along = :diagonals,
    )
    @test_throws ArgumentError scan_continuation!(
        step!,
        state_factory,
        zeros(Float64, 2, 3, 6),
        plane;
        direction = :sideways,
    )
end

@testset "Boolean symbolic boundaries" begin
    plane = ParameterPlane([0.0, 1.0, 2.0], [0.0, 1.0])
    field = Bool[
        false false true
        false false true
    ]
    segments = boolean_contours(plane, field)

    @test length(segments) == 1
    @test segments[1] == ContourSegment(1.5, 0.0, 1.5, 1.0)

    checkerboard_plane = ParameterPlane([0.0, 1.0], [0.0, 1.0])
    checkerboard = Bool[
        true false
        false true
    ]
    checkerboard_segments = boolean_contours(checkerboard_plane, checkerboard)

    @test checkerboard_segments == [
        ContourSegment(0.5, 0.0, 1.0, 0.5),
        ContourSegment(0.5, 1.0, 0.0, 0.5),
    ]
end

jet(values...) = PowerSeriesJet(BigInt[values...])

function denominator_jet(determinant::KneadingDeterminant, coefficient_count)
    coefficients = zeros(BigInt, coefficient_count)
    coefficients[1] = 1
    coefficient_count > 1 && (coefficients[2] = -Int(determinant.epsilon))
    return PowerSeriesJet(coefficients)
end

@testset "Determinant regressions" begin
    # Check both determinant algorithms and both finite-data modes.
    regression_entries = [
        jet(-1, 1, 2) jet(1, 2, -1) jet(0, 1, 0)
        jet(0, 3, 1) jet(-1, -1, 2) jet(1, -2, 1)
    ]
    regression_matrix = KneadingMatrix(
        regression_entries,
        [Increasing, Decreasing, Increasing],
    )

    jet_determinant = kneading_determinant(regression_matrix; deleted_lap = 3)
    berkowitz_determinant = kneading_determinant(
        regression_matrix;
        deleted_lap = 3,
        backend = :berkowitz,
    )
    polynomial_determinant = kneading_determinant(
        regression_matrix;
        deleted_lap = 3,
        mode = :polynomial,
    )

    @test jet_determinant.numerator.coefficients == BigInt[1, -3, -12]
    @test berkowitz_determinant.numerator.coefficients ==
        jet_determinant.numerator.coefficients
    @test polynomial_determinant.numerator.coefficients == BigInt[1, -3, -12, 1, 5]

    # The determinant is independent of the deleted lap.
    deleted_determinants = [
        kneading_determinant(matrix; deleted_lap = lap)
        for lap in axes(matrix, 2)
    ]
    reference = first(deleted_determinants)
    reference_denominator = denominator_jet(reference, 11)

    for determinant in deleted_determinants[2:end]
        determinant_denominator = denominator_jet(determinant, 11)
        @test (reference.numerator * determinant_denominator).coefficients ==
            (determinant.numerator * reference_denominator).coefficients
    end

    # Determinant arithmetic must not overflow Int.
    large = typemax(Int)
    overflow_entries = [
        PowerSeriesJet([
            column == row ? -1 : column == row + 1 ? 1 : 0,
            large,
        ])
        for row in 1:2, column in 1:3
    ]
    overflow_matrix = KneadingMatrix(overflow_entries, fill(Increasing, 3))
    @test kneading_determinant(overflow_matrix).numerator.coefficients ==
        BigInt[1, 3 * BigInt(large)]
end

root_factor(integer) = Polynomial(BigInt[1, -integer])

@testset "Entropy regressions" begin
    near_one = Polynomial(Rational{BigInt}[
        1,
        -(BigInt(1001) // BigInt(1000)),
    ])
    cases = [
        (root_factor(3) * root_factor(2), log(BigFloat(3))),
        (root_factor(2)^2, log(BigFloat(2))),
        (root_factor(1)^2 * root_factor(2), log(BigFloat(2))),
        (root_factor(1000) * root_factor(1001), log(BigFloat(1001))),
        (near_one, log(BigFloat(1001) / BigFloat(1000))),
    ]

    for (polynomial, expected_entropy) in cases
        determinant = KneadingDeterminant(polynomial, Increasing, 1)
        entropy = entropy_estimate(determinant; bits = 96)
        @test entropy.entropy_interval[1] <=
            expected_entropy <=
            entropy.entropy_interval[2]
    end

    determinant = KneadingDeterminant(root_factor(2), Increasing, 1)
    coarse = entropy_estimate(determinant; bits = 64)
    fine = entropy_estimate(determinant; bits = 128)
    @test fine.entropy_interval[2] - fine.entropy_interval[1] <=
        coarse.entropy_interval[2] - coarse.entropy_interval[1]

    no_root = KneadingDeterminant(Polynomial(BigInt[1, 1]), Increasing, 1)
    endpoint_only = KneadingDeterminant(root_factor(1), Increasing, 1)
    @test entropy_estimate(no_root) === nothing
    @test entropy_estimate(endpoint_only) === nothing
end
