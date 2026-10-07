using DynamicalSystemsBase
using Kneading.FlowNormalTangents
using Kneading.RealSaddleInitialization
using LinearAlgebra
using Test

const REAL_SADDLE_TOLERANCES = RealSaddleTolerances(
    1.0e-12,
    1.0e-8,
    1.0e-8,
    1.0e-8,
    1.0e-8,
    1.0e-8,
    1.0e-12,
    1.0e-12,
)

linear_flow(u, p, t) = typeof(u)(p.matrix * (u - p.equilibrium))

function linear_flow!(du, u, p, t)
    mul!(du, p.matrix, u - p.equilibrium)
    return nothing
end

function capture_error(f)
    try
        f()
    catch error
        return error
    end
    return nothing
end

function real_saddle_seed(system; stable_reference = [0.0, -1.0, 0.0])
    return init_real_saddle(
        system,
        [1.2, -1.7, 0.8],
        0.1,
        REAL_SADDLE_TOLERANCES;
        max_root_iterations = 5,
        equilibrium_branch = :primary,
        equilibrium_branch_check = equilibrium ->
            norm(equilibrium - [1.0, -2.0, 0.5]) <= 1.0e-10,
        unstable_reference = [1.0, 0.0, 0.0],
        stable_reference,
        unstable_branch = -1,
    )
end

@testset "Real-saddle tolerances" begin
    promoted = RealSaddleTolerances(1, 2.0, 3.0f0, 4, 5, 6, 7, 8)

    @test promoted.equilibrium_residual === 1.0
    @test promoted.projected_norm === 8.0
    @test_throws ArgumentError RealSaddleTolerances(
        -1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
    )
end

@testset "Real-saddle seed" begin
    parameters = (
        matrix = Matrix(Diagonal([2.0, -1.0, -3.0])),
        equilibrium = [1.0, -2.0, 0.5],
    )
    out_of_place_system = CoupledODEs(linear_flow, zeros(3), parameters)
    in_place_system = CoupledODEs(linear_flow!, zeros(3), parameters)

    for system in (out_of_place_system, in_place_system)
        seed = real_saddle_seed(system)

        @test seed isa RealSaddleSeed
        @test seed.initial_state == [1.2, -1.7, 0.8]
        @test seed.equilibrium ≈ parameters.equilibrium atol = 1.0e-12
        @test seed.jacobian ≈ parameters.matrix atol = 1.0e-12
        @test seed.unstable_eigenvalue ≈ 2.0
        @test seed.leading_stable_eigenvalue ≈ -1.0
        @test seed.unstable_direction ≈ [-1.0, 0.0, 0.0]
        @test seed.leading_stable_direction ≈ [0.0, -1.0, 0.0]
        @test seed.u0 ≈ [0.9, -2.0, 0.5]
        @test vec(seed.Q0) ≈ [0.0, -1.0, 0.0]
        @test norm(seed.Q0) ≈ 1.0
        @test seed.root_iterations == 1
        @test seed.equilibrium_branch == :primary
        @test seed.unstable_branch == -1
        @test seed.equilibrium_residual <= REAL_SADDLE_TOLERANCES.equilibrium_residual
        @test seed.launch_flow_norm ≈ 0.2
        @test seed.projected_tangent_norm ≈ 1.0
        @test length(seed.eigenvalues) == 3
    end
end

@testset "Real-saddle tangent seed options" begin
    parameters = (
        matrix = Matrix(Diagonal([2.0, -1.0, -3.0])),
        equilibrium = [1.0, -2.0, 0.5],
    )
    system = CoupledODEs(linear_flow, zeros(3), parameters)
    seed_with(tangent_seed) = init_real_saddle(
        system,
        [1.2, -1.7, 0.8],
        0.1,
        REAL_SADDLE_TOLERANCES;
        max_root_iterations = 5,
        equilibrium_branch = :primary,
        equilibrium_branch_check = equilibrium -> true,
        unstable_reference = [1.0, 0.0, 0.0],
        stable_reference = [0.0, -1.0, 0.0],
        tangent_seed,
    )

    default_seed = real_saddle_seed(system)
    @test default_seed.tangent_seed === :auto
    @test default_seed.seed_direction ≈ default_seed.leading_stable_direction

    leading = seed_with(:leading_stable)
    @test leading.tangent_seed === :leading_stable
    @test vec(leading.Q0) ≈ vec(default_seed.Q0)

    coordinate = seed_with(3)
    @test coordinate.tangent_seed == 3
    @test coordinate.seed_direction == [0.0, 0.0, 1.0]
    @test vec(coordinate.Q0) ≈ [0.0, 0.0, 1.0]
    @test coordinate.leading_stable_direction ≈ [0.0, -1.0, 0.0]

    # The launch flow is along x, so the x part of a user vector is removed.
    vector = seed_with([1.0, 0.0, 2.0])
    @test vector.tangent_seed == [1.0, 0.0, 2.0]
    @test vector.seed_direction ≈ [1.0, 0.0, 2.0] / sqrt(5.0)
    @test vec(vector.Q0) ≈ [0.0, 0.0, 1.0]
    @test vector.projected_tangent_norm ≈ 2 / sqrt(5.0)

    @test_throws ArgumentError seed_with(:weakest)
    @test_throws ArgumentError seed_with(0)
    @test_throws ArgumentError seed_with(4)
    @test_throws DimensionMismatch seed_with([1.0, 0.0])
    @test_throws ArgumentError seed_with(zeros(3))
    @test_throws ArgumentError seed_with([NaN, 0.0, 1.0])
    # A seed along the launch flow has no flow-normal part.
    @test_throws DomainError seed_with(1)
end

@testset "Real-saddle seed initializes tangent integration" begin
    parameters = (
        matrix = Matrix(Diagonal([2.0, -1.0, -3.0])),
        equilibrium = [1.0, -2.0, 0.5],
    )
    system = CoupledODEs(linear_flow, zeros(3), parameters)
    seed = real_saddle_seed(system)
    tangent_system = init_flow_normal(
        system,
        seed.Q0,
        FlowNormalTolerances(1.0e-12, 1.0e-12, 1.0e-12);
        u0 = seed.u0,
    )

    solve_flow_normal!(tangent_system, 0.1)
    tangent = vec(current_deviations(tangent_system))
    flow = linear_flow(current_state(tangent_system), parameters, 0.1)

    @test norm(tangent) ≈ 1.0
    @test dot(tangent, flow) ≈ 0.0 atol = 1.0e-12
end

@testset "Real-saddle failures" begin
    equilibrium = zeros(3)
    wrong_branch_system = CoupledODEs(
        linear_flow,
        zeros(3),
        (matrix = Matrix(Diagonal([2.0, -1.0, -3.0])), equilibrium),
    )
    wrong_branch_error = capture_error() do
        init_real_saddle(
            wrong_branch_system,
            zeros(3),
            0.1,
            REAL_SADDLE_TOLERANCES;
            max_root_iterations = 2,
            equilibrium_branch = :other,
            equilibrium_branch_check = equilibrium -> false,
            unstable_reference = [1.0, 0.0, 0.0],
            stable_reference = [0.0, 1.0, 0.0],
        )
    end
    @test wrong_branch_error isa InvalidRealSaddleInitialState
    @test wrong_branch_error.stage == :equilibrium_branch

    two_unstable_system = CoupledODEs(
        linear_flow,
        zeros(3),
        (matrix = Matrix(Diagonal([2.0, 1.0, -3.0])), equilibrium),
    )
    two_unstable_error = capture_error() do
        init_real_saddle(
            two_unstable_system,
            zeros(3),
            0.1,
            REAL_SADDLE_TOLERANCES;
            max_root_iterations = 2,
            equilibrium_branch = :primary,
            equilibrium_branch_check = equilibrium -> true,
            unstable_reference = [1.0, 0.0, 0.0],
            stable_reference = [0.0, 0.0, 1.0],
        )
    end
    @test two_unstable_error isa InvalidRealSaddleInitialState
    @test two_unstable_error.stage == :real_saddle_spectrum

    complex_stable_matrix = [
        2.0 0.0  0.0
        0.0 -1.0 -1.0
        0.0 1.0  -1.0
    ]
    complex_stable_system = CoupledODEs(
        linear_flow,
        zeros(3),
        (matrix = complex_stable_matrix, equilibrium),
    )
    complex_stable_error = capture_error() do
        init_real_saddle(
            complex_stable_system,
            zeros(3),
            0.1,
            REAL_SADDLE_TOLERANCES;
            max_root_iterations = 2,
            equilibrium_branch = :primary,
            equilibrium_branch_check = equilibrium -> true,
            unstable_reference = [1.0, 0.0, 0.0],
            stable_reference = [0.0, 1.0, 0.0],
        )
    end
    @test complex_stable_error isa InvalidRealSaddleInitialState
    @test complex_stable_error.stage == :leading_stable_direction

    valid_system = CoupledODEs(
        linear_flow,
        zeros(3),
        (
            matrix = Matrix(Diagonal([2.0, -1.0, -3.0])),
            equilibrium = [1.0, -2.0, 0.5],
        ),
    )
    ambiguous_orientation_error = capture_error() do
        real_saddle_seed(valid_system; stable_reference = [0.0, 0.0, 1.0])
    end
    @test ambiguous_orientation_error isa DomainError

    constant_system = CoupledODEs(
        (u, p, t) -> typeof(u)(ones(2)),
        zeros(2),
    )
    root_error = capture_error() do
        init_real_saddle(
            constant_system,
            zeros(2),
            0.1,
            REAL_SADDLE_TOLERANCES;
            max_root_iterations = 2,
            equilibrium_branch = :none,
            equilibrium_branch_check = equilibrium -> true,
            unstable_reference = [1.0, 0.0],
            stable_reference = [0.0, 1.0],
        )
    end
    @test root_error isa InvalidRealSaddleInitialState
    @test root_error.stage == :root_solve

    @test_throws ArgumentError init_real_saddle(
        wrong_branch_system,
        zeros(3),
        0.0,
        REAL_SADDLE_TOLERANCES;
        max_root_iterations = 2,
        equilibrium_branch = :primary,
        equilibrium_branch_check = equilibrium -> true,
        unstable_reference = [1.0, 0.0, 0.0],
        stable_reference = [0.0, 1.0, 0.0],
    )
end
