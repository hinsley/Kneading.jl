using DynamicalSystemsBase
using Kneading.FlowNormalTangents
using LinearAlgebra
using SciMLBase
using Test

const FLOW_NORMAL_TOLERANCES = FlowNormalTolerances(1.0e-12, 1.0e-12, 1.0e-12)

function constant_flow!(du, u, p, t)
    du[1] = 1
    du[2] = 0
    return nothing
end

function zero_jacobian!(J, u, p, t)
    fill!(J, 0)
    return nothing
end

constant_flow(u, p, t) = SVector(1.0, 0.0)
zero_jacobian(u, p, t) = SMatrix{2, 2}(0.0, 0.0, 0.0, 0.0)

@testset "Flow-normal tolerances" begin
    promoted = FlowNormalTolerances(1, 2.0, 3.0f0)

    @test promoted.flow_norm === 1.0
    @test promoted.projected_norm === 2.0
    @test promoted.unit_norm === 3.0
    @test_throws ArgumentError FlowNormalTolerances(-1.0, 1.0, 1.0)
    @test_throws ArgumentError FlowNormalTolerances(1.0, Inf, 1.0)
end

@testset "Flow-normal projection" begin
    tangent = [2.0, 3.0]
    returned_tangent = project_flow_normal!(
        tangent,
        [1.0, 0.0],
        FLOW_NORMAL_TOLERANCES,
    )

    @test returned_tangent === tangent
    @test tangent == [0.0, 1.0]
    @test dot(tangent, [1.0, 0.0]) == 0
    @test norm(tangent) == 1
    @test_throws DimensionMismatch project_flow_normal!(
        [1.0],
        [1.0, 0.0],
        FLOW_NORMAL_TOLERANCES,
    )
    @test_throws DomainError project_flow_normal!(
        [1.0, 1.0],
        [0.0, 0.0],
        FLOW_NORMAL_TOLERANCES,
    )
    @test_throws DomainError project_flow_normal!(
        [1.0, 0.0],
        [1.0, 0.0],
        FLOW_NORMAL_TOLERANCES,
    )
end

@testset "Out-of-place tangent evolution" begin
    base_system = CoupledODEs(constant_flow, zeros(2))
    tangent_system = init_flow_normal(
        base_system,
        reshape([1.0, 1.0], 2, 1),
        FLOW_NORMAL_TOLERANCES;
        J = zero_jacobian,
    )

    solve_flow_normal!(tangent_system, 0.1)
    tangent = vec(current_deviations(tangent_system))

    @test norm(tangent) ≈ 1
    @test dot(tangent, [1.0, 0.0]) ≈ 0 atol = 1.0e-12
end

@testset "Flow-normal tangent evolution" begin
    base_system = CoupledODEs(
        constant_flow!,
        zeros(2);
        diffeq = (abstol = 1.0e-10, reltol = 1.0e-10),
    )
    Q0 = reshape([1.0, 1.0], 2, 1)
    callback_times = Float64[]
    user_callback = DiscreteCallback(
        (u, t, integrator) -> true,
        integrator -> push!(callback_times, integrator.t);
        save_positions = (false, false),
    )
    tangent_system = init_flow_normal(
        base_system,
        Q0,
        FLOW_NORMAL_TOLERANCES;
        u0 = [2.0, 3.0],
        J = zero_jacobian!,
        J0 = zeros(2, 2),
        callback = user_callback,
    )

    @test tangent_system isa TangentDynamicalSystem
    @test current_state(tangent_system) == [2.0, 3.0]
    @test current_deviations(tangent_system) == Q0

    solution = solve_flow_normal!(tangent_system, 0.2)
    tangent = vec(current_deviations(tangent_system))

    @test solution isa SciMLBase.AbstractODESolution
    @test current_time(tangent_system) == 0.2
    @test norm(tangent) ≈ 1
    @test dot(tangent, [1.0, 0.0]) ≈ 0 atol = 1.0e-12
    @test !isempty(callback_times)
    @test callback_times[end] == 0.2
    @test_throws ArgumentError solve_flow_normal!(tangent_system, 0.1)
    @test_throws ArgumentError init_flow_normal(
        base_system,
        ones(2, 2),
        FLOW_NORMAL_TOLERANCES;
        J = zero_jacobian!,
        J0 = zeros(2, 2),
    )
end
