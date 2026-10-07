using DynamicalSystemsBase
using Kneading.FlowKneading
using LinearAlgebra
using Test

@testset "Flow event capture and word encoding" begin
    rotation(u, p, t) = SVector(-p[1] * u[2], p[1] * u[1])
    function rotation!(du, u, p, t)
        du[1] = -p[1] * u[2]
        du[2] = p[1] * u[1]
        return nothing
    end
    for rule in (rotation, rotation!)
        system = CoupledODEs(rule, [1.0, 0.0], [1.0])
        seed = (u0 = [1.0, 0.0], Q0 = reshape([1.0, 0.0], :, 1))
        state_before = copy(current_state(system))
        time_before = current_time(system)
        for integration in (:adaptive, :rk4)
            problem = FlowKneadingProblem(system; initializer = seed,
                capture = LocalMaximum(1), word_length = 3, maximum_time = 30.0,
                integration, dt = 0.01,
            )
            result = flow_kneading(problem)
            @test result.complete
            @test result.status === :complete
            @test result.raw_word == Int8[1, 1, 1, 1]
            @test result.transition_word == Int8[1, 1, 1]
            @test (result.raw_code, result.raw_length) == (15, 4)
            @test (result.transition_code, result.transition_length) == (7, 3)
            @test length(result.events) == 4
            @test result.return_times ≈ fill(2pi, 3) atol = 1e-5
            for event in result.events
                @test event.value ≈ 1.0 atol = 2e-5
                @test norm(event.tangent) ≈ 1.0 atol = 1e-12
                @test abs(dot(event.tangent, rotation(event.state, [1.0], event.time))) < 1e-12
            end
            flipped = (u0 = seed.u0, Q0 = -seed.Q0)
            opposite = flow_kneading(FlowKneadingProblem(system; initializer = flipped,
                capture = problem.capture, word_length = 3, maximum_time = 30.0,
                integration, dt = 0.01,
            ))
            @test opposite.complete
            @test opposite.raw_word == -result.raw_word
            @test opposite.raw_code == 0
            @test opposite.transition_word == result.transition_word
            @test opposite.transition_code == result.transition_code
        end
        minima = flow_kneading(FlowKneadingProblem(system; initializer = seed,
            capture = LocalMinimum(1), word_length = 2, maximum_time = 20.0,
        ))
        @test minima.complete
        @test minima.raw_word == Int8[-1, -1, -1]
        @test all(event -> event.value ≈ -1.0 && event.rate > 0, minima.events)
        @test current_state(system) == state_before
        @test current_time(system) == time_before
        @test current_parameters(system) == [1.0]
        current_parameters(system)[1] = 2.0
        @test minima.metadata.parameters == [1.0]
    end
end

@testset "Initial events, filters, and incomplete words" begin
    rotation(u, p, t) = SVector(-u[2], u[1])
    system = CoupledODEs(rotation, [1.0, 0.0])
    seed = (u0 = [1.0, 0.0], Q0 = reshape([1.0, 0.0], :, 1))
    initial = flow_kneading(FlowKneadingProblem(system; initializer = seed,
        capture = LocalMaximum(1), include_initial_event = true,
        word_length = 0, maximum_time = 1.0,
    ))
    @test initial.complete
    @test only(initial.events).time == 0.0
    @test initial.raw_word == Int8[1]
    @test initial.transition_word == Int8[]
    @test initial.transition_code == 0
    transient = flow_kneading(FlowKneadingProblem(system; initializer = seed,
        capture = LocalMaximum(1), include_initial_event = true,
        transient_events = 1, word_length = 0, maximum_time = 7.0,
    ))
    @test transient.complete
    @test only(transient.events).time ≈ 2pi atol = 1e-8
    missing = flow_kneading(FlowKneadingProblem(system; initializer = seed,
        capture = LocalMaximum(1), word_length = 2, maximum_time = 7.0,
    ))
    @test !missing.complete
    @test missing.status === :maximum_time
    @test length(missing.raw_word) == 1
    @test missing.raw_code == missing.transition_code == -1
    rejected = flow_kneading(FlowKneadingProblem(system; initializer = seed,
        capture = LocalMaximum(1; accept = (u, p, t) -> u[1] < 0),
        word_length = 2, maximum_time = 7.0,
    ))
    @test !rejected.complete
    @test isempty(rejected.events)
    ambiguous = flow_kneading(FlowKneadingProblem(system; initializer = seed,
        capture = LocalMaximum(1), observable = CoordinateComponent(2),
        word_length = 2, maximum_time = 7.0,
    ))
    @test !ambiguous.complete
    @test ambiguous.status === :ambiguous_sign
    @test only(ambiguous.events).sign == 0
    @test isempty(ambiguous.raw_word)
    direction = flow_kneading(FlowKneadingProblem(system; initializer = seed,
        capture = LocalMaximum(1), observable = (u, p, t) -> [2.0, 0.0],
        word_length = 0, maximum_time = 7.0,
    ))
    @test direction.complete
    @test only(direction.events).component ≈ 2.0
    @test_throws ArgumentError FlowKneadingProblem(system; initializer = seed,
        capture = LocalMaximum(1), word_length = -1)
    @test_throws ArgumentError FlowKneadingProblem(system; initializer = seed,
        capture = LocalMaximum(3))
    @test_throws ArgumentError flow_kneading(FlowKneadingProblem(system;
        initializer = seed, capture = LocalMinimum(1), include_initial_event = true))
end

@testset "Real-saddle flow initializer and Lorenz positive maxima" begin
    linear(u, p, t) = SVector(u[1], -u[2], -2u[3])
    system = CoupledODEs(linear, zeros(3))
    initializer = RealSaddleInitializer(equilibrium_guess = zeros(3))
    result = flow_kneading(FlowKneadingProblem(system; initializer,
        capture = LocalMaximum(1), word_length = 1, maximum_time = 0.2))
    @test result.initialization.u0 ≈ [1e-6, 0.0, 0.0]
    @test vec(result.initialization.Q0) ≈ [0.0, 1.0, 0.0]
    @test !result.complete
    @test result.status === :maximum_time
    stable(u, p, t) = -u
    @test_throws Kneading.RealSaddleInitialization.InvalidRealSaddleInitialState flow_kneading(
        FlowKneadingProblem(CoupledODEs(stable, zeros(3)); initializer,
            capture = LocalMaximum(1), maximum_time = 0.2))
    lorenz_rule(u, p, t) = SVector(p[1] * (u[2] - u[1]),
        u[1] * (p[2] - u[3]) - u[2], u[1] * u[2] - p[3] * u[3])
    lorenz = CoupledODEs(lorenz_rule, zeros(3), [10.0, 28.0, 8 / 3])
    problem = FlowKneadingProblem(lorenz; initializer,
        capture = LocalMaximum(1; accept = (u, p, t) -> u[1] > 0),
        word_length = 6, maximum_time = 200.0,
    )
    positive = flow_kneading(problem)
    @test positive.complete
    @test !positive.metadata.tangent_seed_fallback
    @test length(positive.events) == 7
    @test all(event -> event.value > 0 && event.rate < 0, positive.events)
    @test all(event -> abs(event.state[2] - event.state[1]) < 1e-8, positive.events)
    @test all(positive.return_times .> 0)
    for event in positive.events
        @test norm(event.tangent) ≈ 1.0 atol = 1e-12
        @test abs(dot(event.tangent, lorenz_rule(event.state, positive.metadata.parameters, event.time))) < 1e-9
    end
end

@testset "Degenerate real-saddle tangent seeds" begin
    # Lorenz with a slow linear variable w that does not feed back. Its rate
    # makes e_w the leading stable direction, and span(e_w) is invariant under
    # the variational equation, so the leading-stable seed never reaches z.
    for drive in (0.0, 0.5)
        rule(u, p, t) = SVector(10.0 * (u[2] - u[1]), u[1] * (28.0 - u[3]) - u[2],
            u[1] * u[2] - 8 / 3 * u[3], -p[1] * u[4] + drive * u[1])
        system = CoupledODEs(rule, zeros(4), [0.1])
        solve(tangent_seed; kwargs...) = flow_kneading(FlowKneadingProblem(system;
            initializer = RealSaddleInitializer(; equilibrium_guess = zeros(4), tangent_seed),
            capture = LocalMaximum(3), word_length = 8, maximum_time = 200.0, kwargs...))

        unchecked = solve(:leading_stable; seed_invariance_tolerance = nothing)
        @test unchecked.initialization.leading_stable_eigenvalue ≈ -0.1
        @test unchecked.initialization.seed_direction ≈ [0.0, 0.0, 0.0, 1.0]
        @test abs(first(unchecked.events).component) < 1e-10

        rejected = solve(:leading_stable)
        @test !rejected.complete
        @test rejected.status === :degenerate_tangent_seed
        @test isempty(rejected.events)
        @test !isempty(rejected.metadata.detail)

        automatic = solve(:auto)
        @test automatic.complete
        @test automatic.metadata.tangent_seed_fallback
        @test automatic.initialization.tangent_seed == 3
        @test automatic.initialization.seed_direction == [0.0, 0.0, 1.0, 0.0]
        @test all(event -> abs(event.component) > 0.1, automatic.events)
        for tangent_seed in (1, [1.0, 2.0, -1.0, 0.0])
            other = solve(tangent_seed)
            @test other.complete
            @test !other.metadata.tangent_seed_fallback
            @test other.transition_word == automatic.transition_word
        end

        # A returned seed, as reused by parameter scans, falls back in the same way.
        seed = Kneading.RealSaddleInitialization.init_real_saddle(system, zeros(4), 1e-6,
            Kneading.RealSaddleInitialization.RealSaddleTolerances(
                1e-10, 1e-8, 1e-8, 1e-8, 1e-8, 1e-10, 1e-12, 1e-12);
            max_root_iterations = 30, equilibrium_branch = :selected,
            equilibrium_branch_check = u -> true,
            unstable_reference = ones(4), stable_reference = ones(4))
        @test seed.tangent_seed === :auto
        reused = flow_kneading(FlowKneadingProblem(system; initializer = seed,
            capture = LocalMaximum(3), word_length = 8, maximum_time = 200.0))
        @test reused.metadata.tangent_seed_fallback
        @test reused.transition_word == automatic.transition_word
    end
    @test_throws ArgumentError RealSaddleInitializer(equilibrium_guess = zeros(4), tangent_seed = :weakest)
    @test_throws ArgumentError FlowKneadingProblem(CoupledODEs((u, p, t) -> -u, zeros(3));
        initializer = RealSaddleInitializer(equilibrium_guess = zeros(3)),
        capture = LocalMaximum(1), seed_invariance_tolerance = -1.0)
end
