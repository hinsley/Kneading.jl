using DynamicalSystemsBase
using KernelAbstractions
using Kneading
using Kneading.Diagrams: ParameterPlane
using LinearAlgebra
using Test

function same_results(cpu, device)
    @test cpu.status === device.status
    @test cpu.raw_word == device.raw_word
    @test cpu.transition_word == device.transition_word
    @test cpu.raw_code == device.raw_code
    @test cpu.transition_code == device.transition_code
    @test cpu.metadata.accepted_events == device.metadata.accepted_events
    @test cpu.metadata.terminal_time == device.metadata.terminal_time
    @test cpu.metadata.tangent_seed_fallback == device.metadata.tangent_seed_fallback
    @test length(cpu.events) == length(device.events)
    for (a, b) in zip(cpu.events, device.events)
        @test a.time == b.time
        @test a.state == b.state
        @test a.tangent ≈ b.tangent atol = 1e-12
        @test a.component ≈ b.component atol = 1e-12
        @test a.sign == b.sign
        @test sign(a.rate) == sign(b.rate)
    end
    @test cpu.return_times == device.return_times
end

@testset "Batched Lorenz words match fixed-step RK4" begin
    lorenz(u, p, t) = SVector(p[1] * (u[2] - u[1]), u[1] * (p[2] - u[3]) - u[2], u[1] * u[2] - p[3] * u[3])
    initializer = RealSaddleInitializer(equilibrium_guess = zeros(3))
    for (capture, observable) in ((LocalMaximum(3), CoordinateComponent(3)),
        (LocalMaximum(1; accept = (u, p, t) -> u[1] > 0), CoordinateComponent(1)))
        problems = [FlowKneadingProblem(CoupledODEs(lorenz, zeros(3), [10.0, rho, 8 / 3]);
            initializer, capture, observable, word_length = 16, maximum_time = 300.0,
            integration = :rk4, dt = 0.01) for rho in range(24.5, 60.0; length = 8)]
        cpu = flow_kneading.(problems)
        device = flow_kneading(problems; backend = CPU())
        @test all(result -> result.complete, cpu)
        foreach(same_results, cpu, device)
        @test [result.raw_word for result in flow_kneading(problems)] == [result.raw_word for result in cpu]
        mixed = flow_kneading(problems; backend = CPU(), precision = Float32)
        @test all(result -> result.complete, mixed)
        @test count(k -> mixed[k].raw_word[1:8] == cpu[k].raw_word[1:8], eachindex(cpu)) >= 7
        @test all(result -> result.metadata.precision === Float32, mixed)
    end
end

@testset "Batched Rössler scans match fixed-step RK4" begin
    rossler(u, p, t) = SVector(-u[2] - u[3], u[1] + p[1] * u[2], 0.3 * u[1] + u[3] * (u[1] - p[2]))
    plane = ParameterPlane([5.0, 5.5, 6.0], [0.30, 0.32]; xname = "c", yname = "a")
    settings = (; initializer = SaddleFocusInitializer(equilibrium_guess = zeros(3), critical_kind = :minimum),
        capture = LocalMinimum(2), word_length = 12, maximum_time = 2000.0, integration = :rk4, dt = 0.02)
    builder = (c, a) -> CoupledODEs(rossler, zeros(3), [a, c])
    cpu = scan_flow_kneading(builder, plane; settings..., store_results = true)
    device = scan_flow_kneading(builder, plane; settings..., backend = CPU(), store_results = true)
    @test all(==(:complete), cpu.statuses)
    @test device.statuses == cpu.statuses
    @test device.raw_words == cpu.raw_words
    @test device.raw_codes == cpu.raw_codes
    @test device.transition_codes == cpu.transition_codes
    @test device.critical_states == cpu.critical_states
    foreach(same_results, cpu.results, device.results)
    @test all(result -> result.metadata.initial_event_included, device.results)
    words = scan_flow_kneading(builder, plane; settings..., backend = CPU())
    @test words.raw_words == cpu.raw_words
    @test all(isnothing, words.results)
    mktempdir() do directory
        @test readlines(write_flow_scan(joinpath(directory, "cpu.tsv"), cpu)) ==
            readlines(write_flow_scan(joinpath(directory, "device.tsv"), words))
    end
end

@testset "Batched statuses and incomplete words" begin
    lorenz(u, p, t) = SVector(p[1] * (u[2] - u[1]), u[1] * (p[2] - u[3]) - u[2], u[1] * u[2] - p[3] * u[3])
    system = CoupledODEs(lorenz, zeros(3), [10.0, 28.0, 8 / 3])
    initializer = RealSaddleInitializer(equilibrium_guess = zeros(3))
    for options in ((; maximum_time = 5.0), (; max_state = 30.0), (; transient_events = 3),
        (; maxiters = 300), (; word_length = 70, maximum_time = 300.0),
        (; word_length = 127, maximum_time = 600.0))
        problem = FlowKneadingProblem(system; initializer, capture = LocalMaximum(3),
            integration = :rk4, dt = 0.01, options...)
        cpu = flow_kneading(problem)
        same_results(cpu, only(flow_kneading([problem]; backend = CPU())))
    end
    statuses = [only(flow_kneading([FlowKneadingProblem(system; initializer, capture = LocalMaximum(3),
        integration = :rk4, dt = 0.01, options...)]; backend = CPU())).status
        for options in ((; maximum_time = 5.0), (; max_state = 30.0), (; maxiters = 300))]
    @test statuses == [:maximum_time, :state_limit, :maximum_iterations]
    short = only(flow_kneading([FlowKneadingProblem(system; initializer, capture = LocalMaximum(3),
        integration = :rk4, dt = 0.01, maximum_time = 5.0)]; backend = CPU()))
    @test !short.complete
    @test short.raw_code == short.transition_code == -1
    @test 0 < short.raw_length == length(short.raw_word) < 17

    rotation(u, p, t) = SVector(-u[2], u[1])
    planar = CoupledODEs(rotation, [1.0, 0.0])
    seed = (u0 = [1.0, 0.0], Q0 = reshape([1.0, 0.0], :, 1))
    for options in ((; observable = CoordinateComponent(2)),
        (; observable = (u, p, t) -> SVector(2.0, 0.0), word_length = 0),
        (; observable = (u, p, t) -> SVector(NaN, 0.0)),
        (; include_initial_event = true, transient_events = 1, word_length = 0),
        (; capture = LocalMaximum(1; accept = (u, p, t) -> u[1] < 0)))
        problem = FlowKneadingProblem(planar; initializer = seed, capture = LocalMaximum(1),
            word_length = 2, maximum_time = 7.0, integration = :rk4, dt = 0.01, options...)
        cpu = flow_kneading(problem)
        device = only(flow_kneading([problem]; backend = CPU()))
        same_results(cpu, device)
    end
    ambiguous = only(flow_kneading([FlowKneadingProblem(planar; initializer = seed,
        capture = LocalMaximum(1), observable = CoordinateComponent(2), word_length = 2,
        maximum_time = 7.0, integration = :rk4, dt = 0.01)]; backend = CPU()))
    @test ambiguous.status === :ambiguous_sign
    @test only(ambiguous.events).sign == 0
    failed = only(flow_kneading([FlowKneadingProblem(planar; initializer = seed,
        capture = LocalMaximum(1), observable = (u, p, t) -> SVector(NaN, 0.0),
        maximum_time = 7.0, integration = :rk4, dt = 0.01)]; backend = CPU()))
    @test failed.status === :numerical_failure
    @test !isempty(failed.metadata.detail)
end

@testset "Batched degenerate real-saddle seeds" begin
    rule(u, p, t) = SVector(10.0 * (u[2] - u[1]), u[1] * (28.0 - u[3]) - u[2],
        u[1] * u[2] - 8 / 3 * u[3], -p[1] * u[4] + 0.5 * u[1])
    system = CoupledODEs(rule, zeros(4), [0.1])
    problems = [FlowKneadingProblem(system;
        initializer = RealSaddleInitializer(; equilibrium_guess = zeros(4), tangent_seed),
        capture = LocalMaximum(3), word_length = 8, maximum_time = 200.0, integration = :rk4, dt = 0.01)
        for tangent_seed in (:auto, :leading_stable, 1)]
    cpu = flow_kneading.(problems)
    device = flow_kneading(problems; backend = CPU())
    foreach(same_results, cpu, device)
    @test [result.status for result in device] == [:complete, :degenerate_tangent_seed, :complete]
    @test device[1].metadata.tangent_seed_fallback
    @test device[1].initialization.tangent_seed == 3
end

@testset "Batched backend requirements" begin
    lorenz(u, p, t) = SVector(p[1] * (u[2] - u[1]), u[1] * (p[2] - u[3]) - u[2], u[1] * u[2] - p[3] * u[3])
    function lorenz!(du, u, p, t)
        du .= lorenz(u, p, t)
        return nothing
    end
    initializer = RealSaddleInitializer(equilibrium_guess = zeros(3))
    problem(rule; kwargs...) = FlowKneadingProblem(CoupledODEs(rule, zeros(3), [10.0, 28.0, 8 / 3]);
        initializer, capture = LocalMaximum(3), word_length = 4, kwargs...)
    @test_throws ArgumentError flow_kneading([problem(lorenz)]; backend = CPU())
    @test_throws ArgumentError flow_kneading([problem(lorenz!; integration = :rk4)]; backend = CPU())
    @test_throws ArgumentError flow_kneading([problem(lorenz; integration = :rk4)]; precision = Float32)
    @test_throws ArgumentError flow_kneading([problem(lorenz; integration = :rk4),
        problem(lorenz; integration = :rk4, dt = 0.01)]; backend = CPU())
    @test_throws ArgumentError flow_kneading([problem(lorenz; integration = :rk4)]; backend = :gpu)
    @test isempty(flow_kneading(FlowKneadingProblem[]; backend = CPU()))
    @test only(flow_kneading([problem(lorenz)])).complete
end
