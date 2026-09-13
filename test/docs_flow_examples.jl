using DynamicalSystemsBase
using Kneading
using LinearAlgebra
using Test

@testset "End-user Lorenz flow-kneading examples" begin
    function lorenz_example!(du, u, p, t)
        x, y, z = u
        sigma, rho, beta = p
        du[1] = sigma * (y - x)
        du[2] = x * (rho - z) - y
        du[3] = x * y - beta * z
        return nothing
    end

    parameters = [10.0, 28.0, 8 / 3]
    lorenz = CoupledODEs(lorenz_example!, zeros(3), parameters)
    initializer = RealSaddleInitializer(
        equilibrium_guess=zeros(3),
        launch_distance=1e-6,
        unstable_branch=1,
    )
    maxima_z = flow_kneading(FlowKneadingProblem(
        lorenz;
        initializer,
        capture=LocalMaximum(3),
        observable=CoordinateComponent(3),
        word_length=16,
        maximum_time=1000.0,
    ))

    @test maxima_z.complete
    @test maxima_z.status == :complete
    @test length(maxima_z.raw_word) == 17
    @test length(maxima_z.transition_word) == 16
    @test maxima_z.raw_code >= 0
    @test length(maxima_z.events) == 17
    @test length(maxima_z.return_times) == 16
    @test all(>(0), maxima_z.return_times)
    @test !maxima_z.metadata.initial_event_included
    @test norm(maxima_z.initialization.equilibrium) < 1e-10
    @test norm(maxima_z.initialization.u0) ≈ 1e-6
    @test all(event -> event.rate < 0, maxima_z.events)
    @test all(event -> abs(event.state[1] * event.state[2] - parameters[3] * event.state[3]) < 1e-6, maxima_z.events)

    maxima_positive_x = flow_kneading(FlowKneadingProblem(
        lorenz;
        initializer,
        capture=LocalMaximum(1; accept=(u, p, t) -> u[1] > 0),
        observable=CoordinateComponent(1),
        word_length=16,
        maximum_time=1000.0,
    ))

    @test maxima_positive_x.complete
    @test length(maxima_positive_x.raw_word) == 17
    @test length(maxima_positive_x.transition_word) == 16
    @test all(event -> event.state[1] > 0, maxima_positive_x.events)
    @test all(event -> abs(event.state[1] - event.state[2]) < 1e-6, maxima_positive_x.events)
    @test all(event -> event.rate < 0, maxima_positive_x.events)
    @test maxima_positive_x.return_times ≈ diff([event.time for event in maxima_positive_x.events])

    for result in (maxima_z, maxima_positive_x), event in result.events
        flow = zeros(3)
        lorenz_example!(flow, event.state, parameters, event.time)
        @test norm(event.tangent) ≈ 1.0 atol=1e-12
        @test abs(dot(event.tangent, flow)) <= 1e-10 * norm(flow)
    end
    @test DynamicalSystemsBase.current_state(lorenz) == zeros(3)
    @test DynamicalSystemsBase.current_parameters(lorenz) == parameters

    short = flow_kneading(FlowKneadingProblem(
        lorenz;
        initializer,
        capture=LocalMaximum(3),
        word_length=16,
        maximum_time=0.01,
    ))
    @test !short.complete
    @test short.status == :maximum_time
    @test short.raw_code == -1
    @test short.transition_code == -1
    @test isempty(short.events)
end

@testset "End-user refined Rössler initialization and continuation" begin
    function rossler_example!(du, u, p, t)
        x, y, z = u
        a, c = p
        du[1] = -y - z
        du[2] = x + a * y
        du[3] = 0.3 * x + z * (x - c)
        return nothing
    end

    rossler = CoupledODEs(rossler_example!, zeros(3), [0.3, 5.5])
    capture = LocalMinimum(2)
    seed = init_saddle_focus(
        rossler;
        equilibrium_guess=zeros(3),
        capture,
        critical_kind=:minimum,
    )

    @test seed.diagnostics.converged
    @test seed.diagnostics.refinement_checked
    @test seed.diagnostics.state_error <= 1e-6
    @test seed.diagnostics.tangent_error <= 1e-6
    @test seed.diagnostics.refinement_steps >= 1
    @test seed.event_index >= 4
    @test seed.diagnostics.refinement_event_index > seed.event_index
    @test norm(seed.Q0) ≈ 1.0 atol=1e-12
    @test seed.equilibrium == zeros(3)
    @test seed.parameters == [0.3, 5.5]
    @test norm(seed.seed_direction) ≈ 1.0
    @test norm(seed.equilibrium + exp(seed.rho) * seed.seed_direction) ≈ exp(seed.rho)

    result = flow_kneading(FlowKneadingProblem(
        rossler;
        initializer=seed,
        capture,
        observable=CoordinateComponent(2),
        word_length=7,
        maximum_time=2000.0,
    ))

    @test result.complete
    @test length(result.raw_word) == 8
    @test length(result.transition_word) == 7
    @test result.metadata.initial_event_included
    @test result.events[1].time == 0.0
    @test result.events[1].state == seed.u0
    @test all(event -> event.rate > 0, result.events)

    next_system = CoupledODEs(rossler_example!, zeros(3), [0.301, 5.5])
    next_seed = init_saddle_focus(next_system; previous=seed)

    @test next_seed.diagnostics.converged
    @test next_seed.parameters == [0.301, 5.5]
    @test seed.parameters == [0.3, 5.5]
    @test result.metadata.parameters == [0.3, 5.5]
    @test norm(next_seed.u0 - seed.u0) < 0.02
    @test dot(vec(next_seed.Q0), vec(seed.Q0)) > 0.99
    @test DynamicalSystemsBase.current_state(rossler) == zeros(3)
    @test DynamicalSystemsBase.current_state(next_system) == zeros(3)
    @test_throws ArgumentError flow_kneading(FlowKneadingProblem(
        next_system;
        initializer=seed,
        capture,
        word_length=7,
    ))
end

@testset "Literal documentation examples with a narrow 2×2 scan smoke test" begin
    document_path = joinpath(@__DIR__, "..", "docs", "src", "flow-kneading.md")
    source = read(document_path, String)
    blocks = [match.captures[1] for match in eachmatch(r"(?ms)^```julia\n(.*?)^```", source)]
    @test length(blocks) >= 6
    examples_module = Module(gensym(:FlowKneadingDocumentation))
    mktempdir() do directory
        cd(directory) do
            for (index, block) in enumerate(blocks)
                reduced = replace(block,
                    "range(2.0, 7.0; length = 64)" => "range(5.45, 5.55; length = 2)",
                    "range(0.30, 0.55; length = 64)" => "range(0.30, 0.301; length = 2)",
                    "maximum_time = 2_000.0" => "maximum_time = 100.0",
                )
                @testset "Documentation block $index" begin
                    Base.include_string(examples_module, reduced, document_path)
                    if isdefined(examples_module, :result)
                        result = getfield(examples_module, :result)
                        @test result isa FlowKneadingResult
                        @test result.complete
                    end
                end
            end
            initialization = getfield(examples_module, :initialization)
            next_initialization = getfield(examples_module, :next_initialization)
            diagram = getfield(examples_module, :diagram)
            @test initialization.diagnostics.converged
            @test next_initialization.diagnostics.converged
            @test size(diagram.statuses) == (2, 2)
            @test all(!=(:not_started), diagram.statuses)
            @test any(==(:complete), diagram.statuses)
            @test isfile("output/rossler-flow-kneading.tsv")
            @test length(readlines("output/rossler-flow-kneading.tsv")) == 5
        end
    end
end
