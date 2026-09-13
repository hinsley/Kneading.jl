using Kneading.FlowKneading

struct FlowScanAnchorProbe
    predecessors::Matrix{Vector{Float64}}
end

function Kneading.FlowKneading._scan_initialize(system, probe::FlowScanAnchorProbe, capture, previous)
    parameters = copy(DynamicalSystemsBase.current_parameters(system))
    j, i = Int.(parameters)
    probe.predecessors[i, j] = isnothing(previous) ? Float64[] : copy(previous.parameters)
    yield()
    return (u0=[1.0, 0.0, 0.0], Q0=reshape([1.0, 0.0, 0.0], :, 1), parameters)
end

@testset "Flow parameter scans" begin
    rule = (u, p, t) -> DynamicalSystemsBase.SVector(-u[2], u[1], -u[3])
    builder = (x, y) -> DynamicalSystemsBase.CoupledODEs(rule, [1.0, 0.0, 0.0], [x, y])
    seed = (u0 = [1.0, 0.0, 0.0], Q0 = reshape([1.0, 0.0, 0.0], :, 1))
    plane = ParameterPlane([1.0, 2.0], [3.0, 4.0])
    serial = scan_flow_kneading(builder, plane; initializer = seed,
        capture = LocalMaximum(1), include_initial_event = true,
        word_length = 2, maximum_time = 15.0, threaded = false, store_results = true)
    threaded = scan_flow_kneading(builder, plane; initializer = seed,
        capture = LocalMaximum(1), include_initial_event = true,
        word_length = 2, maximum_time = 15.0, threaded = true)
    @test all(==(:complete), serial.statuses)
    @test serial.raw_codes == threaded.raw_codes == fill(big(7), 2, 2)
    @test serial.transition_codes == threaded.transition_codes == fill(big(3), 2, 2)
    @test serial.results[1, 1].complete
    @test all(isnothing, threaded.results)
    @test all(==(3), serial.raw_lengths)
    @test all(==(2), serial.transition_lengths)
    failed = scan_flow_kneading(builder, plane; initializer=seed,
        capture=LocalMaximum(1), word_length=1, maximum_time=7.0,
        observable=(u, p, t) -> throw(DomainError(NaN, "unresolved observable")))
    @test all(==(:numerical_failure), failed.statuses)
    @test all(message -> occursin("unresolved observable", message), failed.errors)
    mktempdir() do directory
        path = write_flow_scan(joinpath(directory, "scan.tsv"), serial)
        rows = readlines(path)
        @test length(rows) == 5
        @test startswith(first(rows), "x\ty\tstatus\traw_word")
        @test occursin("\t111\t", rows[2])
    end
    @test_throws ArgumentError scan_flow_kneading(builder, plane;
        initializer = seed, capture = LocalMaximum(1), on_error = :ignore)
end

@testset "Independent threaded continuation anchors" begin
    if Threads.nthreads() > 1
        rule = (u, p, t) -> DynamicalSystemsBase.SVector(-u[2], u[1], -u[3])
        builder = (x, y) -> DynamicalSystemsBase.CoupledODEs(rule, [1.0, 0.0, 0.0], [x, y])
        plane = ParameterPlane(collect(1.0:32.0), collect(1.0:32.0))
        probe = FlowScanAnchorProbe([Float64[] for _ in plane.y, _ in plane.x])
        result = scan_flow_kneading(builder, plane; initializer=probe,
            capture=LocalMaximum(1), include_initial_event=true,
            word_length=0, threaded=true)
        @test all(==(:complete), result.statuses)
        @test isempty(probe.predecessors[1, 1])
        @test all(probe.predecessors[1, j] == [j-1, 1] for j in 2:length(plane.x))
        @test all(probe.predecessors[i, j] == [j, i-1] for j in eachindex(plane.x), i in 2:length(plane.y))
    else
        @test_skip Threads.nthreads() > 1
    end
end
