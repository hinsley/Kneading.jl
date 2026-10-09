"""
Results of [`scan_flow_kneading`](@ref), indexed as `[y_index, x_index]`.

`raw_codes` and `transition_codes` are `-1` for incomplete words; the length
matrices retain the available prefix lengths. `statuses` records orbit failures,
and `errors` records available failure details. `critical_states` and
`critical_residuals` retain initialization diagnostics. `results` contains full
pointwise results only when `store_results=true`.
"""
struct FlowKneadingDiagram{P,M}
    plane::P
    raw_codes::Matrix{BigInt}
    transition_codes::Matrix{BigInt}
    raw_lengths::Matrix{Int}
    transition_lengths::Matrix{Int}
    raw_words::Matrix{Vector{Int8}}
    statuses::Matrix{Symbol}
    errors::Matrix{String}
    critical_states::Matrix{Vector{Float64}}
    critical_residuals::Matrix{Float64}
    critical_rhos::Matrix{Float64}
    event_indices::Matrix{Int}
    results::Matrix{Any}
    metadata::M
end

function _scan_initialize(system, initializer::SaddleFocusInitializer, capture, previous)
    options = merge(initializer.options, (; capture, previous))
    if !isnothing(previous)
        options = merge(options, (; equilibrium_guess = previous.equilibrium,
            initial_event_index = previous.event_index))
    end
    return init_saddle_focus(system; options...)
end

_scan_initialize(system, initializer, capture, previous) =
    _initialize_flow(system, initializer, capture)

function _try_scan_initialize(system, initializer, capture, previous, on_error)
    try
        return _scan_initialize(system, initializer, capture, previous), ""
    catch exception
        if on_error == :throw || !(exception isa SaddleFocusInitializationError ||
            exception isa InvalidRealSaddleInitialState || exception isa DomainError)
            rethrow()
        end
        return nothing, sprint(showerror, exception)
    end
end

function _record_scan_seed!(critical_states, critical_residuals, critical_rhos, event_indices, i, j, seed)
    critical_states[i, j] = copy(seed.u0)
    if seed isa SaddleFocusSeed
        critical_residuals[i, j] = seed.diagnostics.residual
        critical_rhos[i, j] = seed.rho
        event_indices[i, j] = seed.event_index
    end
    return nothing
end

function _continuation_sweep(visit, dims, threaded)
    anchors = Vector{Any}(undef, dims[2])
    previous = nothing
    for j in 1:dims[2]
        previous = visit(1, j, previous)
        anchors[j] = previous
    end
    function column(j)
        column_previous = anchors[j]
        for i in 2:dims[1]
            column_previous = visit(i, j, column_previous)
        end
    end
    parallel = threaded && Threads.nthreads() > 1
    if parallel
        Threads.@threads :dynamic for j in 1:dims[2]
            column(j)
        end
    else
        for j in 1:dims[2]
            column(j)
        end
    end
    return parallel
end

"""
    scan_flow_kneading(builder, plane; initializer, capture,
                      observable=CoordinateComponent(capture.index),
                      threaded=true, store_results=false, progress=nothing,
                      backend=nothing, precision=Float64, kwargs...)

Build a fresh `CoupledODEs` with `builder(x,y)` at each parameter point and
calculate its flow-kneading word. The first row is continued serially; each
column then continues along increasing `y`, independently across threads.
The builder and custom event or observable callbacks must not mutate shared state.

Pass a `SaddleFocusInitializer` to reuse equilibrium and critical-point anchors
between neighboring parameters. `kwargs` are forwarded to `FlowKneadingProblem`.
Expected numerical initialization failures are recorded rather than painted as
valid words; unexpected programming errors are rethrown. Set `on_error=:throw`
to stop on any initialization failure. The optional `progress(i,j,result)` runs
under a lock, with `result=nothing` when initialization fails.

With a KernelAbstractions `backend`, every point is first initialized on the CPU
with the same continuation, and all words are then integrated together on the
device, as in the vector method of [`flow_kneading`](@ref). This requires
`integration=:rk4`; `precision` is passed on as well.
"""
function scan_flow_kneading(
    builder,
    plane::ParameterPlane;
    initializer,
    capture::ExtremumCapture,
    observable = CoordinateComponent(capture.index),
    threaded::Bool = true,
    store_results::Bool = false,
    progress = nothing,
    on_error::Symbol = :record,
    backend = nothing,
    precision::Type = Float64,
    kwargs...,
)
    on_error in (:record, :throw) || throw(ArgumentError("on_error must be :record or :throw"))
    isnothing(backend) || return _scan_flow_kneading_batch(builder, plane, initializer, capture,
        observable, threaded, store_results, progress, on_error, backend, precision; kwargs...)
    precision === Float64 || throw(ArgumentError("precision requires a device backend"))
    dims = (length(plane.y), length(plane.x))
    raw_codes = fill(big(-1), dims)
    transition_codes = fill(big(-1), dims)
    raw_lengths = zeros(Int, dims)
    transition_lengths = zeros(Int, dims)
    raw_words = [Int8[] for _ in 1:dims[1], _ in 1:dims[2]]
    statuses = fill(:not_started, dims)
    errors = fill("", dims)
    critical_states = [Float64[] for _ in 1:dims[1], _ in 1:dims[2]]
    critical_residuals = fill(NaN, dims)
    critical_rhos = fill(NaN, dims)
    event_indices = zeros(Int, dims)
    results = Matrix{Any}(nothing, dims...)
    progress_lock = ReentrantLock()

    function visit(i, j, previous)
        system = builder(plane.x[j], plane.y[i])
        seed, message = _try_scan_initialize(system, initializer, capture, previous, on_error)
        if isnothing(seed)
            statuses[i, j] = :initialization_failed
            errors[i, j] = message
            if !isnothing(progress)
                lock(progress_lock) do
                    progress(i, j, nothing)
                end
            end
            return previous
        end
        result = flow_kneading(FlowKneadingProblem(
            system; initializer = seed, capture, observable, kwargs...,
        ))
        raw_codes[i, j] = result.raw_code
        transition_codes[i, j] = result.transition_code
        raw_lengths[i, j] = result.raw_length
        transition_lengths[i, j] = result.transition_length
        raw_words[i, j] = copy(result.raw_word)
        statuses[i, j] = result.status
        errors[i, j] = result.metadata.detail
        _record_scan_seed!(critical_states, critical_residuals, critical_rhos, event_indices, i, j, seed)
        store_results && (results[i, j] = result)
        if !isnothing(progress)
            lock(progress_lock) do
                progress(i, j, result)
            end
        end
        return seed
    end

    parallel = _continuation_sweep(visit, dims, threaded)
    metadata = (; initializer, capture, observable, options = (; kwargs...),
        threaded = parallel, word_encoding = :preservation_is_one)
    return FlowKneadingDiagram(plane, raw_codes, transition_codes, raw_lengths,
        transition_lengths, raw_words, statuses, errors, critical_states, critical_residuals, critical_rhos,
        event_indices, results, metadata)
end

"""
    write_flow_scan(path, diagram)

Write a tab-separated parameter scan with complete codes, available raw-word
prefixes, failure status, and initialization diagnostics. The first two columns
use the parameter-plane axis names. Create parent directories when needed.
"""
function write_flow_scan(path::AbstractString, diagram::FlowKneadingDiagram)
    mkpath(dirname(abspath(path)))
    n = maximum(length, diagram.critical_states)
    state_names = n == 3 ? ["critical_x", "critical_y", "critical_z"] :
        ["critical_u_$i" for i in 1:n]
    header = [diagram.plane.xname, diagram.plane.yname, "status", "raw_word",
        state_names..., "critical_rho", "critical_residual", "critical_event_index",
        "raw_code", "transition_code", "raw_length", "transition_length", "error"]
    open(path, "w") do io
        println(io, join(header, '\t'))
        for i in eachindex(diagram.plane.y), j in eachindex(diagram.plane.x)
            state = diagram.critical_states[i, j]
            values = length(state) == n ? state : fill(NaN, n)
            word = join(sign > 0 ? '1' : '0' for sign in diagram.raw_words[i, j])
            detail = replace(diagram.errors[i, j], '\t' => ' ', '\n' => ' ', '\r' => ' ')
            row = (diagram.plane.x[j], diagram.plane.y[i], diagram.statuses[i, j], word,
                values..., diagram.critical_rhos[i, j], diagram.critical_residuals[i, j],
                diagram.event_indices[i, j], diagram.raw_codes[i, j], diagram.transition_codes[i, j],
                diagram.raw_lengths[i, j], diagram.transition_lengths[i, j], detail)
            println(io, join(row, '\t'))
        end
    end
    return abspath(path)
end
