const _BATCH_STATUSES = (:complete, :maximum_time, :ambiguous_sign, :state_limit,
    :numerical_failure, :nonfinite_state, :maximum_iterations, :degenerate_tangent_seed)
const _BATCH_COMPLETE = Int8(1)
const _BATCH_MAXIMUM_TIME = Int8(2)
const _BATCH_AMBIGUOUS_SIGN = Int8(3)
const _BATCH_STATE_LIMIT = Int8(4)
const _BATCH_NUMERICAL_FAILURE = Int8(5)
const _BATCH_NONFINITE_STATE = Int8(6)
const _BATCH_MAXIMUM_ITERATIONS = Int8(7)
const _BATCH_DEGENERATE_SEED = Int8(8)

struct _BatchClock{T}
    dt::T
    last_step::T
    steps::Int
    maximum_time::Float64
end

function _BatchClock{T}(dt, maximum_time) where {T}
    step = T(dt)
    count = maximum_time / step
    steps = count < typemax(Int) ? ceil(Int, count) : typemax(Int)
    last_step = T(maximum_time - (steps - 1) * Float64(step))
    last_step > 0 || ((steps, last_step) = (steps - 1, step))
    return _BatchClock{T}(step, last_step, steps, Float64(maximum_time))
end

@inline _batch_start(::_BatchClock{Float64}) = 0.0
@inline _batch_start(::_BatchClock) = 0
@inline _batch_running(c::_BatchClock, t::Float64) = t < c.maximum_time
@inline _batch_running(c::_BatchClock, n::Int) = n < c.steps
@inline _batch_step(c::_BatchClock, t::Float64) = min(c.dt, c.maximum_time - t)
@inline _batch_step(c::_BatchClock, n::Int) = n + 1 < c.steps ? c.dt : c.last_step
@inline _batch_tick(t::Float64, step) = t + step
@inline _batch_tick(n::Int, step) = n + 1
@inline _batch_time(c::_BatchClock, t::Float64) = t
@inline _batch_time(c::_BatchClock, n::Int) = n < c.steps ? n * Float64(c.dt) : c.maximum_time
@inline _batch_rule_time(c::_BatchClock, t::Float64) = t
@inline _batch_rule_time(c::_BatchClock{T}, n::Int) where {T} = T(n) * c.dt

struct _BatchOptions{T,A,O}
    capture_index::Int
    maximum::Bool
    accept::A
    observable::O
    word_length::Int
    transient_events::Int
    clock::_BatchClock{T}
    max_state::T
    sign_atol::Float64
    sign_rtol::Float64
    flow_norm::T
    projected_norm::T
    unit_norm::T
    minimum_event_separation::Float64
    maxiters::Int
    seed_invariance_tolerance::Float64
end

function _BatchOptions{T}(capture, observable, o) where {T}
    tolerance = isnothing(o.seed_invariance_tolerance) ? 0.0 : o.seed_invariance_tolerance
    return _BatchOptions{T,typeof(capture.accept),typeof(observable)}(
        capture.index, capture isa LocalMaximum, capture.accept, observable, o.word_length,
        o.transient_events, _BatchClock{T}(o.dt, o.maximum_time), T(o.max_state), o.sign_atol,
        o.sign_rtol, T(o.tolerances.flow_norm), T(o.tolerances.projected_norm),
        T(max(o.tolerances.unit_norm, 64eps(T))), o.minimum_event_separation, o.maxiters, tolerance)
end

struct _BatchTag end

struct _BatchBits128
    high::UInt64
    low::UInt64
end

Base.zero(::Type{_BatchBits128}) = _BatchBits128(0, 0)

@inline _batch_push(bits::UInt64, positive::Bool) = (bits << 1) | UInt64(positive)
@inline _batch_push(bits::_BatchBits128, positive::Bool) =
    _BatchBits128((bits.high << 1) | (bits.low >> 63), (bits.low << 1) | UInt64(positive))

_batch_bit(bits::UInt64, j) = (bits >> j) & 0x01 == 0x01
_batch_bit(bits::_BatchBits128, j) = j < 64 ? _batch_bit(bits.low, j) : _batch_bit(bits.high, j - 64)

@inline _batch_dot(a, b) = sum(a .* b)

@inline _batch_split(u::SVector{N,Float64}, ::Type{T}) where {N,T} =
    SVector{N,T}(u), SVector{N,T}(u - SVector{N,Float64}(SVector{N,T}(u)))

@inline _batch_value(u::SVector{N}, error = zero(u)) where {N} =
    SVector{N,Float64}(u) + SVector{N,Float64}(error)

@inline _batch_sum(u::SVector{N,Float64}, error, increment) where {N} = u + increment, error

@inline function _batch_sum(u, error, increment)
    y = increment + error
    s = u + y
    w = s - u
    return s, (u - (s - w)) + (y - w)
end

@inline function _batch_field(rule, u::SVector{N,S}, p, t, ::Type{T}) where {N,S,T}
    return SVector{N,S}(rule(SVector{N,T}(u), p, T(t)))
end

@inline function _batch_jvp(rule, u::SVector{N,S}, v::SVector{N}, p, t, ::Type{T}) where {N,S,T}
    D = ForwardDiff.Dual{_BatchTag,T,1}
    dual = SVector{N,D}(ntuple(i -> D(T(u[i]), ForwardDiff.Partials((T(v[i]),))), Val(N)))
    image = rule(dual, p, T(t))
    flow = SVector{N,S}(ntuple(i -> S(ForwardDiff.value(image[i])), Val(N)))
    derivative = SVector{N,S}(ntuple(i -> S(ForwardDiff.partials(image[i], 1)), Val(N)))
    return flow, derivative
end

@inline function _batch_rk4(rule, u, u_error, v, v_error, p, t, dt, ::Type{T}) where {T}
    k1u, k1v = _batch_jvp(rule, u, v, p, t, T)
    k2u, k2v = _batch_jvp(rule, u + dt / 2 * k1u, v + dt / 2 * k1v, p, t + dt / 2, T)
    k3u, k3v = _batch_jvp(rule, u + dt / 2 * k2u, v + dt / 2 * k2v, p, t + dt / 2, T)
    k4u, k4v = _batch_jvp(rule, u + dt * k3u, v + dt * k3v, p, t + dt, T)
    next_u, next_u_error = _batch_sum(u, u_error, dt / 6 * (k1u + 2k2u + 2k3u + k4u))
    next_v, next_v_error = _batch_sum(v, v_error, dt / 6 * (k1v + 2k2v + 2k3v + k4v))
    return next_u, next_u_error, next_v, next_v_error
end

@inline function _batch_project(v, flow, o::_BatchOptions)
    flow_length = sqrt(_batch_dot(flow, flow))
    flow_length > o.flow_norm || return v, false
    tangent = v - (_batch_dot(v, flow) / (flow_length * flow_length)) * flow
    projected_length = sqrt(_batch_dot(tangent, tangent))
    projected_length > o.projected_norm || return v, false
    tangent = tangent / projected_length
    abs(sqrt(_batch_dot(tangent, tangent)) - 1) <= o.unit_norm || return v, false
    return tangent, true
end

@inline function _batch_rate(rule, o::_BatchOptions{T}, p, u, t) where {T}
    _, rate = _batch_jvp(rule, u, _batch_field(rule, u, p, t, T), p, t, T)
    return rate[o.capture_index]
end

@inline function _batch_accept(rate, o::_BatchOptions, p, u, t)
    direction = o.maximum ? rate < 0 : rate > 0
    return direction && o.accept(u, p, t)
end

@inline function _batch_invariant(rule, o::_BatchOptions{T}, p, u::SVector{N,Float64}, t, direction) where {T,N}
    _, image = _batch_jvp(rule, u, direction, p, t, T)
    residual = image - _batch_dot(direction, image) * direction
    jacobian_norm2 = 0.0
    for k in 1:N
        basis = SVector{N,Float64}(ntuple(i -> i == k ? 1.0 : 0.0, Val(N)))
        _, column = _batch_jvp(rule, u, basis, p, t, T)
        jacobian_norm2 += _batch_dot(column, column)
    end
    return sqrt(_batch_dot(residual, residual)) <= o.seed_invariance_tolerance * sqrt(jacobian_norm2)
end

@inline _batch_observable(observable::CoordinateComponent, u, v, p, t) = v[observable.index], 1.0, true

@inline function _batch_observable(observable, u, v, p, t)
    direction = typeof(u)(observable(u, p, t))
    all(isfinite, direction) || return 0.0, 0.0, false
    return _batch_dot(v, direction), sqrt(_batch_dot(direction, direction)), true
end

@inline _batch_store!(::Nothing, item, slot, t, u, v, rate, component) = nothing

@inline function _batch_store!(events, item, slot, t, u, v, rate, component)
    events.times[slot, item] = t
    events.states[slot, item] = u
    events.tangents[slot, item] = v
    events.rates[slot, item] = rate
    events.components[slot, item] = component
    return nothing
end

@inline function _batch_record(rule, o::_BatchOptions{T}, p, u, v, t, direction, record, events, item) where {T}
    status, bits, len, seen, last_event, check_seed = record
    t - last_event > o.minimum_event_separation || return false, record
    rate = _batch_rate(rule, o, p, u, t)
    _batch_accept(rate, o, p, u, t) || return false, record
    last_event = t
    if check_seed && t > 0
        check_seed = false
        _batch_invariant(rule, o, p, u, t, direction) &&
            return true, (_BATCH_DEGENERATE_SEED, bits, len, seen, last_event, check_seed)
    end
    seen += 1
    seen > o.transient_events || return false, (status, bits, len, seen, last_event, check_seed)
    projected, projected_ok = _batch_project(v, _batch_field(rule, u, p, t, T), o)
    projected_ok || return true, (_BATCH_NUMERICAL_FAILURE, bits, len, seen, last_event, check_seed)
    component, observable_norm, observable_ok = _batch_observable(o.observable, u, projected, p, t)
    observable_ok || return true, (_BATCH_NUMERICAL_FAILURE, bits, len, seen, last_event, check_seed)
    _batch_store!(events, item, len + 1, t, u, projected, rate, component)
    tolerance = o.sign_atol + o.sign_rtol * sqrt(_batch_dot(projected, projected)) * observable_norm
    abs(component) <= tolerance &&
        return true, (_BATCH_AMBIGUOUS_SIGN, bits, len, seen, last_event, check_seed)
    bits = _batch_push(bits, component > 0)
    len += 1
    len == o.word_length + 1 && return true, (_BATCH_COMPLETE, bits, len, seen, last_event, check_seed)
    return false, (status, bits, len, seen, last_event, check_seed)
end

function _batch_word(rule, o::_BatchOptions{T}, p, u0::SVector{N,Float64}, v0::SVector{N,Float64},
    direction::SVector{N,Float64}, flag::UInt8, ::Type{B}, events, item) where {T,N,B}
    record = (_BATCH_MAXIMUM_TIME, zero(B), 0, 0, -Inf, (flag & 0x02) != 0x00)
    v0, projected_ok = _batch_project(v0, _batch_field(rule, u0, p, 0.0, T), o)
    projected_ok || return zero(B), 0, _BATCH_NUMERICAL_FAILURE, 0, 0.0
    stopped = false
    stop_time = 0.0
    if (flag & 0x01) != 0x00
        stopped, record = _batch_record(rule, o, p, u0, v0, 0.0, direction, record, events, item)
    end
    u, u_error = _batch_split(u0, T)
    v = SVector{N,T}(v0)
    clock = _batch_start(o.clock)
    h_previous = _batch_field(rule, u, p, 0.0, T)[o.capture_index]
    iterations = 0
    while !stopped && _batch_running(o.clock, clock)
        iterations += 1
        if iterations > o.maxiters
            record = Base.setindex(record, _BATCH_MAXIMUM_ITERATIONS, 1)
            break
        end
        step = _batch_step(o.clock, clock)
        rule_time = _batch_rule_time(o.clock, clock)
        next_u, next_u_error, next_v, _ = _batch_rk4(rule, u, u_error, v, zero(v), p, rule_time, step, T)
        next_clock = _batch_tick(clock, step)
        if !(all(isfinite, next_u) && all(isfinite, next_v))
            record = Base.setindex(record, _BATCH_NONFINITE_STATE, 1)
            clock = next_clock
            break
        end
        if maximum(abs, next_u) > o.max_state
            record = Base.setindex(record, _BATCH_STATE_LIMIT, 1)
            clock = next_clock
            break
        end
        next_flow = _batch_field(rule, next_u, p, rule_time + step, T)
        next_v, projected_ok = _batch_project(next_v, next_flow, o)
        if !projected_ok
            record = Base.setindex(record, _BATCH_NUMERICAL_FAILURE, 1)
            clock = next_clock
            break
        end
        h_next = next_flow[o.capture_index]
        crossed = o.maximum ? h_previous > 0 && h_next <= 0 : h_previous < 0 && h_next >= 0
        if crossed
            fraction = clamp(h_previous / (h_previous - h_next), 0, 1)
            t = _batch_time(o.clock, clock)
            event_time = t + fraction * step
            if event_time > o.minimum_event_separation
                a, b = _batch_value(u, u_error), _batch_value(next_u, next_u_error)
                tangent = _batch_value(v)
                stopped, record = _batch_record(rule, o, p, a + fraction * (b - a),
                    tangent + fraction * (_batch_value(next_v) - tangent), event_time, direction, record, events, item)
                stop_time = event_time
            end
        end
        u, u_error, v = next_u, next_u_error, next_v
        h_previous = h_next
        clock = next_clock
    end
    terminal_time = stopped && record[1] != _BATCH_NUMERICAL_FAILURE ? stop_time : _batch_time(o.clock, clock)
    return record[2], record[3], record[1], record[4], terminal_time
end

function _launch_flow_words end

function _batch_extension()
    extension = Base.get_extension(parentmodule(@__MODULE__), :KneadingKernelAbstractionsExt)
    isnothing(extension) && throw(ArgumentError(
        "a device backend requires KernelAbstractions; run `using KernelAbstractions` first"))
    return extension
end

_batch_parameters(p::AbstractVector{<:Real}, T) = SVector{length(p),T}(p)
_batch_parameters(p::Tuple{Vararg{Real}}, T) = SVector{length(p),T}(p)
_batch_parameters(p, T) = p

function _batch_rule(system)
    SciMLBase.isinplace(system) && throw(ArgumentError(
        "a device backend requires an out-of-place rule returning an SVector"))
    return DynamicalSystemsBase.dynamic_rule(system)
end

_batch_raw_word(bits, len) = Int8[_batch_bit(bits, len - k) ? 1 : -1 for k in 1:len]

function _batch_initial_direction(observable, u0, p)
    observable isa CoordinateComponent &&
        return Float64[i == observable.index for i in eachindex(u0)], observable.index
    direction = collect(Float64, observable(u0, p, 0.0))
    return direction, direction
end

function _flow_kneading_batch(problems, seeds, backend, precision, record_events)
    backend isa _batch_extension().KernelAbstractions.Backend || throw(ArgumentError(
        "backend must be a KernelAbstractions backend such as CPU() or CUDABackend()"))
    count = length(problems)
    reference = first(problems)
    o = reference.options
    o.integration === :rk4 || throw(ArgumentError(
        "a device backend integrates with fixed-step RK4; set integration = :rk4 and dt"))
    for problem in problems
        problem.options == o && isequal(problem.capture, reference.capture) &&
            isequal(problem.observable, reference.observable) || throw(ArgumentError(
                "batched problems must share their capture, observable, and options"))
    end
    capture, observable = reference.capture, reference.observable
    precision <: AbstractFloat || throw(ArgumentError("precision must be a floating-point type"))
    o.word_length < 128 || throw(ArgumentError("a device backend supports at most 127 transitions"))
    rule = _batch_rule(reference.system)
    for (name, value) in ((:rule, rule), (:accept, capture.accept), (:observable, observable))
        isbits(value) || throw(ArgumentError(
            "the $name must be isbits to run on a device; avoid capturing arrays or globals"))
    end
    contexts = [_flow_context(problem.system) for problem in problems]
    n = length(first(contexts).u0)
    all(problem -> _batch_rule(problem.system) === rule, problems) ||
        throw(ArgumentError("batched problems must share their rule"))
    host_parameters = [_batch_parameters(ctx.p, Float64) for ctx in contexts]
    device_parameters = [_batch_parameters(ctx.p, precision) for ctx in contexts]
    P = typeof(first(device_parameters))
    isbitstype(P) && all(p -> typeof(p) === P, device_parameters) || throw(ArgumentError(
        "parameters must be real vectors of one length or a shared isbits type"))
    host_options = _BatchOptions{Float64}(capture, observable, o)

    states = Vector{SVector{n,Float64}}(undef, count)
    tangents = Vector{SVector{n,Float64}}(undef, count)
    directions = fill(zero(SVector{n,Float64}), count)
    flags = zeros(UInt8, count)
    for (k, (seed, ctx)) in enumerate(zip(seeds, contexts))
        length(seed.u0) == n || throw(DimensionMismatch("seed state must match the system"))
        size(seed.Q0) == (n, 1) || throw(DimensionMismatch("seed Q0 must contain one tangent column"))
        all(isfinite, seed.u0) && all(isfinite, seed.Q0) ||
            throw(ArgumentError("seed state and tangent must be finite"))
        if hasproperty(seed, :parameters) && !isequal(seed.parameters, ctx.p)
            throw(ArgumentError("seed parameters differ from the system; continue the initializer first"))
        end
        states[k] = SVector{n,Float64}(seed.u0)
        tangents[k] = SVector{n,Float64}(vec(seed.Q0))
        initial = o.include_initial_event === true ||
            (o.include_initial_event === :auto && hasproperty(seed, :rho) && hasproperty(seed, :event_index))
        if initial
            u, p = states[k], host_parameters[k]
            abs(_batch_field(rule, u, p, 0.0, Float64)[capture.index]) <= o.initial_event_tolerance &&
                _batch_accept(_batch_rate(rule, host_options, p, u, 0.0), host_options, p, u, 0.0) ||
                throw(ArgumentError("the initial state is not an accepted capture event"))
            flags[k] |= 0x01
        end
        if seed isa RealSaddleSeed && !isnothing(o.seed_invariance_tolerance)
            directions[k] = SVector{n,Float64}(seed.seed_direction)
            flags[k] |= 0x02
        end
    end

    options = _BatchOptions{precision}(capture, observable, o)
    B = o.word_length < 64 ? UInt64 : _BatchBits128
    launch(selected) = _launch_flow_words(backend, rule, options, device_parameters[selected],
        states[selected], tangents[selected], directions[selected], flags[selected], B,
        record_events ? o.word_length + 1 : 0)
    outputs = launch(1:count)
    current_seeds = collect(Any, seeds)
    fallback = falses(count)
    restart = Int[]
    for k in 1:count
        outputs.statuses[k] == _BATCH_DEGENERATE_SEED && seeds[k] isa RealSaddleSeed &&
            seeds[k].tangent_seed === :auto || continue
        direction, request = _batch_initial_direction(observable, states[k], host_parameters[k])
        reseeded = try
            RealSaddleInitialization._reseed_tangent(seeds[k], request, direction,
                collect(_batch_field(rule, states[k], host_parameters[k], 0.0, Float64)))
        catch error
            error isa DomainError || rethrow()
            continue
        end
        current_seeds[k] = reseeded
        tangents[k] = SVector{n,Float64}(vec(reseeded.Q0))
        directions[k] = SVector{n,Float64}(reseeded.seed_direction)
        fallback[k] = true
        push!(restart, k)
    end
    if !isempty(restart)
        retry = launch(restart)
        for name in keys(outputs)
            selectdim(getfield(outputs, name), ndims(getfield(outputs, name)), restart) .=
                getfield(retry, name)
        end
    end

    return map(1:count) do k
        len = Int(outputs.lengths[k])
        raw = _batch_raw_word(outputs.bits[k], len)
        status = _BATCH_STATUSES[outputs.statuses[k]]
        events = FlowKneadingEvent[]
        if record_events
            recorded = len + (status === :ambiguous_sign)
            for slot in 1:recorded
                state = collect(outputs.states[slot, k])
                push!(events, FlowKneadingEvent(outputs.times[slot, k], state,
                    collect(outputs.tangents[slot, k]), state[capture.index], outputs.rates[slot, k],
                    outputs.components[slot, k], slot <= len ? raw[slot] : Int8(0)))
            end
        end
        transitions = Int8[raw[i] * raw[i + 1] for i in 1:max(0, len - 1)]
        complete = status === :complete
        detail = status === :degenerate_tangent_seed ?
            "the tangent seed direction spans an invariant subspace of the variational equation" :
            status === :numerical_failure ?
            "the flow-normal projection or the observable direction was degenerate or nonfinite" : ""
        times = Float64[events[i + 1].time - events[i].time for i in 1:max(0, length(events) - 1)]
        metadata = (;
            parameters = deepcopy(contexts[k].p), capture, observable, options = o,
            accepted_events = Int(outputs.accepted[k]), initial_event_included = (flags[k] & 0x01) != 0x00,
            terminal_time = Float64(outputs.terminal_times[k]), tangent_gauge = :accepted_step_flow_normal,
            detail, tangent_seed_fallback = fallback[k], backend, precision,
        )
        FlowKneadingResult(raw, transitions, complete ? _word_code(raw) : BigInt(-1),
            complete ? _word_code(transitions) : BigInt(-1), len, length(transitions), events,
            times, status, complete, current_seeds[k], metadata)
    end
end

"""
    flow_kneading(problems::AbstractArray{<:FlowKneadingProblem};
                  backend=nothing, precision=Float64, threaded=true)

Solve an array of problems and return an array of results with the same shape.
Without a `backend`, each problem is solved by `flow_kneading(problem)`, across
threads when `threaded`.

With a KernelAbstractions backend, such as `CPU()` or `CUDABackend()`, the
initializers run on the CPU, except `SaddleFocusInitializer(...; selection=:first)`,
whose seeds are searched for all problems together on the device; a problem
without such a seed throws its `SaddleFocusInitializationError`. All words are
then integrated together on the device by the `integration = :rk4` method; the
problems must select it and share their rule, capture, observable, and options. This requires
`using KernelAbstractions`, an out-of-place rule returning an `SVector`, and real
parameter vectors (or an isbits parameter value). The rule, the capture's
`accept`, and a function observable run on the device with
`u::SVector{N,Float64}`, so they must avoid allocation and captured arrays.
`precision=Float32` integrates the state and tangent in `Float32`, with a
compensated state update and a step count for time, and reports events in
`Float64`.
"""
function flow_kneading(problems::AbstractArray{<:FlowKneadingProblem};
    backend = nothing, precision::Type = Float64, threaded::Bool = true)
    items = vec(collect(problems))
    parallel = threaded && Threads.nthreads() > 1
    if isnothing(backend)
        precision === Float64 || throw(ArgumentError("precision requires a device backend"))
        results = Vector{FlowKneadingResult}(undef, length(items))
        if parallel
            Threads.@threads :dynamic for k in eachindex(items)
                results[k] = flow_kneading(items[k])
            end
        else
            for k in eachindex(items)
                results[k] = flow_kneading(items[k])
            end
        end
        return reshape(results, size(problems))
    end
    isempty(items) && return reshape(FlowKneadingResult[], size(problems))
    seeds = Vector{Any}(undef, length(items))
    searched = filter(k -> _sf_first_selection(items[k].initializer), eachindex(items))
    others = setdiff(eachindex(items), searched)
    initialize(k) = seeds[k] = _initialize_flow(items[k].system, items[k].initializer, items[k].capture)
    if parallel
        Threads.@threads :dynamic for k in others
            initialize(k)
        end
    else
        foreach(initialize, others)
    end
    if !isempty(searched)
        all(k -> items[k].options.integration === :rk4, searched) || throw(ArgumentError(
            "a device backend integrates with fixed-step RK4; set integration = :rk4 and dt"))
        entries = [_sf_scan_entry(items[k].system, items[k].initializer, items[k].capture, nothing) for k in searched]
        for (k, outcome) in zip(searched, _sf_batch_seeds(entries, backend, precision, items[searched[1]].options.dt))
            outcome isa Exception && throw(outcome)
            seeds[k] = outcome
        end
    end
    return reshape(_flow_kneading_batch(items, seeds, backend, precision, true), size(problems))
end

function _scan_flow_kneading_batch(builder, plane, initializer, capture, observable, threaded,
    store_results, progress, on_error, backend, precision; kwargs...)
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
    problems = Matrix{Any}(nothing, dims...)
    seeds = Matrix{Any}(nothing, dims...)
    anchors = Matrix{Any}(nothing, dims...)
    progress_lock = ReentrantLock()
    searched = _sf_first_selection(initializer)
    report(i, j, result) = isnothing(progress) || lock(() -> progress(i, j, result), progress_lock)

    function visit(i, j, previous)
        system = builder(plane.x[j], plane.y[i])
        if searched
            anchor, message = _try_scan_initialize(system, _SaddleFocusAnchor(initializer), capture, previous, on_error)
            if isnothing(anchor)
                statuses[i, j] = :initialization_failed
                errors[i, j] = message
                report(i, j, nothing)
                return previous
            end
            anchors[i, j] = anchor
            return anchor
        end
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
        problems[i, j] = FlowKneadingProblem(system; initializer = seed, capture, observable, kwargs...)
        seeds[i, j] = seed
        _record_scan_seed!(critical_states, critical_residuals, critical_rhos, event_indices, i, j, seed)
        return seed
    end

    parallel = _continuation_sweep(visit, dims, threaded)
    anchored = findall(!isnothing, anchors)
    if !isempty(anchored)
        template = FlowKneadingProblem(anchors[first(anchored)].system; initializer, capture, observable, kwargs...)
        template.options.integration === :rk4 || throw(ArgumentError(
            "a device backend integrates with fixed-step RK4; set integration = :rk4 and dt"))
        outcomes = _sf_batch_seeds([anchors[c] for c in anchored], backend, precision, template.options.dt)
        for (c, seed) in zip(anchored, outcomes)
            if seed isa Exception
                on_error === :throw && throw(seed)
                statuses[c] = :initialization_failed
                errors[c] = sprint(showerror, seed)
                report(c[1], c[2], nothing)
                continue
            end
            problems[c] = FlowKneadingProblem(anchors[c].system; initializer = seed, capture, observable, kwargs...)
            seeds[c] = seed
            _record_scan_seed!(critical_states, critical_residuals, critical_rhos, event_indices, c[1], c[2], seed)
        end
    end
    initialized = findall(!isnothing, seeds)
    if !isempty(initialized)
        batch = _flow_kneading_batch([problems[c] for c in initialized],
            [seeds[c] for c in initialized], backend, precision, store_results)
        for (c, result) in zip(initialized, batch)
            raw_codes[c] = result.raw_code
            transition_codes[c] = result.transition_code
            raw_lengths[c] = result.raw_length
            transition_lengths[c] = result.transition_length
            raw_words[c] = result.raw_word
            statuses[c] = result.status
            errors[c] = result.metadata.detail
            store_results && (results[c] = result)
            isnothing(progress) || progress(c[1], c[2], result)
        end
    end
    metadata = (; initializer, capture, observable, options = (; kwargs...),
        threaded = parallel, word_encoding = :preservation_is_one, backend, precision)
    return FlowKneadingDiagram(plane, raw_codes, transition_codes, raw_lengths,
        transition_lengths, raw_words, statuses, errors, critical_states, critical_residuals,
        critical_rhos, event_indices, results, metadata)
end
