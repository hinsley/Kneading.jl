const _SF_FIRST_POINTS = 8

const _SF_BATCH_MESSAGES = (
    "ok",
    "launch radius is below the resolvable minimum",
    "not enough accepted extrema before the integration limit",
    "the seeded trajectory or sensitivity became unbounded",
    "an event is insufficiently transverse",
    "not enough accepted extrema before the integration limit (MaxIters)",
)

_sf_invalid_point(rho, message) = (; valid=false, rho=Float64(rho), residual=NaN, derivative=NaN, message)

function _sf_point(rho, time, state, tangent, transversality, next_time, next_state, next_tangent,
    next_transversality, capture, options)
    derivative = tangent[capture.index]
    next_derivative = next_tangent[capture.index]
    abs(derivative) > options.denominator_tolerance ||
        return _sf_invalid_point(rho, "return-coordinate derivative denominator is too small")
    residual = next_derivative / derivative
    all(isfinite, tangent) && isfinite(residual) || return _sf_invalid_point(rho, "nonfinite event sensitivity")
    return (; valid=true, rho=Float64(rho), residual, derivative, time, state, tangent, transversality,
        next_time, next_state, next_derivative, next_transversality, message="ok")
end

function _sf_cpu_point(entry, rho, M)
    (; ctx, s, ray) = entry
    exp(rho) > s.options.minimum_radius || return _sf_invalid_point(rho, "launch radius is below the resolvable minimum")
    events, failure = _sf_events(ctx, s.capture, ray, rho, M + 1, s.options)
    isempty(failure) || return _sf_invalid_point(rho, failure)
    current, following = events[M], events[M+1]
    return _sf_point(rho, current.time, current.state, current.tangent, current.transversality,
        following.time, following.state, following.tangent, following.transversality, s.capture, s.options)
end

function _sf_first_entry(s, system=nothing)
    ray = _sf_anchor(s)
    guard = isnothing(s.options.launch_guard_time) ? pi / abs(imag(ray.eigenvalue)) : Float64(s.options.launch_guard_time)
    revolution = 2pi * real(ray.eigenvalue) / abs(imag(ray.eigenvalue))
    start = isnothing(s.initial_rho) ? log(s.initial_radius) - s.initial_event_index * revolution : Float64(s.initial_rho)
    return (; system, ctx=s.ctx, s, ray, guard, revolution, start)
end

mutable struct _SFFirstSearch
    entry::Int
    M::Int
    start::Float64
    step::Float64
    window::Int
    upper::Float64
    next_sample::Int
    previous::Any
    left::Any
    right::Any
    launches::Int
    status::Symbol
    root::Any
    curvature::Float64
    message::String
end

function _SFFirstSearch(entries, k, M)
    entry = entries[k]
    start = entry.start - (M - entry.s.initial_event_index) * entry.revolution
    samples = entry.s.options.rho_samples
    return _SFFirstSearch(k, M, start, entry.revolution / samples, samples, entry.s.options.rho_range[2], 0,
        nothing, nothing, nothing, 0, :scan, nothing, NaN, "")
end

function _sf_first_bracket(a, b, kind)
    a.valid && b.valid || return false
    direction = sign(a.derivative)
    direction != 0 && direction == sign(b.derivative) || return false
    (a.residual > 0 && b.residual <= 0) || (a.residual < 0 && b.residual >= 0) || return false
    curvature = (b.residual - a.residual) * direction
    return kind == :minimum ? curvature > 0 : curvature < 0
end

function _sf_search_rhos(search)
    if search.status === :scan
        samples = search.next_sample:search.next_sample+search.window-1
        return filter(rho -> rho <= search.upper, [search.start + j * search.step for j in samples])
    end
    a, b = search.left, search.right
    rhos = [a.rho + (b.rho - a.rho) * k / _SF_FIRST_POINTS for k in 1:_SF_FIRST_POINTS-1]
    secant = a.rho - a.residual * (b.rho - a.rho) / (b.residual - a.residual)
    isfinite(secant) && a.rho < secant < b.rho && push!(rhos, secant)
    return sort!(rhos)
end

function _sf_search_converge!(search, options)
    a, b = search.left, search.right
    best = abs(a.residual) <= abs(b.residual) ? a : b
    if abs(best.residual) <= options.criticality_tolerance || b.rho - a.rho <= options.criticality_tolerance
        search.status = :found
        search.root = best
        search.curvature = (b.residual - a.residual) / (b.rho - a.rho) / best.derivative
    elseif search.launches >= options.max_bisection_iterations
        search.status = :failed
        search.message = "the bracketed critical point did not converge in $(search.launches) evaluation rounds"
    end
    return search
end

function _sf_search_update!(search, points, entry)
    options = entry.s.options
    kind = options.critical_kind
    search.launches += 1
    if search.status === :scan
        for point in points
            if !isnothing(search.previous) && _sf_first_bracket(search.previous, point, kind)
                search.left, search.right = search.previous, point
                search.status = :refine
                return _sf_search_converge!(search, options)
            end
            search.previous = point
        end
        search.next_sample += search.window
        if search.start + search.next_sample * search.step > search.upper
            search.status = :failed
            search.message = "no $(kind) critical point was found between rho=$(round(search.start; digits=4)) and rho_range[2]=$(search.upper) at event index $(search.M); adjust initial_radius, initial_event_index, or rho_samples"
        end
        return search
    end
    sequence = vcat(search.left, points, search.right)
    index = findfirst(k -> _sf_first_bracket(sequence[k], sequence[k+1], kind), 1:length(sequence)-1)
    if isnothing(index)
        search.status = :failed
        search.message = "the bracket of the first $(kind) critical point was lost at event index $(search.M): $(join(unique(p.message for p in points if !p.valid), "; "))"
        return search
    end
    search.left, search.right = sequence[index], sequence[index+1]
    return _sf_search_converge!(search, options)
end

function _sf_first_compare(entry, root, next_root)
    tangent = _sf_first_tangent(entry, root)
    next_tangent = _sf_first_tangent(entry, next_root)
    dot(tangent, next_tangent) < 0 && (next_tangent .*= -1)
    return norm(collect(next_root.state) - collect(root.state)), norm(next_tangent - tangent)
end

function _sf_on_section(ctx, capture, state, time)
    u = collect(Float64, state)
    for _ in 1:3
        flow = collect(ctx.f(u, ctx.p, time))
        h = flow[capture.index]
        rate = _event_rate(ctx, capture, u, time)
        isfinite(h) && isfinite(rate) && rate != 0 && h != 0 || break
        u .-= (h / rate) .* flow
    end
    return u
end

function _sf_first_event(entry, point)
    u = _sf_on_section(entry.ctx, entry.s.capture, point.state, point.time)
    n = length(u)
    return _SaddleFocusEvent(point.time, u, collect(Float64, point.tangent), fill(NaN, n),
        point.derivative, NaN, point.transversality)
end

_sf_first_tangent(entry, point) = _sf_unit_tangent(entry.ctx, _sf_first_event(entry, point), entry.s.options)

function _sf_first_seed(entry, search, next_search, converged, state_error, tangent_error, refinement_steps,
    root_event_index, launches, integration)
    (; ctx, s, ray) = entry
    root = search.root
    M = search.M
    event = _sf_first_event(entry, root)
    tangent = _sf_unit_tangent(ctx, event, s.options)
    diagnostics = (; converged, root_converged=true, refinement_checked=s.refine,
        initial_event_index=s.initial_event_index, root_event_index,
        root_rho_range=(search.start, search.upper), root_fallback_used=false,
        failed_root_attempts=NamedTuple[], residual=root.residual, curvature=search.curvature,
        derivative_method=:bracket, state_error, tangent_error, target_distance=0.0,
        branch_tolerance=s.branch_tolerance, iterations=launches, refinement_steps,
        refinement_event_index=s.refine ? M + (converged ? 1 : 0) : nothing,
        event_time=root.time, next_state=collect(Float64, root.next_state),
        event_derivatives=(root.derivative, root.next_derivative),
        event_transversality=(root.transversality, root.next_transversality),
        eigenvalue=ray.eigenvalue, selection=:first, integration)
    if s.refine && !converged
        return SaddleFocusInitializationError(:refinement, "state and tangent did not converge before maximum_event_index=$(s.maximum_event_index)", diagnostics)
    end
    return SaddleFocusSeed(event.state, reshape(tangent, :, 1), copy(ray.equilibrium), copy(ray.direction),
        root.rho, M, deepcopy(ctx.p), s.capture, s.critical_kind, diagnostics, s.configuration)
end

function _sf_first_seeds(entries, evaluate; integration=:adaptive)
    count = length(entries)
    outcomes = Vector{Any}(nothing, count)
    current = [entry.s.initial_event_index for entry in entries]
    steps = zeros(Int, count)
    searches = Dict{Tuple{Int,Int},_SFFirstSearch}()
    for k in 1:count
        searches[(k, current[k])] = _SFFirstSearch(entries, k, current[k])
        if entries[k].s.refine && current[k] < entries[k].s.maximum_event_index
            searches[(k, current[k] + 1)] = _SFFirstSearch(entries, k, current[k] + 1)
        end
    end
    rounds = 0
    while true
        active = [search for search in values(searches) if search.status in (:scan, :refine) && isnothing(outcomes[search.entry])]
        if !isempty(active)
            requests = [(search, _sf_search_rhos(search)) for search in active]
            items = [(search.entry, rho, search.M) for (search, rhos) in requests for rho in rhos]
            points = evaluate(items)
            rounds += 1
            offset = 0
            for (search, rhos) in requests
                _sf_search_update!(search, points[offset+1:offset+length(rhos)], entries[search.entry])
                offset += length(rhos)
            end
        end
        for k in 1:count
            isnothing(outcomes[k]) || continue
            entry = entries[k]
            M = current[k]
            search = searches[(k, M)]
            if search.status === :failed
                outcomes[k] = SaddleFocusInitializationError(:root, search.message, (; event_index=M, rho_range=(search.start, search.upper)))
                continue
            end
            search.status === :found || continue
            if !entry.s.refine || M >= entry.s.maximum_event_index
                outcomes[k] = _sf_first_seed(entry, search, nothing, false, NaN, NaN, steps[k],
                    entry.s.initial_event_index, rounds, integration)
                continue
            end
            next_search = searches[(k, M + 1)]
            if next_search.status === :failed
                outcomes[k] = SaddleFocusInitializationError(:root, next_search.message, (; event_index=M + 1))
                continue
            end
            next_search.status === :found || continue
            state_error, tangent_error = try
                _sf_first_compare(entry, search.root, next_search.root)
            catch error
                error isa SaddleFocusInitializationError || rethrow()
                outcomes[k] = error
                continue
            end
            steps[k] += 1
            if state_error <= entry.s.options.state_tolerance && tangent_error <= entry.s.options.tangent_tolerance
                outcomes[k] = _sf_first_seed(entry, search, next_search, true, state_error, tangent_error,
                    steps[k], entry.s.initial_event_index, rounds, integration)
            elseif M + 1 < entry.s.maximum_event_index
                current[k] = M + 1
                searches[(k, M + 2)] = _SFFirstSearch(entries, k, M + 2)
            else
                outcomes[k] = _sf_first_seed(entry, next_search, nothing, false, state_error, tangent_error,
                    steps[k], entry.s.initial_event_index, rounds, integration)
            end
        end
        all(!isnothing, outcomes) && break
    end
    return outcomes
end

function _sf_first_initialize(s)
    entry = _sf_first_entry(s)
    outcome = only(_sf_first_seeds([entry], items -> [_sf_cpu_point(entry, rho, M) for (_, rho, M) in items]))
    outcome isa Exception && throw(outcome)
    return outcome
end

struct _SFBatchOptions{T,A}
    capture_index::Int
    maximum::Bool
    accept::A
    dt::Float64
    max_time::Float64
    max_state::Float64
    maxiters::Int
    event_tolerance::Float64
end

_SFBatchOptions{T}(capture, options, dt) where {T} = _SFBatchOptions{T,typeof(capture.accept)}(
    capture.index, capture isa LocalMaximum, capture.accept, Float64(dt), Float64(options.max_time),
    Float64(options.max_state), Int(options.maxiters), Float64(options.event_tolerance))

struct _SFBatchItem{N,P}
    parameters::P
    equilibrium::SVector{N,Float64}
    direction::SVector{N,Float64}
    rho::Float64
    event_index::Int
    guard::Float64
end

struct _SFBatchEvaluation{N}
    status::Int8
    time::Float64
    state::SVector{N,Float64}
    tangent::SVector{N,Float64}
    transversality::Float64
    next_time::Float64
    next_state::SVector{N,Float64}
    next_tangent::SVector{N,Float64}
    next_transversality::Float64
end

@inline function _sf_batch_failure(::Type{_SFBatchEvaluation{N}}, status) where {N}
    zero_state = zero(SVector{N,Float64})
    return _SFBatchEvaluation{N}(Int8(status), NaN, zero_state, zero_state, NaN, NaN, zero_state, zero_state, NaN)
end

@inline function _sf_batch_rate(rule, o::_SFBatchOptions{T}, p, u, t) where {T}
    flow = _batch_field(rule, u, p, t, T)
    _, rate = _batch_jvp(rule, u, flow, p, t, T)
    return flow, rate[o.capture_index]
end

@inline function _sf_batch_event(rule, o::_SFBatchOptions{T}, p, u, v, t, step, fraction) where {T}
    tau = fraction * step
    for _ in 1:3
        trial, _ = _batch_rk4(rule, u, v, p, t, tau, T)
        flow, rate = _sf_batch_rate(rule, o, p, trial, t + tau)
        correction = flow[o.capture_index] / rate
        isfinite(correction) && (tau = clamp(tau - correction, 0.0, step))
    end
    state, sensitivity = _batch_rk4(rule, u, v, p, t, tau, T)
    flow, rate = _sf_batch_rate(rule, o, p, state, t + tau)
    _, image = _batch_jvp(rule, state, sensitivity, p, t + tau, T)
    tangent = sensitivity - (image[o.capture_index] / rate) * flow
    return t + tau, state, tangent, rate
end

function _batch_sf_evaluate(rule, o::_SFBatchOptions{T}, item::_SFBatchItem{N}) where {T,N}
    E = _SFBatchEvaluation{N}
    radius = exp(item.rho)
    p = item.parameters
    u = item.equilibrium + radius * item.direction
    v = radius * item.direction
    t = 0.0
    h_previous = _batch_field(rule, u, p, t, T)[o.capture_index]
    found = 0
    time = NaN
    state = zero(SVector{N,Float64})
    tangent = state
    transversality = NaN
    iterations = 0
    while t < o.max_time
        iterations += 1
        iterations > o.maxiters && return _sf_batch_failure(E, 6)
        step = min(o.dt, o.max_time - t)
        next_u, next_v = _batch_rk4(rule, u, v, p, t, step, T)
        all(isfinite, next_u) && all(isfinite, next_v) && maximum(abs, next_u) <= o.max_state ||
            return _sf_batch_failure(E, 4)
        h_next = _batch_field(rule, next_u, p, t + step, T)[o.capture_index]
        crossed = o.maximum ? h_previous > 0 && h_next <= 0 : h_previous < 0 && h_next >= 0
        if crossed
            fraction = clamp(h_previous / (h_previous - h_next), 0.0, 1.0)
            event_time = t + fraction * step
            if event_time >= item.guard
                candidate = u + fraction * (next_u - u)
                _, rate = _sf_batch_rate(rule, o, p, candidate, event_time)
                if (o.maximum ? rate < 0 : rate > 0) && o.accept(candidate, p, event_time)
                    abs(rate) > o.event_tolerance || return _sf_batch_failure(E, 5)
                    found += 1
                    if found >= item.event_index
                        event_time, event_state, event_tangent, event_rate =
                            _sf_batch_event(rule, o, p, u, v, t, step, fraction)
                        isfinite(event_rate) && abs(event_rate) > o.event_tolerance || return _sf_batch_failure(E, 5)
                        if found == item.event_index
                            time, state, tangent, transversality = event_time, event_state, event_tangent, event_rate
                        else
                            return E(Int8(1), time, state, tangent, transversality,
                                event_time, event_state, event_tangent, event_rate)
                        end
                    end
                end
            end
        end
        u, v = next_u, next_v
        h_previous = h_next
        t += step
    end
    return _sf_batch_failure(E, 3)
end

function _launch_saddle_focus end

function _sf_batch_point(evaluation, rho, capture, options)
    evaluation.status == 1 || return _sf_invalid_point(rho, _SF_BATCH_MESSAGES[evaluation.status])
    return _sf_point(rho, evaluation.time, evaluation.state, evaluation.tangent, evaluation.transversality,
        evaluation.next_time, evaluation.next_state, evaluation.next_tangent, evaluation.next_transversality,
        capture, options)
end

function _sf_batch_seeds(entries, backend, precision, dt)
    isempty(entries) && return Any[]
    backend isa _batch_extension().KernelAbstractions.Backend || throw(ArgumentError(
        "backend must be a KernelAbstractions backend such as CPU() or CUDABackend()"))
    precision <: AbstractFloat || throw(ArgumentError("precision must be a floating-point type"))
    reference = first(entries)
    rule = _batch_rule(reference.system)
    capture = reference.s.capture
    keys = (:max_time, :max_state, :maxiters, :event_tolerance, :minimum_radius, :denominator_tolerance)
    for entry in entries
        _batch_rule(entry.system) === rule && isequal(entry.s.capture, capture) &&
            all(key -> getfield(entry.s.options, key) == getfield(reference.s.options, key), keys) ||
            throw(ArgumentError("batched saddle-focus seeds must share their rule, capture, and integration limits"))
    end
    for (name, value) in ((:rule, rule), (:accept, capture.accept))
        isbits(value) || throw(ArgumentError(
            "the $name must be isbits to run on a device; avoid capturing arrays or globals"))
    end
    n = length(reference.ctx.u0)
    parameters = [_batch_parameters(entry.ctx.p, precision) for entry in entries]
    P = typeof(first(parameters))
    isbitstype(P) && all(p -> typeof(p) === P, parameters) || throw(ArgumentError(
        "parameters must be real vectors of one length or a shared isbits type"))
    options = _SFBatchOptions{precision}(capture, reference.s.options, dt)
    function evaluate(requests)
        items = [_SFBatchItem{n,P}(parameters[k], SVector{n,Float64}(entries[k].ray.equilibrium),
            SVector{n,Float64}(entries[k].ray.direction), rho, M, entries[k].guard) for (k, rho, M) in requests]
        evaluations = _launch_saddle_focus(backend, rule, options, items, _SFBatchEvaluation{n})
        return [_sf_batch_point(evaluation, rho, capture, entries[k].s.options)
            for (evaluation, (k, rho, _)) in zip(evaluations, requests)]
    end
    return _sf_first_seeds(entries, evaluate; integration=(; method=:rk4, dt, precision))
end

_sf_first_selection(initializer) =
    initializer isa SaddleFocusInitializer && get(initializer.options, :selection, :continuation) === :first

function _sf_scan_entry(system, initializer, capture, guess)
    options = merge(initializer.options, (; capture))
    isnothing(guess) || (options = merge(options, (; equilibrium_guess=guess)))
    return _sf_first_entry(_sf_settings(system; options...), system)
end
