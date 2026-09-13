"""
    RealSaddleInitializer(; equilibrium_guess, launch_distance=1e-6, kwargs...)

Configure initialization from a real saddle with one unstable direction and a
real, simple leading stable direction. The unstable branch is `1` or `-1`.
Explicit `unstable_reference` and `stable_reference` vectors fix orientation;
their defaults are all-ones vectors. Spectral or orientation failures are errors.
"""
struct RealSaddleInitializer{O<:NamedTuple}
    options::O
end

function RealSaddleInitializer(;
    equilibrium_guess,
    launch_distance = 1e-6,
    unstable_branch = 1,
    unstable_reference = nothing,
    stable_reference = nothing,
    equilibrium_branch = :selected,
    equilibrium_branch_check = u -> true,
    max_root_iterations = 30,
    tolerances = RealSaddleTolerances(1e-10, 1e-8, 1e-8, 1e-8, 1e-8, 1e-10, 1e-12, 1e-12),
)
    launch_distance > 0 && isfinite(launch_distance) ||
        throw(ArgumentError("launch_distance must be finite and positive"))
    unstable_branch in (-1, 1) || throw(ArgumentError("unstable_branch must be 1 or -1"))
    max_root_iterations > 0 || throw(ArgumentError("max_root_iterations must be positive"))
    return RealSaddleInitializer((;
        equilibrium_guess = collect(float.(equilibrium_guess)), launch_distance,
        unstable_branch, unstable_reference, stable_reference, equilibrium_branch,
        equilibrium_branch_check, max_root_iterations, tolerances,
    ))
end

function _initialize_flow(system, initializer::RealSaddleInitializer, capture)
    o = initializer.options
    n = length(o.equilibrium_guess)
    unstable = isnothing(o.unstable_reference) ? ones(n) : o.unstable_reference
    stable = isnothing(o.stable_reference) ? ones(n) : o.stable_reference
    return init_real_saddle(
        system, o.equilibrium_guess, o.launch_distance, o.tolerances;
        max_root_iterations = o.max_root_iterations,
        equilibrium_branch = o.equilibrium_branch,
        equilibrium_branch_check = o.equilibrium_branch_check,
        unstable_reference = unstable, stable_reference = stable,
        unstable_branch = o.unstable_branch,
    )
end

function _initialize_flow(system, initializer::SaddleFocusInitializer, capture)
    return init_saddle_focus(system; merge(initializer.options, (; capture))...)
end

function _initialize_flow(system, initializer, capture)
    hasproperty(initializer, :u0) && hasproperty(initializer, :Q0) ||
        throw(ArgumentError("initializer must be an initializer configuration or a seed with u0 and Q0"))
    return deepcopy(initializer)
end

"""
    FlowKneadingProblem(system; initializer, capture, word_length=16, kwargs...)

Capture extrema and signs of a unit flow-normal tangent. `word_length` is the
requested number of adjacent-sign transitions, requiring one additional raw
sign. `observable=CoordinateComponent(capture.index)` selects the component;
a function `(u, p, t) -> direction` can instead supply an observable vector.

The initializer can be `RealSaddleInitializer`, `SaddleFocusInitializer`, a
returned seed, or `(u0=state, Q0=reshape(tangent, :, 1))`. With
`include_initial_event=:auto`, saddle-focus seeds include the critical-point
event at time zero. Explicit `true` requires an accepted extremum at the seed.
`transient_events` skips accepted events before encoding.

Adaptive Tsit5 integration projects the tangent after every accepted step and
at captured events. `integration=:rk4` selects fixed-step RK4 with linearly
interpolated events, provided for reproducing fixed-step scans; set `dt`
explicitly for those scans and check results under step refinement.

`maximum_time` bounds elapsed integration time. A missing or ambiguous symbol
produces an incomplete result, not a padded word. No user system is advanced.
"""
struct FlowKneadingProblem{S,I,C,O,P<:NamedTuple}
    system::S
    initializer::I
    capture::C
    observable::O
    options::P
end

function FlowKneadingProblem(
    system;
    initializer,
    capture,
    observable = CoordinateComponent(capture.index),
    transient_events = 0,
    word_length = 16,
    maximum_time = 1000.0,
    include_initial_event = :auto,
    sign_atol = 1e-10,
    sign_rtol = 1e-10,
    tolerances = FlowNormalTolerances(1e-12, 1e-12, 1e-10),
    max_state = Inf,
    abstol = 1e-10,
    reltol = 1e-9,
    dtmax = 0.05,
    maxiters = 10^7,
    integration = :adaptive,
    dt = 0.02,
    initial_event_tolerance = 1e-7,
    minimum_event_separation = 1e-7,
)
    transient_events isa Integer && transient_events >= 0 ||
        throw(ArgumentError("transient_events must be a nonnegative integer"))
    word_length isa Integer && word_length >= 0 ||
        throw(ArgumentError("word_length must be a nonnegative integer"))
    include_initial_event in (:auto, true, false) ||
        throw(ArgumentError("include_initial_event must be :auto, true, or false"))
    integration in (:adaptive, :rk4) || throw(ArgumentError("integration must be :adaptive or :rk4"))
    for (name, value) in ((:maximum_time, maximum_time), (:abstol, abstol),
        (:reltol, reltol), (:dtmax, dtmax), (:dt, dt))
        isfinite(value) && value > 0 || throw(ArgumentError("$name must be finite and positive"))
    end
    for (name, value) in ((:sign_atol, sign_atol), (:sign_rtol, sign_rtol),
        (:initial_event_tolerance, initial_event_tolerance),
        (:minimum_event_separation, minimum_event_separation))
        isfinite(value) && value >= 0 || throw(ArgumentError("$name must be finite and nonnegative"))
    end
    max_state > 0 || throw(ArgumentError("max_state must be positive"))
    maxiters isa Integer && maxiters > 0 || throw(ArgumentError("maxiters must be a positive integer"))
    n = length(_flow_context(system).u0)
    1 <= capture.index <= n || throw(ArgumentError("capture index is outside the state"))
    if observable isa CoordinateComponent
        1 <= observable.index <= n || throw(ArgumentError("observable index is outside the state"))
    end
    return FlowKneadingProblem(system, initializer, capture, observable, (;
        transient_events = Int(transient_events), word_length = Int(word_length),
        maximum_time = Float64(maximum_time), include_initial_event,
        sign_atol = Float64(sign_atol), sign_rtol = Float64(sign_rtol), tolerances,
        max_state = Float64(max_state), abstol = Float64(abstol), reltol = Float64(reltol),
        dtmax = Float64(dtmax), maxiters = Int(maxiters), integration, dt = Float64(dt),
        initial_event_tolerance = Float64(initial_event_tolerance),
        minimum_event_separation = Float64(minimum_event_separation),
    ))
end

"""
A captured extremum with elapsed `time`, full `state`, unit flow-normal
`tangent`, captured coordinate `value`, event `rate`, observable `component`,
and its `sign`. A zero sign denotes an unresolved component.
"""
struct FlowKneadingEvent
    time::Float64
    state::Vector{Float64}
    tangent::Vector{Float64}
    value::Float64
    rate::Float64
    component::Float64
    sign::Int8
end

"""
Result of [`flow_kneading`](@ref). `raw_word` holds component signs and
`transition_word` their adjacent products. Codes encode positive signs as one,
negative signs as zero, from left to right; lengths preserve leading zeros.
Both codes are `-1` if `complete` is false. Recorded events and the valid word
prefix remain available in incomplete results. `return_times` contains the
elapsed times between successive recorded events. `initialization` retains
the seed and `metadata` records independent parameters and solver settings.
"""
struct FlowKneadingResult{I,M<:NamedTuple}
    raw_word::Vector{Int8}
    transition_word::Vector{Int8}
    raw_code::BigInt
    transition_code::BigInt
    raw_length::Int
    transition_length::Int
    events::Vector{FlowKneadingEvent}
    return_times::Vector{Float64}
    status::Symbol
    complete::Bool
    initialization::I
    metadata::M
end

function _word_code(word)
    code = BigInt(0)
    for sign in word
        code = 2code + (sign > 0)
    end
    return code
end

function _flow_observable(observable::CoordinateComponent, u, v, p, t)
    return Float64(v[observable.index]), 1.0
end

function _flow_observable(observable, u, v, p, t)
    direction = observable(u, p, t)
    length(direction) == length(u) || throw(DimensionMismatch("observable direction must match the state"))
    all(isfinite, direction) || throw(DomainError(direction, "observable direction must be finite"))
    return Float64(dot(v, direction)), Float64(norm(direction))
end

function _word_tangent(ctx, u, v, t, tolerances)
    tangent = collect(Float64, v)
    project_flow_normal!(tangent, ctx.f(u, ctx.p, t), tolerances)
    return tangent
end

function _word_rk4_step(ctx, u, v, t, dt)
    k1u = ctx.f(u, ctx.p, t)
    k1v = ctx.J(u, ctx.p, t) * v
    u2, v2 = u + dt / 2 * k1u, v + dt / 2 * k1v
    k2u = ctx.f(u2, ctx.p, t + dt / 2)
    k2v = ctx.J(u2, ctx.p, t + dt / 2) * v2
    u3, v3 = u + dt / 2 * k2u, v + dt / 2 * k2v
    k3u = ctx.f(u3, ctx.p, t + dt / 2)
    k3v = ctx.J(u3, ctx.p, t + dt / 2) * v3
    u4, v4 = u + dt * k3u, v + dt * k3v
    k4u = ctx.f(u4, ctx.p, t + dt)
    k4v = ctx.J(u4, ctx.p, t + dt) * v4
    return u + dt / 6 * (k1u + 2k2u + 2k3u + k4u),
        v + dt / 6 * (k1v + 2k2v + 2k3v + k4v)
end

"""
    flow_kneading(problem::FlowKneadingProblem)

Initialize a trajectory and its tangent, capture the requested extrema, and
return raw and transition words with completeness and event diagnostics.
"""
function flow_kneading(problem::FlowKneadingProblem)
    ctx = _flow_context(problem.system)
    seed = _initialize_flow(problem.system, problem.initializer, problem.capture)
    u = collect(Float64, seed.u0)
    n = length(u)
    length(ctx.u0) == n || throw(DimensionMismatch("seed state must match the system"))
    size(seed.Q0) == (n, 1) || throw(DimensionMismatch("seed Q0 must contain one tangent column"))
    all(isfinite, u) && all(isfinite, seed.Q0) || throw(ArgumentError("seed state and tangent must be finite"))
    if hasproperty(seed, :parameters) && !isequal(seed.parameters, ctx.p)
        throw(ArgumentError("seed parameters differ from the system; continue the initializer first"))
    end
    o = problem.options
    v = _word_tangent(ctx, u, vec(seed.Q0), 0.0, o.tolerances)
    events = FlowKneadingEvent[]
    raw = Int8[]
    seen = Ref(0)
    last_event = Ref(-Inf)
    status = Ref(:maximum_time)
    detail = Ref("")
    terminal_time = Ref(0.0)

    function record!(state, tangent, time)
        time - last_event[] > o.minimum_event_separation || return false
        _accept_event(ctx, problem.capture, state, time) || return false
        last_event[] = time
        seen[] += 1
        seen[] > o.transient_events || return false
        projected = _word_tangent(ctx, state, tangent, time, o.tolerances)
        component, observable_norm = _flow_observable(problem.observable, state, projected, ctx.p, time)
        tolerance = o.sign_atol + o.sign_rtol * norm(projected) * observable_norm
        sign = abs(component) <= tolerance ? Int8(0) : Int8(component > 0 ? 1 : -1)
        push!(events, FlowKneadingEvent(time, collect(Float64, state), projected,
            state[problem.capture.index], _event_rate(ctx, problem.capture, state, time), component, sign))
        if sign == 0
            status[] = :ambiguous_sign
            return true
        end
        push!(raw, sign)
        if length(raw) == o.word_length + 1
            status[] = :complete
            return true
        end
        return false
    end

    include_initial = o.include_initial_event === true ||
        (o.include_initial_event === :auto && hasproperty(seed, :rho) && hasproperty(seed, :event_index))
    if include_initial
        abs(_event_value(ctx, problem.capture, u, 0.0)) <= o.initial_event_tolerance &&
            _accept_event(ctx, problem.capture, u, 0.0) ||
            throw(ArgumentError("the initial state is not an accepted capture event"))
        record!(u, v, 0.0)
    end

    if status[] ∉ (:complete, :ambiguous_sign)
        try
            if o.integration === :rk4
                t = 0.0
                h_previous = _event_value(ctx, problem.capture, u, t)
                iterations = 0
                while t < o.maximum_time
                    iterations += 1
                    if iterations > o.maxiters
                        status[] = :maximum_iterations
                        break
                    end
                    step = min(o.dt, o.maximum_time - t)
                    next_u, next_v = _word_rk4_step(ctx, u, v, t, step)
                    terminal_time[] = t + step
                    if !all(isfinite, next_u) || !all(isfinite, next_v)
                        status[] = :nonfinite_state
                        break
                    end
                    if maximum(abs, next_u) > o.max_state
                        status[] = :state_limit
                        break
                    end
                    next_v = _word_tangent(ctx, next_u, next_v, t + step, o.tolerances)
                    h_next = _event_value(ctx, problem.capture, next_u, t + step)
                    crossed = problem.capture isa LocalMaximum ? h_previous > 0 && h_next <= 0 :
                        h_previous < 0 && h_next >= 0
                    if crossed
                        fraction = clamp(h_previous / (h_previous - h_next), 0.0, 1.0)
                        event_time = t + fraction * step
                        if event_time > o.minimum_event_separation
                            stopped = record!(u + fraction * (next_u - u), v + fraction * (next_v - v), event_time)
                            if stopped
                                terminal_time[] = event_time
                                break
                            end
                        end
                    end
                    u, v = next_u, next_v
                    h_previous = h_next
                    t += step
                end
            else
                function augmented!(dz, z, p, t)
                    state = @view z[1:n]
                    tangent = @view z[n + 1:2n]
                    dz[1:n] .= ctx.f(state, p, t)
                    dz[n + 1:2n] .= ctx.J(state, p, t) * tangent
                    return nothing
                end
                condition = (z, t, integrator) -> _event_value(ctx, problem.capture, @view(z[1:n]), t)
                function capture!(integrator)
                    t = integrator.t
                    t > o.minimum_event_separation || return nothing
                    stopped = record!(@view(integrator.u[1:n]), @view(integrator.u[n + 1:2n]), t)
                    stopped && SciMLBase.terminate!(integrator)
                    return nothing
                end
                event_callback = problem.capture isa LocalMaximum ?
                    SciMLBase.ContinuousCallback(condition, nothing, capture!; save_positions = (false, false), abstol = 1e-12, reltol = 0.0) :
                    SciMLBase.ContinuousCallback(condition, capture!, nothing; save_positions = (false, false), abstol = 1e-12, reltol = 0.0)
                function project!(integrator)
                    terminal_time[] = integrator.t
                    if !all(isfinite, integrator.u)
                        status[] = :nonfinite_state
                        SciMLBase.terminate!(integrator)
                    elseif maximum(abs, @view(integrator.u[1:n])) > o.max_state
                        status[] = :state_limit
                        SciMLBase.terminate!(integrator)
                    else
                        tangent = @view integrator.u[n + 1:2n]
                        project_flow_normal!(tangent, ctx.f(@view(integrator.u[1:n]), ctx.p, integrator.t), o.tolerances)
                        SciMLBase.derivative_discontinuity!(integrator, true)
                    end
                    return nothing
                end
                step_callback = SciMLBase.DiscreteCallback((z, t, integrator) -> true, project!; save_positions = (false, false))
                ode = SciMLBase.ODEProblem(augmented!, vcat(u, v), (0.0, o.maximum_time), ctx.p)
                solution = SciMLBase.solve(ode, Tsit5();
                    callback = SciMLBase.CallbackSet(event_callback, step_callback),
                    abstol = o.abstol, reltol = o.reltol, dtmax = o.dtmax,
                    maxiters = o.maxiters, save_everystep = false, save_start = false,
                    save_end = true, dense = false,
                )
                !isempty(solution.t) && (terminal_time[] = last(solution.t))
                if !SciMLBase.successful_retcode(solution) && status[] === :maximum_time
                    status[] = :integration_failure
                    detail[] = string(solution.retcode)
                end
            end
        catch error
            error isa DomainError || rethrow()
            status[] = :numerical_failure
            detail[] = sprint(showerror, error)
        end
    end

    transitions = Int8[raw[i] * raw[i + 1] for i in 1:max(0, length(raw) - 1)]
    complete = status[] === :complete
    raw_code = complete ? _word_code(raw) : BigInt(-1)
    transition_code = complete ? _word_code(transitions) : BigInt(-1)
    times = Float64[events[i + 1].time - events[i].time for i in 1:max(0, length(events) - 1)]
    metadata = (;
        parameters = deepcopy(ctx.p), capture = problem.capture,
        observable = problem.observable, options = o, accepted_events = seen[],
        initial_event_included = include_initial, terminal_time = terminal_time[],
        tangent_gauge = :accepted_step_flow_normal, detail = detail[],
    )
    return FlowKneadingResult(raw, transitions, raw_code, transition_code,
        length(raw), length(transitions), events, times, status[], complete, seed, metadata)
end
