"""
    SaddleFocusInitializer(; kwargs...)

Configure saddle-focus critical-point initialization for a flow-kneading scan.
Keywords are passed to [`init_saddle_focus`](@ref).
"""
struct SaddleFocusInitializer{K<:NamedTuple}
    options::K
end

SaddleFocusInitializer(; kwargs...) = SaddleFocusInitializer((; kwargs...))

"""
The critical state `u0`, unit flow-normal tangent `Q0`, launch anchor, and
diagnostics returned by [`init_saddle_focus`](@ref). Parameter values are copied
so changing the source system does not change the recorded initialization.
"""
struct SaddleFocusSeed{P,C,D,K}
    u0::Vector{Float64}
    Q0::Matrix{Float64}
    equilibrium::Vector{Float64}
    seed_direction::Vector{Float64}
    rho::Float64
    event_index::Int
    parameters::P
    capture::C
    critical_kind::Symbol
    diagnostics::D
    configuration::K
end

"""
An initialization failure with its stage and numerical diagnostics.
"""
struct SaddleFocusInitializationError{D} <: Exception
    stage::Symbol
    message::String
    diagnostics::D
end

function Base.showerror(io::IO, error::SaddleFocusInitializationError)
    print(io, "Saddle-focus initialization failed (", error.stage, "): ", error.message)
end

function _sf_equilibrium(ctx, guess; tolerance, max_iterations)
    u = collect(Float64, guess)
    for iteration in 0:max_iterations
        residual = collect(ctx.f(u, ctx.p, 0.0))
        residual_norm = norm(residual)
        isfinite(residual_norm) || break
        residual_norm <= tolerance && return u
        iteration == max_iterations && break
        delta = try
            ctx.J(u, ctx.p, 0.0) \ residual
        catch error
            error isa LinearAlgebra.SingularException || rethrow()
            break
        end
        accepted = false
        for scale in (1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125)
            candidate = u - scale * delta
            if norm(ctx.f(candidate, ctx.p, 0.0)) < residual_norm
                u = candidate
                accepted = true
                break
            end
        end
        accepted || break
    end
    throw(SaddleFocusInitializationError(:equilibrium, "equilibrium correction did not converge", (; guess=collect(guess), state=u)))
end

function _sf_seed_ray(ctx, equilibrium, capture; eigenvalue_tolerance=1e-10)
    J = Matrix(ctx.J(equilibrium, ctx.p, 0.0))
    eig = eigen(J)
    unstable = findall(x -> real(x) > eigenvalue_tolerance, eig.values)
    candidates = filter(i -> imag(eig.values[i]) > eigenvalue_tolerance, unstable)
    if length(unstable) != 2 || length(candidates) != 1 || any(x -> abs(real(x)) <= eigenvalue_tolerance, eig.values)
        throw(SaddleFocusInitializationError(:eigenspace, "expected a hyperbolic saddle-focus with exactly two unstable complex directions", (; eigenvalues=eig.values)))
    end
    length(eig.values) >= 3 || throw(ArgumentError("a saddle-focus requires at least three state variables"))
    index = only(candidates)
    eigenvector = eig.vectors[:, index]
    basis = hcat(real.(eigenvector), imag.(eigenvector))
    Dh = vec(J[capture.index, :])
    row = vec(transpose(Dh) * basis)
    norm(row) > eigenvalue_tolerance || throw(SaddleFocusInitializationError(:seed_ray, "the event tangent plane does not select a seed line", (; event_gradient=Dh)))
    direction = basis * [-row[2], row[1]]
    direction ./= norm(direction)
    acceleration = dot(Dh, J * direction)
    abs(acceleration) > eigenvalue_tolerance || throw(SaddleFocusInitializationError(:seed_ray, "the extremum orientation is degenerate", (; acceleration)))
    desired_sign = capture isa LocalMaximum ? -1.0 : 1.0
    direction .*= desired_sign * sign(acceleration)
    return (; equilibrium, direction, eigenvalue=ComplexF64(eig.values[index]))
end

struct _SaddleFocusEvent
    time::Float64
    state::Vector{Float64}
    tangent::Vector{Float64}
    second_tangent::Vector{Float64}
    derivative::Float64
    second_derivative::Float64
    transversality::Float64
end

function _sf_event(ctx, capture, u, v, s, t, second_order, event_tolerance)
    f = collect(ctx.f(u, ctx.p, t))
    J = Matrix(ctx.J(u, ctx.p, t))
    Dh = vec(J[capture.index, :])
    denominator = dot(Dh, f)
    if !isfinite(denominator) || abs(denominator) <= event_tolerance
        return nothing
    end
    time_derivative = -dot(Dh, v) / denominator
    tangent = v + time_derivative * f
    if second_order
        Hh = ForwardDiff.hessian(x -> ctx.f(x, ctx.p, t)[capture.index], u)
        fixed_second = s + 2time_derivative * (J * v) + time_derivative^2 * (J * f)
        second_time_derivative = -(dot(Dh, fixed_second) + dot(tangent, Hh * tangent)) / denominator
        second_tangent = fixed_second + second_time_derivative * f
    else
        second_tangent = fill(NaN, length(u))
    end
    return _SaddleFocusEvent(t, u, tangent, second_tangent, tangent[capture.index], second_tangent[capture.index], denominator)
end

function _sf_events(ctx, capture, ray, rho, count, options; second_order=false)
    radius = exp(rho)
    seed = ray.equilibrium + radius * ray.direction
    sensitivity = radius * ray.direction
    n = length(seed)
    augmented_state = second_order ? vcat(seed, sensitivity, sensitivity) : vcat(seed, sensitivity)
    events = _SaddleFocusEvent[]
    failure = Ref("")
    guard = isnothing(options.launch_guard_time) ? pi / abs(imag(ray.eigenvalue)) : options.launch_guard_time
    function augmented!(dy, y, p, t)
        u = @view y[1:n]
        v = @view y[n+1:2n]
        J = ctx.J(u, ctx.p, t)
        dy[1:n] .= ctx.f(u, ctx.p, t)
        dy[n+1:2n] .= J * v
        if second_order
            s = @view y[2n+1:3n]
            dy[2n+1:3n] .= J * s
            for i in 1:n
                H = ForwardDiff.hessian(z -> ctx.f(z, ctx.p, t)[i], u)
                dy[2n+i] += dot(v, H * v)
            end
        end
        return nothing
    end
    condition(y, t, integrator) = _event_value(ctx, capture, @view(y[1:n]), t)
    function affect!(integrator)
        integrator.t >= guard || return
        u = collect(@view integrator.u[1:n])
        _accept_event(ctx, capture, u, integrator.t) || return
        v = collect(@view integrator.u[n+1:2n])
        s = second_order ? collect(@view integrator.u[2n+1:3n]) : Float64[]
        event = _sf_event(ctx, capture, u, v, s, Float64(integrator.t), second_order, options.event_tolerance)
        if isnothing(event)
            failure[] = "an event is insufficiently transverse"
            SciMLBase.terminate!(integrator)
            return
        end
        push!(events, event)
        length(events) >= count && SciMLBase.terminate!(integrator)
    end
    callback = capture isa LocalMaximum ?
        SciMLBase.ContinuousCallback(condition, nothing, affect!; save_positions=(false, false), abstol=1e-14, reltol=0.0) :
        SciMLBase.ContinuousCallback(condition, affect!, nothing; save_positions=(false, false), abstol=1e-14, reltol=0.0)
    function unbounded(y, t, integrator)
        return any(x -> !isfinite(x), y) || maximum(abs, @view(y[1:n])) > options.max_state
    end
    function stop_unbounded!(integrator)
        failure[] = "the seeded trajectory or sensitivity became unbounded"
        SciMLBase.terminate!(integrator)
    end
    bound_callback = SciMLBase.DiscreteCallback(unbounded, stop_unbounded!; save_positions=(false, false))
    problem = SciMLBase.ODEProblem(augmented!, augmented_state, (0.0, options.max_time))
    solution = SciMLBase.solve(problem, options.alg;
        callback=SciMLBase.CallbackSet(callback, bound_callback),
        abstol=options.abstol, reltol=options.reltol, dtmax=options.dtmax,
        maxiters=options.maxiters, save_start=false, save_end=false, save_everystep=false,
    )
    if length(events) < count && isempty(failure[])
        failure[] = "not enough accepted extrema before the integration limit ($(solution.retcode))"
    end
    return events, failure[]
end

function _sf_residual(ctx, capture, ray, rho, M, options; second_order=false)
    invalid(message, events=_SaddleFocusEvent[]) = (; valid=false, rho=Float64(rho), residual=NaN, slope=NaN, events, message)
    isfinite(rho) && exp(rho) > options.minimum_radius || return invalid("launch radius is below the resolvable minimum")
    events, failure = _sf_events(ctx, capture, ray, rho, M + 1, options; second_order)
    isempty(failure) || return invalid(failure, events)
    current, following = events[M], events[M+1]
    abs(current.derivative) > options.denominator_tolerance || return invalid("return-coordinate derivative denominator is too small", events)
    residual = following.derivative / current.derivative
    slope = second_order ? (following.second_derivative - residual * current.second_derivative) / current.derivative : NaN
    all(isfinite, current.tangent) && isfinite(residual) || return invalid("nonfinite event sensitivity", events)
    return (; valid=true, rho=Float64(rho), residual, slope, events, message="ok")
end

function _sf_slope(evaluate, result, options)
    options.newton_derivative == :second_order_sensitivity && return result.slope
    h = options.finite_difference_step
    plus, minus = evaluate(result.rho + h), evaluate(result.rho - h)
    if !(plus.valid && minus.valid)
        return NaN
    end
    M = length(result.events) - 1
    if sign(plus.events[M].derivative) != sign(minus.events[M].derivative)
        return NaN
    end
    return (plus.residual - minus.residual) / (2h)
end

function _sf_newton(evaluate, guess, options)
    result = evaluate(guess)
    for iteration in 0:options.max_newton_iterations
        result.valid || return (; result, iterations=iteration, slope=NaN, converged=false)
        slope = _sf_slope(evaluate, result, options)
        if abs(result.residual) <= options.criticality_tolerance
            return (; result, iterations=iteration, slope, converged=isfinite(slope))
        end
        if iteration == options.max_newton_iterations || !isfinite(slope) || abs(slope) <= eps(Float64)
            return (; result, iterations=iteration, slope, converged=false)
        end
        step = clamp(result.residual / slope, -options.max_newton_step, options.max_newton_step)
        next_rho = result.rho - step
        log(options.minimum_radius) < next_rho <= options.rho_range[2] || return (; result, iterations=iteration, slope, converged=false)
        result = evaluate(next_rho)
    end
    error("unreachable Newton iteration state")
end

function _sf_bisect(evaluate, left, right, M, options)
    left.valid && right.valid || return nothing
    sign(left.events[M].derivative) == sign(right.events[M].derivative) || return nothing
    sign(left.residual) != sign(right.residual) || return nothing
    for _ in 1:options.max_bisection_iterations
        middle = evaluate((left.rho + right.rho) / 2)
        middle.valid || return nothing
        sign(middle.events[M].derivative) == sign(left.events[M].derivative) || return nothing
        abs(middle.residual) <= options.criticality_tolerance && return middle
        (middle.rho == left.rho || middle.rho == right.rho) && return nothing
        if sign(middle.residual) == sign(left.residual)
            left = middle
        else
            right = middle
        end
    end
    return nothing
end

function _sf_kind_matches(candidate, M, kind)
    candidate.converged || return false
    curvature = candidate.slope / candidate.result.events[M].derivative
    return isfinite(curvature) && (kind == :minimum ? curvature > 0 : curvature < 0)
end

function _sf_target_distance(candidate, target, M, capture)
    isnothing(target) && return 0.0
    state = candidate.result.events[M].state
    distance = target isa Number ? abs(state[capture.index] - target) : norm(state - target)
    return distance / max(1.0, norm(target))
end

function _sf_find_root(ctx, capture, ray, guess, M, options; target=nothing, scan=true)
    second_order = options.newton_derivative == :second_order_sensitivity
    evaluate(rho) = _sf_residual(ctx, capture, ray, rho, M, options; second_order)
    direct = _sf_newton(evaluate, guess, options)
    direct_matches_kind = _sf_kind_matches(direct, M, options.critical_kind)
    target_mismatch = direct_matches_kind && _sf_target_distance(direct, target, M, capture) > options.branch_tolerance
    if direct_matches_kind && !target_mismatch
        return direct
    end
    if !scan && target_mismatch
        throw(SaddleFocusInitializationError(:branch, "corrected critical point is too far from the target branch; reduce the parameter step or review branch_tolerance", (; event_index=M, target_distance=_sf_target_distance(direct, target, M, capture), branch_tolerance=options.branch_tolerance)))
    end
    scan || throw(SaddleFocusInitializationError(:root, "critical-point correction failed; reduce the parameter step or adjust the seed guess", (; rho=guess, event_index=M, message=direct.result.message)))
    roots = Any[]
    rhos = range(options.rho_range[1], options.rho_range[2]; length=options.rho_samples)
    previous = evaluate(first(rhos))
    for rho in Iterators.drop(rhos, 1)
        current = evaluate(rho)
        if current.valid && abs(current.residual) <= options.criticality_tolerance
            candidate = _sf_newton(evaluate, current.rho, options)
            if _sf_kind_matches(candidate, M, options.critical_kind)
                if _sf_target_distance(candidate, target, M, capture) <= options.branch_tolerance
                    push!(roots, candidate)
                else
                    target_mismatch = true
                end
            end
        end
        root = _sf_bisect(evaluate, previous, current, M, options)
        if !isnothing(root)
            candidate = _sf_newton(evaluate, root.rho, options)
            if _sf_kind_matches(candidate, M, options.critical_kind)
                if _sf_target_distance(candidate, target, M, capture) <= options.branch_tolerance
                    push!(roots, candidate)
                else
                    target_mismatch = true
                end
            end
        end
        previous = current
    end
    if isempty(roots)
        target_mismatch && throw(SaddleFocusInitializationError(:branch, "no corrected critical point is close enough to the target branch; reduce the parameter step or review branch_tolerance", (; event_index=M, branch_tolerance=options.branch_tolerance)))
        throw(SaddleFocusInitializationError(:root, "no $(options.critical_kind) critical point was found; adjust initial_radius, rho_range, or initial_event_index", (; event_index=M, rho_range=options.rho_range, message=direct.result.message)))
    end
    function score(candidate)
        point = candidate.result.events[M]
        if !isnothing(target)
            return target isa Number ? abs(point.state[capture.index] - target) : norm(point.state - target)
        end
        next_value = candidate.result.events[M+1].state[capture.index]
        return options.critical_kind == :minimum ? next_value : -next_value
    end
    return roots[argmin(score.(roots))]
end

function _sf_unit_tangent(ctx, event, options)
    tangent = copy(event.tangent)
    flow = collect(ctx.f(event.state, ctx.p, event.time))
    norm(flow) > options.event_tolerance || throw(SaddleFocusInitializationError(:tangent, "the flow vanishes at the selected critical point", (; state=event.state)))
    tangent .-= (dot(tangent, flow) / dot(flow, flow)) .* flow
    tangent_norm = norm(tangent)
    tangent_norm > options.denominator_tolerance || throw(SaddleFocusInitializationError(:tangent, "the projected tangent vanishes at the selected critical point", (; state=event.state)))
    tangent ./= tangent_norm
    return tangent
end

function _sf_previous_option(previous, name, fallback)
    isnothing(previous) && return fallback
    previous isa SaddleFocusSeed || throw(ArgumentError("previous must be a SaddleFocusSeed"))
    return get(previous.configuration, name, fallback)
end

"""
    init_saddle_focus(system; equilibrium_guess, capture=LocalMinimum(2), kwargs...)

Find a critical point of an extremum return map using the unstable spiral at a
saddle-focus. Return its full state `u0`, one-column unit flow-normal tangent
`Q0`, and launch data in a [`SaddleFocusSeed`](@ref).

`critical_kind` is `:minimum` or `:maximum` of the return map, independently of
the captured extremum type. `initial_radius` and `initial_event_index` supply
the first root guess. If correction fails, scan `rho_range` for candidate roots.
Use `critical_target` to select a particular critical state or return coordinate.
When a target or `previous` result is supplied, require its distance from the
corrected state, divided by `max(1, norm(target))`, to be at most
`branch_tolerance=0.25`. This is a heuristic guard against branch jumps, not a
proof of branch identity. Reduce the parameter step if it rejects a correction.

`newton_derivative=:finite_difference` uses central differences for the Newton
slope; `:second_order_sensitivity` integrates second-order sensitivities and
corrects both orders for event time. Both integrate first-order sensitivities
for the return-map slope itself.

By default, check each candidate against a smaller launch radius and an extra
event. Retain the smaller event index when the critical state and oriented
tangent agree within `state_tolerance` and `tangent_tolerance`; otherwise
advance to the finer candidate and repeat. `refine=false` performs an unchecked fixed-index
initialization. Pass `previous=seed` after changing parameters to continue the
same critical-point branch. Solver settings and refinement mode are inherited
from `previous` unless explicitly overridden. The vector field must be autonomous and support
ForwardDiff for second-order sensitivities.

Control ODE accuracy with `abstol`, `reltol`, and `dtmax`. The root solver uses
`finite_difference_step`, `max_newton_step`, and `max_newton_iterations`;
`rho_range` and `rho_samples` control its fallback search. `event_tolerance`
guards event transversality and `denominator_tolerance` guards the coordinate
derivative and projected tangent. The equilibrium corrector uses
`equilibrium_tolerance` and `max_equilibrium_iterations`.

`max_time`, `max_state`, `minimum_radius`, and `maximum_event_index` bound the
search. An unsuccessful root solve or refinement throws
[`SaddleFocusInitializationError`](@ref), carrying a failure stage and numerical
diagnostics. Inspect `seed.diagnostics` for the accepted residual, refinement
errors, derivative method, and target distance.
"""
function init_saddle_focus(system;
    previous=nothing,
    equilibrium_guess=isnothing(previous) ? nothing : previous.equilibrium,
    capture=isnothing(previous) ? LocalMinimum(2) : previous.capture,
    critical_kind=isnothing(previous) ? :minimum : previous.critical_kind,
    newton_derivative=_sf_previous_option(previous, :newton_derivative, :finite_difference),
    initial_radius=_sf_previous_option(previous, :initial_radius, 1e-3),
    initial_rho=nothing,
    initial_event_index=isnothing(previous) ? 4 : previous.event_index,
    maximum_event_index=_sf_previous_option(previous, :maximum_event_index, 12),
    criticality_tolerance=_sf_previous_option(previous, :criticality_tolerance, 1e-8),
    state_tolerance=_sf_previous_option(previous, :state_tolerance, 1e-6),
    tangent_tolerance=_sf_previous_option(previous, :tangent_tolerance, 1e-6),
    equilibrium_tolerance=_sf_previous_option(previous, :equilibrium_tolerance, 1e-12),
    max_equilibrium_iterations=_sf_previous_option(previous, :max_equilibrium_iterations, 30),
    event_tolerance=_sf_previous_option(previous, :event_tolerance, 1e-10),
    denominator_tolerance=_sf_previous_option(previous, :denominator_tolerance, 1e-9),
    eigenvalue_tolerance=_sf_previous_option(previous, :eigenvalue_tolerance, 1e-10),
    finite_difference_step=_sf_previous_option(previous, :finite_difference_step, 1e-4),
    max_newton_step=_sf_previous_option(previous, :max_newton_step, 0.5),
    max_newton_iterations=_sf_previous_option(previous, :max_newton_iterations, 16),
    max_bisection_iterations=_sf_previous_option(previous, :max_bisection_iterations, 60),
    rho_range=_sf_previous_option(previous, :rho_range, (-24.0, -1.0)),
    rho_samples=_sf_previous_option(previous, :rho_samples, 65),
    critical_target=isnothing(previous) ? nothing : previous.u0,
    branch_tolerance=_sf_previous_option(previous, :branch_tolerance, 0.25),
    minimum_radius=_sf_previous_option(previous, :minimum_radius, 1e-14),
    launch_guard_time=_sf_previous_option(previous, :launch_guard_time, nothing),
    max_time=_sf_previous_option(previous, :max_time, 650.0),
    max_state=_sf_previous_option(previous, :max_state, 1e6),
    abstol=_sf_previous_option(previous, :abstol, 1e-11),
    reltol=_sf_previous_option(previous, :reltol, 1e-11),
    dtmax=_sf_previous_option(previous, :dtmax, 0.25),
    maxiters=_sf_previous_option(previous, :maxiters, 10_000_000),
    alg=_sf_previous_option(previous, :alg, Tsit5()),
    refine=_sf_previous_option(previous, :refine, true),
)
    isnothing(previous) || previous isa SaddleFocusSeed || throw(ArgumentError("previous must be a SaddleFocusSeed"))
    isnothing(equilibrium_guess) && throw(ArgumentError("equilibrium_guess is required for the first initialization"))
    capture isa Union{LocalMaximum,LocalMinimum} || throw(ArgumentError("capture must be LocalMaximum or LocalMinimum"))
    critical_kind in (:minimum, :maximum) || throw(ArgumentError("critical_kind must be :minimum or :maximum"))
    newton_derivative in (:finite_difference, :second_order_sensitivity) || throw(ArgumentError("newton_derivative must be :finite_difference or :second_order_sensitivity"))
    1 <= initial_event_index <= maximum_event_index || throw(ArgumentError("require 1 <= initial_event_index <= maximum_event_index"))
    initial_event_index isa Integer && maximum_event_index isa Integer || throw(ArgumentError("event indices must be integers"))
    rho_samples isa Integer && rho_samples >= 2 || throw(ArgumentError("rho_samples must be an integer of at least two"))
    for (name, value) in pairs((; max_equilibrium_iterations, max_newton_iterations))
        value isa Integer && value >= 0 || throw(ArgumentError("$name must be a nonnegative integer"))
    end
    for (name, value) in pairs((; max_bisection_iterations, maxiters))
        value isa Integer && value > 0 || throw(ArgumentError("$name must be a positive integer"))
    end
    if !isnothing(launch_guard_time)
        isfinite(launch_guard_time) && launch_guard_time >= 0 || throw(ArgumentError("launch_guard_time must be nonnegative and finite"))
    end
    isnothing(initial_rho) || isfinite(initial_rho) || throw(ArgumentError("initial_rho must be finite"))
    isfinite(initial_radius) && initial_radius > 0 || throw(ArgumentError("initial_radius must be positive and finite"))
    all(isfinite, rho_range) && rho_range[1] < rho_range[2] || throw(ArgumentError("rho_range must have finite increasing endpoints"))
    for (name, value) in pairs((; criticality_tolerance, state_tolerance, tangent_tolerance, equilibrium_tolerance, event_tolerance, denominator_tolerance, eigenvalue_tolerance, finite_difference_step, max_newton_step, branch_tolerance, minimum_radius, max_time, max_state, abstol, reltol, dtmax))
        isfinite(value) && value > 0 || throw(ArgumentError("$name must be positive and finite"))
    end
    ctx = _flow_context(system)
    1 <= capture.index <= length(ctx.u0) || throw(ArgumentError("capture index is outside the state"))
    length(equilibrium_guess) == length(ctx.u0) || throw(DimensionMismatch("equilibrium_guess has the wrong state dimension"))
    if !isnothing(critical_target)
        if critical_target isa Number
            isfinite(critical_target) || throw(ArgumentError("critical_target must be finite"))
        else
            length(critical_target) == length(ctx.u0) || throw(DimensionMismatch("critical_target has the wrong state dimension"))
            all(isfinite, critical_target) || throw(ArgumentError("critical_target must be finite"))
        end
    end
    equilibrium = _sf_equilibrium(ctx, equilibrium_guess; tolerance=equilibrium_tolerance, max_iterations=max_equilibrium_iterations)
    ray = _sf_seed_ray(ctx, equilibrium, capture; eigenvalue_tolerance)
    options = (; critical_kind, newton_derivative, criticality_tolerance, state_tolerance, tangent_tolerance,
        event_tolerance, denominator_tolerance, finite_difference_step, max_newton_step,
        max_newton_iterations, max_bisection_iterations, rho_range, rho_samples, branch_tolerance, minimum_radius,
        launch_guard_time, max_time, max_state, abstol, reltol, dtmax, maxiters, alg)
    rho = isnothing(initial_rho) ? (isnothing(previous) ? log(initial_radius) : previous.rho) : Float64(initial_rho)
    M = Int(initial_event_index)
    corrected = _sf_find_root(ctx, capture, ray, rho, M, options; target=critical_target)
    event = corrected.result.events[M]
    tangent = _sf_unit_tangent(ctx, event, options)
    if !isnothing(previous) && dot(tangent, vec(previous.Q0)) < 0
        tangent .*= -1
    end
    state_error = NaN
    tangent_error = NaN
    refinement_steps = 0
    converged = false
    total_iterations = corrected.iterations
    while refine && M < maximum_event_index
        next_M = M + 1
        guess = corrected.result.rho - 2pi * real(ray.eigenvalue) / abs(imag(ray.eigenvalue))
        next_corrected = _sf_find_root(ctx, capture, ray, guess, next_M, options; target=event.state, scan=false)
        next_event = next_corrected.result.events[next_M]
        next_tangent = _sf_unit_tangent(ctx, next_event, options)
        dot(tangent, next_tangent) < 0 && (next_tangent .*= -1)
        state_error = norm(next_event.state - event.state)
        tangent_error = norm(next_tangent - tangent)
        refinement_steps += 1
        total_iterations += next_corrected.iterations
        if state_error <= state_tolerance && tangent_error <= tangent_tolerance
            converged = true
            break
        end
        M, corrected, event, tangent = next_M, next_corrected, next_event, next_tangent
    end
    diagnostics = (; converged, root_converged=true, refinement_checked=refine,
        residual=corrected.result.residual,
        curvature=corrected.slope / event.derivative,
        derivative_method=newton_derivative, state_error, tangent_error,
        target_distance=_sf_target_distance(corrected, critical_target, M, capture),
        branch_tolerance,
        iterations=total_iterations, refinement_steps,
        refinement_event_index=refine ? M + (converged ? 1 : 0) : nothing,
        event_time=event.time,
        next_state=copy(corrected.result.events[M+1].state),
        event_derivatives=(event.derivative, corrected.result.events[M+1].derivative),
        event_transversality=(event.transversality, corrected.result.events[M+1].transversality),
        eigenvalue=ray.eigenvalue)
    if refine && !converged
        throw(SaddleFocusInitializationError(:refinement, "state and tangent did not converge before maximum_event_index=$maximum_event_index", diagnostics))
    end
    configuration = (; newton_derivative, initial_radius, maximum_event_index, criticality_tolerance,
        state_tolerance, tangent_tolerance, equilibrium_tolerance, max_equilibrium_iterations,
        event_tolerance, denominator_tolerance, eigenvalue_tolerance, finite_difference_step,
        max_newton_step, max_newton_iterations, max_bisection_iterations, rho_range, rho_samples, branch_tolerance,
        minimum_radius, launch_guard_time, max_time, max_state, abstol, reltol, dtmax, maxiters, alg, refine)
    return SaddleFocusSeed(copy(event.state), reshape(tangent, :, 1), copy(equilibrium),
        copy(ray.direction), corrected.result.rho, M, deepcopy(ctx.p), capture,
        critical_kind, diagnostics, configuration)
end

(initializer::SaddleFocusInitializer)(system; kwargs...) = init_saddle_focus(system; initializer.options..., kwargs...)
