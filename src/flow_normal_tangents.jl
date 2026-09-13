"""
Flow-normal tangent integration for continuous dynamical systems.
"""
module FlowNormalTangents

import DynamicalSystemsBase
import SciMLBase

using LinearAlgebra: dot, norm

export FlowNormalTolerances,
    init_flow_normal,
    project_flow_normal!,
    project_integrator_tangent!,
    projection_callback,
    solve_flow_normal!

"""
    FlowNormalTolerances(flow_norm, projected_norm, unit_norm)

Store the tolerances for flow-normal tangent projection.

`flow_norm` rejects a small flow before division. `projected_norm` rejects a
small transverse component before normalization. `unit_norm` bounds the final
unit-length error.
"""
struct FlowNormalTolerances{T<:Real}
    flow_norm::T
    projected_norm::T
    unit_norm::T

    function FlowNormalTolerances{T}(
        flow_norm::T,
        projected_norm::T,
        unit_norm::T,
    ) where {T<:Real}
        values = (flow_norm, projected_norm, unit_norm)
        all(isfinite, values) || throw(ArgumentError(
            "flow-normal tolerances must be finite",
        ))
        all(value -> value >= 0, values) || throw(ArgumentError(
            "flow-normal tolerances must be nonnegative",
        ))
        return new{T}(values...)
    end
end

function FlowNormalTolerances(
    flow_norm::Real,
    projected_norm::Real,
    unit_norm::Real,
)
    values = promote(float(flow_norm), float(projected_norm), float(unit_norm))
    return FlowNormalTolerances{typeof(values[1])}(values...)
end

"""
    project_flow_normal!(tangent, flow, tolerances)

Remove the component of `tangent` parallel to `flow`, normalize the remaining
component in place, and return the updated tangent.
"""
function project_flow_normal!(
    tangent,
    flow,
    tolerances::FlowNormalTolerances,
)
    length(tangent) == length(flow) || throw(DimensionMismatch(
        "the tangent and flow must have the same length",
    ))

    flow_length = norm(flow)
    flow_length > tolerances.flow_norm || throw(DomainError(
        flow_length,
        "the flow norm must be greater than $(tolerances.flow_norm)",
    ))

    parallel_scale = dot(tangent, flow) / (flow_length * flow_length)
    tangent .-= parallel_scale .* flow

    projected_length = norm(tangent)
    projected_length > tolerances.projected_norm || throw(DomainError(
        projected_length,
        "the projected tangent norm must be greater than $(tolerances.projected_norm)",
    ))
    tangent ./= projected_length

    unit_error = abs(norm(tangent) - one(projected_length))
    unit_error <= tolerances.unit_norm || throw(DomainError(
        unit_error,
        "the normalized tangent error must not exceed $(tolerances.unit_norm)",
    ))
    return tangent
end

"""
    project_integrator_tangent!(integrator, tolerances)

Read the current derivative from `SciMLBase.get_du(integrator)`. Project and
normalize the single tangent in the augmented integrator state without another
vector-field evaluation. Store and return the updated tangent.
"""
function project_integrator_tangent!(
    integrator,
    tolerances::FlowNormalTolerances,
)
    state = integrator.u
    ndims(state) == 2 && size(state, 2) == 2 || throw(ArgumentError(
        "the tangent integrator must contain one deviation vector",
    ))

    derivative = SciMLBase.get_du(integrator)
    size(derivative) == size(state) || throw(DimensionMismatch(
        "the integrator state and derivative must have the same size",
    ))
    flow = @view derivative[:, 1]

    if Base.ismutable(state)
        tangent = @view state[:, 2]
        project_flow_normal!(tangent, flow, tolerances)
        SciMLBase.set_u!(integrator, state)
        return tangent
    end

    tangent = collect(@view state[:, 2])
    project_flow_normal!(tangent, flow, tolerances)
    updated_state = typeof(state)(hcat(state[:, 1], tangent))
    SciMLBase.set_u!(integrator, updated_state)
    return updated_state[:, 2]
end

"""
    projection_callback(tolerances)

Build a callback that projects the tangent after every accepted solver step.
"""
function projection_callback(tolerances::FlowNormalTolerances)
    condition = (u, t, integrator) -> true
    affect! = integrator -> project_integrator_tangent!(integrator, tolerances)
    return SciMLBase.DiscreteCallback(
        condition,
        affect!;
        save_positions = (false, false),
    )
end

"""
    init_flow_normal(system, Q0, tolerances;
                     u0=nothing, J=nothing, J0=nothing, callback=nothing)

Build a `TangentDynamicalSystem` with one deviation vector. Preserve callbacks
already attached to `system`, add `callback` when supplied, and add the
accepted-step projection callback. `J` and `J0` are passed to
`TangentDynamicalSystem` when supplied.
"""
function init_flow_normal(
    system,
    Q0,
    tolerances::FlowNormalTolerances;
    u0 = nothing,
    J = nothing,
    J0 = nothing,
    callback = nothing,
)
    system isa DynamicalSystemsBase.CoupledODEs || throw(ArgumentError(
        "system must be a CoupledODEs",
    ))
    Q0 isa AbstractMatrix || throw(ArgumentError(
        "Q0 must be a matrix with one deviation vector",
    ))
    size(Q0, 2) == 1 || throw(ArgumentError(
        "Q0 must contain one deviation vector",
    ))
    size(Q0, 1) == length(DynamicalSystemsBase.current_state(system)) ||
        throw(DimensionMismatch(
            "Q0 must have one row for each system dimension",
        ))
    all(isfinite, Q0) || throw(ArgumentError("Q0 must be finite"))
    if !isnothing(u0)
        length(u0) == size(Q0, 1) || throw(DimensionMismatch(
            "u0 must have one entry for each system dimension",
        ))
        all(isfinite, u0) || throw(ArgumentError("u0 must be finite"))
    end

    existing_callback = haskey(system.diffeq, :callback) ?
        system.diffeq.callback : nothing
    callback_values = filter(
        !isnothing,
        (existing_callback, callback, projection_callback(tolerances)),
    )
    callbacks = SciMLBase.CallbackSet(callback_values...)
    configured_system = DynamicalSystemsBase.CoupledODEs(
        system,
        (callback = callbacks,),
    )

    tangent_options = (Q0 = Q0,)
    isnothing(u0) || (tangent_options = merge(tangent_options, (u0 = u0,)))
    isnothing(J) || (tangent_options = merge(tangent_options, (J = J,)))
    isnothing(J0) || (tangent_options = merge(tangent_options, (J0 = J0,)))
    return DynamicalSystemsBase.TangentDynamicalSystem(
        configured_system;
        tangent_options...,
    )
end

"""
    solve_flow_normal!(system, final_time)

Advance a continuous `TangentDynamicalSystem` to `final_time` and return its
DifferentialEquations.jl solution. Leave the system at the requested time.
"""
function solve_flow_normal!(system, final_time)
    system isa DynamicalSystemsBase.TangentDynamicalSystem || throw(ArgumentError(
        "system must be a TangentDynamicalSystem",
    ))
    system.ds isa DynamicalSystemsBase.CoupledODEs || throw(ArgumentError(
        "system must use continuous dynamics",
    ))
    final_time isa Real && isfinite(final_time) || throw(ArgumentError(
        "final_time must be a finite real number",
    ))

    current_time = DynamicalSystemsBase.current_time(system)
    final_time >= current_time || throw(ArgumentError(
        "final_time must not be less than the current time",
    ))
    if final_time > current_time
        SciMLBase.step!(system, final_time - current_time, true)
    end

    DynamicalSystemsBase.successful_step(system) || error(
        "the tangent integration did not complete successfully",
    )
    DynamicalSystemsBase.current_time(system) == final_time || error(
        "the tangent integration stopped before final_time",
    )
    return system.ds.integ.sol
end

end
