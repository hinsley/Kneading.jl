"""
Critical-orbit initialization, event capture, and symbolic parameter scans for ODEs.
"""
module FlowKneading

import DynamicalSystemsBase
import ForwardDiff
import SciMLBase
using OrdinaryDiffEqTsit5: Tsit5
using LinearAlgebra
using ..FlowNormalTangents
using ..RealSaddleInitialization
using ..Diagrams: ParameterPlane

export LocalMaximum, LocalMinimum, CoordinateComponent,
    RealSaddleInitializer, SaddleFocusInitializer, SaddleFocusSeed,
    SaddleFocusInitializationError, init_saddle_focus,
    FlowKneadingProblem, FlowKneadingEvent, FlowKneadingResult, flow_kneading,
    FlowKneadingDiagram, scan_flow_kneading, write_flow_scan

abstract type ExtremumCapture end

"""
    LocalMaximum(index; accept=(u, p, t) -> true)

Capture strict local maxima of state coordinate `index`. The optional predicate
filters events before consecutive returns are paired.
"""
struct LocalMaximum{F} <: ExtremumCapture
    index::Int
    accept::F
    function LocalMaximum(index::Integer; accept = (u, p, t) -> true)
        index > 0 || throw(ArgumentError("the event coordinate must be positive"))
        return new{typeof(accept)}(Int(index), accept)
    end
end

"""
    LocalMinimum(index; accept=(u, p, t) -> true)

Capture strict local minima of state coordinate `index`. The optional predicate
filters events before consecutive returns are paired.
"""
struct LocalMinimum{F} <: ExtremumCapture
    index::Int
    accept::F
    function LocalMinimum(index::Integer; accept = (u, p, t) -> true)
        index > 0 || throw(ArgumentError("the event coordinate must be positive"))
        return new{typeof(accept)}(Int(index), accept)
    end
end

"""
    CoordinateComponent(index)

Encode the sign of coordinate `index` of the unit flow-normal tangent.
"""
struct CoordinateComponent
    index::Int
    function CoordinateComponent(index::Integer)
        index > 0 || throw(ArgumentError("the observable coordinate must be positive"))
        return new(Int(index))
    end
end

function _flow_context(system)
    system isa DynamicalSystemsBase.CoupledODEs ||
        throw(ArgumentError("flow kneading requires a CoupledODEs system"))
    rule = DynamicalSystemsBase.dynamic_rule(system)
    p = deepcopy(DynamicalSystemsBase.current_parameters(system))
    u0 = collect(Float64, DynamicalSystemsBase.current_state(system))
    all(isfinite, u0) || throw(ArgumentError("the system state must be finite"))
    f = if SciMLBase.isinplace(system)
        (u, p, t) -> begin
            du = similar(u)
            rule(du, u, p, t)
            du
        end
    else
        (u, p, t) -> collect(rule(u, p, t))
    end
    J = (u, p, t) -> ForwardDiff.jacobian(z -> f(z, p, t), u)
    return (; f, p, u0, J)
end

_event_value(ctx, capture::ExtremumCapture, u, t) = ctx.f(u, ctx.p, t)[capture.index]

function _event_rate(ctx, capture::ExtremumCapture, u, t)
    gradient = ForwardDiff.gradient(z -> _event_value(ctx, capture, z, t), u)
    return dot(gradient, ctx.f(u, ctx.p, t))
end

function _accept_event(ctx, capture::ExtremumCapture, u, t)
    rate = _event_rate(ctx, capture, u, t)
    direction = capture isa LocalMaximum ? rate < 0 : rate > 0
    return direction && capture.accept(u, ctx.p, t)
end

include("saddle_focus_initialization.jl")
include("flow_words.jl")
include("flow_kneading_scan.jl")

end
