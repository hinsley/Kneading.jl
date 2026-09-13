"""
Initialization from the leading stable direction of a real saddle.
"""
module RealSaddleInitialization

import DynamicalSystemsBase
import SciMLBase

using LinearAlgebra: dot, eigen, norm

export InvalidRealSaddleInitialState,
    RealSaddleSeed,
    RealSaddleTolerances,
    init_real_saddle

"""
    RealSaddleTolerances(
        equilibrium_residual,
        eigenvalue_real_part,
        eigenvalue_imaginary_part,
        eigenvalue_separation,
        stable_real_part_gap,
        orientation,
        flow_norm,
        projected_norm,
    )

Store the tolerances for equilibrium, spectrum, orientation, launch-flow, and
tangent-projection validation.
"""
struct RealSaddleTolerances{T<:Real}
    equilibrium_residual::T
    eigenvalue_real_part::T
    eigenvalue_imaginary_part::T
    eigenvalue_separation::T
    stable_real_part_gap::T
    orientation::T
    flow_norm::T
    projected_norm::T

    function RealSaddleTolerances{T}(values::Vararg{T, 8}) where {T<:Real}
        all(isfinite, values) || throw(ArgumentError(
            "real-saddle tolerances must be finite",
        ))
        all(value -> value >= 0, values) || throw(ArgumentError(
            "real-saddle tolerances must be nonnegative",
        ))
        return new{T}(values...)
    end
end

function RealSaddleTolerances(values::Vararg{Real, 8})
    promoted = promote(map(float, values)...)
    return RealSaddleTolerances{typeof(promoted[1])}(promoted...)
end

"""
    InvalidRealSaddleInitialState(initial_state, stage, detail)

Report that an initial state did not produce a valid real-saddle seed. `stage`
identifies equilibrium solving, branch validation, saddle-spectrum validation,
or leading-stable-direction selection.
"""
struct InvalidRealSaddleInitialState{U} <: Exception
    initial_state::U
    stage::Symbol
    detail::String
end

function Base.showerror(io::IO, error::InvalidRealSaddleInitialState)
    print(
        io,
        "invalid initial state for real-saddle initialization at ",
        error.stage,
        ": ",
        error.detail,
    )
end

struct _RealSaddleProblem{S, U, T, B, C, RU, RS}
    system::S
    initial_state::U
    launch_distance::T
    max_root_iterations::Int
    equilibrium_branch::B
    equilibrium_branch_check::C
    unstable_reference::RU
    stable_reference::RS
    unstable_branch::Int8
end

"""
Store the real-saddle equilibrium, spectrum, launch state, initial tangent, and
all validation diagnostics returned by [`init_real_saddle`](@ref).
"""
struct RealSaddleSeed{
    UI,
    UE,
    J,
    UL,
    Q,
    EU,
    ES,
    LU,
    LS,
    E,
    T,
    B,
    RU,
    RS,
    C,
}
    initial_state::UI
    equilibrium::UE
    jacobian::J
    u0::UL
    Q0::Q
    unstable_direction::EU
    leading_stable_direction::ES
    unstable_eigenvalue::LU
    leading_stable_eigenvalue::LS
    eigenvalues::E
    root_iterations::Int
    max_root_iterations::Int
    equilibrium_residual::T
    unstable_orientation_dot::T
    stable_orientation_dot::T
    launch_flow_norm::T
    projected_tangent_norm::T
    launch_distance::T
    equilibrium_branch::B
    unstable_reference::RU
    stable_reference::RS
    unstable_branch::Int8
    tolerances::C
end

function _problem(
    system,
    initial_state::AbstractVector{<:Real},
    launch_distance::Real;
    max_root_iterations::Integer,
    equilibrium_branch,
    equilibrium_branch_check,
    unstable_reference::AbstractVector{<:Real},
    stable_reference::AbstractVector{<:Real},
    unstable_branch::Integer,
)
    system isa DynamicalSystemsBase.CoupledODEs || throw(ArgumentError(
        "system must be a CoupledODEs",
    ))
    max_root_iterations > 0 || throw(ArgumentError(
        "max_root_iterations must be positive",
    ))
    unstable_branch in (-1, 1) || throw(ArgumentError(
        "unstable_branch must be 1 or -1",
    ))
    isnothing(equilibrium_branch) && throw(ArgumentError(
        "equilibrium_branch must identify the selected branch",
    ))
    isfinite(launch_distance) && launch_distance > 0 || throw(ArgumentError(
        "launch_distance must be finite and positive",
    ))

    dimension = length(DynamicalSystemsBase.current_state(system))
    length(initial_state) == dimension || throw(DimensionMismatch(
        "initial_state must have one entry for each system dimension",
    ))
    length(unstable_reference) == dimension || throw(DimensionMismatch(
        "unstable_reference must have one entry for each system dimension",
    ))
    length(stable_reference) == dimension || throw(DimensionMismatch(
        "stable_reference must have one entry for each system dimension",
    ))
    all(isfinite, initial_state) || throw(ArgumentError(
        "initial_state must be finite",
    ))
    all(isfinite, unstable_reference) || throw(ArgumentError(
        "unstable_reference must be finite",
    ))
    all(isfinite, stable_reference) || throw(ArgumentError(
        "stable_reference must be finite",
    ))
    norm(unstable_reference) > 0 || throw(ArgumentError(
        "unstable_reference must have nonzero norm",
    ))
    norm(stable_reference) > 0 || throw(ArgumentError(
        "stable_reference must have nonzero norm",
    ))
    applicable(equilibrium_branch_check, initial_state) || throw(ArgumentError(
        "equilibrium_branch_check must accept an equilibrium state",
    ))

    state = collect(float.(initial_state))
    unstable = collect(float.(unstable_reference))
    stable = collect(float.(stable_reference))
    distance = float(launch_distance)
    return _RealSaddleProblem(
        system,
        state,
        distance,
        Int(max_root_iterations),
        equilibrium_branch,
        equilibrium_branch_check,
        unstable,
        stable,
        Int8(unstable_branch),
    )
end

function _system_flow(problem::_RealSaddleProblem, state)
    system = problem.system
    rule = DynamicalSystemsBase.dynamic_rule(system)
    parameters = DynamicalSystemsBase.current_parameters(system)
    time = DynamicalSystemsBase.current_time(system)
    if SciMLBase.isinplace(system)
        flow = similar(state)
        rule(flow, state, parameters, time)
    else
        flow = collect(rule(state, parameters, time))
    end
    length(flow) == length(state) || throw(DimensionMismatch(
        "the vector field must return one derivative for each state entry",
    ))
    all(isfinite, flow) || throw(DomainError(
        flow,
        "the vector field returned a non-finite derivative",
    ))
    return flow
end

function _system_jacobian(problem::_RealSaddleProblem, state)
    system = problem.system
    rule = DynamicalSystemsBase.jacobian(system)
    parameters = DynamicalSystemsBase.current_parameters(system)
    time = DynamicalSystemsBase.current_time(system)
    dimension = length(state)
    if SciMLBase.isinplace(system)
        jacobian = Matrix{eltype(state)}(undef, dimension, dimension)
        rule(jacobian, state, parameters, time)
    else
        jacobian = Matrix(rule(state, parameters, time))
    end
    size(jacobian) == (dimension, dimension) || throw(DimensionMismatch(
        "the state Jacobian must be square with the system dimension",
    ))
    all(isfinite, jacobian) || throw(DomainError(
        jacobian,
        "the state Jacobian contains a non-finite value",
    ))
    return jacobian
end

function _root_failure(problem, detail)
    return InvalidRealSaddleInitialState(
        copy(problem.initial_state),
        :root_solve,
        detail,
    )
end

function _solve_equilibrium(
    problem::_RealSaddleProblem,
    tolerances::RealSaddleTolerances,
)
    equilibrium = copy(problem.initial_state)
    for iteration in 0:problem.max_root_iterations
        residual_vector = try
            _system_flow(problem, equilibrium)
        catch error
            throw(_root_failure(problem, sprint(showerror, error)))
        end
        residual = norm(residual_vector)
        residual <= tolerances.equilibrium_residual &&
            return equilibrium, residual, iteration
        iteration == problem.max_root_iterations && break

        jacobian = try
            _system_jacobian(problem, equilibrium)
        catch error
            throw(_root_failure(problem, sprint(showerror, error)))
        end
        update = try
            jacobian \ residual_vector
        catch error
            throw(_root_failure(problem, sprint(showerror, error)))
        end
        all(isfinite, update) || throw(_root_failure(
            problem,
            "the Newton update contains a non-finite value",
        ))
        equilibrium .-= update
        all(isfinite, equilibrium) || throw(_root_failure(
            problem,
            "the Newton state contains a non-finite value",
        ))
    end
    throw(_root_failure(
        problem,
        "the equilibrium solve did not converge within $(problem.max_root_iterations) iterations",
    ))
end

function _validate_equilibrium_branch(problem, equilibrium)
    accepted = try
        problem.equilibrium_branch_check(copy(equilibrium))
    catch error
        throw(InvalidRealSaddleInitialState(
            copy(problem.initial_state),
            :equilibrium_branch,
            sprint(showerror, error),
        ))
    end
    accepted isa Bool || throw(ArgumentError(
        "equilibrium_branch_check must return Bool",
    ))
    accepted || throw(InvalidRealSaddleInitialState(
        copy(problem.initial_state),
        :equilibrium_branch,
        "the equilibrium does not belong to the selected branch",
    ))
    return nothing
end

function _real_eigenvector(vector, tolerance, stage, problem)
    imaginary_norm = norm(imag.(vector))
    imaginary_norm <= tolerance || throw(InvalidRealSaddleInitialState(
        copy(problem.initial_state),
        stage,
        "the selected real eigenvalue does not have a numerically real eigendirection",
    ))
    direction = collect(real.(vector))
    direction_norm = norm(direction)
    direction_norm > 0 || throw(InvalidRealSaddleInitialState(
        copy(problem.initial_state),
        stage,
        "the selected eigendirection has zero norm",
    ))
    direction ./= direction_norm
    return direction
end

function _validate_real_saddle_spectrum(
    problem::_RealSaddleProblem,
    jacobian,
    tolerances::RealSaddleTolerances,
)
    decomposition = eigen(jacobian)
    eigenvalues = collect(decomposition.values)
    real_parts = real.(eigenvalues)
    neutral = findall(value -> abs(value) <= tolerances.eigenvalue_real_part, real_parts)
    isempty(neutral) || throw(InvalidRealSaddleInitialState(
        copy(problem.initial_state),
        :real_saddle_spectrum,
        "the spectrum contains a neutral eigenvalue",
    ))

    unstable_indices = findall(value -> value > tolerances.eigenvalue_real_part, real_parts)
    length(unstable_indices) == 1 || throw(InvalidRealSaddleInitialState(
        copy(problem.initial_state),
        :real_saddle_spectrum,
        "the spectrum must contain exactly one unstable eigenvalue",
    ))
    unstable_index = only(unstable_indices)
    unstable_value = eigenvalues[unstable_index]
    abs(imag(unstable_value)) <= tolerances.eigenvalue_imaginary_part ||
        throw(InvalidRealSaddleInitialState(
            copy(problem.initial_state),
            :real_saddle_spectrum,
            "the unstable eigenvalue must be real",
        ))
    for index in eachindex(eigenvalues)
        index == unstable_index && continue
        abs(unstable_value - eigenvalues[index]) >=
            tolerances.eigenvalue_separation ||
            throw(InvalidRealSaddleInitialState(
                copy(problem.initial_state),
                :real_saddle_spectrum,
                "the unstable eigenvalue must be simple",
            ))
    end

    stable_indices = filter(!=(unstable_index), eachindex(eigenvalues))
    isempty(stable_indices) && throw(InvalidRealSaddleInitialState(
        copy(problem.initial_state),
        :real_saddle_spectrum,
        "the spectrum must contain at least one stable eigenvalue",
    ))
    all(index -> real_parts[index] < -tolerances.eigenvalue_real_part, stable_indices) ||
        throw(InvalidRealSaddleInitialState(
            copy(problem.initial_state),
            :real_saddle_spectrum,
            "every non-unstable eigenvalue must have negative real part",
        ))
    unstable_direction = _real_eigenvector(
        decomposition.vectors[:, unstable_index],
        tolerances.eigenvalue_imaginary_part,
        :real_saddle_spectrum,
        problem,
    )
    return (
        eigenvalues = eigenvalues,
        eigenvectors = decomposition.vectors,
        unstable_eigenvalue = real(unstable_value),
        unstable_direction = unstable_direction,
        stable_indices = stable_indices,
    )
end

function _select_leading_stable_direction(
    problem::_RealSaddleProblem,
    spectrum,
    tolerances::RealSaddleTolerances,
)
    real_parts = real.(spectrum.eigenvalues)
    candidate_position = argmax(real_parts[spectrum.stable_indices])
    candidate_index = spectrum.stable_indices[candidate_position]
    candidate_value = spectrum.eigenvalues[candidate_index]
    abs(imag(candidate_value)) <= tolerances.eigenvalue_imaginary_part ||
        throw(InvalidRealSaddleInitialState(
            copy(problem.initial_state),
            :leading_stable_direction,
            "the leading stable eigenvalue must be real",
        ))

    for index in eachindex(spectrum.eigenvalues)
        index == candidate_index && continue
        abs(candidate_value - spectrum.eigenvalues[index]) >=
            tolerances.eigenvalue_separation ||
            throw(InvalidRealSaddleInitialState(
                copy(problem.initial_state),
                :leading_stable_direction,
                "the leading stable eigenvalue must be simple",
            ))
    end
    for index in spectrum.stable_indices
        index == candidate_index && continue
        real(candidate_value) - real_parts[index] >=
            tolerances.stable_real_part_gap ||
            throw(InvalidRealSaddleInitialState(
                copy(problem.initial_state),
                :leading_stable_direction,
                "the leading stable real part is not sufficiently separated",
            ))
    end

    direction = _real_eigenvector(
        spectrum.eigenvectors[:, candidate_index],
        tolerances.eigenvalue_imaginary_part,
        :leading_stable_direction,
        problem,
    )
    return real(candidate_value), direction
end

function _orient_direction(direction, reference, tolerance)
    reference_norm = norm(reference)
    reference_norm > 0 || throw(ArgumentError(
        "orientation reference vectors must have nonzero norm",
    ))
    unit_reference = reference ./ reference_norm
    orientation_dot = dot(direction, unit_reference)
    abs(orientation_dot) > tolerance || throw(DomainError(
        orientation_dot,
        "the eigendirection orientation is ambiguous",
    ))
    oriented = orientation_dot > 0 ? copy(direction) : -direction
    return oriented, abs(orientation_dot), unit_reference
end

function _initial_tangent(leading_stable_direction, launch_flow, tolerances)
    flow_norm = norm(launch_flow)
    flow_norm > tolerances.flow_norm || throw(DomainError(
        flow_norm,
        "the launch flow norm is too small",
    ))
    tangent = leading_stable_direction .-
        (dot(leading_stable_direction, launch_flow) / dot(launch_flow, launch_flow)) .*
        launch_flow
    projected_norm = norm(tangent)
    projected_norm > tolerances.projected_norm || throw(DomainError(
        projected_norm,
        "the projected leading stable direction is too small",
    ))
    tangent ./= projected_norm
    return reshape(tangent, :, 1), flow_norm, projected_norm
end

"""
    init_real_saddle(
        system,
        initial_state,
        launch_distance,
        tolerances;
        max_root_iterations,
        equilibrium_branch,
        equilibrium_branch_check,
        unstable_reference,
        stable_reference,
        unstable_branch=1,
    )

Find and validate a real saddle from `initial_state`. Select deterministic
unstable and leading stable eigendirections, launch on the selected unstable
ray, and return the displaced state and one-column flow-normal tangent seed.
"""
function init_real_saddle(
    system,
    initial_state::AbstractVector{<:Real},
    launch_distance::Real,
    tolerances::RealSaddleTolerances;
    max_root_iterations::Integer,
    equilibrium_branch,
    equilibrium_branch_check,
    unstable_reference::AbstractVector{<:Real},
    stable_reference::AbstractVector{<:Real},
    unstable_branch::Integer = 1,
)
    problem = _problem(
        system,
        initial_state,
        launch_distance;
        max_root_iterations,
        equilibrium_branch,
        equilibrium_branch_check,
        unstable_reference,
        stable_reference,
        unstable_branch,
    )
    equilibrium, equilibrium_residual, root_iterations =
        _solve_equilibrium(problem, tolerances)
    _validate_equilibrium_branch(problem, equilibrium)
    jacobian = try
        _system_jacobian(problem, equilibrium)
    catch error
        throw(InvalidRealSaddleInitialState(
            copy(problem.initial_state),
            :real_saddle_spectrum,
            sprint(showerror, error),
        ))
    end
    spectrum = _validate_real_saddle_spectrum(problem, jacobian, tolerances)
    leading_stable_eigenvalue, leading_stable_direction =
        _select_leading_stable_direction(problem, spectrum, tolerances)

    oriented_unstable, unstable_orientation_dot, normalized_unstable_reference =
        _orient_direction(
            spectrum.unstable_direction,
            problem.unstable_reference,
            tolerances.orientation,
        )
    oriented_stable, stable_orientation_dot, normalized_stable_reference =
        _orient_direction(
            leading_stable_direction,
            problem.stable_reference,
            tolerances.orientation,
        )
    selected_unstable = problem.unstable_branch .* oriented_unstable
    u0 = equilibrium .+ problem.launch_distance .* selected_unstable
    launch_flow = _system_flow(problem, u0)
    Q0, launch_flow_norm, projected_tangent_norm = _initial_tangent(
        oriented_stable,
        launch_flow,
        tolerances,
    )
    diagnostic_values = promote(
        equilibrium_residual,
        unstable_orientation_dot,
        stable_orientation_dot,
        launch_flow_norm,
        projected_tangent_norm,
        problem.launch_distance,
    )

    return RealSaddleSeed(
        copy(problem.initial_state),
        equilibrium,
        jacobian,
        u0,
        Q0,
        selected_unstable,
        oriented_stable,
        spectrum.unstable_eigenvalue,
        leading_stable_eigenvalue,
        spectrum.eigenvalues,
        root_iterations,
        problem.max_root_iterations,
        diagnostic_values...,
        problem.equilibrium_branch,
        normalized_unstable_reference,
        normalized_stable_reference,
        problem.unstable_branch,
        tolerances,
    )
end

end
