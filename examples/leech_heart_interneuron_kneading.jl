using DynamicalSystemsBase
using ForwardDiff
using Kneading
using Kneading.Diagrams: ParameterPlane
using LinearAlgebra
using StaticArrays

boltzmann(a, b, V) = 1 / (1 + exp(a * (b + V)))

function leech_heart_interneuron(u, p, t)
    V, h, m = u
    shift, current = p
    sodium = 200 * boltzmann(-150, 0.0305, V)^3 * h * (V - 0.045)
    potassium = 30 * m^2 * (V + 0.07)
    leak = 8 * (V + 0.046)
    return SVector(2 * (current - sodium - potassium - leak),
        24.69 * (boltzmann(500, 0.0333, V) - h),
        4 * (boltzmann(-83, 0.018 + shift, V) - m))
end

leech(shift, current) = CoupledODEs(leech_heart_interneuron, SVector(-0.03, 0.5, 0.1), [shift, current])

rest_state(V, shift) = [V, boltzmann(500, 0.0333, V), boltzmann(-83, 0.018 + shift, V)]

holding_current(V, shift) = -leech_heart_interneuron(SVector(rest_state(V, shift)...), (shift, 0.0), 0.0)[1] / 2

function equilibrium_voltages(shift, current; voltages = range(-0.08, 0.05; length = 2601))
    g(V) = holding_current(V, shift) - current
    roots = Float64[]
    for (a, b) in zip(voltages[1:end-1], voltages[2:end])
        g(a) * g(b) < 0 || continue
        for _ in 1:60
            c = (a + b) / 2
            g(a) * g(c) <= 0 ? (b = c) : (a = c)
        end
        push!(roots, (a + b) / 2)
    end
    return roots
end

const settings = (
    capture = LocalMinimum(3),
    word_length = 63,
    maximum_time = 80.0,
    include_initial_event = false,
    integration = :rk4,
    dt = 2e-4,
)

function leech_plane(resolution)
    return ParameterPlane(
        collect(range(parse.(Float64, split(get(ENV, "LEECH_SHIFT", "-0.028,-0.008"), ','))...; length = resolution)),
        collect(range(parse.(Float64, split(get(ENV, "LEECH_CURRENT", "-0.03,0.05"), ','))...; length = resolution));
        xname = "shift", yname = "current",
    )
end

function saddle_focus_diagram(plane; backend = nothing, progress = nothing)
    initializer = SaddleFocusInitializer(
        equilibrium_guess = rest_state(-0.0272, plane.x[1]),
        critical_kind = :minimum,
        initial_event_index = 8,
        maximum_event_index = 30,
        max_time = 60.0,
    )
    return scan_flow_kneading(leech, plane; initializer, settings..., store_results = true, backend, progress)
end

function real_saddle_results(plane; backend = nothing)
    cells = CartesianIndex{2}[]
    problems = FlowKneadingProblem[]
    for (j, shift) in enumerate(plane.x), (i, current) in enumerate(plane.y)
        voltages = equilibrium_voltages(shift, current)
        length(voltages) == 3 || continue
        initializer = RealSaddleInitializer(
            equilibrium_guess = rest_state(voltages[2], shift),
            unstable_reference = [1.0, 0.0, 0.0],
        )
        push!(cells, CartesianIndex(i, j))
        push!(problems, FlowKneadingProblem(leech(shift, current); initializer, settings...))
    end
    results = isempty(problems) ? [] :
        isnothing(backend) ? flow_kneading(identity.(problems)) : flow_kneading(identity.(problems); backend)
    return cells, results
end

function write_orbit_table(path, plane, rows)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, join(("shift", "current", "seed", "status", "events", "longest_return", "last_event_time", "transition_word"), '\t'))
        for (cell, seed, result) in rows
            isnothing(result) && continue
            times = result.return_times
            println(io, join((plane.x[cell[2]], plane.y[cell[1]], seed, result.status, length(result.events),
                isempty(times) ? NaN : maximum(times), isempty(result.events) ? NaN : result.events[end].time,
                join(map(s -> s > 0 ? '1' : '0', result.transition_word))), '\t'))
        end
    end
end

function equilibrium_type(V, shift, current)
    jacobian = ForwardDiff.jacobian(u -> leech_heart_interneuron(u, (shift, current), 0.0), SVector(rest_state(V, shift)...))
    eigenvalues = eigvals(Matrix(jacobian))
    unstable = count(λ -> real(λ) > 0, eigenvalues)
    return string(unstable, any(λ -> abs(imag(λ)) > 0, eigenvalues) ? "f" : "")
end

function write_equilibrium_table(path, plane; refinement = 4)
    shifts = range(first(plane.x), last(plane.x); length = refinement * length(plane.x))
    currents = range(first(plane.y), last(plane.y); length = refinement * length(plane.y))
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "shift\tcurrent\tvoltages\ttypes")
        for shift in shifts, current in currents
            voltages = equilibrium_voltages(shift, current)
            types = [equilibrium_type(V, shift, current) for V in voltages]
            println(io, shift, '\t', current, '\t', join(voltages, ','), '\t', join(types, ','))
        end
    end
end

function main(; backend = nothing)
    resolution = parse(Int, get(ENV, "LEECH_RESOLUTION", "32"))
    output = get(ENV, "LEECH_OUTPUT", joinpath(@__DIR__, "..", "output", "leech-heart-interneuron"))
    plane = leech_plane(resolution)
    start = time()
    initialized = Ref(NaN)
    progress = (i, j, result) -> isnan(initialized[]) && !isnothing(result) && (initialized[] = time() - start)
    focus_time = @elapsed focus = saddle_focus_diagram(plane; backend, progress)
    saddle_time = @elapsed cells, saddle = real_saddle_results(plane; backend)
    write_flow_scan(joinpath(output, "saddle-focus.tsv"), focus)
    write_equilibrium_table(joinpath(output, "equilibria.tsv"), plane)
    rows = vcat(vec([(c, "saddle_focus", focus.results[c]) for c in CartesianIndices(focus.results)]),
        [(c, "real_saddle", r) for (c, r) in zip(cells, saddle)])
    write_orbit_table(joinpath(output, "orbits.tsv"), plane, rows)
    println("Saddle-focus scan: ", round(focus_time; digits = 1), " s (first word after ",
        round(initialized[]; digits = 1), " s), complete words ",
        count(==(:complete), focus.statuses), "/", length(focus.statuses),
        ", initialization failures ", count(==(:initialization_failed), focus.statuses))
    println("Real-saddle orbits: ", round(saddle_time; digits = 1), " s, ", length(saddle), " points")
    println("Saved ", abspath(output))
    return focus, cells, saddle
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
