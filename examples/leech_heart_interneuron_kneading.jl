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
    word_length = 100,
    maximum_time = 200.0,
    include_initial_event = false,
    integration = :rk4,
    dt = 2e-4,
)

function leech_plane(resolution)
    return ParameterPlane(
        collect(range(parse.(Float64, split(get(ENV, "LEECH_SHIFT", "-0.032,-0.008"), ','))...; length = resolution)),
        collect(range(parse.(Float64, split(get(ENV, "LEECH_CURRENT", "-0.03,0.035"), ','))...; length = resolution));
        xname = "shift", yname = "current",
    )
end

const focus_options = (critical_kind = :minimum, initial_event_index = 8, maximum_event_index = 30, max_time = 60.0)

function saddle_focus_diagram(plane; backend = nothing)
    initializer = SaddleFocusInitializer(; equilibrium_guess = rest_state(-0.0272, plane.x[1]), focus_options...)
    return scan_flow_kneading(leech, plane; initializer, settings..., store_results = true, backend)
end

function fresh_saddle_focus_results(plane, focus; backend = nothing)
    cells = [c for c in CartesianIndices(focus.statuses) if focus.statuses[c] == :initialization_failed &&
        last(equilibrium_types(plane.x[c[2]], plane.y[c[1]])) == "2f"]
    seeds = Vector{Any}(nothing, length(cells))
    Threads.@threads :dynamic for k in eachindex(cells)
        shift, current = plane.x[cells[k][2]], plane.y[cells[k][1]]
        seeds[k] = try
            init_saddle_focus(leech(shift, current); capture = settings.capture, focus_options...,
                initial_event_index = 16, event_index_fallback = false,
                equilibrium_guess = rest_state(last(equilibrium_voltages(shift, current)), shift))
        catch exception
            exception isa Union{SaddleFocusInitializationError,DomainError} || rethrow()
            nothing
        end
    end
    found = findall(!isnothing, seeds)
    problems = [FlowKneadingProblem(leech(plane.x[cells[k][2]], plane.y[cells[k][1]]); initializer = seeds[k],
        settings...) for k in found]
    return cells[found], isempty(problems) ? [] : flow_kneading(identity.(problems); backend)
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
    return cells, isempty(problems) ? [] : flow_kneading(identity.(problems); backend)
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

equilibrium_types(shift, current) =
    [equilibrium_type(V, shift, current) for V in equilibrium_voltages(shift, current)]

function write_equilibrium_table(path, plane; refinement = 2)
    shifts = range(first(plane.x), last(plane.x); length = refinement * length(plane.x))
    currents = range(first(plane.y), last(plane.y); length = refinement * length(plane.y))
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "shift\tcurrent\tvoltages\ttypes")
        for shift in shifts, current in currents
            println(io, shift, '\t', current, '\t', join(equilibrium_voltages(shift, current), ','), '\t',
                join(equilibrium_types(shift, current), ','))
        end
    end
end

function main(; backend = nothing)
    resolution = parse(Int, get(ENV, "LEECH_RESOLUTION", "32"))
    output = get(ENV, "LEECH_OUTPUT", joinpath(@__DIR__, "..", "output", "leech-heart-interneuron"))
    plane = leech_plane(resolution)
    focus_time = @elapsed focus = saddle_focus_diagram(plane; backend)
    write_flow_scan(joinpath(output, "saddle-focus.tsv"), focus)
    fresh_time = @elapsed fresh_cells, fresh = fresh_saddle_focus_results(plane, focus; backend)
    saddle_time = @elapsed cells, saddle = real_saddle_results(plane; backend)
    results = copy(focus.results)
    foreach((c, r) -> results[c] = r, fresh_cells, fresh)
    rows = vcat(vec([(c, "saddle_focus", results[c]) for c in CartesianIndices(results)]),
        [(c, "real_saddle", r) for (c, r) in zip(cells, saddle)])
    write_orbit_table(joinpath(output, "orbits.tsv"), plane, rows)
    write_equilibrium_table(joinpath(output, "equilibria.tsv"), plane)
    statuses = [isnothing(r) ? :no_seed : r.status for r in results]
    println("Saddle-focus scan: ", round(focus_time; digits = 1), " s; fresh starts: ", round(fresh_time; digits = 1),
        " s for ", length(fresh), " points; real saddle: ", round(saddle_time; digits = 1), " s for ",
        length(saddle), " points")
    println("Saddle-focus statuses: ", sort!(collect(pairs(Dict(s => count(==(s), statuses) for s in unique(statuses)))); by = last, rev = true))
    println("Saved ", abspath(output))
    return focus, results, cells, saddle
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
