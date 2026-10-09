include(joinpath(@__DIR__, "leech_heart_interneuron_kneading.jl"))

const trace_points = [
    ("A", -0.0285, 0.005, :attractor),
    ("B", -0.021, 0.0, :attractor),
    ("C", -0.016, 0.0, :attractor),
    ("D", -0.011, 0.015, :attractor),
    ("E", -0.020, -0.020, :separatrix),
]

function trace(shift, current, kind; duration = 6.0, transient = 40.0)
    if kind == :separatrix
        voltages = equilibrium_voltages(shift, current)
        initializer = RealSaddleInitializer(equilibrium_guess = rest_state(voltages[2], shift),
            unstable_reference = [1.0, 0.0, 0.0])
        result = flow_kneading(FlowKneadingProblem(leech(shift, current); initializer, settings...))
        start, transient = result.initialization.u0, 0.0
    else
        start = [-0.04, 0.5, 0.1]
    end
    system = CoupledODEs(leech_heart_interneuron, SVector(start...), [shift, current])
    path, times = trajectory(system, duration; Ttr = transient, Δt = 1e-3)
    return times .- first(times), path[:, 1]
end

function write_traces(path)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "label\tshift\tcurrent\ttime\tvoltage")
        for (label, shift, current, kind) in trace_points
            times, voltages = trace(shift, current, kind)
            for (t, V) in zip(times, voltages)
                println(io, label, '\t', shift, '\t', current, '\t', t, '\t', V)
            end
        end
    end
    return path
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) &&
    println("Saved ", write_traces(isempty(ARGS) ? joinpath(@__DIR__, "..", "output", "leech-heart-interneuron", "traces.tsv") : ARGS[1]))
