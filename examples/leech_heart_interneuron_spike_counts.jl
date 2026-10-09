include(joinpath(@__DIR__, "leech_heart_interneuron_kneading.jl"))

const count_settings = (bursts = 14, gap = 0.4, threshold = 0.0, dt = 1e-4, maximum_time = 60.0)

function rk4_step(u, p, h)
    f(x) = leech_heart_interneuron(x, p, 0.0)
    k1 = f(u)
    k2 = f(u + h / 2 * k1)
    k3 = f(u + h / 2 * k2)
    k4 = f(u + h * k3)
    return u + h / 6 * (k1 + 2k2 + 2k3 + k4)
end

function spike_counts(u0, shift, current; transient = 0.0, bursts = count_settings.bursts, gap = count_settings.gap,
    threshold = count_settings.threshold, dt = count_settings.dt, maximum_time = count_settings.maximum_time)
    p = SVector(shift, current)
    u = SVector{3,Float64}(u0...)
    rate = leech_heart_interneuron(u, p, 0.0)[1]
    counts = Int[]
    spikes = 0
    last_spike = -Inf
    started = transient == 0
    steps = round(Int, (transient + maximum_time) / dt)
    for step in 1:steps
        t = step * dt
        next = rk4_step(u, p, dt)
        next_rate = leech_heart_interneuron(next, p, 0.0)[1]
        if rate > 0 && next_rate <= 0 && u[1] > threshold
            if t - last_spike > gap
                started && spikes > 0 && push!(counts, spikes)
                length(counts) == bursts && return counts, :complete
                started = t > transient
                spikes = 0
            end
            spikes += 1
            last_spike = t
        end
        u, rate = next, next_rate
        if step % 10_000 == 0 && t - last_spike > 1 && norm(leech_heart_interneuron(u, p, 0.0)) < 1e-9
            started && spikes > 0 && push!(counts, spikes)
            return counts, :rest
        end
    end
    started && spikes > 0 && transient + maximum_time - last_spike > gap && push!(counts, spikes)
    return counts, isempty(counts) && spikes == 0 ? :no_spikes : :maximum_time
end

function read_seeds(path)
    lines = readlines(path)
    header = split(first(lines), '\t')
    column(name) = findfirst(==(name), header)
    return [(parse(Float64, f[column("shift")]), parse(Float64, f[column("current")]), f[column("seed")],
        parse.(Float64, f[column.(["u1", "u2", "u3"])]))
        for f in (split(line, '\t') for line in lines[2:end] if !isempty(line))]
end

function count_orbits(seeds)
    focus = filter(seed -> seed[3] == "saddle_focus", seeds)
    rows = Vector{Any}(undef, length(focus))
    Threads.@threads :dynamic for k in eachindex(focus)
        shift, current, _, u0 = focus[k]
        counts, status = spike_counts(u0, shift, current)
        rows[k] = (shift, current, status, counts)
    end
    return rows
end

function write_counts(path, rows)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "shift\tcurrent\tstatus\tcounts")
        for (shift, current, status, counts) in rows
            println(io, shift, '\t', current, '\t', status, '\t', join(counts, ','))
        end
    end
    return path
end

function main()
    directory = get(ENV, "LEECH_OUTPUT", joinpath(@__DIR__, "..", "output", "leech-heart-interneuron"))
    elapsed = @elapsed rows = count_orbits(read_seeds(joinpath(directory, "orbits.tsv")))
    println("Counted ", length(rows), " orbits in ", round(elapsed; digits = 1), " s")
    println("Saved ", write_counts(joinpath(directory, "spike-counts.tsv"), rows))
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
