using DynamicalSystemsBase
using Kneading
using Kneading.Diagrams: ParameterPlane

function rossler!(du, u, p, t)
    x, y, z = u
    a, c = p
    du[1] = -y - z
    du[2] = x + a * y
    du[3] = 0.3 * x + z * (x - c)
    return nothing
end

function rossler_diagram(; resolution = 256, progress = nothing)
    plane = ParameterPlane(
        collect(range(2.0, 7.0; length = resolution)),
        collect(range(0.30, 0.55; length = resolution));
        xname = "c", yname = "a",
    )
    initializer = SaddleFocusInitializer(
        equilibrium_guess = zeros(3),
        critical_kind = :minimum,
        initial_event_index = 4,
        refine = false,
        newton_derivative = :finite_difference,
        criticality_tolerance = 1e-6,
        abstol = 1e-9, reltol = 1e-9,
        rho_samples = 45,
    )
    return scan_flow_kneading(
        (c, a) -> CoupledODEs(rossler!, zeros(3), [a, c]), plane;
        initializer,
        capture = LocalMinimum(2),
        observable = CoordinateComponent(2),
        word_length = 7,
        maximum_time = 2000.0,
        integration = :rk4,
        dt = 0.05,
        minimum_event_separation = 0.025,
        sign_atol = 0.0, sign_rtol = 0.0,
        max_state = 1e6,
        progress,
    )
end

function main()
    resolution = parse(Int, get(ENV, "ROSSLER_RESOLUTION", "256"))
    destination = get(ENV, "ROSSLER_OUTPUT", joinpath(@__DIR__, "..", "output", "rossler-flow-kneading.tsv"))
    visited = Ref(0)
    last_report = Ref(time())
    progress = (i, j, result) -> begin
        visited[] += 1
        if time() - last_report[] > 30 || visited[] == resolution^2
            println("Completed $(visited[])/$(resolution^2) parameter points")
            flush(stdout)
            last_report[] = time()
        end
    end
    diagram = rossler_diagram(; resolution, progress)
    write_flow_scan(destination, diagram)
    println("Complete words: ", count(==(:complete), diagram.statuses))
    println("Initialization failures: ", count(==(:initialization_failed), diagram.statuses))
    println("Saved ", abspath(destination))
    return diagram
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
