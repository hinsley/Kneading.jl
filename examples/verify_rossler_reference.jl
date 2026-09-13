using Kneading
using DynamicalSystemsBase
using LinearAlgebra
using Printf
using Test

function reference_rows(path)
    lines = readlines(path)
    names = Symbol.(split(first(lines), '\t'))
    return [Dict(zip(names, split(line, '\t'; keepempty=true))) for line in lines[2:end]]
end

reference_number(row, name) = parse(Float64, row[name])
reference_vector(row, names) = [reference_number(row, name) for name in names]
binary_word(signs) = join(sign > 0 ? '1' : sign < 0 ? '0' : '.' for sign in signs)

function rossler_reference_system(a, c)
    function rule!(du, u, p, t)
        x, y, z = u
        a, c = p
        du[1] = -y - z
        du[2] = x + a * y
        du[3] = 0.3 * x + z * (x - c)
        return nothing
    end
    return CoupledODEs(rule!, zeros(3), [a, c])
end

function reference_orbit(system, seed)
    return flow_kneading(FlowKneadingProblem(
        system;
        initializer=seed,
        capture=LocalMinimum(2),
        observable=CoordinateComponent(2),
        word_length=7,
        maximum_time=2000.0,
        include_initial_event=true,
        integration=:rk4,
        dt=0.05,
        max_state=1e6,
        minimum_event_separation=0.025,
        sign_atol=0.0,
        sign_rtol=0.0,
    ))
end

function verify_rossler_reference(path=joinpath(@__DIR__, "..", "test", "fixtures", "rossler_reference.tsv"))
    rows = reference_rows(path)
    largest_state_error = 0.0
    largest_tangent_error = 0.0
    largest_derivative_state_error = 0.0
    words_matched = 0
    masks_matched = 0
    saved_words_matched = 0
    derivative_words_matched = 0
    derivative_rows = Set((1, 7, 13, 19, 25))

    @testset "Rössler reference: independent initialization" begin
        for (index, row) in enumerate(rows)
            a, c = reference_number(row, :a), reference_number(row, :c)
            system = rossler_reference_system(a, c)
            seed = init_saddle_focus(
                system;
                equilibrium_guess=zeros(3),
                capture=LocalMinimum(2),
                critical_kind=:minimum,
                initial_event_index=4,
                newton_derivative=:finite_difference,
                criticality_tolerance=1e-6,
                rho_range=(-24.0, -1.0),
                rho_samples=45,
                abstol=1e-9,
                reltol=1e-9,
                refine=false,
            )
            state = reference_vector(row, (:critical_x, :critical_y, :critical_z))
            tangent = reference_vector(row, (:tangent_x, :tangent_y, :tangent_z))
            reference_seed = (u0=state, Q0=reshape(tangent, :, 1))
            saved_orbit = reference_orbit(system, reference_seed)
            flow = zeros(3)
            DynamicalSystemsBase.dynamic_rule(system)(flow, state, [a, c], 0.0)
            tangent -= dot(tangent, flow) / dot(flow, flow) * flow
            tangent /= norm(tangent)
            state_error = norm(seed.u0 - state)
            tangent_error = norm(vec(seed.Q0) - tangent)
            largest_state_error = max(largest_state_error, state_error)
            largest_tangent_error = max(largest_tangent_error, tangent_error)
            result = reference_orbit(system, seed)
            expected_complete = row[:status] == "ok"
            expected_word = row[:raw_word]
            words_matched += binary_word(result.raw_word) == expected_word
            masks_matched += result.complete == expected_complete
            saved_words_matched += binary_word(saved_orbit.raw_word) == expected_word

            @testset "a=$a, c=$c" begin
                @test seed.event_index == 4
                @test seed.diagnostics.root_converged
                @test !seed.diagnostics.refinement_checked
                @test abs(seed.diagnostics.residual) <= 1e-6
                @test state_error < 1e-4
                @test tangent_error < 1e-4
                @test result.complete == expected_complete
                @test binary_word(result.raw_word) == expected_word
                @test binary_word(saved_orbit.raw_word) == expected_word
                @test length(result.events) == parse(Int, row[:events])
                @test result.events[1].time == 0.0
                @test result.transition_word == result.raw_word[1:end-1] .* result.raw_word[2:end]
                @test DynamicalSystemsBase.current_state(system) == zeros(3)
            end

            if index in derivative_rows
                second_order = init_saddle_focus(
                    system;
                    equilibrium_guess=zeros(3),
                    capture=LocalMinimum(2),
                    critical_kind=:minimum,
                    initial_event_index=4,
                    initial_rho=seed.rho + 0.02,
                    newton_derivative=:second_order_sensitivity,
                    criticality_tolerance=1e-6,
                    abstol=1e-9,
                    reltol=1e-9,
                    refine=false,
                )
                second_result = reference_orbit(system, second_order)
                difference = norm(second_order.u0 - seed.u0)
                largest_derivative_state_error = max(largest_derivative_state_error, difference)
                derivative_words_matched += second_result.raw_word == result.raw_word
                @testset "second-order sensitivity a=$a, c=$c" begin
                    @test difference < 1e-4
                    @test norm(second_order.Q0 - seed.Q0) < 1e-4
                    @test second_order.diagnostics.derivative_method == :second_order_sensitivity
                    @test second_order.diagnostics.iterations > 0
                    @test second_result.raw_word == result.raw_word
                    @test second_result.complete == result.complete
                end
            end
            @printf("Reference %d/%d: a=%.6f c=%.6f state_error=%.3g word=%s expected=%s status=%s\n",
                index, length(rows), a, c, state_error, binary_word(result.raw_word), expected_word, result.status)
            flush(stdout)
        end
    end

    @printf("Independent words: %d/%d; completeness masks: %d/%d; saved-state orbit words: %d/%d\n",
        words_matched, length(rows), masks_matched, length(rows), saved_words_matched, length(rows))
    @printf("Maximum critical-state error %.6g; tangent error %.6g; derivative-method state difference %.6g\n",
        largest_state_error, largest_tangent_error, largest_derivative_state_error)
    @printf("Derivative-method words: %d/%d\n", derivative_words_matched, length(derivative_rows))
    return (; words_matched, masks_matched, saved_words_matched, derivative_words_matched,
        largest_state_error, largest_tangent_error, largest_derivative_state_error)
end

function compare_rossler_scans(candidate_path, reference_path)
    candidates = reference_rows(candidate_path)
    references = reference_rows(reference_path)
    parameter_key(row) = (round(reference_number(row, :a); digits=9), round(reference_number(row, :c); digits=9))
    candidate_lookup = Dict(parameter_key(row) => row for row in candidates)
    @test length(candidate_lookup) == length(candidates)
    matched = 0
    masks_matched = 0
    complete_reference = 0
    complete_candidate = 0
    common_complete = 0
    complete_words_matched = 0
    prefix_words_matched = 0
    states_compared = 0
    maximum_state_error = 0.0
    squared_state_error = 0.0
    disagreements = zeros(Int, 8)
    for row in references
        key = parameter_key(row)
        haskey(candidate_lookup, key) || continue
        result = candidate_lookup[key]
        matched += 1
        expected_complete = row[:status] == "ok"
        observed_word = result[:raw_word]
        expected_word = get(row, :word, get(row, :raw_word, ""))
        observed_complete = length(observed_word) == 8 && result[:status] in ("ok", "complete")
        complete_reference += expected_complete
        complete_candidate += observed_complete
        masks_matched += observed_complete == expected_complete
        prefix_words_matched += observed_word == expected_word
        if observed_complete && expected_complete
            common_complete += 1
            complete_words_matched += observed_word == expected_word
            for bit in 1:8
                disagreements[bit] += observed_word[bit] != expected_word[bit]
            end
        end
        expected_state = reference_vector(row, (:critical_x, :critical_y, :critical_z))
        observed_state = reference_vector(result, (:critical_x, :critical_y, :critical_z))
        if all(isfinite, expected_state) && all(isfinite, observed_state)
            error = norm(observed_state - expected_state)
            maximum_state_error = max(maximum_state_error, error)
            squared_state_error += error^2
            states_compared += 1
        end
    end
    state_rmse = states_compared == 0 ? NaN : sqrt(squared_state_error / states_compared)
    println("Reference rows: $(length(references)); candidate rows: $(length(candidates)); matching grid coordinates: $matched")
    println("Complete reference: $complete_reference; complete candidate: $complete_candidate")
    println("Completeness masks matching: $masks_matched/$matched")
    println("All complete words and incomplete prefixes matching: $prefix_words_matched/$matched")
    println("Complete raw words matching: $complete_words_matched/$common_complete")
    println("Raw-bit disagreements among common complete rows: $disagreements")
    @printf("Critical states compared: %d; maximum error %.9g; RMS error %.9g\n",
        states_compared, maximum_state_error, state_rmse)
    @test matched == length(references)
    @test matched == length(candidates)
    return (; matched, masks_matched, complete_reference, complete_candidate,
        complete_words_matched, common_complete, prefix_words_matched,
        states_compared, maximum_state_error, state_rmse, disagreements)
end

if abspath(PROGRAM_FILE) == @__FILE__
    if !isempty(ARGS) && first(ARGS) == "--full"
        length(ARGS) == 3 || throw(ArgumentError("use --full candidate.tsv reference.tsv"))
        compare_rossler_scans(ARGS[2], ARGS[3])
    else
        verify_rossler_reference(isempty(ARGS) ? joinpath(@__DIR__, "..", "test", "fixtures", "rossler_reference.tsv") : only(ARGS))
    end
end
