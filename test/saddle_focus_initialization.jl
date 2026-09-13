using Test
using LinearAlgebra
using DynamicalSystemsBase
using Kneading.FlowKneading

function saddle_test_rossler!(du, u, p, t)
    x, y, z = u
    a, c = p
    du[1] = -y - z
    du[2] = x + a * y
    du[3] = 0.3 * x + z * (x - c)
    return nothing
end

function saddle_test_rossler(u, p, t)
    x, y, z = u
    a, c = p
    return SVector(-y - z, x + a * y, 0.3 * x + z * (x - c))
end

@testset "Saddle-focus initialization" begin
    FK = Kneading.FlowKneading
    system = CoupledODEs(saddle_test_rossler!, zeros(3), [0.3, 5.5])
    fixed = init_saddle_focus(system;
        equilibrium_guess=zeros(3), initial_rho=-1.47, refine=false)
    @test fixed.event_index == 4
    @test fixed.rho ≈ -1.46729502895 atol=1e-6
    @test fixed.u0 ≈ [1.48188585, -4.9396195, 0.037785522] atol=3e-7
    @test abs(fixed.diagnostics.residual) < 1e-8
    @test fixed.diagnostics.root_converged
    @test !fixed.diagnostics.converged
    @test !fixed.diagnostics.refinement_checked
    @test size(fixed.Q0) == (3, 1)
    @test norm(fixed.Q0) ≈ 1.0 atol=1e-14
    ctx = FK._flow_context(system)
    @test abs(dot(vec(fixed.Q0), ctx.f(fixed.u0, ctx.p, 0.0))) < 1e-12
    @test FK._event_rate(ctx, fixed.capture, fixed.u0, 0.0) > 0
    @test abs(FK._event_value(ctx, fixed.capture, fixed.u0, 0.0)) < 1e-12

    second = init_saddle_focus(system;
        equilibrium_guess=zeros(3), initial_rho=-1.47, refine=false,
        newton_derivative=:second_order_sensitivity)
    @test second.u0 ≈ fixed.u0 atol=2e-8
    @test second.rho ≈ fixed.rho atol=1e-8
    @test second.diagnostics.curvature ≈ fixed.diagnostics.curvature rtol=2e-6

    ray = FK._sf_seed_ray(ctx, fixed.equilibrium, fixed.capture)
    rho = fixed.rho + 0.07
    options = merge(fixed.configuration, (; critical_kind=:minimum))
    result = FK._sf_residual(ctx, fixed.capture, ray, rho, 4, options; second_order=true)
    h = 1e-4
    plus = FK._sf_residual(ctx, fixed.capture, ray, rho+h, 4, options)
    minus = FK._sf_residual(ctx, fixed.capture, ray, rho-h, 4, options)
    @test result.valid && plus.valid && minus.valid
    @test result.slope ≈ (plus.residual-minus.residual)/(2h) rtol=3e-6
    for k in (4, 5)
        event = result.events[k]
        @test event.tangent ≈ (plus.events[k].state-minus.events[k].state)/(2h) rtol=3e-6 atol=1e-7
        @test event.second_tangent ≈ (plus.events[k].tangent-minus.events[k].tangent)/(2h) rtol=3e-6 atol=1e-6
        gradient = ctx.J(event.state, ctx.p, event.time)[2, :]
        @test abs(dot(gradient, event.tangent)) < 1e-11
    end
    nonlinear_capture = LocalMaximum(3)
    nonlinear_events, nonlinear_failure = FK._sf_events(ctx, nonlinear_capture, ray, rho, 3, options; second_order=true)
    nonlinear_plus, _ = FK._sf_events(ctx, nonlinear_capture, ray, rho+h, 3, options)
    nonlinear_minus, _ = FK._sf_events(ctx, nonlinear_capture, ray, rho-h, 3, options)
    @test isempty(nonlinear_failure)
    @test nonlinear_events[3].second_tangent ≈ (nonlinear_plus[3].tangent-nonlinear_minus[3].tangent)/(2h) rtol=3e-6 atol=1e-6

    refined = init_saddle_focus(system; equilibrium_guess=zeros(3), initial_rho=-1.47)
    @test refined.event_index == fixed.event_index
    @test refined.rho ≈ fixed.rho atol=1e-9
    @test refined.diagnostics.refinement_event_index == refined.event_index + 1
    @test refined.diagnostics.converged
    @test refined.diagnostics.refinement_checked
    @test refined.diagnostics.state_error <= 1e-6
    @test refined.diagnostics.tangent_error <= 1e-6
    @test refined.u0 ≈ fixed.u0 atol=1e-6
    set_parameter!(system, 2, 5.51)
    continued = init_saddle_focus(system; previous=fixed, refine=false)
    @test continued.event_index == fixed.event_index
    @test continued.parameters == [0.3, 5.51]
    @test fixed.parameters == [0.3, 5.5]
    @test norm(continued.u0-fixed.u0) < 0.1
    @test dot(vec(continued.Q0), vec(fixed.Q0)) > 0
    inherited = init_saddle_focus(system; previous=second)
    @test inherited.diagnostics.derivative_method == :second_order_sensitivity
    @test !inherited.diagnostics.refinement_checked
    @test inherited.configuration.abstol == second.configuration.abstol
    overridden = init_saddle_focus(system; previous=second, newton_derivative=:finite_difference,
        criticality_tolerance=2e-8)
    @test overridden.diagnostics.derivative_method == :finite_difference
    @test overridden.configuration.criticality_tolerance == 2e-8
    @test continued.diagnostics.target_distance < continued.diagnostics.branch_tolerance
    rejected = try
        init_saddle_focus(system; previous=fixed, refine=false, branch_tolerance=1e-8, rho_samples=2)
    catch error
        error
    end
    @test rejected isa SaddleFocusInitializationError
    @test rejected.stage == :branch

    oop = CoupledODEs(saddle_test_rossler, zeros(3), [0.3, 5.5])
    oop_seed = init_saddle_focus(oop; equilibrium_guess=zeros(3), initial_rho=-1.47, refine=false)
    @test oop_seed.u0 ≈ fixed.u0 atol=1e-9
    targeted = init_saddle_focus(oop; equilibrium_guess=zeros(3), initial_rho=-1.47,
        critical_target=fixed.u0[2], refine=false)
    @test targeted.u0 ≈ fixed.u0 atol=1e-9
    wrong_target = try
        init_saddle_focus(oop; equilibrium_guess=zeros(3), initial_rho=-1.47,
            critical_target=fixed.u0 .+ 100, rho_samples=2, refine=false)
    catch error
        error
    end
    @test wrong_target isa SaddleFocusInitializationError
    @test wrong_target.stage == :branch
    configured = SaddleFocusInitializer(equilibrium_guess=zeros(3), initial_rho=-1.47, refine=false)
    @test configured(oop).rho ≈ fixed.rho atol=1e-9

    maximum_ray = FK._sf_seed_ray(ctx, fixed.equilibrium, LocalMaximum(2))
    @test maximum_ray.direction ≈ -ray.direction atol=1e-14
    maximum_events, failure = FK._sf_events(ctx, LocalMaximum(2), maximum_ray, -1.5, 5, options)
    @test isempty(failure)
    @test length(maximum_events) == 5
    @test all(event -> event.transversality < 0, maximum_events)
    @test all(event -> abs(FK._event_value(ctx, LocalMaximum(2), event.state, event.time)) < 1e-10, maximum_events)
    maximum_seed = init_saddle_focus(oop; equilibrium_guess=zeros(3), capture=LocalMaximum(2),
        critical_kind=:maximum, initial_rho=-1.86, refine=false)
    @test maximum_seed.diagnostics.curvature < 0
    @test maximum_seed.diagnostics.event_transversality[1] < 0
    @test abs(maximum_seed.diagnostics.residual) <= 1e-8

    @test_throws ArgumentError init_saddle_focus(system)
    @test_throws ArgumentError init_saddle_focus(system; equilibrium_guess=zeros(3), initial_radius=0.0)
    @test_throws ArgumentError init_saddle_focus(system; equilibrium_guess=zeros(3), newton_derivative=:unknown)
    @test_throws ArgumentError init_saddle_focus(system; equilibrium_guess=zeros(3), critical_kind=:unknown)
    @test_throws SaddleFocusInitializationError init_saddle_focus(system;
        equilibrium_guess=zeros(3), initial_rho=-1.47, maximum_event_index=4)
    @test_throws SaddleFocusInitializationError init_saddle_focus(system;
        equilibrium_guess=zeros(3), initial_rho=-1.47, max_time=0.1, rho_samples=2, refine=false)
    stable(u, p, t) = -u
    stable_system = CoupledODEs(stable, zeros(3), nothing)
    @test_throws SaddleFocusInitializationError init_saddle_focus(stable_system; equilibrium_guess=zeros(3))
end
