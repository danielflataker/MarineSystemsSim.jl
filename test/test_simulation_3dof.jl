# test/test_simulation_3dof.jl

using Test
using MarineSystemsSim
using StaticArrays
using DifferentialEquations

@testset "3DOF end-to-end: constant surge thrust" begin
    # --- 1. Set up a simple 3DOF model without damping / added mass ---

    m   = 10.0
    Iz  = 20.0
    xG  = 0.0

    rb = RigidBody3DOF(m, Iz, xG)

    # All hydrodynamic derivatives = 0 → no added mass, no damping
    hydro = hydroparams_fossen3dof(
        # Added mass derivatives
        0.0, 0.0, 0.0, 0.0,   # Xudot, Yvdot, Yrdot, Nrdot
        # Linear damping derivatives (Xu, Yv, Yr, Nv, Nr)
        0.0, 0.0, 0.0, 0.0, 0.0,
        # Quadratic damping derivatives (Xuu, Yvv, Nrr)
        0.0, 0.0, 0.0
    )

    params = VesselParams3DOF(rb, hydro)
    model  = build_cached_vessel(params)

    # --- 2. Constant thrust in surge ---

    τu = 1.0                    # [N]
    τ  = @SVector [τu, 0.0, 0.0]

    # --- 3. Initial state ---

    # X = [x, y, ψ, u, v, r]
    X0_s = @SVector [0.0, 0.0, 0.0,   # start at origin, ψ = 0
                     0.0, 0.0, 0.0]   # no initial velocity
    # ODEProblem expects a standard Vector, so we copy once
    X0 = collect(X0_s)

    # --- 4. RHS wrapper using vessel_dynamics ---

    function rhs!(dX, X, p, t)
        model, τ = p
        Xs  = SVector{6}(X)
        dXs = vessel_dynamics(Xs, model, τ)
        @inbounds for i in 1:6
            dX[i] = dXs[i]
        end
        return nothing
    end

    p     = (model, τ)
    tspan = (0.0, 10.0)

    prob = ODEProblem(rhs!, X0, tspan, p)
    sol  = solve(prob, Tsit5(); reltol=1e-10, abstol=1e-10)

    # --- 5. Compare analytic solution with numerical solution ---

    t_end = sol.t[end]
    X_end = sol.u[end]

    x_end = X_end[1]
    u_end = X_end[4]
    v_end = X_end[5]
    r_end = X_end[6]

    a = τu / m          # expected constant surge acceleration

    # u(t) = a t
    @test isapprox(u_end, a * t_end; rtol=1e-3, atol=1e-6)

    # x(t) = 0.5 a t^2
    @test isapprox(x_end, 0.5 * a * t_end^2; rtol=1e-3, atol=1e-6)

    # v(t), r(t) should remain ~0
    @test isapprox(v_end, 0.0; atol=1e-6)
    @test isapprox(r_end, 0.0; atol=1e-6)
end
