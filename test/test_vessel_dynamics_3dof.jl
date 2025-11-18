# test/test_vessel_dynamics_3dof.jl

using Test
using MarineSystemsSim
using StaticArrays

@testset "3DOF vessel dynamics: shape and basic kinematics" begin
    rb = RigidBody3DOF(10.0, 20.0, 0.0)

    # No added mass, no damping
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

    # State X = [x, y, ψ, u, v, r]
    X  = @SVector [0.0, 0.0, 0.0,   1.0, 0.0, 0.0]  # ψ = 0, pure surge
    τ  = @SVector [0.0, 0.0, 0.0]

    Ẋ = vessel_dynamics(X, model, τ)

    @test length(Ẋ) == 6

    # For ψ = 0, we should have η̇ = ν
    @test Ẋ[1:3] ≈ @SVector [1.0, 0.0, 0.0]
end
