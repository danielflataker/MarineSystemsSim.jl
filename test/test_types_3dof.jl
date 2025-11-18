# test/test_vessel_dynamics_ad_3dof.jl

using Test
using MarineSystemsSim
using StaticArrays
using ForwardDiff

@testset "3DOF dynamics: works with Dual numbers" begin
    rb = RigidBody3DOF(10.0, 20.0, 0.0)

    # Hydrodynamic parameters given as Fossen-style derivatives:
    # - added mass > 0
    # - linear damping derivatives < 0 (drag)
    # - quadratic damping derivatives < 0 (drag)
    hydro = hydroparams_fossen3dof(
        # Added mass derivatives
        0.1, 0.2, 0.0, 0.3,     # Xudot, Yvdot, Yrdot, Nrdot
        # Linear damping derivatives (Xu, Yv, Yr, Nv, Nr)
        -1.0, -1.0, -1.0, -1.0, -1.0,
        # Quadratic damping derivatives (Xuu, Yvv, Nrr)
        -0.1, -0.1, -0.1
    )

    params = VesselParams3DOF(rb, hydro)
    model  = build_cached_vessel(params)

    # ν with real components; ForwardDiff will replace these with Dual numbers
    ν0 = @SVector [1.0, 0.5, -0.2]
    τ0 = @SVector [0.0, 0.0, 0.0]

    # Wrap body_dynamics in a function that takes a standard Vector{T},
    # so ForwardDiff can compute the Jacobian w.r.t. ν.
    f(ν_vec) = begin
        ν = SVector{3}(ν_vec)
        body_dynamics(ν, model, τ0)
    end

    J = ForwardDiff.jacobian(f, collect(ν0))

    @test size(J) == (3, 3)
end
