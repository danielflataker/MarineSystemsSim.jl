# test/test_vessel_dynamics_ad_3dof.jl

using Test
using MarineSystemsSim
using StaticArrays
using ForwardDiff

# Helper to build a reasonable 3-DOF vessel model (non-singular, physically sound signs)
function _make_model_for_ad()
    rb = RigidBody3DOF(10.0, 25.0, 1.5)

    # Hydrodynamic parameters given as Fossen-style derivatives:
    # - added mass derivatives > 0
    # - linear damping derivatives < 0 (drag)
    # - quadratic damping derivatives < 0 (drag)
    hydro = hydroparams_fossen3dof(
        # Added mass
        1.0, 2.0, 0.3, 0.7,           # Xudot, Yvdot, Yrdot, Nrdot
        # Linear damping derivatives (Xu, Yv, Yr, Nv, Nr)
        -1.0, -0.5, -0.2, -0.3, -0.7,
        # Quadratic damping derivatives (Xuu, Yvv, Nrr)
        -0.1, -0.2, -0.3;
        # Optional cross-derivatives
        Xv = -0.4,
        Xr = -0.6,
        Yu = -0.8,
        Nu = -0.9,
    )

    params = VesselParams3DOF(rb, hydro)
    return Vessel3DOF(params)
end

@testset "3DOF vessel_dynamics works with Dual numbers" begin
    model = _make_model_for_ad()

    # State: [x, y, ψ, u, v, r]
    X0 = @SVector [0.0, 0.0, 0.1,
                   1.0, 0.5, -0.2]

    τ0 = @SVector [0.1, 0.0, 0.05]

    # Do not type-annotate X_vec; ForwardDiff will pass Vector{Dual} here.
    f(X_vec) = begin
        X = SVector{6}(X_vec)             # works for Float64 and Dual
        vessel_dynamics(X, model, τ0)     # returns SVector{6, S}
    end

    J = ForwardDiff.jacobian(f, collect(X0))

    @test size(J) == (6, 6)
    @test all(isfinite, J)
end
