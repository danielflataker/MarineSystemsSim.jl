# test/test_cached_vessel_3dof.jl

using Test
using MarineSystemsSim
using LinearAlgebra   # for I

@testset "CachedVessel3DOF assembly" begin
    # Simple, but non-trivial rigid-body parameters
    rb = RigidBody3DOF(10.0, 25.0, 1.5)

    # Hydrodynamic parameters given as Fossen-style derivatives.
    # - Added mass derivatives: positive
    # - Linear damping derivatives: typically ≤ 0 (drag)
    # - Quadratic damping derivatives: typically ≤ 0 (drag)
    hydro = hydroparams_fossen3dof(
        # Added mass
        1.0, 2.0, 0.5, 0.7,      # Xudot, Yvdot, Yrdot, Nrdot
        # Linear damping
        -1.0, -0.5, -0.2, -0.3, -0.7,  # Xu, Yv, Yr, Nv, Nr
        # Quadratic damping
        -0.1, -0.2, -0.3;        # Xuu, Yvv, Nrr
        # Optional cross-derivatives
        Xv = -0.4,
        Xr = -0.6,
        Yu = -0.8,
        Nu = -0.9,
    )

    params = VesselParams3DOF(rb, hydro)
    model  = build_cached_vessel(params)

    # 1) M must match the standalone mass_matrix
    @test model.M == mass_matrix(params)

    # 2) Minv should be a true inverse (within numerical tolerance)
    @test model.Minv * model.M ≈ I(3)
    @test model.M * model.Minv ≈ I(3)

    # 3) D_lin must match the linear damping matrix returned by the internal helper
    D_ref = MarineSystemsSim._D_lin(hydro)
    @test model.D_lin == D_ref
end
