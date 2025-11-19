using Test
using MarineSystemsSim
using StaticArrays

@testset "3DOF damping: zero velocity" begin
    rb = RigidBody3DOF(10.0, 20.0, 0.0)

    # Fossen-style derivatives:
    # - no added mass
    # - linear damping derivatives (Xu, Yv, Yr, Nv, Nr) negative on the diagonal (drag)
    # - quadratic damping derivatives negative (drag)
    hydro = hydroparams_fossen3dof(
        # Added mass
        0.0, 0.0, 0.0, 0.0,          # Xudot, Yvdot, Yrdot, Nrdot
        # Linear damping derivatives (Xu, Yv, Yr, Nv, Nr)
        -1.0, -2.0, -3.0, -4.0, -5.0,
        # Quadratic damping derivatives (Xuu, Yvv, Nrr)
        -0.1, -0.2, -0.3;
        # Optional cross-derivatives
        Xv = -0.5,
        Xr = -0.6,
        Yu = -0.7,
        Nu = -0.8,
    )

    params = VesselParams3DOF(rb, hydro)
    model  = Vessel3DOF(params)

    ν  = @SVector [0.0, 0.0, 0.0]
    Dν = damping_forces(ν, model)

    @test Dν == @SVector [0.0, 0.0, 0.0]
end


@testset "3DOF damping: pure surge" begin
    rb = RigidBody3DOF(10.0, 20.0, 0.0)

    # Only surge derivatives active:
    # Xu and Xuu < 0 (Fossen-style), no other contributions
    Xu  = -1.0
    Xuu = -0.5

    hydro = hydroparams_fossen3dof(
        # Added mass
        0.0, 0.0, 0.0, 0.0,     # Xudot, Yvdot, Yrdot, Nrdot
        # Linear damping derivatives
        Xu,  0.0, 0.0, 0.0, 0.0,  # Xu, Yv, Yr, Nv, Nr
        # Quadratic damping derivatives
        Xuu, 0.0, 0.0            # Xuu, Yvv, Nrr
    )

    params = VesselParams3DOF(rb, hydro)
    model  = Vessel3DOF(params)

    u = 2.0
    ν = @SVector [u, 0.0, 0.0]

    Dν = damping_forces(ν, model)

    # Internal representation:
    # D_lin(1,1) = -Xu
    # D_nl(1,1)  = -Xuu*|u|
    # => D_total(1,1) = -Xu - Xuu*|u|
    # => (Dν)_u = (-Xu - Xuu*|u|)*u
    Du_expected = (-Xu - Xuu * abs(u)) * u

    @test Dν ≈ @SVector [Du_expected, 0.0, 0.0]
end


@testset "3DOF damping: dissipative" begin
    rb = RigidBody3DOF(10.0, 20.0, 0.0)

    # Damping configuration:
    # - diagonal Fossen derivatives (Xu, Yv, Nr) negative (drag)
    # - quadratic derivatives negative (drag)
    hydro = hydroparams_fossen3dof(
        # Added mass
        0.0, 0.0, 0.0, 0.0,        # Xudot, Yvdot, Yrdot, Nrdot
        # Linear damping derivatives (Xu, Yv, Yr, Nv, Nr)
        -1.0, -2.0, 0.0, 0.0, -3.0,
        # Quadratic damping derivatives (Xuu, Yvv, Nrr)
        -0.1, -0.2, -0.3
    )

    params = VesselParams3DOF(rb, hydro)
    model  = Vessel3DOF(params)

    ν  = @SVector [1.0, -0.5, 0.3]
    Dν = damping_forces(ν, model)

    # With our sign convention, mechanical power from damping is
    #   P_damp = -νᵀ D ν,
    # which should be ≤ 0 for dissipative damping.
    # Here we test νᵀ D ν ≥ 0.
    P = ν[1] * Dν[1] + ν[2] * Dν[2] + ν[3] * Dν[3]

    @test P ≥ -1e-10  # should be non-negative up to numerical roundoff
end
