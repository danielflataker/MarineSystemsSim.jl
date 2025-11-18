# test/test_damping_3dof.jl

using Test
using MarineSystemsSim
using StaticArrays

@testset "3DOF damping: zero velocity" begin
    rb = RigidBody3DOF(10.0, 20.0, 0.0)

    # Fossen-style derivater:
    # - ingen added mass
    # - lineær dampingderivater (Xu, Yv, Yr, Nv, Nr) negative på diagonalen
    # - kvadratisk dampingderivater negative
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
    model  = build_cached_vessel(params)

    ν  = @SVector [0.0, 0.0, 0.0]
    Dν = damping_forces(ν, model)

    @test Dν == @SVector [0.0, 0.0, 0.0]
end


@testset "3DOF damping: pure surge" begin
    rb = RigidBody3DOF(10.0, 20.0, 0.0)

    # Kun surge-derivater aktiv:
    # Xu og Xuu < 0 (Fossen-style), ingen andre bidrag
    Xu  = -1.0
    Xuu = -0.5

    hydro = hydroparams_fossen3dof(
        # Added mass
        0.0, 0.0, 0.0, 0.0,   # Xudot, Yvdot, Yrdot, Nrdot
        # Linear damping derivatives
        Xu,  0.0, 0.0, 0.0, 0.0,  # Xu, Yv, Yr, Nv, Nr
        # Quadratic damping derivatives
        Xuu, 0.0, 0.0           # Xuu, Yvv, Nrr
    )

    params = VesselParams3DOF(rb, hydro)
    model  = build_cached_vessel(params)

    u = 2.0
    ν = @SVector [u, 0.0, 0.0]

    Dν = damping_forces(ν, model)

    # I intern representasjon:
    # D_lin(1,1) = -Xu
    # D_nl(1,1)  = -Xuu*|u|
    # => D_total(1,1) = -Xu - Xuu*|u|
    # => (Dν)_u = (-Xu - Xuu*|u|)*u
    Du_expected = (-Xu - Xuu*abs(u)) * u

    @test Dν ≈ @SVector [Du_expected, 0.0, 0.0]
end


@testset "3DOF damping: dissipative" begin
    rb = RigidBody3DOF(10.0, 20.0, 0.0)

    # Dempende system:
    # - diagonal Fossen-derivater (Xu, Yv, Nr) negative
    # - kvadratisk derivater negative
    hydro = hydroparams_fossen3dof(
        # Added mass
        0.0, 0.0, 0.0, 0.0,     # Xudot, Yvdot, Yrdot, Nrdot
        # Linear damping derivatives (Xu, Yv, Yr, Nv, Nr)
        -1.0, -2.0, 0.0, 0.0, -3.0,
        # Quadratic damping derivatives (Xuu, Yvv, Nrr)
        -0.1, -0.2, -0.3
    )

    params = VesselParams3DOF(rb, hydro)
    model  = build_cached_vessel(params)

    ν  = @SVector [1.0, -0.5, 0.3]
    Dν = damping_forces(ν, model)

    # Med vår signkonvensjon er mekanisk effekt fra damping
    # P_damp = -νᵀ D ν, som skal være ≤ 0 for dissipativ damping.
    # Her tester vi νᵀ D ν ≥ 0.
    P = ν[1]*Dν[1] + ν[2]*Dν[2] + ν[3]*Dν[3]

    @test P ≥ -1e-10  # skal være ikke-negativ (opp til numerisk støv)
end
