using Test
using MarineSystemsSim
using StaticArrays

# Helper to build a 3DOF model and a test velocity ν
function _test_coriolis_setup()
    m   = 12.0
    Iz  = 30.0
    xG  = 1.2

    # Fossen-style added-mass derivatives (typically > 0)
    Xudot = 0.8
    Yvdot = 1.1
    Yrdot = 0.4
    Nrdot = 0.9

    rb = RigidBody3DOF(m, Iz, xG)

    # No linear or quadratic damping in this test
    hydro = hydroparams_fossen3dof(
        # Added mass
        Xudot, Yvdot, Yrdot, Nrdot,
        # Linear damping derivatives (Xu, Yv, Yr, Nv, Nr)
        0.0, 0.0, 0.0, 0.0, 0.0,
        # Quadratic damping derivatives (Xuu, Yvv, Nrr)
        0.0, 0.0, 0.0,
    )

    params = VesselParams3DOF(rb, hydro)
    model  = Vessel3DOF(params)

    ν = @SVector [0.7, -0.3, 0.5]   # [u, v, r]

    return model, ν, Xudot, Yvdot, Yrdot, m, Iz, xG
end


@testset "3DOF Coriolis force: matches Fossen closed form" begin
    model, ν, Xudot, Yvdot, Yrdot, m, Iz, xG = _test_coriolis_setup()

    Cν = coriolis_forces(ν, model)

    u, v, r = ν

    # Fossen closed form: C_RB(ν)ν + C_A(ν)ν
    cRB1 = -m * xG * r^2 - m * v * r
    cRB2 =  m * u * r
    cRB3 =  m * xG * u * r

    cA1 =  Yvdot * v * r + Yrdot * r^2
    cA2 = -Xudot * u * r
    cA3 = (Xudot - Yvdot) * u * v - Yrdot * u * r

    Cν_ref = @SVector [
        cRB1 + cA1,
        cRB2 + cA2,
        cRB3 + cA3,
    ]

    @test isapprox(Cν, Cν_ref; atol = 1e-12, rtol = 1e-12)
end


@testset "3DOF Coriolis term: νᵀ C(ν) ν = 0 (power conservation)" begin
    model, ν, _, _, _, _, _, _ = _test_coriolis_setup()

    Cν = coriolis_forces(ν, model)

    scalar = ν[1] * Cν[1] + ν[2] * Cν[2] + ν[3] * Cν[3]
    @test isapprox(scalar, 0.0; atol = 1e-10)
end
