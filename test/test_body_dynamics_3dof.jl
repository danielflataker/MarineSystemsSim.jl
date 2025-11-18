using Test
using MarineSystemsSim
using StaticArrays

@testset "3DOF body dynamics: rest equilibrium" begin
    rb = RigidBody3DOF(10.0, 20.0, 0.0)

    # Fysisk fornuftige Fossen-derivater:
    # - ingen added mass
    # - negativ lineær damping (Xu, Yv, Nr) -> D_lin med positiv diag
    # - negativ kvadratisk dampingderivater
    h = hydroparams_fossen3dof(
        # Added mass derivatives (none)
        0.0, 0.0, 0.0, 0.0,          # Xudot, Yvdot, Yrdot, Nrdot
        # Linear damping derivatives (Fossen-style, typisk ≤ 0)
        -1.0, -1.0, 0.0, 0.0, -1.0, # Xu,   Yv,   Yr,  Nv,  Nr
        # Quadratic damping derivatives (Fossen-style, typisk ≤ 0)
        -0.1, -0.1, -0.1            # Xuu,  Yvv,  Nrr
    )

    params = VesselParams3DOF(rb, h)
    model  = build_cached_vessel(params)

    ν    = @SVector [0.0, 0.0, 0.0]
    τ    = @SVector [0.0, 0.0, 0.0]
    νdot = body_dynamics(ν, model, τ)

    @test νdot ≈ @SVector [0.0, 0.0, 0.0]
end


@testset "3DOF body dynamics: M⁻¹ τ" begin
    m, Iz, xG = 10.0, 20.0, 0.0
    rb = RigidBody3DOF(m, Iz, xG)

    # Ingen added mass, ingen damping
    h = hydroparams_fossen3dof(
        # Added mass
        0.0, 0.0, 0.0, 0.0,  # Xudot, Yvdot, Yrdot, Nrdot
        # Linear damping
        0.0, 0.0, 0.0, 0.0, 0.0,  # Xu, Yv, Yr, Nv, Nr
        # Quadratic damping
        0.0, 0.0, 0.0       # Xuu, Yvv, Nrr
    )

    params = VesselParams3DOF(rb, h)
    model  = build_cached_vessel(params)

    ν    = @SVector [0.0, 0.0, 0.0]
    τ    = @SVector [1.0, 0.0, 0.0]   # 1 N i surge

    νdot = body_dynamics(ν, model, τ)

    @test νdot ≈ @SVector [τ[1]/m, 0.0, 0.0]
end
