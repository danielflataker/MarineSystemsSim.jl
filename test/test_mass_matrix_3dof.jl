# test/test_mass_matrix_3dof.jl

using Test
using MarineSystemsSim
using StaticArrays

@testset "3DOF mass matrix: rigid-body + added mass (Fossen form)" begin
    m   = 10.0
    Iz  = 25.0
    xG  = 1.5

    # Fossen-style added-mass derivatives (simplified symmetric 3DOF form)
    Xudot = 1.0
    Yvdot = 2.0
    Yrdot = 0.3
    Nrdot = 0.7

    rb = RigidBody3DOF(m, Iz, xG)

    hydro = hydroparams_fossen3dof(
        # Added mass derivatives (Fossen 3DOF form: symmetric, no separate Nvdot)
        Xudot, Yvdot, Yrdot, Nrdot,
        # Linear damping derivatives (Xu, Yv, Yr, Nv, Nr) – all zero in this test
        0.0,   0.0,  0.0,  0.0,  0.0,
        # Quadratic damping derivatives (Xuu, Yvv, Nrr)
        0.0,   0.0,  0.0,
    )

    params = VesselParams3DOF(rb, hydro)

    M = mass_matrix(params)

    # Reference mass matrix:
    # M = M_RB + M_A, with
    # M_RB =
    # [ m    0       0;
    #   0    m      m xG;
    #   0   m xG    Iz  ]
    #
    # M_A =
    # [ -Xudot      0          0;
    #    0        -Yvdot    -Yrdot;
    #    0        -Yrdot    -Nrdot ]
    #
    # => M =
    # [ m - Xudot              0                 0;
    #   0                      m - Yvdot         m xG - Yrdot;
    #   0                      m xG - Yrdot      Iz - Nrdot ]
    M_ref = @SMatrix [
        m - Xudot           0.0               0.0;
        0.0                 m - Yvdot         m * xG - Yrdot;
        0.0                 m * xG - Yrdot    Iz - Nrdot
    ]

    @test M == M_ref
end
