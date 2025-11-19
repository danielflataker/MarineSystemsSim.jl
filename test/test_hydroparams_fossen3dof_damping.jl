# test/test_hydroparams_fossen3dof_damping.jl

using Test
using MarineSystemsSim

@testset "3DOF linear damping matrix assembly" begin
    # Distinct Fossen-style linear damping derivatives (typically ≤ 0 for drag)
    # so we can see if any entries are misplaced.
    Xu  = -1.0
    Xv  = -2.0
    Xr  = -3.0
    Yu  = -4.0
    Yv  = -5.0
    Yr  = -6.0
    Nu  = -7.0
    Nv  = -8.0
    Nr  = -9.0

    hydro = hydroparams_fossen3dof(
        # Added mass (not used here, but must be provided)
        0.0, 0.0, 0.0, 0.0,      # Xudot, Yvdot, Yrdot, Nrdot
        # Linear damping derivatives (Xu, Yv, Yr, Nv, Nr)
        Xu,  Yv,  Yr,  Nv,  Nr,
        # Quadratic damping derivatives (Xuu, Yvv, Nrr) – not used here
        0.0, 0.0, 0.0;
        # Cross-derivatives (Xv, Xr, Yu, Nu)
        Xv = Xv,
        Xr = Xr,
        Yu = Yu,
        Nu = Nu,
    )

    D = MarineSystemsSim._D_lin(hydro)

    @test size(D) == (3, 3)

    # D_lin is defined as:
    # D_lin = -[Xu  Xv  Xr;
    #           Yu  Yv  Yr;
    #           Nu  Nv  Nr]
    @test D[1, 1] == -Xu
    @test D[1, 2] == -Xv
    @test D[1, 3] == -Xr

    @test D[2, 1] == -Yu
    @test D[2, 2] == -Yv
    @test D[2, 3] == -Yr

    @test D[3, 1] == -Nu
    @test D[3, 2] == -Nv
    @test D[3, 3] == -Nr
end
