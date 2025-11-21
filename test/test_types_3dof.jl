# test/test_types_3dof.jl

using Test
using MarineSystemsSim
using StaticArrays

import MarineSystemsSim: AbstractVesselModel

@testset "3DOF parameter and model types" begin
    # Use Float32 to verify that the element type propagates through all structs.
    m   = 10.0f0
    Iz  = 25.0f0
    xG  = 1.5f0

    rb = RigidBody3DOF(m, Iz, xG)
    @test rb isa RigidBody3DOF{Float32}

    # Distinct Fossen-style derivatives so any permutation bugs would show up.
    Xudot = 1.0f0
    Yvdot = 2.0f0
    Yrdot = 0.3f0
    Nrdot = 0.7f0

    Xu  = -1.0f0
    Yv  = -0.5f0
    Yr  = -0.2f0
    Nv  = -0.3f0
    Nr  = -0.7f0

    Xuu = -0.1f0
    Yvv = -0.2f0
    Nrr = -0.3f0

    hydro = hydroparams_fossen3dof(
        Xudot, Yvdot, Yrdot, Nrdot,
        Xu,    Yv,    Yr,    Nv,    Nr,
        Xuu,   Yvv,   Nrr;
        Xv = -0.4f0,
        Xr = -0.6f0,
        Yu = -0.8f0,
        Nu = -0.9f0,
    )

    @test hydro isa HydroParams3DOF{Float32}
    @test hydro.M_A isa SMatrix{3,3,Float32}
    @test hydro.D_lin isa SMatrix{3,3,Float32}

    params = VesselParams3DOF(rb, hydro)
    @test params isa VesselParams3DOF{Float32}

    model = Vessel3DOF(params)
    @test model isa Vessel3DOF{Float32}
    @test model isa AbstractVesselModel{3,Float32}

    @test dofs(model) == 3

    # Mass matrix from the model should have the same scalar type.
    M = mass_matrix(model)
    @test M isa SMatrix{3,3,Float32}
end


@testset "3DOF core dynamics type stability" begin
    # Use a Float64 model here; we only care that types are concrete and stable.
    rb = RigidBody3DOF(10.0, 25.0, 1.5)

    hydro = hydroparams_fossen3dof(
        # Added mass
        1.0, 2.0, 0.3, 0.7,
        # Linear damping derivatives
        -1.0, -0.5, -0.2, -0.3, -0.7,
        # Quadratic damping derivatives
        -0.1, -0.2, -0.3;
        Xv = -0.4,
        Xr = -0.6,
        Yu = -0.8,
        Nu = -0.9,
    )

    params = VesselParams3DOF(rb, hydro)
    model  = Vessel3DOF(params)

    # Body dynamics: ν̇ should be an SVector{3,Float64}.
    ν   = @SVector [1.0, 0.5, -0.2]
    τ   = @SVector [0.1, 0.0, 0.05]

    νdot = @inferred body_dynamics(ν, model, τ)
    @test νdot isa SVector{3,Float64}

    # Full vessel dynamics: Ẋ should be an SVector{6,Float64}.
    X = @SVector [0.0, 0.0, 0.1,
                  1.0, 0.5, -0.2]

    Xdot = @inferred vessel_dynamics(X, model, τ)
    @test Xdot isa SVector{6,Float64}
end
