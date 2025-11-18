using Test
using MarineSystemsSim
using StaticArrays
using ForwardDiff
using DifferentialEquations
using Unitful
using Documenter

@testset "MarineSystemsSim 3DOF" begin
    include("test_kinematics_3dof.jl")
    include("test_mass_matrix_3dof.jl")
    include("test_coriolis_3dof.jl")
    include("test_hydroparams_fossen3dof_damping.jl")
    include("test_damping_3dof.jl")
    include("test_cached_vessel_3dof.jl")
    include("test_body_dynamics_3dof.jl")
    include("test_vessel_dynamics_3dof.jl")
    include("test_vessel_dynamics_ad_3dof.jl")
    include("test_simulation_3dof.jl")
    include("test_types_3dof.jl")
    include("test_unitful_rigidbody_3dof.jl")
    include("test_physical_warnings_3dof.jl")

    # Run all doctests (docstrings + docs/src/*.md).
    # This behaves like a nested @testset.
    doctest(MarineSystemsSim)
end
