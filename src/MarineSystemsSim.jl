__precompile__(false)

module MarineSystemsSim

    using StaticArrays
    using DocStringExtensions
    using Unitful
    using LinearAlgebra: cond, dot

    include("params_3dof.jl")
    include("frames_3dof.jl")
    include("assemble_3dof.jl")
    include("dynamics_3dof.jl")
    include("checks.jl")

    export
        RigidBody3DOF,
        HydroParams3DOF,
        VesselParams3DOF,
        kinematics,
        CachedVessel3DOF,
        build_cached_vessel,
        mass_matrix,
        damping_forces,
        body_dynamics,
        vessel_dynamics,
        hydroparams_fossen3dof
end # module MarineSystemsSim
