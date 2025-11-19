# src/MarineSystemsSim.jl

__precompile__(false)

module MarineSystemsSim

    using StaticArrays
    using DocStringExtensions
    using Unitful
    using LinearAlgebra: cond, dot

    include("core.jl")
    include("params_3dof.jl")
    include("frames_3dof.jl")
    include("matrices_3dof.jl")
    include("model_3dof.jl")
    include("checks.jl")

    export
        # types
        RigidBody3DOF,
        HydroParams3DOF,
        VesselParams3DOF,
        Vessel3DOF,
        # core API
        kinematics,
        mass_matrix,
        damping_forces,
        body_dynamics,
        vessel_dynamics,
        vessel_rhs!,
        dofs,
        # 3DOF Fossen helpers
        hydroparams_fossen3dof,
        vesselparams_fossen3dof,
        build_vessel3dof_fossen

end