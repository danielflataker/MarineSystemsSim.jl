# MarineSystemsSim.jl

```@meta
CurrentModule = MarineSystemsSim
DocTestSetup  = quote
    using MarineSystemsSim
    using StaticArrays
end
```

`MarineSystemsSim.jl` is a small, Fossen-inspired library for 3-DOF
surface vessel models in Julia. The goals are:

* a clean 3-DOF surge–sway–yaw model in standard Fossen form,
* high performance (StaticArrays, cached mass matrix and inverse),
* easy integration with `DifferentialEquations.jl` and automatic differentiation.

## What the package provides

**Parameter types**

* [`RigidBody3DOF`](@ref): mass, yaw inertia, and CG position.
* [`QuadraticDamping3DOF`](@ref): quadratic manoeuvring derivatives.
* [`HydroParams3DOF`](@ref): added mass and damping matrices.
* [`VesselParams3DOF`](@ref): rigid body + hydrodynamics in one bundle.
* [`Vessel3DOF`](@ref): cached 3-DOF vessel model (implements [`AbstractVesselModel`](@ref)).

**Kinematics and frames**

* [`rotation_body_to_earth`](@ref): body → Earth rotation matrix in the horizontal plane.
* [`kinematics`](@ref): computes $\dot{\eta}$ from $\eta$ and $\nu$.

**Core dynamics and interface**

* [`AbstractVesselModel`](@ref): DOF-agnostic vessel model interface.
* [`mass_matrix`](@ref): builds $\mathbf{M} = \mathbf{M}_{RB} + \mathbf{M}_A$.
* [`damping_forces`](@ref): evaluates $\mathbf{D}(\nu)\,\nu$.
* [`coriolis_forces`](@ref): evaluates $\mathbf{C}(\nu)\,\nu$.
* [`body_dynamics`](@ref): computes $\dot{\nu}$ from $\nu$, model, and $\tau$.
* [`vessel_dynamics`](@ref): computes the full $[\,\dot{\eta}\;\dot{\nu}\,]$ state derivative.
* [`vessel_rhs!`](@ref): in-place ODE right-hand side for use with `DifferentialEquations.jl`.

**Fossen-style helpers**

* [`hydroparams_fossen3dof`](@ref): build `HydroParams3DOF` from manoeuvring derivatives.
* [`vesselparams_fossen3dof`](@ref): assemble `VesselParams3DOF` from rigid-body + manoeuvring data.
* [`build_vessel3dof_fossen`](@ref): convenience constructor for a cached 3-DOF model.

See:

* [Getting started](@ref) for a complete end-to-end example with an ODE solver.
* [3-DOF model](@ref) for the equations and how they map to code.
* [API reference](@ref) for a grouped list of types and functions.

## A tiny taste of the API

Here is a minimal sketch of how the pieces fit together. For a full runnable
example (including an ODE solver), see [Getting started](@ref).

```julia
using MarineSystemsSim
using StaticArrays

# 1. Rigid body
rb = RigidBody3DOF(m, Iz, xG)

# 2. Hydrodynamics (Fossen-style derivatives)
hydro = hydroparams_fossen3dof(
    Xudot, Yvdot, Yrdot, Nrdot,   # added mass
    Xu,    Yv,    Yr,    Nv,  Nr, # linear damping
    Xuu,   Yvv,   Nrr;            # quadratic damping
    Xv = Xv, Xr = Xr, Yu = Yu, Nu = Nu,
)

# 3. Bundle parameters and build cached model
params = VesselParams3DOF(rb, hydro)
model  = Vessel3DOF(params)  # cached M, M⁻¹ and D_lin

# 4. Evaluate state derivative once
X = @SVector [x, y, ψ, u, v, r]   # 6×1 state: [η; ν]
τ = @SVector [τx, τy, τψ]         # 3×1 generalised forces

Xdot = vessel_dynamics(X, model, τ)
```

This is the same structure you use inside an ODE right-hand side; the
[Getting started](@ref) page shows the full workflow with `OrdinaryDiffEq.jl`
and concrete numbers.