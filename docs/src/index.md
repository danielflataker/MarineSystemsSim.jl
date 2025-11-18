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
* [`CachedVessel3DOF`](@ref): preassembled inertia, inverse, and linear damping.

**Kinematics and frames**

* [`rotation`](@ref): body → Earth rotation matrix in the horizontal plane.
* [`kinematics`](@ref): computes (\dot{\eta}) from (\eta) and (\nu).

**Core dynamics**

* [`mass_matrix`](@ref): builds (\mathbf{M} = \mathbf{M}_{RB} + \mathbf{M}_A).
* [`build_cached_vessel`](@ref): caches (\mathbf{M}), (\mathbf{M}^{-1}) and (\mathbf{D}_\text{lin}).
* [`damping_forces`](@ref): evaluates (\mathbf{D}(\nu),\nu).
* [`coriolis_forces`](@ref): evaluates (\mathbf{C}(\nu),\nu).
* [`body_dynamics`](@ref): computes (\dot{\nu}).
* [`vessel_dynamics`](@ref): computes the full 6-state derivative.

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
    Xu,    Yv,    Yr,    Nv, Nr,  # linear damping
    Xuu,   Yvv,   Nrr,            # quadratic damping
)

# 3. Bundle parameters and build cached model
params = VesselParams3DOF(rb, hydro)
model  = build_cached_vessel(params)

# 4. Evaluate state derivative once
X = @SVector [x, y, ψ, u, v, r]   # 6×1 state
τ = @SVector [τx, τy, τψ]         # 3×1 generalized forces

Xdot = vessel_dynamics(X, model, τ)
```

This is the same structure you use inside an ODE right-hand side; the
[Getting started](@ref) page shows the full workflow with `OrdinaryDiffEq.jl`
and concrete numbers.