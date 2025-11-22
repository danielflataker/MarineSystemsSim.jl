# Getting started

```@meta
CurrentModule = MarineSystemsSim
DocTestSetup  = quote
    using MarineSystemsSim
    using StaticArrays
    using Unitful
    using ForwardDiff
end
```

This page walks through two typical workflows:

1. **Fossen-style quick start** (recommended): use
   [`build_vessel3dof_fossen`](@ref) or
   [`hydroparams_fossen3dof`](@ref) with manoeuvring derivatives.
2. **Lower-level setup**: build `HydroParams3DOF` from your own matrices.

We assume you will also use an ODE solver, e.g.:

```julia
using DifferentialEquations   # or OrdinaryDiffEq
```

---

## 1. Quick start: Fossen-style derivatives

### 1.1 Define parameters

The 3-DOF model uses:

* [`RigidBody3DOF`](@ref) for mass, inertia, and CG position,
* [`HydroParams3DOF`](@ref) for added mass and damping,
* [`VesselParams3DOF`](@ref) as a container,
* [`Vessel3DOF`](@ref) as the cached 3-DOF model,
* [`build_vessel3dof_fossen`](@ref) and
  [`vesselparams_fossen3dof`](@ref) as convenience constructors from
  Fossen-style derivatives.

The recommended “one-shot” constructor from Fossen derivatives is
[`build_vessel3dof_fossen`](@ref):

```julia
using MarineSystemsSim
using StaticArrays

# Rigid-body parameters
m  = 10.0
Iz = 25.0
xG = 1.5

# Hydrodynamic coefficients (Fossen-style)
Xudot =  1.0
Yvdot =  2.0
Yrdot =  0.3
Nrdot =  0.7

Xu  = -1.0
Yv  = -0.5
Yr  = -0.2
Nv  = -0.3
Nr  = -0.7

Xuu = -0.1
Yvv = -0.2
Nrr = -0.3

model = build_vessel3dof_fossen(
    m, Iz, xG,
    Xudot, Yvdot, Yrdot, Nrdot,
    Xu,    Yv,    Yr,    Nv,    Nr,
    Xuu,   Yvv,   Nrr;
    # Optional cross-derivatives:
    Xv = 0.0,
    Xr = 0.0,
    Yu = 0.0,
    Nu = 0.0,
)
```

Internally this constructs `RigidBody3DOF`, `HydroParams3DOF`,
`VesselParams3DOF`, and finally a cached [`Vessel3DOF`](@ref).

If you prefer to build the pieces explicitly, you can do:

```julia
rb = RigidBody3DOF(m, Iz, xG)

hydro = hydroparams_fossen3dof(
    Xudot, Yvdot, Yrdot, Nrdot,
    Xu,    Yv,    Yr,    Nv,    Nr,
    Xuu,   Yvv,   Nrr;
    Xv = 0.0,
    Xr = 0.0,
    Yu = 0.0,
    Nu = 0.0,
)

params = VesselParams3DOF(rb, hydro)
model  = Vessel3DOF(params)
```

Both approaches produce the same `Vessel3DOF` model.

### 1.2 State and input

The 3-DOF state is

```math
\mathbf{X} = [x, y, \psi, u, v, r]^T,
```

where:

* `x, y` – horizontal position (Earth-fixed frame),
* `ψ` – yaw angle,
* `u, v` – surge and sway velocities (body frame),
* `r` – yaw rate.

Example initial condition and input:

```julia
X0 = @SVector [0.0, 0.0, 0.0,   # x, y, ψ
               0.0, 0.0, 0.0]   # u, v, r

τ_const  = @SVector [1.0, 0.0, 0.0]   # constant generalized forces in surge, sway, yaw
```

We will treat the generalized forces as a function
`τfun(::SVector{6}, t)` that returns a `SVector{3}`. For a constant
thrust in surge:

```julia
τfun(X, t) = τ_const
```

### 1.3 Right-hand side using `vessel_rhs!`

The high-level function [`vessel_dynamics`](@ref) computes `Ẋ` from
state, model, and input for any [`AbstractVesselModel`](@ref). For ODE
solvers it is often convenient to use the in-place wrapper
[`vessel_rhs!`](@ref):

```julia
function rhs!(dX, X, p, t)
    model, τfun = p
    vessel_rhs!(dX, X, model, τfun, t)
end
```

This function:

* treats `X` as a length-`6` state vector `[η; ν]`,
* converts it to an `SVector` internally,
* calls [`vessel_dynamics`](@ref),
* writes the result into `dX`.

### 1.4 Solve the ODE

Here is an example using `DifferentialEquations.jl`:

```julia
using DifferentialEquations

tspan = (0.0, 100.0)
p     = (model, τfun)         # parameters: model + input function

prob = ODEProblem(
    rhs!,          # in-place RHS
    collect(X0),   # ODEProblem likes a plain Vector
    tspan,
    p,
)

sol = solve(prob, Tsit5(); reltol = 1e-8, abstol = 1e-8)
```

You can now inspect the time series:

```julia
x = sol[1, :]
y = sol[2, :]
ψ = sol[3, :]
u = sol[4, :]
v = sol[5, :]
r = sol[6, :]
```

and plot them with your preferred plotting library.

---

## 2. Lower-level usage: build `HydroParams3DOF` by hand

Sometimes you already have the hydrodynamic matrices in **physical form**
(e.g. from CFD, your own identification code, or another model that
does not use Fossen derivatives directly). In that case you can construct
[`HydroParams3DOF`](@ref) manually.

### 2.1 Rigid-body parameters

Same as before:

```julia
rb = RigidBody3DOF(10.0, 25.0, 1.5)
```

### 2.2 Added mass and damping matrices

You provide:

* `M_A`: the `3×3` added-mass matrix (physical sign),
* `D_lin`: the `3×3` linear damping matrix such that `D_lin * ν` is the
  physical linear damping force/moment,
* `QuadraticDamping3DOF`: coefficients used to build the quadratic
  damping matrix internally.

```julia
using StaticArrays

# Example: diagonal added mass
M_A = @SMatrix [
    1.0   0.0   0.0;
    0.0   2.0   0.0;
    0.0   0.0   0.5,
]

# Example: diagonal linear damping (positive diagonal → dissipative)
D_lin = @SMatrix [
    1.0   0.0   0.0;
    0.0   2.0   0.0;
    0.0   0.0   3.0,
]

# Quadratic damping coefficients.
# When constructed via hydroparams_fossen3dof these are Fossen-style
# derivatives Xuu, Yvv, Nrr (typically ≤ 0).
D_quad = QuadraticDamping3DOF(
    -0.1,   # Xuu
    -0.2,   # Yvv
    -0.3,   # Nrr
)

hydro = HydroParams3DOF(M_A, D_lin, D_quad)
params = VesselParams3DOF(rb, hydro)
model  = Vessel3DOF(params)
```

From here you can use the same high-level functions as before:

```julia
ν   = @SVector [1.0, 0.2, 0.01]   # body-fixed velocity [u, v, r]
τ   = @SVector [0.0, 0.0, 0.0]

Dν  = damping_forces(ν, model)    # D(ν) ν
Cν  = coriolis_forces(ν, model)   # C(ν) ν
ν̇   = body_dynamics(ν, model, τ) # ν̇ given τ

ψ    = 0.3
η    = @SVector [0.0, 0.0, ψ]
η̇   = kinematics(η, ν)           # [ẋ, ẏ, ψ̇]
```

This “by hand” workflow gives you full control over the matrices while
reusing the kinematics and dynamics code from `MarineSystemsSim.jl`.

## 3. Units and automatic differentiation

### Units with `Unitful.jl`

`RigidBody3DOF` has a constructor that accepts Unitful quantities and converts
them to SI units internally:

```jldoctest units
julia> using Unitful

julia> rb = RigidBody3DOF(100u"lb", 1.0e6u"lb*ft^2", 10u"ft");

julia> rb.m # mass is stored as kg internally
45.359237

julia> rb.xG # position is stored as meters internally
3.048
```

This allows you to work in whatever mass/length units you like at the API
boundary while keeping the internal state dimensionless and fast.

### Automatic differentiation with `ForwardDiff.jl`

The dynamics functions are written to be compatible with dual numbers, so you
can differentiate them with `ForwardDiff.jl`:

```jldoctest ad
julia> using ForwardDiff, StaticArrays

julia> rb = RigidBody3DOF(10.0, 20.0, 0.0);

julia> hydro = hydroparams_fossen3dof(
           0.1, 0.2, 0.0, 0.3,
           -1.0, -1.0, -1.0, -1.0, -1.0,
           -0.1, -0.1, -0.1,
       );

julia> params = VesselParams3DOF(rb, hydro);

julia> model = Vessel3DOF(params);

julia> ν0 = @SVector [1.0, 0.5, -0.2];

julia> τ0 = @SVector [0.0, 0.0, 0.0];

julia> f(ν_vec) = body_dynamics(SVector{3}(ν_vec), model, τ0);

julia> J = ForwardDiff.jacobian(f, collect(ν0));

julia> size(J)
(3, 3)
```

Here `J` is the Jacobian $\partial \dot{\nu} / \partial \nu$, which you can
use for linearization, MPC, or sensitivity analysis.
