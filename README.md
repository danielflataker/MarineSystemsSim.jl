# MarineSystemsSim.jl

Fossen-style 3-DOF surface vessel dynamics (surge–sway–yaw) in Julia.

`MarineSystemsSim.jl` is a small, lightweight port of the core 3-DOF model
from Thor I. Fossen’s Marine Systems Simulator (MSS) to the Julia ecosystem.
It focuses on clean equations of motion and efficient numerics, and delegates
all heavy lifting in time integration and numerics to existing Julia packages.

> Note: This package is under active development and is not yet registered in the Julia General registry. Expect breaking changes.

## Documentation

The latest docs for the latest commit in main branch (experimental) is available at  
https://danielflataker.github.io/MarineSystemsSim.jl/dev/

## Features

- **3-DOF Fossen model**  
  Classical surge–sway–yaw equations:

  $$ \mathbf{M}\dot{\nu} + \mathbf{C}(\nu) \\, \nu + \mathbf{D}(\nu) \\, \nu = \tau $$

  with rigid-body + added mass, Coriolis/centripetal, and linear + quadratic damping.

- **Lightweight, “MSS-inspired” core**  
  - Implements a subset of MSS: the planar 3-DOF dynamics and kinematics.
  - Not a full simulator, scenario engine, or toolchain like MSS (and likely
    never will be in that sense).

- **Fossen-style hydrodynamic derivatives**  
  - Helper `hydroparams_fossen3dof` maps added-mass and damping derivatives
    directly to internal matrices.

- **Fast small systems**  
  - Uses `StaticArrays` for 3×3 and 3×1 quantities.
  - Caches mass matrix and its inverse in `CachedVessel3DOF` for tight
    simulation loops.

- **Unit-aware rigid body parameters**  
  - `RigidBody3DOF` has a `Unitful.jl` constructor: pass `kg`, `kg*m^2`, `m`, etc.,
    and they are converted to SI internally.

- **AD-friendly**  
  - Dynamics are written to work with `ForwardDiff.Dual` numbers for
    linearization, sensitivity analysis, or MPC.

## Relation to Fossen’s Marine Systems Simulator (MSS)

This package is **not** a replacement for MSS. It is:

- a small, Julia-based reimplementation of the **3-DOF Fossen model**, and
- heavily inspired by the notation, sign conventions, and structure in MSS.

For full models, extensive examples, and original mathematical derivations,
you should refer to Fossen’s own resources:

- Marine Systems Simulator (MATLAB/Simulink):  
  <https://github.com/cybergalactic/MSS>

and his textbooks and papers. `MarineSystemsSim.jl` is intended as a Julia-flavoured,
minimal core you can script around, not a full simulation environment.

## Design philosophy

- **Do one thing well:** implement the 3-DOF equations cleanly and transparently.
- **Reuse the Julia ecosystem:** leave ODE solving, stiffness handling, and
  numerics to packages like `OrdinaryDiffEq.jl` instead of reinventing them.
- **Be a good AD citizen:** keep the dynamics as generic and type-stable as
  possible so tools like `ForwardDiff.jl` and `Optimization.jl` can be used
  on top.

## Installation

This package is **not yet registered**. If you want to try it out, you can still install it locally.

### Option 1: Using Julia’s package manager (recommended)

In Julia, press `]` to enter the Pkg REPL and run:

```julia
pkg> develop https://github.com/danielflataker/MarineSystemsSim.jl
````

This clones the repo into your `.julia/dev` folder and makes it available as a package:

```julia
julia> using MarineSystemsSim
```

### Option 2: Download ZIP manually

1. On GitHub, click **Code → Download ZIP**.
2. Extract the ZIP somewhere on your computer (for example `C:\Users\you\Projects\MarineSystemsSim.jl` or `~/Projects/MarineSystemsSim.jl`).
3. In Julia, enter the Pkg REPL (`]`) and run:

```julia
pkg> develop /full/path/to/MarineSystemsSim.jl
```

Now you can use the package as normal:

```julia
julia> using MarineSystemsSim
```

> Note: The package is under active development; APIs may change without notice.

## Quick start

Typical workflow:

1. Define rigid-body parameters.
2. Define hydrodynamic parameters (Fossen-style derivatives).
3. Build a cached model.
4. Use `vessel_dynamics` inside an ODE solver.

```julia
using MarineSystemsSim
using StaticArrays
using OrdinaryDiffEq

# 1. Rigid body: m [kg], Iz [kg m^2], xG [m]
rb = RigidBody3DOF(10.0, 20.0, 0.0)

# 2. Hydrodynamics from Fossen-style derivatives
#    (Signs follow Fossen: damping derivatives are typically ≤ 0.)
hydro = hydroparams_fossen3dof(
    1.0, 1.5, 0.2, 0.8,          # Xudot, Yvdot, Yrdot, Nrdot
   -1.0, -1.2, -0.5, -0.6, -1.0, # Xu, Yv, Yr, Nv, Nr
   -0.1, -0.2, -0.3              # Xuu, Yvv, Nrr
)

params = VesselParams3DOF(rb, hydro)

# 3. Build cached model (precomputes M, M⁻¹, D_lin, and runs some sanity checks)
model = build_cached_vessel(params)

# 4. State and input
# X = [x, y, ψ, u, v, r]
X0 = @SVector [0.0, 0.0, 0.0,
               0.0, 0.0, 0.0]   # start at rest, ψ = 0

τ  = @SVector [1.0, 0.0, 0.0]   # 1 N thrust in surge

# Right-hand side for ODEProblem (in-place form)
function rhs!(dX, X, p, t)
    Xs  = SVector{6}(X)
    dXs = vessel_dynamics(Xs, p, τ)
    @inbounds for i in 1:6
        dX[i] = dXs[i]
    end
    return nothing
end

prob = ODEProblem(rhs!, collect(X0), (0.0, 50.0), model)
sol  = solve(prob, Tsit5())

# Example: extract x(t), y(t) trajectory
x = sol[1, :]
y = sol[2, :]
```

For lower-level access you can also work directly with:

* `mass_matrix(params)` – build ($M = M_{RB} + M_A$)
* `damping_forces(ν, model)` – compute ($D(\nu)\,\nu$)
* `body_dynamics(ν, model, τ)` – compute (\dot{ν})
* `kinematics(η, ν)` – compute ($\dot{\eta}$) from pose and body velocity

## Status

This is an **early-stage** package:

* API focused on a clean 3-DOF core that mirrors Fossen’s notation.
* Tested with recent Julia 1.x versions.
* Breaking changes may still occur while the API settles.

Bug reports, suggestions, and small PRs are very welcome.

## Credits

* Model structure, notation and sign conventions closely follow the work of
  **Thor I. Fossen** and his *Marine Systems Simulator (MSS)* project.
* For full derivations, extended models, and MATLAB/Simulink implementations,
  please refer to his textbooks and the MSS repository:
  [https://github.com/cybergalactic/MSS](https://github.com/cybergalactic/MSS).
* Any mistakes or deviations in this Julia port are the author’s, not Fossen’s.

## Documentation

The package is documented using `Documenter.jl`. Once hosted, the docs link
will be added here. In the meantime, see:

* docstrings in the source, and
* the examples under `docs/` and `test/`.

## License

MarineSystemsSim.jl is intended to be released under a permissive open-source
license (e.g. MIT). See the `LICENSE` file in the repository once added.
