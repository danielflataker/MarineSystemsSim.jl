# 3-DOF model

```@meta
CurrentModule = MarineSystemsSim
```

This page summarizes the 3-DOF horizontal-plane model implemented in
`MarineSystemsSim.jl` and how it maps to the main API functions.

---

## Equation form

We follow the standard Fossen 3-DOF surface craft model in surge–sway–yaw:

```math
\mathbf{M} \dot{\boldsymbol{\nu}} +
\mathbf{C}(\boldsymbol{\nu}) \boldsymbol{\nu} +
\mathbf{D}(\boldsymbol{\nu}) \boldsymbol{\nu}
= \boldsymbol{\tau},
```

with

* body-fixed velocity vector
  (\boldsymbol{\nu} = [u, v, r]^T),
* generalized forces/moments
  (\boldsymbol{\tau} = [\tau_x, \tau_y, \tau_\psi]^T),
* constant inertia matrix (\mathbf{M} = \mathbf{M}_{RB} + \mathbf{M}_A),
* Coriolis/centripetal term (\mathbf{C}(\boldsymbol{\nu}) \boldsymbol{\nu}),
* linear + quadratic hydrodynamic damping
  (\mathbf{D}(\boldsymbol{\nu}) \boldsymbol{\nu}).

The kinematics relate body-fixed velocities to Earth-fixed pose
(\boldsymbol{\eta} = [x, y, \psi]^T) via

```math
\dot{\boldsymbol{\eta}} = \mathbf{R}(\psi)\, \boldsymbol{\nu}.
```

The full state is

```math
\mathbf{X} =
[x, y, \psi, u, v, r]^T.
```

---

## Parameters and mass matrix

The rigid-body parameters and hydrodynamic coefficients are stored in
[`RigidBody3DOF`](@ref), [`QuadraticDamping3DOF`](@ref),
[`HydroParams3DOF`](@ref) and [`VesselParams3DOF`](@ref). The assembled
inertia and damping are cached in [`CachedVessel3DOF`](@ref), and the
inertia matrix is constructed by [`mass_matrix`](@ref) and cached via
[`build_cached_vessel`](@ref).

### Rigid-body mass

For a surface vessel with mass `m`, yaw inertia `Iz` and centre of
gravity at `(x_G, 0)` in the body-fixed frame, the rigid-body mass matrix
is

```math
\mathbf{M}_{RB} =
\begin{bmatrix}
 m &   0 &   0 \\
 0 &   m & m x_G \\
 0 & m x_G & I_z
\end{bmatrix}.
```

This is assembled internally from [`RigidBody3DOF`](@ref) by a private
helper and used by [`mass_matrix`](@ref).

### Added mass

The added-mass matrix (\mathbf{M}_A) is stored directly inside
[`HydroParams3DOF`](@ref). When constructed via
[`hydroparams_fossen3dof`](@ref) it has the simplified Fossen 3-DOF
structure

```math
\mathbf{M}_A =
\begin{bmatrix}
- X_{\dot{u}} & 0              & 0 \\
0             & - Y_{\dot{v}}  & - Y_{\dot{r}} \\
0             & - Y_{\dot{r}}  & - N_{\dot{r}}
\end{bmatrix},
```

where (X_{\dot{u}}, Y_{\dot{v}}, Y_{\dot{r}}, N_{\dot{r}}) are the usual
Fossen added-mass derivatives.

The total inertia matrix is

```math
\mathbf{M} = \mathbf{M}_{RB} + \mathbf{M}_A,
```

and is constructed by [`mass_matrix`](@ref) and cached in
[`CachedVessel3DOF`](@ref).

---

## Coriolis and centripetal terms

The Coriolis/centripetal contribution is modeled as

```math
\mathbf{C}(\boldsymbol{\nu}) \boldsymbol{\nu}
=
\bigl(
\mathbf{C}_{RB}(\boldsymbol{\nu})
+
\mathbf{C}_A(\boldsymbol{\nu})
\bigr)\, \boldsymbol{\nu},
```

with rigid-body and added-mass parts following Fossen’s 3-DOF expressions.
Both matrices are constructed internally from the parameters and satisfy

```math
\boldsymbol{\nu}^T \mathbf{C}_{RB}(\boldsymbol{\nu}) \boldsymbol{\nu}
= 0,
\qquad
\boldsymbol{\nu}^T \mathbf{C}_A(\boldsymbol{\nu}) \boldsymbol{\nu}
= 0,
```

so that they do not create or dissipate energy.

The product (\mathbf{C}(\boldsymbol{\nu}) \boldsymbol{\nu}) is provided by [`coriolis_forces`](@ref).

---

## Hydrodynamic damping

Hydrodynamic damping is represented as a linear part plus a diagonal
quadratic part:

```math
\mathbf{D}(\boldsymbol{\nu}) \boldsymbol{\nu}
=
\mathbf{D}_{\text{lin}} \boldsymbol{\nu}
+
\mathbf{D}_{\text{n}}(\boldsymbol{\nu}) \boldsymbol{\nu}.
```

### Linear damping

The physical linear damping matrix (\mathbf{D}_{\text{lin}}) is stored
directly in [`HydroParams3DOF`](@ref) and cached in
[`CachedVessel3DOF`](@ref). It is used exactly as it appears in the
equations:

```math
\mathbf{D}_{\text{lin}} \boldsymbol{\nu}.
```

When [`hydroparams_fossen3dof`](@ref) is used, (\mathbf{D}_{\text{lin}})
is constructed from the Fossen derivatives (X_u, X_v, X_r, \dots) as

```math
\mathbf{D}_{\text{lin}} =
- \begin{bmatrix}
X_u & X_v & X_r \\
Y_u & Y_v & Y_r \\
N_u & N_v & N_r
\end{bmatrix}.
```

### Quadratic damping

The quadratic damping coefficients are stored in
[`QuadraticDamping3DOF`](@ref), holding the entries
(X_{uu}, Y_{vv}, N_{rr}). For a given velocity
(\boldsymbol{\nu} = [u, v, r]^T), the diagonal quadratic matrix is

```math
\mathbf{D}_{\text{n}}(\boldsymbol{\nu})
=
- \operatorname{diag}\bigl(
X_{uu} |u|,
Y_{vv} |v|,
N_{rr} |r|
\bigr),
```

so that the total damping becomes

```math
\mathbf{D}(\boldsymbol{\nu}) =
\mathbf{D}_{\text{lin}} + \mathbf{D}_{\text{n}}(\boldsymbol{\nu}).
```

The product (\mathbf{D}(\boldsymbol{\nu}) \boldsymbol{\nu}) is computed
by [`damping_forces`](@ref).

---

## Kinematics

The 3-DOF kinematics in the horizontal plane are

```math
\dot{\boldsymbol{\eta}} = \mathbf{R}(\psi)\, \boldsymbol{\nu},
```

where (\mathbf{R}(\psi)) is the rotation from body-fixed to Earth-fixed
frame.

[`rotation`](@ref) returns the \(3 \times 3\) rotation matrix
\(\mathbf{R}(\psi)\), and [`kinematics`](@ref) computes
\(\dot{\boldsymbol{\eta}}\) for a given pose and velocity.

* `rotation(ψ)` returns the `3×3` rotation matrix (\mathbf{R}(\psi)),
* `kinematics(η, ν)` computes (\dot{\boldsymbol{\eta}}) for a given
  pose and velocity.

---

## Full state dynamics

Combining kinematics and dynamics, the full state derivative is

```math
\dot{\mathbf{X}} =
[\dot{x}, \dot{y}, \dot{\psi}, \dot{u}, \dot{v}, \dot{r}]^T,
```

where

```math
\dot{\boldsymbol{\eta}} =
\mathbf{R}(\psi)\, \boldsymbol{\nu},
\qquad
\dot{\boldsymbol{\nu}} =
\mathbf{M}^{-1} \bigl(
\boldsymbol{\tau}
- \mathbf{C}(\boldsymbol{\nu}) \boldsymbol{\nu}
- \mathbf{D}(\boldsymbol{\nu}) \boldsymbol{\nu}
\bigr).
```

The body dynamics are computed by [`body_dynamics`](@ref), and the full
state derivative by [`vessel_dynamics`](@ref), which splits the state into
\(\eta\) and \(\nu\), applies the kinematics and dynamics, and returns the
combined 6-element derivative.

* [`body_dynamics`](@ref) computes (\dot{\boldsymbol{\nu}}),
* [`vessel_dynamics`](@ref) splits `X` into `(η, ν)`, computes `η̇` and
  `ν̇`, and returns the full 6-element state derivative.

These functions are designed to be used directly with differential
equation solvers and automatic differentiation.