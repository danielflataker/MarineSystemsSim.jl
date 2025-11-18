# src/params_3dof.jl

using Unitful: AbstractQuantity, ustrip


"""
$(TYPEDEF)

Rigid-body parameters for 3-DOF motion in the horizontal plane.

# Fields
- `m::T`: Mass of the vessel [``\\mathrm{kg}``].
- `Iz::T`: Yaw moment of inertia about the z-axis [``\\mathrm{kg \\, m^2}``].
- `xG::T`: x-position of the center of gravity in the body-fixed frame [``\\mathrm{m}``].

# Notes

The 3-DOF rigid-body mass matrix is

```math
\\mathbf{M}_{RB} =
\\begin{bmatrix}
 m & 0 & 0 \\\\
 0 & m & m x_G \\\\
 0 & m x_G & I_z
\\end{bmatrix}.
```
"""
struct RigidBody3DOF{T<:Real}
    m  :: T
    Iz :: T
    xG :: T

    # Inner constructor: enforces physical constraints
    function RigidBody3DOF(m::T, Iz::T, xG::T) where {T<:Real}
        if m ≤ zero(T)
            throw(ArgumentError("RigidBody3DOF: mass m must be > 0 (got $m)"))
        end
        if Iz ≤ zero(T)
            throw(ArgumentError("RigidBody3DOF: yaw inertia Iz must be > 0 (got $Iz)"))
        end
        new{T}(m, Iz, xG)
    end
end

# Generic numeric convenience constructor
# Disabled in order to force same types for efficiency etc.
# RigidBody3DOF(m::Real, Iz::Real, xG::Real) =
#     RigidBody3DOF(promote(m, Iz, xG)...)

# Unitful constructor for user-friendly API
RigidBody3DOF(
    m::AbstractQuantity,
    Iz::AbstractQuantity,
    xG::AbstractQuantity,
) = RigidBody3DOF(promote(
    ustrip(u"kg",      m),
    ustrip(u"kg*m^2",  Iz),
    ustrip(u"m",       xG),
)...)


"""
$(TYPEDEF)

Quadratic (velocity-dependent) damping parameters for the 3-DOF model.

This type stores Fossen-style manoeuvring derivatives for the diagonal
approximation to the nonlinear damping matrix. They appear in the
velocity-dependent damping term

```math
\\mathbf{D}_\\text{n}(\\boldsymbol{\\nu})\\,\\boldsymbol{\\nu}
=
-\\operatorname{diag}(
  X_{uu}\\lvert u\\rvert,\\,
  Y_{vv}\\lvert v\\rvert,\\,
  N_{rr}\\lvert r\\rvert
)\\,\\boldsymbol{\\nu}.
```

Sign convention:

* `Xuu`, `Yvv`, `Nrr` correspond to the hydrodynamic derivatives
  `X_{uu}, Y_{vv}, N_{rr}` in Fossen's notation.
* In many identification schemes these are non-positive so that the
  effective damping is dissipative, but this is not enforced here.
"""
struct QuadraticDamping3DOF{T<:Real}
    Xuu :: T
    Yvv :: T
    Nrr :: T
    # senere: Yvr, Nvr, osv.
end


"""
$(TYPEDEF)

Hydrodynamic parameters for 3-DOF surge–sway–yaw motion.

The 3-DOF vessel model is written as

```math
\\mathbf{M} \\dot{\\boldsymbol{\\nu}} +
\\mathbf{C}(\\boldsymbol{\\nu})\\,\\boldsymbol{\\nu} +
\\mathbf{D}(\\boldsymbol{\\nu})\\,\\boldsymbol{\\nu}
= \\boldsymbol{\\tau},
```

where `\\mathbf{M} = \\mathbf{M}_{RB} + \\mathbf{M}_A` and
`\\mathbf{D} = \\mathbf{D}_\\text{lin} + \\mathbf{D}_\\text{n}(\\boldsymbol{\\nu})`.

`HydroParams3DOF` stores the hydrodynamic contribution in matrix form:

* `M_A`   — the 3×3 added-mass matrix.
* `D_lin` — the 3×3 linear manoeuvring damping matrix.
* `D_quad` — quadratic damping coefficients used to build
  the velocity-dependent matrix `\\mathbf{D}_\\text{n}(\\boldsymbol{\\nu})`.

In the standard Fossen 3-DOF surface craft model these matrices have
a particular structure. That structure is not required by the core
routines in this package, but it is enforced if you construct
`HydroParams3DOF` via [`hydroparams_fossen3dof`](@ref).
"""
struct HydroParams3DOF{T<:Real}
    M_A   :: SMatrix{3,3,T}              # added mass matrix
    D_lin :: SMatrix{3,3,T}              # linear damping matrix
    D_quad:: QuadraticDamping3DOF{T}     # quadratic coeffs
end


"""
$(TYPEDSIGNATURES)

Construct hydrodynamic parameters for the 3-DOF model from Fossen-style
manoeuvring derivatives.

This helper maps the usual surge–sway–yaw derivatives

* Added mass:
  `Xudot`, `Yvdot`, `Yrdot`, `Nrdot`
* Linear damping:
  `Xu`, `Yv`, `Yr`, `Nv`, `Nr` plus optional cross-terms `Xv`, `Xr`,
  `Yu`, `Nu`
* Quadratic damping:
  `Xuu`, `Yvv`, `Nrr`

to the internal matrices `M_A` and `D_lin` and the quadratic
coefficients `D_quad`.

The resulting added-mass matrix is

```math
\\mathbf{M}_A =
\\begin{bmatrix}
- X_{\\dot{u}} & 0              & 0 \\\\
0              & - Y_{\\dot{v}} & - Y_{\\dot{r}} \\\\
0              & - Y_{\\dot{r}} & - N_{\\dot{r}}
\\end{bmatrix},
```
assuming the simplified symmetric 3-DOF Fossen form,
where the off-diagonal terms are given by a single derivative ``Y_{\\dot{r}}``, so that
``\\mathbf{M}_{A,23} = \\mathbf{M}_{A,32} = -Y_{\\dot{r}}``.

The linear damping matrix is

```math
\\mathbf{D}_\\text{lin} =
-\\begin{bmatrix}
X_u & X_v & X_r \\\\
Y_u & Y_v & Y_r \\\\
N_u & N_v & N_r
\\end{bmatrix}.
```

The quadratic derivatives `Xuu`, `Yvv`, `Nrr` are stored in
[`QuadraticDamping3DOF`](@ref) and used to build the diagonal
approximation to `\\mathbf{D}_\\text{n}(\\boldsymbol{\\nu})`.

All derivatives are given with Fossen's sign convention; no sign
flips are applied to the inputs.
"""
function hydroparams_fossen3dof(
    Xudot::T, Yvdot::T, Yrdot::T, Nrdot::T,
    Xu::T, Yv::T, Yr::T, Nv::T, Nr::T,
    Xuu::T, Yvv::T, Nrr::T;
    Xv::T = zero(T), Xr::T = zero(T),
    Yu::T = zero(T), Nu::T = zero(T),
) where {T<:Real}
    z = zero(T)

    # Added mass in standard 3-DOF Fossen form
    M_A = @SMatrix [
        -Xudot  z      z;
         z     -Yvdot -Yrdot;
         z     -Yrdot -Nrdot;
    ]

    # Linear damping matrix D such that:
    # M ν̇ + C(ν)ν + D ν = τ
    # and Xu, Yv, Yr, Nv, Nr are Fossen hydrodynamic derivatives
    D_lin = @SMatrix [
        -Xu   -Xv   -Xr;
        -Yu   -Yv   -Yr;
        -Nu   -Nv   -Nr;
    ]

    # Quadratic terms: store Fossen derivatives Xuu, Yvv, Nrr (typically ≤ 0)
    D_quad = QuadraticDamping3DOF{T}(Xuu, Yvv, Nrr)

    return HydroParams3DOF{T}(M_A, D_lin, D_quad)
end


function vesselparams_fossen3dof(
    m::T, Iz::T, xG::T,
    Xudot::T, Yvdot::T, Yrdot::T, Nrdot::T,
    Xu::T,   Yv::T,   Yr::T,   Nv::T,   Nr::T,
    Xuu::T,  Yvv::T,  Nrr::T;
    Xv::T = zero(T), Xr::T = zero(T),
    Yu::T = zero(T), Nu::T = zero(T),
) where {T<:Real}
    rb = RigidBody3DOF(m, Iz, xG)

    hydro = hydroparams_fossen3dof(
        Xudot, Yvdot, Yrdot, Nrdot,
        Xu,    Yv,    Yr,    Nv,    Nr,
        Xuu,   Yvv,   Nrr;
        Xv = Xv, Xr = Xr,
        Yu = Yu, Nu = Nu,
    )

    return VesselParams3DOF(rb, hydro)
end


"""
$(TYPEDSIGNATURES)

Convenience constructor for a cached 3-DOF vessel model from
Fossen-style rigid-body and manoeuvring derivatives.

Wraps [`vesselparams_fossen3dof`](@ref) and
[`build_cached_vessel`](@ref) into a single call.

All derivatives are given with Fossen's sign convention; no sign flips
are applied.

Keyword:
- `check_physical::Bool = true`: run heuristic checks and emit
  warnings if the resulting model looks suspicious.
"""
function build_cached_vessel_fossen3dof(
    args...;
    check_physical::Bool = true,
    kwargs...
)
    params = vesselparams_fossen3dof(args...; kwargs...)
    return build_cached_vessel(params; check_physical=check_physical)
end



"""
$(TYPEDEF)

Collects all 3-DOF vessel parameters: rigid-body and hydrodynamics.

This is the main input data structure for the rest of the 3-DOF MSS.
"""
struct VesselParams3DOF{T<:Real}
    rb    :: RigidBody3DOF{T}
    hydro :: HydroParams3DOF{T}
end
