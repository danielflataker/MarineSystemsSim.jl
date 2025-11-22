# src/params_3dof.jl

using Unitful: AbstractQuantity, ustrip


"""
$(TYPEDEF)

Rigid-body parameters for 3-DOF motion in the horizontal plane.

# Fields
- `m::T`: Mass of the vessel [kg].
- `Iz::T`: Yaw moment of inertia about the z-axis [kg*m^2].
- `xG::T`: x-position of the center of gravity in the body-fixed frame [m].

# Notes

The 3-DOF rigid-body mass matrix is

    [ m      0      0
      0      m    m*xG
      0    m*xG   Iz ]
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

    D_n(ν) ν = -diag(X_uu * |u|, Y_vv * |v|, N_rr * |r|) ν.

Sign convention:

* `Xuu`, `Yvv`, `Nrr` correspond to the hydrodynamic derivatives
  `X_uu, Y_vv, N_rr` in Fossen's notation.
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

    M ν̇ + C(ν) ν + D(ν) ν = τ,

where `M = M_RB + M_A` and
`D = D_lin + D_n(ν)`.

`HydroParams3DOF` stores the hydrodynamic contribution in matrix form:

* `M_A`   — the 3×3 added-mass matrix.
* `D_lin` — the 3×3 linear manoeuvring damping matrix.
* `D_quad` — quadratic damping coefficients used to build
  the velocity-dependent matrix `D_n(ν)`.

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
$(TYPEDEF)

Collects all 3-DOF vessel parameters: rigid-body and hydrodynamics.

This is the main input data structure for the rest of the 3-DOF MSS.
"""
struct VesselParams3DOF{T<:Real}
    rb    :: RigidBody3DOF{T}
    hydro :: HydroParams3DOF{T}
end



"""
$(TYPEDSIGNATURES)

Construct hydrodynamic parameters for the 3-DOF model from Fossen-style
manoeuvring derivatives.

This helper maps the usual surge–sway–yaw derivatives

* Added mass:      `Xudot`, `Yvdot`, `Yrdot`, `Nrdot`
* Linear damping:  `Xu`, `Yv`, `Yr`, `Nv`, `Nr` plus optional cross-terms
                   `Xv`, `Xr`, `Yu`, `Nu`
* Quadratic damping: `Xuu`, `Yvv`, `Nrr`

into the internal representation `HydroParams3DOF`, i.e. the added-mass
matrix `M_A`, the linear damping matrix `D_lin`, and the quadratic
coefficients `D_quad`.

All derivatives are given with Fossen's sign convention; no sign flips
are applied.
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


"""
$(TYPEDSIGNATURES)

Convenience constructor for 3-DOF vessel parameters from Fossen-style
manoeuvring derivatives.

This wraps [`hydroparams_fossen3dof`](@ref) and [`VesselParams3DOF`](@ref)
into a single call:

- builds [`HydroParams3DOF`](@ref) from added-mass and damping derivatives, and
- returns `VesselParams3DOF(rb, hydro)`.

All derivatives use Fossen's sign convention (linear and quadratic drag
derivatives are typically ≤ 0).
"""
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
