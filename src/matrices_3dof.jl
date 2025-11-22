# src/matrices_3dof.jl

"""
$(TYPEDSIGNATURES)

Rigid-body mass matrix `M_RB` for 3-DOF surge–sway–yaw motion.

Implements the planar rigid-body mass matrix for a surface craft with
mass `m`, yaw inertia `Iz` and centre of gravity at `(x_G, 0)` in the
body-fixed frame:

    M_RB = [  m      0      0
              0      m    m*x_G
              0    m*x_G   Iz   ]

"""
@inline function _mass_rb(rb::RigidBody3DOF{T}) where {T}
    m  = rb.m
    Iz = rb.Iz
    xG = rb.xG
    z = zero(T)

    return @SMatrix [
        m   z   z;
        z   m   m*xG;
        z   m*xG Iz
    ]
end


"""
$(TYPEDSIGNATURES)

Return the 3×3 added-mass matrix `M_A` stored in
`HydroParams3DOF`.

The core routines treat `M_A` as a general symmetric 3×3 matrix. When
constructed via [`hydroparams_fossen3dof`](@ref) it has the standard
Fossen 3-DOF structure.
"""
@inline _mass_added(h::HydroParams3DOF{T}) where {T} = h.M_A


"""
$(TYPEDSIGNATURES)

Build the 3×3 rigid-body + added-mass inertia matrix in surge–sway–yaw:

    M = M_RB + M_A

where `M_RB` is built from `RigidBody3DOF` and `M_A` from `HydroParams3DOF`.
"""
function mass_matrix(params::VesselParams3DOF{T}) where {T<:Real}
    M = _mass_rb(params.rb) + _mass_added(params.hydro)

    if cond(M) > 1e12
        @warn "Mass matrix is near singular (cond(M) = $(cond(M)))."
    end

    return M
end


"""
$(TYPEDSIGNATURES)

Return the 3×3 linear manoeuvring damping matrix `D_lin`.

In the equations of motion the linear damping term appears as

    D_lin ν,

so `D_lin` is stored with the physical sign. When constructed via
[`hydroparams_fossen3dof`](@ref) this matrix is

    D_lin =
        -[ X_u  X_v  X_r
           Y_u  Y_v  Y_r
           N_u  N_v  N_r ],

where `X_u, Y_v, ...` are Fossen-style hydrodynamic derivatives.
"""
@inline _D_lin(h::HydroParams3DOF{T}) where {T} = h.D_lin


"""
$(TYPEDSIGNATURES)

Construct the nonlinear (quadratic) damping matrix Dₙ(ν) for the 3-DOF model.

This front-end method extracts the quadratic coefficients from
`HydroParams3DOF` and delegates to the specialised method on
[`QuadraticDamping3DOF`](@ref).
"""
# Front-end: takes full HydroParams3DOF, forwards to coefficient form
@inline function _D_nl(
    ν::SVector{3,S},
    h::HydroParams3DOF{T},
) where {S<:Real, T<:Real}
    return _D_nl(ν, h.D_quad)
end


"""
$(TYPEDSIGNATURES)

Construct the diagonal approximation to the nonlinear damping matrix
D_n(ν) from Fossen-style quadratic derivatives.

Given

    ν = [u, v, r]ᵀ,

this method returns

    D_n(ν) = -diag(X_uu * |u|, Y_vv * |v|, N_rr * |r|).

The total damping matrix used in the dynamics is then

    D(ν) = D_lin + D_n(ν).
"""
# Actual 3×3 matrix, based on Fossen-style derivatives in QuadraticDamping3DOF
@inline function _D_nl(
    ν::SVector{3,S},
    q::QuadraticDamping3DOF{T},
) where {S<:Real, T<:Real}
    u, v, r = ν

    return @SMatrix [
        -q.Xuu * abs(u)    0                    0;
         0                -q.Yvv * abs(v)       0;
         0                 0                   -q.Nrr * abs(r);
    ]
end


"""
$(TYPEDSIGNATURES)

Rigid-body Coriolis/centripetal matrix `C_RB(ν)` for the 3-DOF model.

Implements the standard Fossen 3-DOF surface craft expression based on
the rigid-body parameters `m` and `xG` from [`RigidBody3DOF`](@ref).
The returned matrix satisfies

    νᵀ C_RB(ν) ν = 0

for all velocities ν, so that the rigid-body
inertial forces do not create or dissipate energy.
"""
@inline function _coriolis_rb(
    rb::RigidBody3DOF{T},
    nu::SVector{3,S},
) where {T<:Real, S<:Real}
    m  = rb.m
    xG = rb.xG
    u, v, r = nu
    z = zero(S)
    
    return @SMatrix [
        z               z        -m*(xG*r + v);
        z               z         m*u;
        m*(xG*r + v)     -m*u         z;
    ]
end


"""
$(TYPEDSIGNATURES)

Added-mass Coriolis/centripetal matrix `C_A(ν)`
for the 3-DOF model.

This implementation reconstructs the added-mass derivatives
`X_udot, Y_vdot, Y_rdot` from the stored
added-mass matrix `M_A` in [`HydroParams3DOF`](@ref) assuming the
standard Fossen 3-DOF structure

    M_A =
        [ -X_udot      0          0
           0         -Y_vdot   -Y_rdot
           0         -Y_rdot   -N_rdot ].

The returned matrix `C_A(ν)` satisfies

    νᵀ C_A(ν) ν = 0

for all velocities ν.
If `M_A` does not have the above structure the interpretation of the
entries in `C_A` will no longer match Fossen's added-mass derivatives.
"""
# Added coriolis. NB: Assumes M_A is simplified Fossen-form. Might be worth revisiting later.
@inline function _coriolis_added(
    h::HydroParams3DOF{T},
    ν::SVector{3,S},
) where {T<:Real, S<:Real}
    u, v, r = ν
    z = zero(S)

    M_A = h.M_A
    Xudot = -M_A[1,1]
    Yvdot = -M_A[2,2]
    Yrdot = -M_A[2,3]
    # Nrdot = -M_A[3,3]  # not used in 3-DOF C_A

    return @SMatrix [
        z                 z           Yvdot*v + Yrdot*r;
        z                 z          -Xudot*u;
       -Yvdot*v - Yrdot*r  Xudot*u    z;
    ]
end