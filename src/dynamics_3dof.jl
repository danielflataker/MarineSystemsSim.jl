# src/dynamics_3dof.jl

"""
$(TYPEDEF)

Cached representation of a 3-DOF vessel model.

This type stores the physical parameters together with the assembled
inertia and linear damping matrices to avoid recomputing them inside
tight simulation loops.

# Fields

- `params::VesselParams3DOF{T}`:
  Rigid-body and hydrodynamic parameters.
- `M::SMatrix{3,3,T}`:
  Constant 3×3 inertia matrix ``\\mathbf{M} = \\mathbf{M}_{RB} + \\mathbf{M}_A``.
- `Minv::SMatrix{3,3,T}`:
  Precomputed inverse of `M`.
- `D_lin::SMatrix{3,3,T}`:
  Copy of the linear damping matrix `params.hydro.D_lin`, used as the
  velocity-independent part of ``\\mathbf{D}(\\boldsymbol{\\nu})``.
"""
struct CachedVessel3DOF{T<:Real}
    params :: VesselParams3DOF{T}
    M      :: SMatrix{3,3,T}
    Minv   :: SMatrix{3,3,T}
    D_lin  :: SMatrix{3,3,T}
end



"""
$(TYPEDSIGNATURES)

Build the ``3\\times 3`` rigid-body + added-mass inertia matrix in surge–sway–yaw:

```math
\\mathbf{M} = \\mathbf{M}_{RB} + \\mathbf{M}_A
```

where ``M_{RB}`` is built from `RigidBody3DOF` and ``M_A`` from `HydroParams3DOF`.
"""
function mass_matrix(params::VesselParams3DOF{T}) where {T<:Real}
    M = _mass_rb(params.rb) + _mass_added(params.hydro)

    if cond(M) > 1e12
        @warn "Mass matrix is near singular (cond(M) = $(cond(M)))."
    end

    return M
end

# Should we have this, so it works for non-cached vessels too?
#@inline _D_total(ν, h) = _D_lin(h) + _D_nl(ν, h)

@inline function _D_total(
    ν::SVector{3,S},
    model::CachedVessel3DOF{T},
) where {S<:Real, T<:Real}
    D_lin = model.D_lin
    D_nl  = _D_nl(ν, model.params.hydro)
    return D_lin + D_nl
end

"""
$(TYPEDSIGNATURES)

Construct a cached 3-DOF model from physical parameters.
"""
function build_cached_vessel(
    params::VesselParams3DOF{T};
    check_physical::Bool = true,
)::CachedVessel3DOF{T} where {T<:Real}
    M     = mass_matrix(params)
    Minv  = inv(M)
    D_lin = _D_lin(params.hydro)

    model = CachedVessel3DOF{T}(params, M, Minv, D_lin)

    if check_physical
        _warn_if_suspicious(model)
    end

    return model
end


"""
$(TYPEDSIGNATURES)

Compute hydrodynamic damping forces/moments in surge, sway and yaw.

Given body-fixed velocity

```math
\\boldsymbol{\\nu} = [u, v, r]^T,
```

this function evaluates the total damping matrix

```math
\\mathbf{D}(\\boldsymbol{\\nu}) =
\\mathbf{D}_\\text{lin} + \\mathbf{D}_\\text{n}(\\boldsymbol{\\nu}),
```

and returns the product

```math
\\mathbf{D}(\\boldsymbol{\\nu})\\,\\boldsymbol{\\nu}.
```

Here `D_lin` is the cached linear damping matrix from the model and
`\\mathbf{D}_\\text{n}(\\boldsymbol{\\nu})` is constructed from the
quadratic coefficients stored in `model.params.hydro.D_quad`.

Returns an `SVector{3}` containing the damping forces/moment
`[X_d, Y_d, N_d]^T`.
"""
@inline function damping_forces(
    ν::SVector{3,S},
    model::CachedVessel3DOF{T},
) where {S<:Real, T<:Real}
    return _D_total(ν, model) * ν
end


"""
$(TYPEDSIGNATURES)

Compute the Coriolis/centripetal term ``\\mathbf{C}(\\boldsymbol{\\nu})\\,\\boldsymbol{\\nu}``
for the 3-DOF model.

The implementation constructs both the rigid-body part
``\\mathbf{C}_{RB}(\\boldsymbol{\\nu})`` and the added-mass part
``\\mathbf{C}_A(\\boldsymbol{\\nu})`` as ``3\\times 3`` matrices following
the standard Fossen 3-DOF surface craft model, and returns

```math
\\mathbf{C}(\\boldsymbol{\\nu})\\,\\boldsymbol{\\nu}
=
\\bigl( \\mathbf{C}_{RB}(\\boldsymbol{\\nu})
      + \\mathbf{C}_A(\\boldsymbol{\\nu}) \\bigr) \\, \\boldsymbol{\\nu}.
```

Returns an `SVector{3}` representing ``\\mathbf{C}(\\boldsymbol{\\nu}) \\, \\boldsymbol{\\nu}``.
"""
function coriolis_forces(
    nu    :: SVector{3,S},
    model:: CachedVessel3DOF{T},
)::SVector{3,S} where {S<:Real, T<:Real}
    rb = model.params.rb
    h  = model.params.hydro

    C_RB = _coriolis_rb(rb, nu)
    C_A  = _coriolis_added(h, nu)
    C    = C_RB + C_A

    return C * nu
end


"""
$(TYPEDSIGNATURES)

Compute the body-fixed acceleration ``\\dot{\\nu}`` for a 3-DOF vessel model given
the generalized forces ``\\tau``.

The dynamics are

```math
\\mathbf{M} \\dot{\\boldsymbol{\\nu}} +
\\mathbf{C}(\\boldsymbol{\\nu})\\,\\boldsymbol{\\nu} +
\\mathbf{D}(\\boldsymbol{\\nu})\\,\\boldsymbol{\\nu} = \\boldsymbol{\\tau},
```

so that

```math
\\dot{\\boldsymbol{\\nu}} =
\\mathbf{M}^{-1} (\\boldsymbol{\\tau} -
\\mathbf{C}(\\boldsymbol{\\nu})\\,\\boldsymbol{\\nu} -
\\mathbf{D}(\\boldsymbol{\\nu})\\,\\boldsymbol{\\nu}).
```

This version neglects environmental loads (they must be included in τ
by the caller if needed).
"""
function body_dynamics(
    nu::SVector{3,S},
    model::CachedVessel3DOF{T},
    tau::SVector{3,R}
)::SVector{3,S} where {S<:Real, R<:Real, T<:Real}
    Dnu = damping_forces(nu, model)
    Cnu = coriolis_forces(nu, model)
    return model.Minv * (tau - Cnu - Dnu)
end


"""
$(TYPEDSIGNATURES)

Computes the full 3-DOF state derivative for the vessel model.

The state is ordered as

```math
\\mathbf{X} =
[x, y, \\psi, u, v, r]^T,
```

where

``x``, ``y`` are the Earth-fixed position coordinates,

``\\psi`` is the yaw angle,

``u``, ``v``, ``r`` are body-fixed surge, sway, and yaw rates.

We split the state into

```math
\\boldsymbol{\\eta} = [x, y, \\psi]^T, \\qquad \\boldsymbol{\\nu} = [u, v, r]^T.
```

The kinematics are

```math
\\dot{\\boldsymbol{\\eta}} = \\mathbf{R}(\\psi) \\, \\boldsymbol{\\nu},
```

and the dynamics are

```math
\\mathbf{M} \\dot{\\boldsymbol{\\nu}} +
\\mathbf{C}(\\boldsymbol{\\nu})\\,\\boldsymbol{\\nu} +
\\mathbf{D}(\\boldsymbol{\\nu})\\,\\boldsymbol{\\nu} = \\boldsymbol{\\tau},
```

so that

```math
\\dot{\\boldsymbol{\\nu}} =
\\mathbf{M}^{-1} (\\boldsymbol{\\tau} -
\\mathbf{C}(\\boldsymbol{\\nu})\\,\\boldsymbol{\\nu} -
\\mathbf{D}(\\boldsymbol{\\nu})\\,\\boldsymbol{\\nu}).
```

The input ``\\tau`` is a ``3\\times 1`` vector of generalized forces/moments in
surge, sway, and yaw.

The function returns the state derivative

```math
\\dot{\\mathbf{X}} =
[\\dot{x}, \\dot{y}, \\dot{\\psi}, \\dot{u}, \\dot{v}, \\dot{r}]^T,
```

as an `SVector{6}`.
"""
function vessel_dynamics(
    X::SVector{6,S},
    model::CachedVessel3DOF{T},
    tau::SVector{3,R},
)::SVector{6,S} where {S<:Real, R<:Real, T<:Real}
    eta = @SVector [X[1], X[2], X[3]]
    nu = @SVector [X[4], X[5], X[6]]

    etadot = kinematics(eta, nu)
    nudot = body_dynamics(nu, model, tau)

    return @SVector [
        etadot[1], etadot[2], etadot[3],
        nudot[1], nudot[2], nudot[3],
    ]
end


function vessel_rhs_3dof!(
    dX,
    X,
    model::CachedVessel3DOF,
    τfun,
    t
)
    Xs = SVector{6}(X)
    τ  = τfun(Xs, t)
    dXs = vessel_dynamics(Xs, model, τ)
    @. dX = dXs
end