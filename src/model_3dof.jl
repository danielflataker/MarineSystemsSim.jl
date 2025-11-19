# src/model_3dof.jl

"""
$(TYPEDEF)

3-DOF vessel model with cached mass matrices.

This type stores the physical parameters together with the assembled
inertia and linear damping matrices to avoid recomputing them inside
tight simulation loops.
"""
struct Vessel3DOF{T<:Real} <: AbstractVesselModel{3,T}
    params :: VesselParams3DOF{T}
    M      :: SMatrix{3,3,T}
    Minv   :: SMatrix{3,3,T}
    D_lin  :: SMatrix{3,3,T}
end


"""
$(TYPEDSIGNATURES)

Construct a cached 3-DOF model from physical parameters.
"""
function Vessel3DOF(
    params::VesselParams3DOF{T};
    check_physical::Bool = true,
)::Vessel3DOF{T} where {T<:Real}
    M     = mass_matrix(params)
    Minv  = inv(M)
    D_lin = _D_lin(params.hydro)

    model = Vessel3DOF{T}(params, M, Minv, D_lin)

    if check_physical
        _warn_if_suspicious(model)
    end

    return model
end


@inline function _D_total(
    ν::SVector{3,S},
    model::Vessel3DOF{T},
) where {S<:Real, T<:Real}
    D_lin = model.D_lin
    D_nl  = _D_nl(ν, model.params.hydro)
    return D_lin + D_nl
end


"""
$(TYPEDSIGNATURES)

Return the constant 3×3 inertia matrix `M` for the cached 3-DOF model.
"""
mass_matrix(model::Vessel3DOF) = model.M


mass_solve(
    model::Vessel3DOF{T},
    rhs::SVector{3,S},
) where {T<:Real,S<:Real} = model.Minv * rhs


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
    nu::SVector{3,S},
    model::Vessel3DOF{T},
) where {S<:Real, T<:Real}
    return _D_total(nu, model) * nu
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
    model:: Vessel3DOF{T},
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

3-DOF kinematics `η̇ = R(ψ) ν` for the cached vessel model.

This method dispatches to the 3-DOF kinematic relation defined in
[`kinematics(η, ν)`](@ref) and exists to satisfy the interface defined
for [`AbstractVesselModel`](@ref).
"""
kinematics(eta::SVector{3}, nu::SVector{3}, ::Vessel3DOF) =
    kinematics(eta, nu)  # reuse frames_3dof version


"""
$(TYPEDSIGNATURES)

Convenience constructor for a cached 3-DOF vessel model from
Fossen-style rigid-body and manoeuvring derivatives.

Wraps [`vesselparams_fossen3dof`](@ref) and
[`Vessel3DOF`](@ref) into a single call.

All derivatives are given with Fossen's sign convention; no sign flips
are applied.

Keyword:
- `check_physical::Bool = true`: run heuristic checks and emit
  warnings if the resulting model looks suspicious.
"""
function build_vessel3dof_fossen(
    args...;
    check_physical::Bool = true,
    kwargs...,
)
    params = vesselparams_fossen3dof(args...; kwargs...)
    return Vessel3DOF(params; check_physical = check_physical)
end
