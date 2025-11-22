# src/core.jl

"""
$(TYPEDEF)

Abstract supertype for vessel models with `N` degrees of freedom and
scalar type `T`.

Concrete subtypes represent a specific dynamic vessel model (for
example a cached 3-DOF or 6-DOF model) and are expected to implement
at least the following interface:

- [`mass_matrix(model::AbstractVesselModel)`](@ref): return the generalized mass matrix `M`.
- [`damping_forces(nu, model::AbstractVesselModel)`](@ref): return the damping term `D(ν) ν`.
- [`coriolis_forces(nu, model::AbstractVesselModel)`](@ref): return the Coriolis/centripetal term `C(ν) ν`.
- [`kinematics(eta, nu, model::AbstractVesselModel)`](@ref): return the kinematic relation `η̇`.

The high-level routines [`body_dynamics`](@ref),
[`vessel_dynamics`](@ref), and [`vessel_rhs!`](@ref) use these
functions to assemble the equations of motion in a DOF-agnostic way.
"""
abstract type AbstractVesselModel{N,T} end


"""
Return the number of degrees of freedom for a vessel model.
"""
@inline dofs(::AbstractVesselModel{N}) where {N} = N


@noinline function _interface_error(f, model)
    error("$(f) is not implemented for $(typeof(model)). " *
          "Define an appropriate method $(f)(..., model::$(typeof(model))) " *
          "for this vessel model type.")
end


"""
Return the generalized mass matrix for the given vessel model `model`.

Concrete subtypes of [`AbstractVesselModel`](@ref) must implement this
method.
"""
mass_matrix(model::AbstractVesselModel{N,T}) where {N,T} =
    _interface_error(mass_matrix, model)


"""
Return the damping term D(ν) ν for the given vessel model and
body-fixed velocity ν.

Concrete subtypes of [`AbstractVesselModel`](@ref) must implement this
method.
"""
damping_forces(
    ν::SVector{N,S},
    model::AbstractVesselModel{N,T},
) where {N,S<:Real,T} =
    _interface_error(damping_forces, model)


"""
Return the Coriolis/centripetal term C(ν) ν for the given vessel
model and body-fixed velocity ν.

Concrete subtypes of [`AbstractVesselModel`](@ref) must implement this
method.
"""
coriolis_forces(
    ν::SVector{N,S},
    model::AbstractVesselModel{N,T},
) where {N,S<:Real,T} =
    _interface_error(coriolis_forces, model)


"""
Return the kinematic relation η̇ for the given vessel model, based on
position/orientation η and body-fixed velocity ν.

Concrete subtypes of [`AbstractVesselModel`](@ref) must implement this
method.
"""
kinematics(
    η::SVector{N,S},
    ν::SVector{N,S},
    model::AbstractVesselModel{N,T},
) where {N,S<:Real,T} =
    _interface_error(kinematics, model)


"""
Solve M ν̇ = rhs for ν̇, using the mass matrix associated with the
vessel model.

The default implementation calls `mass_matrix(model) \\ rhs`. Concrete
vessel models can override this method to use cached inverses or
factorisations for efficiency.
"""
function mass_solve(
    model::AbstractVesselModel{N,T},
    rhs::SVector{N,S},
)::SVector{N,promote_type(T,S)} where {N,T,S<:Real}
    M = mass_matrix(model)
    return M \ rhs
end


"""
Compute the body-fixed acceleration ν̇ for a vessel model given the
generalised forces τ.

The default implementation uses the interface functions
[`damping_forces`](@ref), [`coriolis_forces`](@ref) and
[`mass_solve`](@ref):

    M ν̇ + C(ν) ν + D(ν) ν = τ
    ν̇ = M⁻¹ (τ - C(ν) ν - D(ν) ν).
"""
function body_dynamics(ν, model::AbstractVesselModel, τ)
    Dν  = damping_forces(ν, model)
    Cν  = coriolis_forces(ν, model)
    rhs = τ - Cν - Dν
    return mass_solve(model, rhs)
end

# Type-stable specialisation for static vectors with matching scalar type
function body_dynamics(
    ν::SVector{N,S},
    model::AbstractVesselModel{N,T},
    τ::SVector{N,S},
)::SVector{N,S} where {N,S<:Real,T<:Real}
    Dν  = damping_forces(ν, model)
    Cν  = coriolis_forces(ν, model)
    rhs = τ - Cν - Dν
    return mass_solve(model, rhs)
end


"""
Compute the full state derivative `[η̇; ν̇]` for a vessel model.

The state vector is ordered as

    X = [η; ν],

where both `η` and `ν` have length `N = dofs(model)`.

The default implementation splits `X`, calls [`kinematics`](@ref) and
[`body_dynamics`](@ref), and concatenates the results. Concrete vessel
models can provide more specialised methods if needed.
"""
function vessel_dynamics(
    X::SVector{M,SX},
    model::AbstractVesselModel{N,T},
    τ::SVector{N,ST},
) where {M,N,SX<:Real,ST<:Real,T}
    @assert M == 2N "State vector X must have length 2N for N-DOF vessel models."
    
    # Promote element types so X and τ live in the same numeric world
    S = promote_type(SX, ST)
    
    η = SVector{N,S}(ntuple(i -> convert(S, X[i]), N))
    ν = SVector{N,S}(ntuple(i -> convert(S, X[N + i]), N))

    τp = SVector{N,S}(τ)

    ηdot = kinematics(η, ν, model)
    νdot = body_dynamics(ν, model, τ)
    
    return vcat(ηdot, νdot)
end


"""
In-place ODE right-hand side for a vessel model.

Evaluates the generalized forces `τ = τfun(X, t)` and computes the
state derivative `dX` using [`vessel_dynamics`](@ref). This function
is intended for use with time integrators such as DifferentialEquations.jl.
"""
function vessel_rhs!(
    dX,
    X,
    model::AbstractVesselModel{N,T},
    τfun,
    t,
) where {N,T}
    @assert length(X) == 2N "State vector X must have length 2N for N-DOF vessel models."

    # Convert mutable state to SVector for fast math / AD
    Xs = SVector{2N, eltype(X)}(Tuple(X))

    τ   = τfun(Xs, t)
    dXs = vessel_dynamics(Xs, model, τ)

    @. dX = dXs
    return nothing
end

