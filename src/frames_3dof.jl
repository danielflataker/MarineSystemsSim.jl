# src/frames_3dof.jl

"""
$(TYPEDSIGNATURES)

Rotation matrix from body-fixed to Earth-fixed frame in the horizontal plane.

The mapping is

```math
\\dot{\\boldsymbol{\\eta}} =
R(\\psi)\\,\\boldsymbol{\\nu}, \\quad
\\boldsymbol{\\eta} = [x, y, \\psi]^T,\\;
\\boldsymbol{\\nu} = [u, v, r]^T.
```
"""
rotation(ψ::T) where {T<:Real} = begin
    c = cos(ψ)
    s = sin(ψ)
    z = zero(T)
    o = one(T)
    @SMatrix [c  -s  z;
              s   c  z;
              z   z  o]
end


"""
$(TYPEDSIGNATURES)

Compute the kinematic relation

```math
\\dot{\\eta} = R(\\psi)\\,\\nu
```
for 3-DOF motion in the horizontal plane.

* ``\\eta = [x, y, \\psi]^T`` is the position and yaw angle,
* ``\\nu = [u, v, r]^T`` is the body-fixed velocity vector.

Returns an `SVector{3}` representing ``\\dot{\\eta}``.
"""
function kinematics(
    η::SVector{3,S},
    ν::SVector{3,S},
)::SVector{3,S} where {S<:Real}
    ψ = η[3]
    R = rotation(ψ)
    return R * ν
end