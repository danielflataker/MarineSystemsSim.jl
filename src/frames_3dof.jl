# src/frames_3dof.jl
# 3DOF kinematic relations in the horizontal plane

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
rotation_body_to_earth(ψ::T) where {T<:Real} = begin
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
    eta::SVector{3,S},
    nu::SVector{3,S},
)::SVector{3,S} where {S<:Real}
    psi = eta[3]
    R = rotation_body_to_earth(psi)
    return R * nu
end