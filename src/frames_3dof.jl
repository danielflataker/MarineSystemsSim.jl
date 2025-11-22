# src/frames_3dof.jl
# 3DOF kinematic relations in the horizontal plane

"""
$(TYPEDSIGNATURES)

Rotation matrix from body-fixed to Earth-fixed frame in the horizontal plane.

The mapping is

    η̇ = R(ψ) ν
    η = [x, y, ψ]ᵀ
    ν = [u, v, r]ᵀ.
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

    η̇ = R(ψ) ν

for 3-DOF motion in the horizontal plane.

* η = [x, y, ψ]ᵀ is the position and yaw angle,
* ν = [u, v, r]ᵀ is the body-fixed velocity vector.

Returns an `SVector{3}` representing η̇.
"""
function kinematics(
    eta::SVector{3,S},
    nu::SVector{3,S},
)::SVector{3,S} where {S<:Real}
    psi = eta[3]
    R = rotation_body_to_earth(psi)
    return R * nu
end