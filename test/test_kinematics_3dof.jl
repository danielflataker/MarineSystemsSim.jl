using Test
using MarineSystemsSim
using StaticArrays

@testset "3DOF kinematics: ψ = 0" begin
    ν = @SVector [1.2, -0.5, 0.3]
    η = @SVector [0.0, 0.0, 0.0]
    ηdot = kinematics(η, ν)
    @test ηdot ≈ ν
end

@testset "3DOF rotation matrix orthonormal" begin
    # Test several angles to verify orthonormality
    for ψ in (-π, -0.7, 0.0, 0.5, π/2, π)
        R = MarineSystemsSim.rotation_body_to_earth(ψ)

        I3 = @SMatrix [1.0 0.0 0.0;
                       0.0 1.0 0.0;
                       0.0 0.0 1.0]

        @test R' * R ≈ I3 atol = 1e-12 rtol = 1e-12
        @test R * R' ≈ I3 atol = 1e-12 rtol = 1e-12
    end
end

@testset "3DOF kinematics: pure surge rotated" begin
    u = 1.0
    ν = @SVector [u, 0.0, 0.0]
    ψ = 0.5
    η = @SVector [0.0, 0.0, ψ]
    ηdot = kinematics(η, ν)
    @test ηdot[1] ≈ u * cos(ψ)
    @test ηdot[2] ≈ u * sin(ψ)
    @test ηdot[3] ≈ 0.0 atol = 1e-14
end
