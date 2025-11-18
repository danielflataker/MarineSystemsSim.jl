# test/test_unitful_rigidbody_3dof.jl

using Test
using MarineSystemsSim
using Unitful

@testset "RigidBody3DOF Unitful constructor" begin
    # Base case: plain SI units
    m1  = 10_000.0u"kg"
    Iz1 = 2.0e6u"kg*m^2"
    xG1 = 1.5u"m"

    rb1 = RigidBody3DOF(m1, Iz1, xG1)

    @test rb1.m  ≈ 10_000.0
    @test rb1.Iz ≈ 2.0e6
    @test rb1.xG ≈ 1.5

    # Same physical values, different units
    m2  = 10.0u"Mg"        # 10_000 kg
    Iz2 = 2.0u"kg*km^2"    # 2 * 1e6 = 2e6 kg*m^2
    xG2 = 150.0u"cm"       # 1.5 m

    rb2 = RigidBody3DOF(m2, Iz2, xG2)

    @test rb2.m  ≈ rb1.m
    @test rb2.Iz ≈ rb1.Iz
    @test rb2.xG ≈ rb1.xG

    # Wrong dimension should throw
    @test_throws Unitful.DimensionError RigidBody3DOF(10u"m", Iz1, xG1)
end

@testset "RigidBody3DOF Unitful: integer-backed quantity" begin
    # Underlying Int in the quantity – ensure constructor handles it correctly
    m  = 10_000u"kg"          # Int64 * kg
    Iz = 2.0e6u"kg*m^2"
    xG = 2.0u"m"

    rb = RigidBody3DOF(m, Iz, xG)

    @test rb.m ≈ 10_000.0
end

@testset "RigidBody3DOF Unitful vs plain constructor" begin
    rb_plain = RigidBody3DOF(10_000.0, 2.0e6, 1.5)
    rb_unit  = RigidBody3DOF(10_000.0u"kg", 2.0e6u"kg*m^2", 1.5u"m")

    @test rb_plain.m  ≈ rb_unit.m
    @test rb_plain.Iz ≈ rb_unit.Iz
    @test rb_plain.xG ≈ rb_unit.xG
end
