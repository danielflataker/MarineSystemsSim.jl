# test/test_physical_warnings_3dof.jl

using Test
using MarineSystemsSim
using StaticArrays
using Logging

# Simple logger that collects all log messages at warn level or above
mutable struct CollectingLogger <: AbstractLogger
    records::Vector{NamedTuple}
end

Logging.min_enabled_level(::CollectingLogger) = Logging.Warn
Logging.shouldlog(::CollectingLogger, level, _mod, _group, _id) =
    level ≥ Logging.Warn

function Logging.handle_message(
    logger::CollectingLogger,
    level,
    message,
    _mod,
    _group,
    _id,
    _file,
    _line;
    kwargs...
)
    push!(logger.records, (level = level, message = string(message)))
    return
end

Logging.catch_exceptions(::CollectingLogger) = false


@testset "3DOF physical warnings" begin
    # 1) Model with "normal" parameters – should not produce any warnings
    rb_ok = RigidBody3DOF(10.0, 25.0, 1.0)

    # Fossen-style derivatives:
    # - added mass > 0
    # - linear damping derivatives (Xu, Yv, Nr) negative → D_lin diagonal ≥ 0
    # - quadratic damping derivatives negative → D_nl diagonal ≥ 0
    hydro_ok = hydroparams_fossen3dof(
        # Added mass
        1.0, 2.0, 0.3, 0.7,        # Xudot, Yvdot, Yrdot, Nrdot
        # Linear damping derivatives (Xu, Yv, Yr, Nv, Nr)
        -1.0, -1.0, 0.0, 0.0, -1.0,
        # Quadratic damping derivatives (Xuu, Yvv, Nrr)
        -0.1, -0.2, -0.3
    )

    params_ok = VesselParams3DOF(rb_ok, hydro_ok)
    model_ok  = build_cached_vessel(params_ok; check_physical = false)

    logger_ok = CollectingLogger(NamedTuple[])
    with_logger(logger_ok) do
        MarineSystemsSim._warn_if_suspicious(model_ok)
    end

    # For a physically reasonable model we expect no warnings
    @test isempty(logger_ok.records)


    # 2) Model with clear "anti-damping" – we expect warnings
    rb_bad = RigidBody3DOF(10.0, 25.0, 1.0)

    # "Bad" Fossen-style derivatives:
    # - linear damping derivatives positive → D_lin diagonal < 0 (anti-damping)
    # - quadratic derivatives positive → D_nl diagonal < 0 for |ν| > 0
    # This should trigger both the diagonal checks and the energy check.
    hydro_bad = hydroparams_fossen3dof(
        # Added mass (same as ok; not critical here)
        1.0, 2.0, 0.3, 0.7,        # Xudot, Yvdot, Yrdot, Nrdot
        # Linear damping derivatives (Xu, Yv, Yr, Nv, Nr)
        1.0, 1.0, 0.0, 0.0, 1.0,
        # Quadratic damping derivatives (Xuu, Yvv, Nrr)
        0.1, 0.2, 0.3
    )

    params_bad = VesselParams3DOF(rb_bad, hydro_bad)
    model_bad  = build_cached_vessel(params_bad; check_physical = false)

    logger_bad = CollectingLogger(NamedTuple[])
    with_logger(logger_bad) do
        MarineSystemsSim._warn_if_suspicious(model_bad)
    end

    @test !isempty(logger_bad.records)

    # We expect at least one warning about linear damping (negative diagonal of D_lin)
    @test any(occursin("Linear damping", rec.message) for rec in logger_bad.records)

    # We expect at least one warning about quadratic damping derivatives
    @test any(occursin("Quadratic damping derivative", rec.message)
              for rec in logger_bad.records)

    # And at least one warning about the energy check (damping_forces injecting energy)
    @test any(occursin("damping_forces appears to inject energy", rec.message)
              for rec in logger_bad.records)
end
