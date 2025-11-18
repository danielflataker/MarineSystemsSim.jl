# checks.jl

"""
$(TYPEDSIGNATURES)

Emit warnings if the hydrodynamic parameters in `model` look suspicious.

This helper performs three classes of checks:

1. **Linear damping diagonal**:
   The diagonal entries of `D_lin` are expected to be non-negative
   (dissipative behaviour) with the sign convention used in this
   package. Negative values may indicate that the Fossen derivatives
   were given with the wrong sign.

2. **Quadratic damping derivatives**:
   The Fossen-style derivatives `Xuu`, `Yvv`, `Nrr` stored in
   `D_quad` are often non-positive so that the resulting quadratic
   drag opposes the motion. Large positive values may indicate a
   sign error.

3. **Energy-based check**:
   For a small set of test velocities the power
   ``\\boldsymbol{\\nu}^T \\mathbf{D}(\\boldsymbol{\\nu})\\,\\boldsymbol{\\nu}``
   is evaluated. If this is significantly negative for any of the
   test points, a warning is issued because the damping appears to
   inject energy.

These checks are heuristic and do not enforce any hard constraints,
but they can help catch sign mistakes and obviously unphysical
parameter sets early.
"""
function _warn_if_suspicious(model::CachedVessel3DOF{T}) where {T<:Real}
    h = model.params.hydro

    # === Linear damping: diagonal entries should usually be >= 0 ===
    D = h.D_lin
    diag_lin = [
        ("D_lin[1,1]", D[1,1]),
        ("D_lin[2,2]", D[2,2]),
        ("D_lin[3,3]", D[3,3]),
    ]

    for (name, val) in diag_lin
        if val < zero(T)
            @warn "Linear damping $name = $val is negative. " *
                  "With the current sign convention, the diagonal of D_lin " *
                  "is normally ≥ 0 for dissipative damping."
        end
    end

    # === Quadratic damping: Fossen-style hydrodynamic derivatives (usually <= 0) ===
    q = h.D_quad
    diag_quad = [
        (:Xuu, q.Xuu),
        (:Yvv, q.Yvv),
        (:Nrr, q.Nrr),
    ]

    for (name, val) in diag_quad
        if val > zero(T)
            @warn "Quadratic damping derivative $name = $val is positive. " *
                  "With the current Fossen-style sign convention, these are " *
                  "normally ≤ 0 for dissipative drag."
        end
    end

    # === Energy-based check: does D ever inject energy? ===
    # We want νᵀ D_total(ν) ν >= 0 for typical velocities.
    test_vals = (T(-1), T(-0.3), T(0.3), T(1))
    for u in test_vals, v in test_vals, r in test_vals
        ν = SVector{3,T}(u, v, r)
        Dν = damping_forces(ν, model)  # this is D_total(ν) ν
        power = dot(ν, Dν)             # νᵀ D ν

        if power < -sqrt(eps(T))  # allow tiny negative from roundoff
            @warn "damping_forces appears to inject energy: νᵀ D ν = $power " *
                  "for ν = $ν. Check sign conventions and magnitudes in D_lin " *
                  "and D_quad."
            break
        end
    end

    return nothing
end
