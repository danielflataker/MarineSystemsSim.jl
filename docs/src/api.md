# API reference

```@meta
CurrentModule = MarineSystemsSim
```

This page lists the main public types and functions, grouped by purpose.
See the other pages for motivation, theory, and examples.

---

## Types

```@docs
RigidBody3DOF
QuadraticDamping3DOF
HydroParams3DOF
VesselParams3DOF
CachedVessel3DOF
```

---

## Fossen-style helpers

```@docs
hydroparams_fossen3dof
```

---

## Kinematics and frames

```@docs
rotation
kinematics
```

---

## Core dynamics

```@docs
mass_matrix
damping_forces
coriolis_forces
body_dynamics
vessel_dynamics
build_cached_vessel
```

If additional modules (e.g. for thrusters, controllers, or 6-DOF
extensions) are added later, their public API can be listed in separate
sections here.
