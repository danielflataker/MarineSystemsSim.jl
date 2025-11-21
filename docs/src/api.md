# API reference

```@meta
CurrentModule = MarineSystemsSim
```

This page lists the main public types and functions, grouped by purpose.
See the other pages for motivation, theory, and examples.

---

## Types

```@docs
AbstractVesselModel
RigidBody3DOF
QuadraticDamping3DOF
HydroParams3DOF
VesselParams3DOF
Vessel3DOF
```

---

## Fossen-style helpers

```@docs
hydroparams_fossen3dof
vesselparams_fossen3dof
build_vessel3dof_fossen
```

---

## Kinematics and frames

```@docs
rotation_body_to_earth
kinematics
```

---

## Core dynamics and ODE helpers

```@docs
dofs
mass_matrix
damping_forces
coriolis_forces
body_dynamics
vessel_dynamics
vessel_rhs!
```

If additional modules (e.g. for thrusters, controllers, or 6-DOF
extensions) are added later, their public API can be listed in separate
sections here.
