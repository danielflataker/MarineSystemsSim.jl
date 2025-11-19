"""
$(TYPEDSIGNATURES)

Computes the full 3-DOF state derivative for the vessel model.

The state is ordered as

```math
\\mathbf{X} =
[x, y, \\psi, u, v, r]^T,
```

where

``x``, ``y`` are the Earth-fixed position coordinates,

``\\psi`` is the yaw angle,

``u``, ``v``, ``r`` are body-fixed surge, sway, and yaw rates.

We split the state into

```math
\\boldsymbol{\\eta} = [x, y, \\psi]^T, \\qquad \\boldsymbol{\\nu} = [u, v, r]^T.
```

The kinematics are

```math
\\dot{\\boldsymbol{\\eta}} = \\mathbf{R}(\\psi) \\, \\boldsymbol{\\nu},
```

and the dynamics are

```math
\\mathbf{M} \\dot{\\boldsymbol{\\nu}} +
\\mathbf{C}(\\boldsymbol{\\nu})\\,\\boldsymbol{\\nu} +
\\mathbf{D}(\\boldsymbol{\\nu})\\,\\boldsymbol{\\nu} = \\boldsymbol{\\tau},
```

so that

```math
\\dot{\\boldsymbol{\\nu}} =
\\mathbf{M}^{-1} (\\boldsymbol{\\tau} -
\\mathbf{C}(\\boldsymbol{\\nu})\\,\\boldsymbol{\\nu} -
\\mathbf{D}(\\boldsymbol{\\nu})\\,\\boldsymbol{\\nu}).
```

The input ``\\tau`` is a ``3\\times 1`` vector of generalized forces/moments in
surge, sway, and yaw.

The function returns the state derivative

```math
\\dot{\\mathbf{X}} =
[\\dot{x}, \\dot{y}, \\dot{\\psi}, \\dot{u}, \\dot{v}, \\dot{r}]^T,
```

as an `SVector{6}`.
"""