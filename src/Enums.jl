# ------------------------------------------------------------------------------
# Enumerated types
# ------------------------------------------------------------------------------

module Kernel

"""
Kernel shape functions

# Shape:

  - `CubicSpline`:   cubic spline shape function
  - `QuinticSpline`: quintic spline shape function

# Example:

```@jldoctest
Kernel.Shape

# output
Enum Pulvis.Kernel.Shape:
CubicSpline = 0
QuinticSpline = 1
```
"""
@enum Shape CubicSpline QuinticSpline

end # submodule


module EquationOfState

"""
Equations of state

# Type:

  - `GenTait`:  general tait
  - `IdealGas`: ideal gas law

# Example:

```@jldoctest
EquationOfState.Type

# output
Enum Pulvis.EquationOfState.Type:
GenTait = 0
IdealGas = 1
```
"""
@enum Type GenTait IdealGas

end # submodule


module MomentumFormulation

"""
Momentum formulations

# Type:

  - `Monagahn`: Monagahn
  - `Adami`:    Adami

# Example:

```@jldoctest
MomentumFormulation.Type

# output
Enum Pulvis.MomentumFormulation.Type:
Monagahn = 0
Adami = 1
```
"""
@enum Type Monagahn Adami

end # submodule
