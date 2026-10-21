# ------------------------------------------------------------------------------
# enumerated types
# ------------------------------------------------------------------------------

module Kernel

"""
Kernel shape functions

# Shape:

  - `cubicspline`:   cubic spline shape function
  - `quinticspline`: quintic spline shape function

# Example:

```@jldoctest
Kernel.Shape

# output
Enum Pulvis.Kernel.Shape:
cubicspline = 0
quinticspline = 1
```
"""
@enum Shape cubicspline quinticspline

end # submodule


module EquationOfState

"""
Equations of state

# Type:

  - `gentait`:  general tait
  - `idealgas`: ideal gas law

# Example:

```@jldoctest
EquationOfState.Type

# output
Enum Pulvis.EquationOfState.Type:
gentait = 0
idealgas = 1
```
"""
@enum Type gentait idealgas

end # submodule
