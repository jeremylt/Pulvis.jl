# ------------------------------------------------------------------------------
# Equations of state
# ------------------------------------------------------------------------------

# source files
include("EquationsOfState/GenTait.jl")
include("EquationsOfState/IdealGas.jl")

# default kernel - GenTait
densitytopressure = densitytopressure_GenTait;

"""
```julia
setequationsofstate(equationstype)
```

Set the equations of state

# Arguments:

  - `type::EquationOfState.Type`: type of equations of state to use

# Returns:

  - none

# Example:

```jldoctest
# set equations of state
setequationsofstate(EquationOfState.GenTait)

# verify
@assert densitytopressure == densitytopressure_GenTait

# set equations of state
setequationsofstate(EquationOfState.IdealGas)

# verify
@assert densitytopressure == densitytopressure_IdealGas

# output

```
"""
function setequationsofstate(equationstype::EquationOfState.Type)
    if equationstype == EquationOfState.GenTait
        global densitytopressure = densitytopressure_GenTait
    elseif equationstype == EquationOfState.IdealGas
        global densitytopressure = densitytopressure_IdealGas
    end
    return nothing
end
