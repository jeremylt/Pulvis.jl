# ------------------------------------------------------------------------------
# Equations of state
# ------------------------------------------------------------------------------

# source files
include("EquationsOfState/GenTait.jl")
include("EquationsOfState/IdealGas.jl")

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
        global pressuretodensity = pressuretodensity_GenTait
        global densitytoenergy = densitytoenergy_GenTait
    elseif equationstype == EquationOfState.IdealGas
        global densitytopressure = densitytopressure_IdealGas
        global pressuretodensity = pressuretodensity_IdealGas
        global densitytoenergy = densitytoenergy_IdealGas
    end
    return nothing
end

# default kernel - GenTait
setequationsofstate(EquationOfState.GenTait)
