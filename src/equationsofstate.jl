# ------------------------------------------------------------------------------
# equations of state
# ------------------------------------------------------------------------------

# source files
include("equationsofstate/gentait.jl")
include("equationsofstate/idealgas.jl")

# default kernel - gentait
densitytopressure = densitytopressure_gentait;

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
setequationsofstate(EquationOfState.gentait)

# verify
@assert densitytopressure == densitytopressure_gentait

# set equations of state
setequationsofstate(EquationOfState.idealgas)

# verify
@assert densitytopressure == densitytopressure_idealgas

# output

```
"""
function setequationsofstate(equationstype::EquationOfState.type)
    if equationstype == EquationOfState.gentait
        densitytopressure = densitytopressure_gentait
    elseif equationstype == EquationOfState.idealgas
        densitytopressure = densitytopressure_idealgas
    end
end
