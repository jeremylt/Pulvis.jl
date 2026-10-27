# ------------------------------------------------------------------------------
# Momentum formulations
# ------------------------------------------------------------------------------

# source files
include("MomentumFormulations/Adami.jl")
include("MomentumFormulations/Monaghan.jl")

"""
```julia
setmomentumformulation(type)
```

Set the momentum formulation

# Arguments:

  - `type::MomentumFormulation.Type`: momentum formulation to use

# Returns:

  - none

# Example:

```jldoctest
# set momentum formulation
setmomentumformulation(MomentumFormulation.Adami)

# verify
@assert specificcoefficient == specificcoefficient_Adami;

# set shape function
setkernelshape(MomentumFormulation.Monaghan)

# verify
@assert specificcoefficient = specificcoefficient_Monaghan;

# output

```
"""
function setmomentumformulation(formulationtype::MomentumFormulation.Type)
    if formulationtype == MomentumFormulation.Adami
        global specificcoefficient = specificcoefficient_Adami
        global pressuregradient = pressuregradient_Adami
        global shearforces = shearforces_Adami
    elseif formulationtype == MomentumFormulation.Monaghan
        global specificcoefficient = specificcoefficient_Monaghan
        global pressuregradient = pressuregradient_Monaghan
        global shearforces = shearforces_Monaghan
    end
    return nothing
end

# default kernel - Monaghan
setmomentumformulation(MomentumFormulation.Monaghan)
