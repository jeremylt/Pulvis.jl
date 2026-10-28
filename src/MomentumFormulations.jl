# ------------------------------------------------------------------------------
# Momentum formulations
# ------------------------------------------------------------------------------

# source files
include("MomentumFormulations/Adami.jl")
include("MomentumFormulations/Monagahn.jl")

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
setkernelshape(MomentumFormulation.Monagahn)

# verify
@assert specificcoefficient = specificcoefficient_Monagahn;

# output

```
"""
function setmomentumformulation(formulationtype::MomentumFormulation.Type)
    if formulationtype == MomentumFormulation.Adami
        global specificcoefficient = specificcoefficient_Adami
        global pressuregradient! = pressuregradient_Adami!
        global shearforces! = shearforces_Adami!
        global standardbackgroundpressure! = standardbackgroundpressure_Adami!
        global generalizedbackgroundpressure! = generalizedbackgroundpressure_Adami!
        global modifiedvelocitycontribution! = modifiedvelocitycontribution_Adami!
    elseif formulationtype == MomentumFormulation.Monagahn
        global specificcoefficient = specificcoefficient_Monagahn
        global pressuregradient! = pressuregradient_Monagahn!
        global shearforces! = shearforces_Monagahn!
        global standardbackgroundpressure! = standardbackgroundpressure_Monagahn!
        global generalizedbackgroundpressure! = generalizedbackgroundpressure_Monagahn!
        global modifiedvelocitycontribution! = modifiedvelocitycontribution_Monagahn!
    end
    return nothing
end

# default kernel - Monagahn
setmomentumformulation(MomentumFormulation.Monagahn)
