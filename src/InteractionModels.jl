# ------------------------------------------------------------------------------
# Interaction physics
# ------------------------------------------------------------------------------

# source files
include("InteractionModels/SmoothedParticleHydrodynamics.jl")

"""
```julia
setinteractionmodel(interactiontype)
```

Set the interaction model

# Arguments:

  - `interactiontype::InteractionModel.Type`: interaction model to use

# Returns:

  - nothing

# Example:

```jldoctest
# set interaction model
Pulvis.setinteractionmodel(InteractionModel.SPH)

# verify
@assert Pulvis.evaluatetimestep! == Pulvis.evaluatetimestep_SPH!;

# output

```
"""
function setinteractionmodel(interactiontype::InteractionModel.Type)
    if interactiontype == InteractionModel.SPH
        global preevaluatetimestep! = preevaluatetimestep_SPH!
        global evaluatetimestep! = evaluatetimestep_SPH!
        global postevaluatetimestep! = postevaluatetimestep_SPH!
    elseif interactiontype == InteractionModel.PD
        error("Peridynamics not yet implemented")
    elseif interactiontype == InteractionModel.DEM
        error("Discrete element method not yet implemented")
    end
    return nothing
end

# Default interaction model - SPH
setinteractionmodel(InteractionModel.SPH)
