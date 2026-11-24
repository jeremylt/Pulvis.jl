# ------------------------------------------------------------------------------
# Simulation setup
# ------------------------------------------------------------------------------

"""
```julia
timeloop!(
    particleids,
    typeids,
    neighborindexes,
    dimension,
    bininfo,
    typeinfo,
    stateinfo,
    statedata,
)
```

Evaluate time step; compute density, convert density to pressure, then various forces, acceleration, etc from density and pressure

# Arguments:

  - `particleids::Array{Int}`:     input/output array of local particle ids by type
  - `typeids::Array{Int}`:         input/output array of particle type ids
  - `neighborindexes::Array{Int}`: input/output array of neighbor particle global indexes
  - `dimension::Int`:              spatial dimension of problem
  - `bininfo::BinInfo`:            info about background bin layout
  - `timeinfo::Array{TimeInfo}`:   input/output array with time object
  - `typeinfo::Array{TypeInfo}`:   input/output array of particle type information
  - `stateinfo::Array{StateInfo}`: input array of particle state information
  - `statedata::Dict`:             input/output dictionary of particle state data arrays

# Returns:

  - nothing
"""
function setup!(
    particleids::Array{Int},
    typeids::Array{Int},
    neighborindexes::Array{Int},
    bininfo::BinInfo,
    timeinfo::Array{TimeInfo},
    typeinfo::Array{TypeInfo},
    stateinfo::Array{StateInfo},
    statedata::Dict,
)
    # Setup particles and materials
    if isrestart
        # loadfromrestart(timeinfo, typeinfo, stateinfo, statedata) # TODO
    else
        generateinitialparticles!(timeinfo, typeinfo, stateinfo, statedata)
    end

    # Setup time integrator
    setuptimeintegrator!()

    return nothing
end
