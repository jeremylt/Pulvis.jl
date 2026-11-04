# ------------------------------------------------------------------------------
# Particle data
# ------------------------------------------------------------------------------

typenames = []
typecounts = []
statedata = Dict()

"""
```julia
addparticletype!(typename, numparticles, particletypes, particlecounts)
```

Add a particle type name and number of particles

# Arguments:

  - `typename::String`:  name of new particle type
  - `numparticles::Int`: number of particles of type
  - `typenames::Array`:  input/output array of particle names
  - `typecounts::Array`: input/output array of particle counts

# Returns:

  - none

# Example:

```jldoctest
# setup
typenames = []
typecounts = []

# add new type
addparticletype!("foo", 2, typenames, typecounts)

# verify
@assert typenames[1] == "foo"
@assert typecounts[1] == 2

# add second type
addparticletype!("bar", 3, typenames, typecounts)

# verify
@assert typenames[2] == "bar"
@assert typecounts[2] == 3

# output

```
"""
function addparticletype!(
    typename::String,
    numparticles::Int,
    typenames::Array,
    typecounts::Array,
)
    @assert length(typenames) == length(typecounts)
    @assert !(typename in typenames)
    push!(typenames, typename)
    push!(typecounts, numparticles)
    return nothing
end

"""
```julia
addparticlestatedata!(
    typename,
    typenames,
    typecounts,
    statename,
    statedim,
    particledata,
    initialize,
)
```

Add new particle state for particle type

# Arguments:

  - `typename::String`:  name of particle type
  - `typenames::Array`:  input array of particle names
  - `typecounts::Array`: input array of particle counts
  - `statename::String`: name of particle state
  - `statedim::Int`:     dimension of components for particle state (typically 1 or spatial dimension)
  - `statedata::Dict`:   input/output dictionary of particle state data arrays
  - `initialize::Bool`:  flag to initialize state data array (default True)

# Returns:

  - none

# Example:

```jldoctest
# setup
typenames = []
typecounts = []
statedata = Dict()

# add new type
addparticletype!("foo", 2, typenames, typecounts)

# add second type
addparticletype!("bar", 3, typenames, typecounts)

# add particle state data
# -- first type, allocate memory
addparticlestatedata!("foo", typenames, typecounts, "baz", 1, statedata, true)

# verify
@assert length(statedata["baz"]) == 1
@assert length(statedata["baz"][1]) == 2

# -- second type, don't allocate memory
addparticlestatedata!("bar", typenames, typecounts, "baz", 1, statedata, false)

# verify
@assert length(statedata["baz"]) == 2
@assert length(statedata["baz"][1]) == 2
@assert length(statedata["baz"][2]) == 0

# -- new state, second type, allocate memory using default flag
addparticlestatedata!("bar", typenames, typecounts, "quux", 2, statedata)

# verify
@assert length(statedata["quux"]) == 2
# -- note array created for first type but no memory allocated
@assert length(statedata["quux"][1]) == 0
@assert length(statedata["quux"][2]) == 6

# output

```
"""
function addparticlestatedata!(
    typename::String,
    typenames::Array,
    typecounts::Array,
    statename::String,
    statedim::Int,
    statedata::Dict,
    initialize::Bool = true,
)
    typeindex = findfirst(typenames .== typename)

    hasstate = haskey(statedata, statename)
    statelength = hasstate ? length(statedata[statename]) : 0
    if typeindex > statelength
        newstatedata = []
        for i = 1:statelength
            push!(newstatedata, statedata[statename][i])
        end
        for i = (statelength+1):typeindex
            push!(newstatedata, PulvisVector{Float64}())
        end
        statedata[statename] =
            adapt(PulvisVector, StateData{typeindex,PulvisVector{Float64}}(newstatedata))
    end
    if initialize
        resize!(statedata[statename][typeindex], typecounts[typeindex] * statedim)
    end
    return nothing
end
