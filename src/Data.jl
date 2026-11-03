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
addparticletype!(\"foo\", 2, typenames, typecounts)

# verify
@assert typenames[1] == \"foo\"
@assert typeof(typecounts[1]) == 2

# add second type
addparticletype!(\"baz\", 3, typenames, typecounts)

# verify
@assert typenames[2] == \"baz\"
@assert typeof(typecounts[2]) == 3

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
particleid = 1
neighborids = [3, 5, 2, 3]
dimension = 2
r = [1.0, 1.0, 1.1, 1.2, 1.3, 1.4]
radius = [1.0, 1.1, 1, 2]
m = [42.0, 13.0, 14.0]
ρ = [1.0, 0.0, 0.0]

# compute density
sumweightedmassneighbors!(particleid, neighborids, dimension, r, radius, m, ρ)

# verify
@assert ρ[particleid] ≈ 210.73260272804137
@assert ρ[neighborids[3]] ≈ 0.0
@assert ρ[neighborids[4]] ≈ 0.0

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
    initialize::Bool,
)
    typeindex = findfirst(typenames .== typename)

    hasstate = haskey(statedata, statename)
    statelength = hasstate ? length(statedata[statename]) : 0
    if typeindex > statelength
        newstatedata = []
        for i = 1:statelength
            newstatedata[i] = statedata[statename][i]
        end
        for i = (statelength+1):typeindex
            newstatedata[i] = CuVector{Float64}()
        end
        statedata[statename] = cu(SVector{typeindex,CuVector{Float64}}(newstatedata))
    end
    if initialize
        resize!(statedata[statename][typeindex], typecounts[typeindex] * statedim)
    end
    return nothing
end
