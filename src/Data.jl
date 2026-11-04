# ------------------------------------------------------------------------------
# Particle data
# ------------------------------------------------------------------------------

typeinfo = []
stateinfo = []
statedata = Dict()

"""
```julia
addparticletype!(name, typeinfo)
```

Add a particle type name

# Arguments:

  - `name::String`:    name of new particle type
  - `count::Int`:      number of particles by type
  - `typeinfo::Array`: input/output array of particle type information

# Returns:

  - none

# Example:

```jldoctest
# setup
typeinfo = []

# add new type
addparticletype!("foo", 3, typeinfo)

# verify
@assert typeinfo[1].name == "foo"
@assert typeinfo[1].count == 3

# add second type
addparticletype!("bar", 4, typeinfo)

# verify
@assert typeinfo[2].name == "bar"
@assert typeinfo[2].count == 4

# output

```
"""
function addparticletype!(name::String, count::Int, typeinfo::Array)
    @assert findfirst(t -> t.name == name, typeinfo) == nothing "Types must have unique names"
    push!(typeinfo, TypeInfo(name, count))
    return nothing
end

"""
```julia
addparticlestate!(name, dimension, stateinfo)
```

Add a particle state name and dimension

# Arguments:

  - `name::String`:     name of new particle state
  - `dimension::Int`:   dimension of new particle state
  - `stateinfo::Array`: input/output array of particle state information

# Returns:

  - none

# Example:

```jldoctest
# setup
stateinfo = []

# add new state
addparticlestate!("baz", 2, stateinfo)

# verify
@assert stateinfo[1].name == "baz"
@assert stateinfo[1].dimension == 2

# add second state
addparticlestate!("quux", 3, stateinfo)

# verify
@assert stateinfo[2].name == "quux"
@assert stateinfo[2].dimension == 3

# output

```
"""
function addparticlestate!(name::String, dimension::Int, stateinfo::Array)
    @assert findfirst(s -> s.name == name, stateinfo) == nothing "States must have unique names"
    push!(stateinfo, StateInfo(name, dimension))
    return nothing
end

"""
```julia
addparticlestatedata!(typename, typeinfo, statenames, stateinfo, particledata, initialize)
```

Add new particle state arrays for particle type

# Arguments:

  - `typename::String`:  name of particle type to add state data for
  - `typeinfo::Array`:   input array of particle type information
  - `statenames::Array`: input array of names of particle states to add
  - `stateinfo::Array`:  input array of particle state information
  - `statedata::Dict`:   input/output dictionary of particle state data arrays
  - `initialize::Bool`:  flag to initialize state data array (default True)

# Returns:

  - none

# Example:

```jldoctest
# setup
typeinfo = []
stateinfo = []
statedata = Dict()

# add new types
addparticletype!("foo", 3, typeinfo)
addparticletype!("bar", 4, typeinfo)

# add new states
addparticlestate!("baz", 2, stateinfo)
addparticlestate!("quux", 3, stateinfo)

# add particle state data
# -- add state for second type, allocate memory
addparticlestatedata!("bar", typeinfo, ["baz"], stateinfo, statedata, true)

# verify
@assert length(statedata["baz"]) == 2
# -- note array created for first type but no memory allocated
@assert length(statedata["baz"][1]) == 0
@assert length(statedata["baz"][2]) == 8

# -- add states for first type, don't allocate memory
addparticlestatedata!("foo", typeinfo, ["baz", "quux"], stateinfo, statedata, false)

# verify
@assert length(statedata["baz"]) == 2
@assert length(statedata["baz"][1]) == 0
@assert length(statedata["quux"]) == 2
@assert length(statedata["quux"][1]) == 0
@assert length(statedata["quux"][2]) == 0

# output

```
"""
function addparticlestatedata!(
    typename::String,
    typeinfo::Array,
    statenames::Array,
    stateinfo::Array,
    statedata::Dict,
    initialize::Bool = true,
)
    typelength = length(typeinfo)
    typeindex = findfirst(t -> t.name == typename, typeinfo)
    @assert typeindex != nothing "Type not known"


    for statename in statenames
        hasstate = haskey(statedata, statename)
        statelength = hasstate ? length(statedata[statename]) : 0
        stateindex = findfirst(s -> s.name == statename, stateinfo)
        @assert stateindex != nothing "State not known"

        if typelength > statelength
            newstatedata = []
            for i = 1:statelength
                push!(newstatedata, statedata[statename][i])
            end
            for i = (statelength+1):typelength
                push!(newstatedata, PulvisVector{Float64}())
            end
            statedata[statename] = adapt(
                PulvisVector,
                StateData{typelength,PulvisVector{Float64}}(newstatedata),
            )
        end
        if initialize
            resize!(
                statedata[statename][typeindex],
                typeinfo[typeindex].count * stateinfo[stateindex].dimension,
            )
        end
    end
    return nothing
end

"""
```julia
squareparticlestatedata!(typeinfo, particledata)
```

Ensure all state data arrays have the same length

# Arguments:

  - `typeinfo::Array`:  input array of particle type information
  - `stateinfo::Array`: input array of particle state information
  - `statedata::Dict`:  input/output dictionary of particle state data arrays

# Returns:

  - none

# Example:

```jldoctest
# setup
typeinfo = []
stateinfo = []
statedata = Dict()

# add first type
addparticletype!("foo", 3, typeinfo)

# add new states
addparticlestate!("baz", 2, stateinfo)
addparticlestate!("quux", 3, stateinfo)

# add particle state data
# -- add states for first type, don't allocate memory
addparticlestatedata!("foo", typeinfo, ["baz", "quux"], stateinfo, statedata, false)

# add new type
addparticletype!("bar", 4, typeinfo)

# add more particle state data
# -- add state for second type, allocate memory
addparticlestatedata!("bar", typeinfo, ["baz"], stateinfo, statedata, true)

# verify
@assert length(statedata["baz"]) == 2
@assert length(statedata["quux"]) == 1

# square off state arrays
squareparticlestatedata!(typeinfo, stateinfo, statedata)

# verify
@assert length(statedata["baz"]) == 2
@assert length(statedata["quux"]) == 2

# output

```
"""
function squareparticlestatedata!(typeinfo::Array, stateinfo::Array, statedata::Dict)
    typelength = length(typeinfo)

    for i = 1:length(stateinfo)
        hasstate = haskey(statedata, stateinfo[i].name)
        statelength = hasstate ? length(statedata[stateinfo[i].name]) : 0

        if typelength > statelength
            newstatedata = []
            for i = 1:statelength
                push!(newstatedata, statedata[stateinfo[i].name][i])
            end
            for i = (statelength+1):typelength
                push!(newstatedata, PulvisVector{Float64}())
            end
            statedata[stateinfo[i].name] = adapt(
                PulvisVector,
                StateData{typelength,PulvisVector{Float64}}(newstatedata),
            )
        end
    end
    return nothing
end

"""
```julia
removeparticles!(typename, particleids, typeinfo, stateinfo, particledata)
```

Add new particle state arrays for particle type

# Arguments:

  - `typename::String`:   name of particle type to remove particles for
  - `particleids::Array`: input array of particle local id numbers for particles to remove
  - `typeinfo::Array`:    input array of particle type information
  - `stateinfo::Array`:   input array of particle state information
  - `statedata::Dict`:    input/output dictionary of particle state data arrays

# Returns:

  - none

# Example:

```jldoctest
# setup
typeinfo = []
stateinfo = []
statedata = Dict()

# add new types
addparticletype!("foo", 3, typeinfo)
addparticletype!("bar", 4, typeinfo)

# add new states
addparticlestate!("baz", 2, stateinfo)
addparticlestate!("quux", 3, stateinfo)

# add particle state data
addparticlestatedata!("bar", typeinfo, ["baz", "quux"], stateinfo, statedata, true)
addparticlestatedata!("foo", typeinfo, ["baz", "quux"], stateinfo, statedata, false)

# remove particles
removeparticles!("foo", [1, 3], typeinfo, stateinfo, statedata)

# verify
@assert typeinfo[1].count == 1
@assert length(statedata[1]["baz"]) == 1
@assert length(statedata[1]["quux"]) == 1

# remove more particles
removeparticles!("bar", [2, 3], typeinfo, stateinfo, statedata)

# verify
@assert typeinfo[1].count == 2
@assert length(statedata[1]["baz"]) == 2
@assert length(statedata[1]["quux"]) == 2

# output

```
"""
function addparticlestatedata!(
    typename::String,
    particleids::Array,
    typeinfo::Array,
    stateinfo::Array,
    statedata::Dict,
)
    return nothing
end
