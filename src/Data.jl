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

    # Loop over state names for the type
    for statename in statenames
        hasstate = haskey(statedata, statename)
        statelength = hasstate ? length(statedata[statename]) : 0
        stateindex = findfirst(s -> s.name == statename, stateinfo)
        @assert stateindex != nothing "State not known"

        # Expand data for specific state for current type, as needed
        if typelength > statelength
            newstatedata = []
            resize!(newstatedata, typelength)

            # Copy current pointers
            for j = 1:statelength
                newstatedata[j] = statedata[statename][j]
            end
            # Add new empty vectors
            for j = (statelength+1):typelength
                newstatedata[j] = PulvisVector{Float64}()
            end
            # And copy to execution space
            statedata[statename] = adapt(
                PulvisVector,
                StateData{typelength,PulvisVector{Float64}}(newstatedata),
            )
        end
        # Allocate vector only if requested
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

    # Check all states
    for i = 1:length(stateinfo)
        hasstate = haskey(statedata, stateinfo[i].name)
        statelength = hasstate ? length(statedata[stateinfo[i].name]) : 0

        # Square off if the state is missing some types
        if typelength > statelength
            newstatedata = []
            resize!(newstatedata, typelength)

            # Copy current pointers
            for j = 1:statelength
                newstatedata[j] = statedata[stateinfo[i].name][j]
            end
            # Add new empty vectors
            for j = (statelength+1):typelength
                newstatedata[j] = PulvisVector{Float64}()
            end
            # And copy to execution space
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

Remove particles by local id numbers

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

# add particle types
addparticletype!("foo", 3, typeinfo)
addparticletype!("bar", 4, typeinfo)

# add particle states
addparticlestate!("baz", 2, stateinfo)
addparticlestate!("quux", 3, stateinfo)

# add particle state data
addparticlestatedata!("foo", typeinfo, ["baz", "quux"], stateinfo, statedata, true)
addparticlestatedata!("bar", typeinfo, ["baz"], stateinfo, statedata, true)

# square off state arrays
squareparticlestatedata!(typeinfo, stateinfo, statedata)

# verify
@assert typeinfo[1].count == 3
# note that allocated lengths should not decrease
@assert length(statedata["quux"][1]) == 9
@assert length(statedata["baz"][1]) == 6

# remove particles
removeparticles!("foo", [1, 3], typeinfo, stateinfo, statedata)

# verify
@assert typeinfo[1].count == 1
# note that allocated lengths should not decrease
@assert length(statedata["quux"][1]) == 9
@assert length(statedata["baz"][1]) == 6

# verify
@assert typeinfo[2].count == 4
# note that allocated lengths should not decrease
@assert length(statedata["quux"][2]) == 0
@assert length(statedata["baz"][2]) == 8

# remove more particles
removeparticles!("bar", [2, 3], typeinfo, stateinfo, statedata)

# verify
@assert typeinfo[2].count == 2
# note that allocated lengths should not decrease
@assert length(statedata["quux"][2]) == 0
@assert length(statedata["baz"][2]) == 8

# output

```
"""
function removeparticles!(
    typename::String,
    particleids::Array,
    typeinfo::Array,
    stateinfo::Array,
    statedata::Dict,
)
    typelength = length(typeinfo)
    typeindex = findfirst(t -> t.name == typename, typeinfo)
    @assert typeindex != nothing "Type not known"
    sort(particleids, rev = true)

    # Loop over all states
    for i = 1:length(stateinfo)
        newstatedata = []
        resize!(newstatedata, typelength)

        # Copy current pointers
        for j = 1:(typeindex-1)
            newstatedata[j] = statedata[stateinfo[i].name][j]
        end
        # Drop the target particles in the target type
        data = adapt(Vector{Float64}, statedata[stateinfo[i].name][typeindex])
        numparticles = typeinfo[typeindex].count
        for id in particleids
            if id < length(data)
                for j =
                        (id*stateinfo[i].dimension):((numparticles-1)*stateinfo[i].dimension)

                    data[j] = data[j+stateinfo[i].dimension]
                end
            end
            numparticles -= 1
        end
        newstatedata[typeindex] = adapt(PulvisVector{Float64}, data)
        # Copy the rest of the current pointers
        for j = (typeindex+1):typelength
            newstatedata[j] = statedata[stateinfo[i].name][j]
        end
        # And copy to execution space
        statedata[stateinfo[i].name] =
            adapt(PulvisVector, StateData{typelength,PulvisVector{Float64}}(newstatedata))
    end
    # Finally decrease particle count
    typeinfo[typeindex] =
        TypeInfo(typeinfo[typeindex].name, typeinfo[typeindex].count - length(particleids))
    return nothing
end

"""
```julia
addparticles!(typename, numparticles, newdata, typeinfo, stateinfo, particledata)
```

Add new particle state arrays for particle type

# Arguments:

  - `typename::String`:  name of particle type to remove particles for
  - `numparticles::Int`: number of new particles to add
  - `newdata::Dict`:     input dictionary of new particle state data arrays
  - `typeinfo::Array`:   input array of particle type information
  - `stateinfo::Array`:  input array of particle state information
  - `statedata::Dict`:   input/output dictionary of particle state data arrays

# Returns:

  - none

# Example:

```jldoctest
# setup
typeinfo = []
stateinfo = []
statedata = Dict()

# add particle types
addparticletype!("foo", 3, typeinfo)
addparticletype!("bar", 4, typeinfo)

# add particle states
addparticlestate!("baz", 2, stateinfo)
addparticlestate!("quux", 3, stateinfo)

# add particle state data
addparticlestatedata!("foo", typeinfo, ["baz", "quux"], stateinfo, statedata, true)
addparticlestatedata!("bar", typeinfo, ["baz"], stateinfo, statedata, true)

# square off state arrays
squareparticlestatedata!(typeinfo, stateinfo, statedata)

# verify
@assert typeinfo[1].count == 3
@assert length(statedata["quux"][1]) == 9
@assert length(statedata["baz"][1]) == 6

# add particles
newdata = Dict()
newdata["baz"] = [1.0, 2.0]
newdata["quux"] = [1.0, 2.0, 3.0]
addparticles!("foo", 1, newdata, typeinfo, stateinfo, statedata)

# verify
@assert typeinfo[1].count == 4
@assert length(statedata["quux"][1]) == 12
@assert Vector(statedata["quux"][1])[(end-2):end] ≈ newdata["quux"]
@assert length(statedata["baz"][1]) == 8
@assert Vector(statedata["baz"][1])[(end-1):end] ≈ newdata["baz"]

# verify
@assert typeinfo[2].count == 4
@assert length(statedata["quux"][2]) == 0
@assert length(statedata["baz"][2]) == 8

# add more particles
newdata = Dict()
newdata["baz"] = [1.0, 2.0, 3.0, 4.0]
addparticles!("bar", 2, newdata, typeinfo, stateinfo, statedata)

# verify
@assert typeinfo[2].count == 6
@assert length(statedata["quux"][2]) == 0
@assert length(statedata["baz"][2]) == 12
@assert Vector(statedata["baz"][2])[(end-3):end] ≈ newdata["baz"]

# output

```
"""
function addparticles!(
    typename::String,
    numparticles::Int,
    newdata::Dict,
    typeinfo::Array,
    stateinfo::Array,
    statedata::Dict,
)
    typelength = length(typeinfo)
    typeindex = findfirst(t -> t.name == typename, typeinfo)
    @assert typeindex != nothing "Type not known"

    # Loop over all states
    for i = 1:length(stateinfo)
        # Sanity check
        @assert !haskey(newdata, stateinfo[i].name) || (
            length(newdata[stateinfo[i].name]) == numparticles * stateinfo[i].dimension
        )

        # Skip if state isn't in data to add
        if !haskey(newdata, stateinfo[i].name)
            continue
        end

        # Check if enough space is already allocated
        if length(newdata[stateinfo[i].name]) > 0
            newstatedata = []
            resize!(newstatedata, typelength)

            # Copy current pointers
            for j = 1:(typeindex-1)
                newstatedata[j] = statedata[stateinfo[i].name][j]
            end
            # Add the target particles in the target type
            data = adapt(Vector{Float64}, statedata[stateinfo[i].name][typeindex])
            if length(statedata[stateinfo[i].name]) <
               ((typeinfo[typeindex].count + numparticles) * stateinfo[i].dimension)
                resize!(
                    data,
                    (typeinfo[typeindex].count + numparticles) * stateinfo[i].dimension,
                )
            end
            for j = 1:length(newdata[stateinfo[i].name])
                data[typeinfo[typeindex].count*stateinfo[i].dimension+j] =
                    newdata[stateinfo[i].name][j]
            end
            newstatedata[typeindex] = adapt(PulvisVector{Float64}, data)
            # Copy the rest of the current pointers
            for j = (typeindex+1):typelength
                newstatedata[j] = statedata[stateinfo[i].name][j]
            end
            # And copy to execution space
            statedata[stateinfo[i].name] = adapt(
                PulvisVector,
                StateData{typelength,PulvisVector{Float64}}(newstatedata),
            )
        end
    end
    # Finally increase particle count
    typeinfo[typeindex] =
        TypeInfo(typeinfo[typeindex].name, typeinfo[typeindex].count + numparticles)
    return nothing
end
