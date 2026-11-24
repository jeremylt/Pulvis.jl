# ------------------------------------------------------------------------------
# Particle data
# ------------------------------------------------------------------------------

typeinfo::Array{TypeInfo} = []
stateinfo::Array{StateInfo} = []
statedata::Dict = Dict()
particleids::AbstractArray{Int} = []
typeids::AbstractArray{Int} = []

"""
```julia
addparticletype!(name, typeinfo)
```

Add a particle type name

# Arguments:

  - `name::String`:              name of new particle type
  - `numparticles::Int`:         number of particles by type
  - `typeinfo::Array{TypeInfo}`: input/output array of particle type information

# Returns:

  - nothing

# Example:

```jldoctest
# setup
typeinfo::Array{TypeInfo} = []

# add new type
Pulvis.addparticletype!("foo", 3, typeinfo)

# verify
@assert typeinfo[1].name == "foo"
@assert typeinfo[1].numparticles == 3

# add second type
Pulvis.addparticletype!("bar", 4, typeinfo)

# verify
@assert typeinfo[2].name == "bar"
@assert typeinfo[2].numparticles == 4

# output

```
"""
function addparticletype!(name::String, numparticles::Int, typeinfo::Array{TypeInfo})
    @assert findfirst(t -> t.name == name, typeinfo) == nothing "Types must have unique names"
    push!(typeinfo, TypeInfo(name, numparticles))
    return nothing
end

"""
```julia
addparticlestate!(name, dimension, stateinfo)
```

Add a particle state name and dimension

# Arguments:

  - `name::String`:                name of new particle state
  - `dimension::Int`:              dimension of new particle state
  - `stateinfo::Array{StateInfo}`: input/output array of particle state information

# Returns:

  - nothing

# Example:

```jldoctest
# setup
stateinfo::Array{StateInfo} = []

# add new state
Pulvis.addparticlestate!("baz", 2, stateinfo)

# verify
@assert stateinfo[1].name == "baz"
@assert stateinfo[1].dimension == 2

# add second state
Pulvis.addparticlestate!("quux", 3, stateinfo)

# verify
@assert stateinfo[2].name == "quux"
@assert stateinfo[2].dimension == 3

# output

```
"""
function addparticlestate!(name::String, dimension::Int, stateinfo::Array{StateInfo})
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

  - `typename::String`:            name of particle type to add state data for
  - `typeinfo::Array{TypeInfo}`:   input array of particle type information
  - `statenames::Array{String}`:   input array of names of particle states to add
  - `stateinfo::Array{StateInfo}`: input array of particle state information
  - `statedata::Dict`:             input/output dictionary of particle state data arrays
  - `initialize::Bool`:            flag to initialize state data array (default True)

# Returns:

  - nothing

# Example:

```jldoctest
# setup
typeinfo::Array{TypeInfo} = []
stateinfo::Array{StateInfo} = []
statedata = Dict()

# add new types
Pulvis.addparticletype!("foo", 3, typeinfo)
Pulvis.addparticletype!("bar", 4, typeinfo)

# add new states
Pulvis.addparticlestate!("baz", 2, stateinfo)
Pulvis.addparticlestate!("quux", 3, stateinfo)

# add particle state data
# -- add state for second type, allocate memory
Pulvis.addparticlestatedata!("bar", typeinfo, ["baz"], stateinfo, statedata, true)

# verify
@assert length(statedata["baz"]) == 2
# -- note array created for first type but no memory allocated
@assert length(statedata["baz"][1]) == 0
@assert length(statedata["baz"][2]) == 8

# -- add states for first type, don't allocate memory
Pulvis.addparticlestatedata!("foo", typeinfo, ["baz", "quux"], stateinfo, statedata, false)

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
    typeinfo::Array{TypeInfo},
    statenames::Array{String},
    stateinfo::Array{StateInfo},
    statedata::Dict,
    initialize::Bool = true,
)
    typelength = length(typeinfo)
    typeindex = findfirst(t -> t.name == typename, typeinfo)
    @assert typeindex != nothing "Type not known"

    # Loop over state names for the type
    for statename ∈ statenames
        hasstate = haskey(statedata, statename)
        statelength = hasstate ? length(statedata[statename]) : 0
        stateindex = findfirst(s -> s.name == statename, stateinfo)
        @assert stateindex != nothing "State not known"

        # Expand data for specific state for current type, as needed
        if typelength > statelength
            newstatedata = []
            resize!(newstatedata, typelength)

            # Copy current pointers
            if statelength > 0
                newstatedata[1:statelength] .= statedata[statename]
            end
            # Add new empty vectors
            for j ∈ (statelength+1):typelength
                newstatedata[j] = PulvisVector{PulvisScalar}()
            end
            # And copy to execution space
            statedata[statename] = adapt(
                PulvisVector,
                StateData{typelength,PulvisVector{PulvisScalar}}(newstatedata),
            )
        end
        # Allocate vector only if requested
        if initialize
            resize!(
                statedata[statename][typeindex],
                typeinfo[typeindex].numparticles * stateinfo[stateindex].dimension,
            )
        end
    end
    return nothing
end

"""
```julia
getparticlestatedata!(particledata, statename)
```

Get the state data for a specific state

# Arguments:

  - `statedata::Dict`:   input dictionary of particle state data arrays
  - `statename::String`: name of the state data to retrieve

# Returns:

  - `StateData` for the requested state

# Example:

```jldoctest
# setup
typeinfo::Array{TypeInfo} = []
stateinfo::Array{StateInfo} = []
statedata = Dict()

# add first type
Pulvis.addparticletype!("foo", 3, typeinfo)

# add new states
Pulvis.addparticlestate!("baz", 2, stateinfo)
Pulvis.addparticlestate!("quux", 3, stateinfo)

# add particle state data
# -- add states for first type, don't allocate memory
Pulvis.addparticlestatedata!("foo", typeinfo, ["baz", "quux"], stateinfo, statedata, false)

# add new type
Pulvis.addparticletype!("bar", 4, typeinfo)

# add more particle state data
# -- add state for second type, allocate memory
Pulvis.addparticlestatedata!("bar", typeinfo, ["baz"], stateinfo, statedata, true)

# verify
@assert length(statedata["baz"]) == 2
@assert length(statedata["quux"]) == 1

# get state data
baz = Pulvis.getparticlestatedata(statedata, "baz")
quux = Pulvis.getparticlestatedata(statedata, "quux")
nix = Pulvis.getparticlestatedata(statedata, "nix")

# verify
@assert length(baz) == 2
@assert length(baz[1]) == 0
@assert length(baz[2]) == 8
@assert length(quux) == 1
@assert length(quux[1]) == 0
@assert length(nix) == 0

# output

```
"""
function getparticlestatedata(statedata::Dict, statename::String)
    if haskey(statedata, statename)
        return statedata[statename]
    else
        return StateData{0,PulvisVector{PulvisScalar}}()
    end
end

"""
```julia
squareparticlestatedata!(typeinfo, particledata)
```

Ensure all state data arrays have the same length

# Arguments:

  - `typeinfo::Array{TypeInfo}`:   input array of particle type information
  - `stateinfo::Array{StateInfo}`: input array of particle state information
  - `statedata::Dict`:             input/output dictionary of particle state data arrays

# Returns:

  - nothing

# Example:

```jldoctest
# setup
typeinfo::Array{TypeInfo} = []
stateinfo::Array{StateInfo} = []
statedata = Dict()

# add first type
Pulvis.addparticletype!("foo", 3, typeinfo)

# add new states
Pulvis.addparticlestate!("baz", 2, stateinfo)
Pulvis.addparticlestate!("quux", 3, stateinfo)

# add particle state data
# -- add states for first type, don't allocate memory
Pulvis.addparticlestatedata!("foo", typeinfo, ["baz", "quux"], stateinfo, statedata, false)

# add new type
Pulvis.addparticletype!("bar", 4, typeinfo)

# add more particle state data
# -- add state for second type, allocate memory
Pulvis.addparticlestatedata!("bar", typeinfo, ["baz"], stateinfo, statedata, true)

# verify
@assert length(statedata["baz"]) == 2
@assert length(statedata["quux"]) == 1

# square off state arrays
Pulvis.squareparticlestatedata!(typeinfo, stateinfo, statedata)

# verify
@assert length(statedata["baz"]) == 2
@assert length(statedata["quux"]) == 2

# output

```
"""
function squareparticlestatedata!(
    typeinfo::Array{TypeInfo},
    stateinfo::Array{StateInfo},
    statedata::Dict,
)
    typelength = length(typeinfo)

    # Check all states
    for i ∈ 1:length(stateinfo)
        statename = stateinfo[i].name
        hasstate = haskey(statedata, statename)
        statelength = hasstate ? length(statedata[statename]) : 0

        # Square off if the state is missing some types
        if typelength > statelength
            newstatedata = []
            resize!(newstatedata, typelength)

            # Copy current pointers
            newstatedata[1:statelength] .= statedata[statename]
            # Add new empty vectors
            for j ∈ (statelength+1):typelength
                newstatedata[j] = PulvisVector{PulvisScalar}()
            end
            # And copy to execution space
            statedata[statename] = adapt(
                PulvisVector,
                StateData{typelength,PulvisVector{PulvisScalar}}(newstatedata),
            )
        end
    end
    return nothing
end

"""
```julia
addparticles!(typename, numparticles, newdata, typeinfo, stateinfo, particledata)
```

Add new particle state arrays for particle type

# Arguments:

  - `typename::String`:            name of particle type to remove particles for
  - `numparticles::Int`:           number of new particles to add
  - `newdata::Dict`:               input dictionary of new particle state data arrays
  - `typeinfo::Array{TypeInfo}`:   input array of particle type information
  - `stateinfo::Array{StateInfo}`: input array of particle state information
  - `statedata::Dict`:             input/output dictionary of particle state data arrays

# Returns:

  - nothing

# Example:

```jldoctest
# setup
typeinfo::Array{TypeInfo} = []
stateinfo::Array{StateInfo} = []
statedata = Dict()

# add particle types
Pulvis.addparticletype!("foo", 3, typeinfo)
Pulvis.addparticletype!("bar", 4, typeinfo)

# add particle states
Pulvis.addparticlestate!("baz", 2, stateinfo)
Pulvis.addparticlestate!("quux", 3, stateinfo)

# add particle state data
Pulvis.addparticlestatedata!("foo", typeinfo, ["baz", "quux"], stateinfo, statedata, true)
Pulvis.addparticlestatedata!("bar", typeinfo, ["baz"], stateinfo, statedata, true)

# square off state arrays
Pulvis.squareparticlestatedata!(typeinfo, stateinfo, statedata)

# verify
@assert typeinfo[1].numparticles == 3
@assert length(statedata["baz"][1]) == 6
@assert length(statedata["quux"][1]) == 9

# add particles
newdata = Dict()
newdata["baz"] = [1.0, 2.0]
newdata["quux"] = [1.0, 2.0, 3.0]
Pulvis.addparticles!("foo", 1, newdata, typeinfo, stateinfo, statedata)

# verify
@assert typeinfo[1].numparticles == 4
@assert length(statedata["baz"][1]) == 8
@assert Vector(statedata["baz"][1])[(end-1):end] ≈ newdata["baz"]
@assert length(statedata["quux"][1]) == 12
@assert Vector(statedata["quux"][1])[(end-2):end] ≈ newdata["quux"]

# verify
@assert typeinfo[2].numparticles == 4
@assert length(statedata["baz"][2]) == 8
@assert length(statedata["quux"][2]) == 0

# add more particles
newdata = Dict()
newdata["baz"] = [1.0, 2.0, 3.0, 4.0]
Pulvis.addparticles!("bar", 2, newdata, typeinfo, stateinfo, statedata)

# verify
@assert typeinfo[2].numparticles == 6
@assert length(statedata["baz"][2]) == 12
@assert Vector(statedata["baz"][2])[(end-3):end] ≈ newdata["baz"]
@assert length(statedata["quux"][2]) == 0

# output

```
"""
function addparticles!(
    typename::String,
    numparticles::Int,
    newdata::Dict,
    typeinfo::Array{TypeInfo},
    stateinfo::Array{StateInfo},
    statedata::Dict,
)
    typelength = length(typeinfo)
    typeindex = findfirst(t -> t.name == typename, typeinfo)
    @assert typeindex != nothing "Type not known"

    # Loop over all states
    for i ∈ 1:length(stateinfo)
        # Sanity check
        statename = stateinfo[i].name
        statedimension = stateinfo[i].dimension
        @assert !haskey(newdata, statename) ||
                (length(newdata[statename]) == numparticles * statedimension)

        # Skip if state isn't in data to add
        if !haskey(newdata, statename)
            continue
        end

        # Check if enough space is already allocated
        if length(newdata[statename]) > 0
            newstatedata = []
            resize!(newstatedata, typelength)

            # Copy current pointers
            newstatedata[1:(typeindex-1)] .= statedata[statename][1:(typeindex-1)]
            # Add the target particles in the target type
            data = adapt(Vector{PulvisScalar}, statedata[statename][typeindex])
            if length(statedata[statename]) <
               ((typeinfo[typeindex].numparticles + numparticles) * statedimension)
                resize!(
                    data,
                    (typeinfo[typeindex].numparticles + numparticles) * statedimension,
                )
            end
            data[(typeinfo[typeindex].numparticles*statedimension+1):((typeinfo[typeindex].numparticles+numparticles)*statedimension)] .=
                newdata[statename]
            newstatedata[typeindex] = adapt(PulvisVector{PulvisScalar}, data)
            # Copy the rest of the current pointers
            newstatedata[(typeindex+1):typelength] .=
                statedata[statename][(typeindex+1):typelength]
            # And copy to execution space
            statedata[statename] = adapt(
                PulvisVector,
                StateData{typelength,PulvisVector{PulvisScalar}}(newstatedata),
            )
        end
    end
    # Finally increase particle count
    typeinfo[typeindex] =
        TypeInfo(typename, typeinfo[typeindex].numparticles + numparticles)
    return nothing
end

"""
```julia
removeparticles!(typename, particleids, typeinfo, stateinfo, particledata)
```

Remove particles by local id numbers

# Arguments:

  - `typename::String`:            name of particle type to remove particles for
  - `particleids::Array{Int}`:     input array of particle local id numbers for particles to remove
  - `typeinfo::Array{TypeInfo}`:   input array of particle type information
  - `stateinfo::Array{StateInfo}`: input array of particle state information
  - `statedata::Dict`:             input/output dictionary of particle state data arrays

# Returns:

  - nothing

# Example:

```jldoctest
# setup
typeinfo::Array{TypeInfo} = []
stateinfo::Array{StateInfo} = []
statedata = Dict()

# add particle types
Pulvis.addparticletype!("foo", 3, typeinfo)
Pulvis.addparticletype!("bar", 4, typeinfo)

# add particle states
Pulvis.addparticlestate!("baz", 2, stateinfo)
Pulvis.addparticlestate!("quux", 3, stateinfo)

# add particle state data
Pulvis.addparticlestatedata!("foo", typeinfo, ["baz", "quux"], stateinfo, statedata, true)
Pulvis.addparticlestatedata!("bar", typeinfo, ["baz"], stateinfo, statedata, true)

# square off state arrays
Pulvis.squareparticlestatedata!(typeinfo, stateinfo, statedata)

# verify
@assert typeinfo[1].numparticles == 3
# note that allocated lengths should not decrease
@assert length(statedata["baz"][1]) == 6
@assert length(statedata["quux"][1]) == 9

# remove particles
Pulvis.removeparticles!("foo", [1, 3], typeinfo, stateinfo, statedata)

# verify
@assert typeinfo[1].numparticles == 1
# note that allocated lengths should not decrease
@assert length(statedata["baz"][1]) == 6
@assert length(statedata["quux"][1]) == 9

# verify
@assert typeinfo[2].numparticles == 4
# note that allocated lengths should not decrease
@assert length(statedata["baz"][2]) == 8
@assert length(statedata["quux"][2]) == 0

# remove more particles
Pulvis.removeparticles!("bar", [2, 3], typeinfo, stateinfo, statedata)

# verify
@assert typeinfo[2].numparticles == 2
# note that allocated lengths should not decrease
@assert length(statedata["baz"][2]) == 8
@assert length(statedata["quux"][2]) == 0

# output

```
"""
function removeparticles!(
    typename::String,
    particleids::Array{Int},
    typeinfo::Array{TypeInfo},
    stateinfo::Array{StateInfo},
    statedata::Dict,
)
    typelength = length(typeinfo)
    typeindex = findfirst(t -> t.name == typename, typeinfo)
    @assert typeindex != nothing "Type not known"
    sort(particleids, rev = true)

    # Loop over all states
    for i ∈ 1:length(stateinfo)
        statename = stateinfo[i].name
        statedimension = stateinfo[i].dimension
        newstatedata = []
        resize!(newstatedata, typelength)

        # Copy current pointers
        newstatedata[1:(typeindex-1)] .= statedata[statename][1:(typeindex-1)]
        # Drop the target particles in the target type
        data = adapt(Vector{PulvisScalar}, statedata[statename][typeindex])
        numparticles = typeinfo[typeindex].numparticles
        for id ∈ particleids
            if id < length(data)
                data[(id*statedimension):((numparticles-1)*statedimension)] .=
                    data[((id+1)*statedimension):(numparticles*statedimension)]
            end
            numparticles -= 1
        end
        newstatedata[typeindex] = adapt(PulvisVector{PulvisScalar}, data)
        # Copy the rest of the current pointers
        newstatedata[(typeindex+1):typelength] .=
            statedata[statename][(typeindex+1):typelength]
        # And copy to execution space
        statedata[statename] = adapt(
            PulvisVector,
            StateData{typelength,PulvisVector{PulvisScalar}}(newstatedata),
        )
    end
    # Finally decrease particle count
    typeinfo[typeindex] = TypeInfo(
        typeinfo[typeindex].name,
        typeinfo[typeindex].numparticles - length(particleids),
    )
    return nothing
end

"""
```julia
copyparticles!(
    typenamesource,
    typenamedestination,
    particleidstocopy,
    particleidstoremove,
    typeinfo,
    stateinfo,
    particledata,
)
```

Copy particles from one type to another by local id numbers in the source type

# Arguments:

  - `typenamesource::String`:          name of particle type to move particles from
  - `typenamedestination::String`:     name of particle type to move particles to
  - `particleidstocopy::Array{Int}`:   input array of particle source local id numbers for particles to copy
  - `particleidstoremove::Array{Int}`: input/output array of particle local id numbers for particles to remove
  - `typeinfo::Array{TypeInfo}`:       input array of particle type information
  - `stateinfo::Array{StateInfo}`:     input array of particle state information
  - `statedata::Dict`:                 input/output dictionary of particle state data arrays

# Returns:

  - nothing

# Example:

```jldoctest
# setup
typeinfo::Array{TypeInfo} = []
stateinfo::Array{StateInfo} = []
statedata = Dict()

# add particle types
Pulvis.addparticletype!("foo", 3, typeinfo)
Pulvis.addparticletype!("bar", 4, typeinfo)

# add particle states
Pulvis.addparticlestate!("baz", 2, stateinfo)
Pulvis.addparticlestate!("quux", 3, stateinfo)

# add particle state data
Pulvis.addparticlestatedata!("foo", typeinfo, ["baz", "quux"], stateinfo, statedata, true)
Pulvis.addparticlestatedata!("bar", typeinfo, ["baz", "quux"], stateinfo, statedata, true)

# verify
@assert typeinfo[1].numparticles == 3
# note that allocated lengths should not decrease
@assert length(statedata["baz"][1]) == 6
@assert length(statedata["quux"][1]) == 9

# verify
@assert typeinfo[2].numparticles == 4
# note that allocated lengths should not decrease
@assert length(statedata["baz"][2]) == 8
@assert length(statedata["quux"][2]) == 12

# set state values that are easy to track
statedata["baz"][1] .= 1:6
statedata["quux"][1] .= 11:19

# move particles
particlestomove::Vector{Int} = [1, 3]
particlestoremove::Vector{Int} = []
Pulvis.copyparticles!(
    "foo",
    "bar",
    particlestomove,
    particlestoremove,
    typeinfo,
    stateinfo,
    statedata,
)

# verify
@assert typeinfo[1].numparticles == 3
# note that allocated lengths should not decrease
@assert length(statedata["baz"][1]) == 6
@assert length(statedata["quux"][1]) == 9

# verify
@assert typeinfo[2].numparticles == 6
# note that allocated lengths should increase
@assert length(statedata["baz"][2]) == 12
@assert Vector(statedata["baz"][2])[9:end] == [1, 2, 5, 6]
@assert length(statedata["quux"][2]) == 18
@assert Vector(statedata["quux"][2])[13:end] == [11, 12, 13, 17, 18, 19]

# remove particles
Pulvis.removeparticles!("foo", particlestoremove, typeinfo, stateinfo, statedata)

# verify
@assert typeinfo[1].numparticles == 1
# note that allocated lengths should not decrease
@assert length(statedata["baz"][1]) == 6
@assert length(statedata["quux"][1]) == 9

# verify
@assert typeinfo[2].numparticles == 6
# note that allocated lengths should not decrease
@assert length(statedata["baz"][2]) == 12
@assert length(statedata["quux"][2]) == 18

# output

```
"""
function copyparticles!(
    typenamesource::String,
    typenamedestination::String,
    particleidstocopy::Array{Int},
    particleidstoremove::Array{Int},
    typeinfo::Array{TypeInfo},
    stateinfo::Array{StateInfo},
    statedata::Dict,
)
    typelength = length(typeinfo)
    typeindexsource = findfirst(t -> t.name == typenamesource, typeinfo)
    @assert typeindexsource != nothing "Source type not known"
    typeindexdestination = findfirst(t -> t.name == typenamedestination, typeinfo)
    @assert typeindexdestination != nothing "Destination type not known"
    sort(particleidstocopy)

    # Loop over all states
    numparticlestomove = length(particleidstocopy)
    for i ∈ 1:length(stateinfo)
        statename = stateinfo[i].name
        statedimension = stateinfo[i].dimension
        newstatedata = []
        resize!(newstatedata, typelength)

        # Copy current pointers
        newstatedata[1:(typeindexdestination-1)] .=
            statedata[statename][1:(typeindexdestination-1)]
        # Move the target particles in the target type
        datasource = adapt(Vector{PulvisScalar}, statedata[statename][typeindexsource])
        datadestination =
            adapt(Vector{PulvisScalar}, statedata[statename][typeindexdestination])
        # Resize if needed
        numparticlesdestination = typeinfo[typeindexdestination].numparticles
        if length(datadestination) <
           ((numparticlesdestination + numparticlestomove) * statedimension)
            resize!(
                datadestination,
                (numparticlesdestination + numparticlestomove) * statedimension,
            )
        end
        for id ∈ particleidstocopy
            datadestination[(numparticlesdestination*statedimension+1):((numparticlesdestination+1)*statedimension)] .=
                datasource[((id-1)*statedimension+1):(id*statedimension)]
            numparticlesdestination += 1
        end
        newstatedata[typeindexdestination] =
            adapt(PulvisVector{PulvisScalar}, datadestination)
        # Copy the rest of the current pointers
        newstatedata[(typeindexdestination+1):typelength] .=
            statedata[statename][(typeindexdestination+1):typelength]
        # And copy to execution space
        statedata[statename] = adapt(
            PulvisVector,
            StateData{typelength,PulvisVector{PulvisScalar}}(newstatedata),
        )
    end
    # Mark particles to remove
    append!(particleidstoremove, particleidstocopy)
    # Finally increase particle count
    typeinfo[typeindexdestination] = TypeInfo(
        typeinfo[typeindexdestination].name,
        typeinfo[typeindexdestination].numparticles + length(particleidstocopy),
    )
    return nothing
end
