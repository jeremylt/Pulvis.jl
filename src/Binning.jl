# ------------------------------------------------------------------------------
# Binning particles
# ------------------------------------------------------------------------------

bins::Vector{Vector{Int}} = [[]]
binids::Vector{Int} = []

"""
```julia
binidtodimids(bininfo, dimension, binid)
```

Convert bin id to x/y/z ids

# Arguments:

  - `bininfo::BinInfo`: info about background bin layout
  - `dimension::Int`:   spatial dimension of problem
  - `binid::Int`:       id of bin to get neighbors for

# Returns:

  - nothing

# Example:

```jldoctest
# setup
lower = [0.0, 0.0, 0.0]
upper = [3.0, 3.0, 3.0]
width = 1.0
counts = Array{Int}(ceil.((upper - lower) / width))
bininfo = BinInfo(lower, upper, width, counts)

# verify
@assert binidtodimids(bininfo, 3, 1) == (1, 1, 1)
@assert binidtodimids(bininfo, 3, 14) == (2, 2, 2)
@assert binidtodimids(bininfo, 3, 27) == (3, 3, 3)

# output

```
"""
function binidtodimids(bininfo::BinInfo, dimension::Int, binid::Int)
    xid = mod1(binid, bininfo.counts[1])
    yid = 1
    if dimension > 1
        yid = ceil(mod1(binid, bininfo.counts[1] * bininfo.counts[2]) / bininfo.counts[1])
    end
    zid = 1
    if dimension > 2
        zid = ceil(binid / (bininfo.counts[1] * bininfo.counts[2]))
    end
    return (xid, yid, zid)
end

"""
```julia
getneighborbinids!(bininfo, dimension, binid, neighborbinids)
```

Get ids of neighboring bin based upon binid

# Arguments:

  - `bininfo::BinInfo`:            info about background bin layout
  - `dimension::Int`:              spatial dimension of problem
  - `binid::Int`:                  id of bin to get neighbors for
  - `neighborbinids::Vector{Int}`: input/output array of neighbor bin ids

# Returns:

  - nothing

# Example:

```jldoctest
# setup
lower = [0.0, 0.0, 0.0]
upper = [3.0, 3.0, 3.0]
width = 1.0
counts = Array{Int}(ceil.((upper - lower) / width))
bininfo = BinInfo(lower, upper, width, counts)
neighborbinids::Array{Int} = []
resize!(neighborbinids, 28)

# middle
getneighborbinids!(bininfo, 3, 14, neighborbinids)

# verify
@assert neighborbinids[1:27] == 1:27

# low corner
getneighborbinids!(bininfo, 3, 1, neighborbinids)

# verify
@assert neighborbinids[1:9] == [1, 2, 4, 5, 10, 11, 13, 14, -1]

# high corner
getneighborbinids!(bininfo, 3, 27, neighborbinids)

# verify
@assert neighborbinids[1:9] == [14, 15, 17, 18, 23, 24, 26, 27, -1]

# output

```
"""
function getneighborbinids!(
    bininfo::BinInfo,
    dimension::Int,
    binid::Int,
    neighborbinids::Vector{Int},
)
    # Reconstruct x/y/z ids
    (xid, yid, zid) = binidtodimids(bininfo, dimension, binid)

    # Build neighbor ids
    xrange = max(xid-1, 1):min(xid+1, bininfo.counts[1])
    yrange = 1:1
    if dimension > 1
        yrange = max(yid-1, 1):min(yid+1, bininfo.counts[2])
    end
    zrange = 1:1
    if dimension > 1
        zrange = max(zid-1, 1):min(zid+1, bininfo.counts[3])
    end

    count = 1
    for k in zrange
        for j in yrange
            for i in xrange
                neighborbinids[count] =
                    i +
                    (j - 1) * bininfo.counts[1] +
                    (k - 1) * bininfo.counts[1] * bininfo.counts[2]
                count += 1
            end
        end
    end
    neighborbinids[count] = -1
    return nothing
end

"""
```julia
particlepositiontobinid(bininfo, dimension, r)
```

Convert particle coordinates to bin id

# Arguments:

  - `bininfo::BinInfo`:       info about background bin layout
  - `dimension::Int`:         spatial dimension of problem
  - `r::Array{PulvisScalar}`: input array of particle position

# Returns:

  - bin id

# Example:

```jldoctest
# setup
lower = [0.0, 1.0, 2.0]
upper = [1.0, 3.0, 5.0]
width = 0.3
counts = Array{Int}(ceil.((upper - lower) / width))
bininfo = BinInfo(lower, upper, width, counts)

# particle location
r = [0.75, 3.0, 4.7]

# verify
@assert particlepositiontobinid(bininfo, 1, r) == 3
@assert particlepositiontobinid(bininfo, 2, r) == 27
@assert particlepositiontobinid(bininfo, 3, r) == 279

# output

```
"""
function particlepositiontobinid(bininfo::BinInfo, dimension::Int, r::Array{PulvisScalar})
    binid = max(ceil(Int, (r[1] - bininfo.lower[1]) / bininfo.width), 1)
    if dimension > 1
        binid +=
            max((ceil(Int, (r[2] - bininfo.lower[2]) / bininfo.width) - 1), 0) *
            bininfo.counts[1]
    end
    if dimension > 2
        binid +=
            max((ceil(Int, (r[3] - bininfo.lower[3]) / bininfo.width) - 1), 0) *
            bininfo.counts[1] *
            bininfo.counts[2]
    end
    return binid
end

"""
```julia
particlestobins(bininfo, dimension, r)
```

Assign particles to bins

# Arguments:

  - `bininfo::BinInfo`:          info about background bin layout
  - `dimension::Int`:            spatial dimension of problem
  - `r::StateData`:              input arrays of particle positions by particle type and id
  - `bins::Vector{Vector{Int}}`: input/output arrays of global particle indexes included in each bin
  - `particleids::Array{Int}`:   input/output array of bin ids for each particle

# Returns:

  - none

# Example:

```jldoctest
# setup
lower = [0.0, 1.0, 2.0]
upper = [1.0, 3.0, 5.0]
width = 1.0
counts = Array{Int}(ceil.((upper - lower) / width))
bininfo = BinInfo(lower, upper, width, counts)
bins::Vector{Vector{Int}} = [[]]
binids = []

# particle locations
r = StateData{2}([[0.0, 1.5, 2.5, 0.0, 1.0, 2.0], [1.0, 3.0, 5.0]])

# assign particles to bins
particlestobins!(bininfo, 3, r, bins, binids)

# verify
@assert bins == [[1, 2], [], [], [], [], [3]]
@assert binids == [1, 1, 6]

# new particle locations
r = StateData{2}([[0.0, 2.5, 2.5, 0.0, 1.0, 2.0, 0.5, 2.0, 3.5], [1.0, 3.0, 5.0]])

# assign particles to new bins
particlestobins!(bininfo, 3, r, bins, binids)

# verify
@assert bins == [[2], [1], [3], [], [], [4]]
@assert binids == [2, 1, 3, 6]

# output

```
"""
function particlestobins!(
    bininfo::BinInfo,
    dimension::Int,
    r::StateData,
    bins::Vector{Vector{Int}},
    binids::Array{Int},
)
    # Allocate bin vectors, if needed
    totalbins = bininfo.counts[1]
    if dimension > 1
        totalbins *= bininfo.counts[2]
    end
    if dimension > 2
        totalbins *= bininfo.counts[3]
    end
    if length(bins) < totalbins
        for i = (length(bins)+1):totalbins
            push!(bins, Vector{Int}([]))
        end
    end

    # Clear old data
    for i = 1:totalbins
        empty!(bins[i])
    end

    # Populate bins
    count = 0
    for t = 1:length(r)
        numparticles = Int(length(r[t]) / dimension)
        if length(binids) < length(r[t]) / dimension + count
            resize!(binids, numparticles + count)
        end
        @threads for i = 1:numparticles
            binid = particlepositiontobinid(
                bininfo,
                dimension,
                r[t][Int((i-1)*dimension+1):Int(i*dimension)],
            )
            push!(bins[binid], i + count)
            binids[i+count] = binid
        end
        count += numparticles
    end
    return nothing
end

"""
```julia
typeinfotoids!(typeinfo, typeids, particleids)
```

Convert `TypeInfo` array into global arrays of local particle ids by type and particle type ids

# Arguments:

  - `typeinfo::Array{TypeInfo}`: input array of particle type information
  - `particleids::Array{Int}`:   output array of local particle ids by type
  - `typeids::Array{Int}`:       output array of particle type ids

# Returns:

  - none

# Example:

```jldoctest
# setup
typeinfo::Array{TypeInfo} = []
particleids::AbstractArray{Int} = []
typeids::AbstractArray{Int} = []

# add types
addparticletype!("foo", 3, typeinfo)
addparticletype!("bar", 4, typeinfo)

# convert to id arrays
typeinfotoids!(typeinfo, particleids, typeids)

# verify
@assert particleids == [1, 2, 3, 1, 2, 3, 4]
@assert typeids == [1, 1, 1, 2, 2, 2, 2]

# output

```
"""
function typeinfotoids!(
    typeinfo::Array{TypeInfo},
    particleids::Array{Int},
    typeids::Array{Int},
)
    # Count total particles
    numparticles = 0
    for info in typeinfo
        numparticles += info.count
    end

    # Resize if needed
    if length(particleids) < numparticles
        resize!(particleids, numparticles)
        resize!(typeids, numparticles)
    end

    # Map local and type ids to global arrays
    current = 0
    for i = 1:length(typeinfo)
        @threads for j = 1:typeinfo[i].count
            particleids[j+current] = j
            typeids[j+current] = i
        end
        current += typeinfo[i].count
    end
    return nothing
end

"""
```julia
buildpairlist!(
    typeinfo,
    bininfo,
    dimension,
    r,
    radius,
    particleids,
    typeids,
    nighborindexes,
)
```

Populate list of particle pairs

# Arguments:

  - `bininfo::BinInfo`:            info about background bin layout
  - `typeinfo::Array{TypeInfo}`:   input array of particle type information
  - `typeinfo::Array{TypeInfo}`:   input array of particle type information
  - `dimension::Int`:              spatial dimension of problem
  - `r::StateData`:                input arrays of particle positions by particle type and id
  - `radius::StateData`:           input arrays of particle radi by particle type and id
  - `neighborindexes::Array{Int}`: output array of neighbor particle global indexes
  - `particleids::Array{Int}`:     output array of local particle ids by type
  - `typeids::Array{Int}`:         output array of particle type ids

# Returns:

  - none

# Example:

Bin layout:

* * *

## |   index 2 |   index 1 |

## |   index 3 |           |

## |           |   index 4 |

```jldoctest
# setup
typeinfo::Array{TypeInfo} = []
addparticletype!("foo", 3, typeinfo)
addparticletype!("bar", 1, typeinfo)

# bins
lower = [0.0, 1.0, 2.0]
upper = [1.0, 3.0, 5.0]
width = 1.0
counts = Array{Int}(ceil.((upper - lower) / width))
bininfo = BinInfo(lower, upper, width, counts)

# particle data
r = StateData{2}([[0.0, 2.5, 2.5, 0.0, 1.0, 2.0, 0.5, 2.0, 3.5], [1.0, 3.0, 5.0]])
radius = StateData{2}([[2.0, 2.0, 1.0], [2.0]])

# build particle pair lists
particleids::AbstractArray{Int} = []
typeids::AbstractArray{Int} = []
neighborindexes::AbstractArray{Int} = []
buildpairlist!(typeinfo, bininfo, 3, r, radius, particleids, typeids, neighborindexes)

# verify
@assert particleids == [1, 2, 3, 1]
@assert typeids == [1, 1, 1, 2]
#                           --- offsets ---  -- neighbors --
#                                           |  1  |  2  | 4|
@assert neighborindexes == [6, 8, 10, 10, 11, 2, 3, 1, 3, 3]

# output

```
"""
function buildpairlist!(
    typeinfo::Array{TypeInfo},
    bininfo::BinInfo,
    dimension::Int,
    r::StateData,
    radius::StateData,
    particleids::Array{Int},
    typeids::Array{Int},
    nighborindexes::Array{Int},
)
    # Update particle and type ids
    typeinfotoids!(typeinfo, particleids, typeids)
    numparticles = length(particleids)

    # Map particles to bins
    particlestobins!(bininfo, dimension, r, bins, binids)

    # Compute neighbors
    empty!(neighborindexes)
    resize!(neighborindexes, numparticles + 1)
    neighborbinids::Array{Int} = []
    resize!(neighborbinids, 3^dimension + 1)

    # Loop over global particle ids
    current = numparticles + 2
    for particleindex = 1:numparticles
        # Mark starting index
        neighborindexes[particleindex] = current

        # Particle location
        typeid = typeids[particleindex]
        particleid = particleids[particleindex]
        r_i = view(r[typeid], ((particleid-1)*dimension+1):(particleid*dimension))

        # Interaction radius
        radius_i = radius[typeid][particleid]

        # Check each neighbor bin
        binid = binids[particleindex]
        getneighborbinids!(bininfo, dimension, binid, neighborbinids)
        println(neighborbinids)
        for neighborbin in neighborbinids
            if neighborbin == -1
                break
            end
            # Check each neighbor in bin
            for neighborindex in bins[neighborbin]
                if neighborindex == particleindex
                    continue
                end

                # Neighbor location
                neighbortypeid = typeids[neighborindex]
                neighborid = particleids[neighborindex]
                r_j = view(
                    r[neighbortypeid],
                    ((neighborid-1)*dimension+1):(neighborid*dimension),
                )
                if norm(r_i - r_j) <= radius_i
                    push!(nighborindexes, neighborindex)
                    current += 1
                end
            end
        end
    end
    neighborindexes[numparticles+1] = current
    return nothing
end
