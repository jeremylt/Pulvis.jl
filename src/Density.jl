# ------------------------------------------------------------------------------
# Density computation at a single point
# ------------------------------------------------------------------------------

"""
```julia
sumweightedmassself!(particleid, radius, m, ρ)
```

Compute the density of the particle from the shape function, radius, and mass

# Arguments:

  - `particleid::Int`: local id number of current particle
  - `dimension::Int`:  spatial dimension of problem
  - `radius::Array`:   input array of particle radi
  - `m::Array`:        input array of particle masses
  - `ρ::Array`:        output array of particle densities

# Returns:

  - none

# Example:

```jldoctest
# setup
particleid = 1
dimension = 2
radius = [1.0]
m = [42.0]
ρ = [0.0]

# compute density
sumweightedmassself!(particleid, dimension, radius, m, ρ)

# verify
@assert ρ[particleid] ≈ 116.29365122090275

# output

```
"""
function sumweightedmassself!(
    particleid::Int,
    dimension::Int,
    radius::Array,
    m::Array,
    ρ::Array,
)
    ρ[particleid] = w0(radius[particleid], dimension) * m[particleid]
    return nothing
end

"""
```julia
sumweightedmassneighbors!(
    particleindex,
    particleids,
    typeids,
    neighborindexes,
    dimension,
    r,
    radius,
    m,
    ρ,
)
```

Sum the density contributions to the particle from neighbor particles

# Arguments:

  - `particleindex::Int`:     global index number of current particle
  - `particleids::Array`:     input array of particle ids by type
  - `typeids::Array`:         input array of particle type ids
  - `neighborindexes::Array`: input array of neighbor particle global indexes
  - `dimension::Int`:         spatial dimension of problem
  - `r::StateData`:           input arrays of particle positions by particle type and id
  - `radius::StateData`:      input arrays of particle radi by particle type and id
  - `m::StateData`:           input arrays of particle masses by particle type and id
  - `ρ::StateData`:           output arrays of particle densities by particle type and id

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
function sumweightedmassneighbors!(
    particleindex::Int,
    particleids::Array,
    typeids::Array,
    neighborindexes::Array,
    dimension::Int,
    r::StateData,
    radius::StateData,
    m::StateData,
    ρ::StateData,
)
    particleid = particleids[particleindex]
    typeid = typeids[particleindex]

    for i = neighborindexes[particleindex]:(neighborindexes[particleindex+1]-1)
        neighborindex = neighborindexes[i]
        neighborid = particleids[neighborindex]
        neighbortypeid = typeids[neighborindex]

        absdistance = norm(
            view(r[typeid], ((particleid-1)*dimension+1):(particleid*dimension)) -
            view(r[neighbortypeid], ((neighborid-1)*dimension+1):(neighborid*dimension)),
        )

        ρ[typeid][particleid] +=
            w(absdistance, radius[typeid][particleid], dimension) * m[typeid][particleid]
    end
    return nothing
end
