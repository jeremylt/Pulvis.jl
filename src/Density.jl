# ------------------------------------------------------------------------------
# Density computation at a single point
# ------------------------------------------------------------------------------

"""
```julia
sumweightedmassself!(particleid, radius, m, ρ)
```

Compute the density of the particle from the shape function, radius, and mass

# Arguments:

  - `particleid::Int`: index number of current particle
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
sumweightedmassneighbors!(particleid, neighborids, dimension, r, radius, m, ρ)
```

Sum the density contributions to the particle from neighbor particles

# Arguments:

  - `particleid::Int`:    index number of current particle
  - `neighborids::Array`: input array of neighbor particle indices
  - `dimension::Int`:     spatial dimension of problem
  - `r::Array`:           input array of particle positions
  - `radius::Array`:      input array of particle radi
  - `m::Array`:           input array of particle masses
  - `ρ::Array`:           output array of particle densities

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
    particleid::Int,
    neighborids::Array,
    dimension::Int,
    r::Array,
    radius::Array,
    m::Array,
    ρ::Array,
)
    for i = neighborids[particleid]:(neighborids[particleid+1]-1)
        neighborid = neighborids[i]
        absdistance = norm(
            view(r, ((particleid-1)*dimension+1):(particleid*dimension)) -
            view(r, ((neighborid-1)*dimension+1):(neighborid*dimension)),
        )

        ρ[particleid] += w(absdistance, radius[particleid], dimension) * m[particleid]
    end
    return nothing
end
