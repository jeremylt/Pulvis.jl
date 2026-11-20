# ------------------------------------------------------------------------------
# Density computation at a single point
# ------------------------------------------------------------------------------

"""
```julia
sumweightedmassself!(particleid, radius, m, ρ)
```

Compute the density of the particle from the shape function, radius, and mass

# Arguments:

  - `particleid::Int`:             local id number of current particle
  - `dimension::Int`:              spatial dimension of problem
  - `radius::Array{PulvisScalar}`: input array of particle radi
  - `m::Array{PulvisScalar}`:      input array of particle masses
  - `ρ::Array{PulvisScalar}`:      output array of particle densities

# Returns:

  - none

# Example:

```jldoctest
# setup
particleid = 1
dimension = 2
radius = StateData{1}([[1.0]])
m = StateData{1}([[42.0]])
ρ = StateData{1}([[0.0]])

# set shape function
setkernelshape(Kernel.CubicSpline)

# compute density
Pulvis.sumweightedmassself!(particleid, dimension, radius[1], m[1], ρ[1])

# verify
@assert ρ[1][particleid] ≈ 76.39437268410974

# output

```
"""
function sumweightedmassself!(
    particleid::Int,
    dimension::Int,
    radius::Array{PulvisScalar},
    m::Array{PulvisScalar},
    ρ::Array{PulvisScalar},
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

  - `particleindex::Int`:          global index number of current particle
  - `particleids::Array{Int}`:     input array of particle ids by type
  - `typeids::Array{Int}`:         input array of particle type ids
  - `neighborindexes::Array{Int}`: input array of neighbor particle global indexes
  - `dimension::Int`:              spatial dimension of problem
  - `r::StateData`:                input arrays of particle positions by particle type and id
  - `radius::StateData`:           input arrays of particle radi by particle type and id
  - `m::StateData`:                input arrays of particle masses by particle type and id
  - `ρ::StateData`:                output arrays of particle densities by particle type and id

# Returns:

  - none

# Example:

```jldoctest
# setup
particleindex = 1
particleids = [1, 2, 1]
typeids = [1, 1, 2]
neighborindexes = [3, 5, 2, 3]
dimension = 2
r = StateData{2}([[1.0, 1.0, 1.1, 1.2], [1.3, 1.4]])
radius = StateData{2}([[1.0, 1.1], [1.2]])
m = StateData{2}([[42.0, 13.0], [14.0]])
ρ = StateData{2}([[1.0, 0.0], [0.0]])

# set shape function
setkernelshape(Kernel.CubicSpline)

# compute density
Pulvis.sumweightedmassneighbors!(
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

# verify
@assert ρ[typeids[particleindex]][particleids[particleindex]] ≈ 78.69934436250794
neighborindex = neighborindexes[3]
@assert ρ[typeids[neighborindex]][particleids[neighborindex]] ≈ 0.0
neighborindex = neighborindexes[4]
@assert ρ[typeids[neighborindex]][particleids[neighborindex]] ≈ 0.0

# output

```
"""
function sumweightedmassneighbors!(
    particleindex::Int,
    particleids::Array{Int},
    typeids::Array{Int},
    neighborindexes::Array{Int},
    dimension::Int,
    r::StateData,
    radius::StateData,
    m::StateData,
    ρ::StateData,
)
    particleid = particleids[particleindex]
    typeid = typeids[particleindex]

    # Loop over pairs for target particle
    for i ∈ neighborindexes[particleindex]:(neighborindexes[particleindex+1]-1)
        neighborindex = neighborindexes[i]
        neighborid = particleids[neighborindex]
        neighbortypeid = typeids[neighborindex]

        absdistance = norm(
            view(r[typeid], ((particleid-1)*dimension+1):(particleid*dimension)) -
            view(r[neighbortypeid], ((neighborid-1)*dimension+1):(neighborid*dimension)),
        )

        # Add contribution to target particle from each pair
        ρ[typeid][particleid] +=
            w(absdistance, radius[typeid][particleid], dimension) * m[typeid][particleid]
    end
    return nothing
end
