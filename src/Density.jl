# ------------------------------------------------------------------------------
# Density computation at a single point
# ------------------------------------------------------------------------------

"""
```julia
sumweightedmassself!(particle, radius, mass, density)
```

Compute the density of the particle from the shape function, radius, and mass

# Arguments:

  - `particle::Int`:  index number of current particle
  - `radius::Array`:  input array of particle radi
  - `mass::Array`:    input array of particle masses
  - `density::Array`: output array of particle densities

# Returns:

  - none

# Example:

```jldoctest
# setup
particle = 1
radius = [1.0]
mass = [42.0]
density = [0.0]

# compute density
sumweightedmassself!(particle, radius, mass, density)

# verify
@assert density[particle] ≈ 42.0

# output

```
"""
function sumweightedmassself!(particle::Int, radius::Array, mass::Array, density::Array)
    density[particle] = w0(radius[particle]) * mass[particle]
end

"""
```julia
sumweightedmassneighbors!(
    particle,
    numneighbors,
    neighbors,
    dimension,
    position,
    radius,
    mass,
    density,
)
```

Sum the density contributions to the particle from neighbor particles

# Arguments:

  - `particle::Int`:     index number of current particle
  - `numneighbors::Int`: number of neighboring particles
  - `neighbors::Array`:  input array of neighbor particle index numbers
  - `dimension::Int`:    spatial dimension of problem
  - `position::Array`:   input array of particle positions
  - `radius::Array`:     input array of particle radi
  - `mass::Array`:       input array of particle masses
  - `density::Array`:    output array of particle densities

# Returns:

  - none

# Example:

```jldoctest
# setup
particle = 1
numneighbors = 2
neighbors = [2, 3]
dimension = 2
position = [1.0, 1.0, 1.1, 1.2, 1.3, 1.4]
radius = [1.0, 1.1, 1, 2]
mass = [42.0, 13.0, 14.0]
density = [1.0, 0.0, 0.0]

# compute density
sumweightedmassneighbors!(
    particle,
    numneighbors,
    neighbors,
    dimension,
    position,
    radius,
    mass,
    density,
)

# verify
@assert density[particle] ≈ 42.0
@assert density[neighbors[1]] = 0.0
@assert density[neighbors[2]] = 0.0

# output

```
"""
function sumweightedmassneighbors!(
    particle::Int,
    numneighbors::Int,
    neighbors::Array,
    dimension::Int,
    position::Array,
    radius::Array,
    mass::Array,
    density::Array,
)
    for neighbor = 1:numneighbors
        neighbor = neighbors[neighborid]
        absdistance = norm(
            view(position, (particle*dimension):((particle+1)*dimension)) -
            view(position, (neighbor*dimension):((neighbor+1)*dimension)),
        )

        density[particle] += w(absdistance, radius[particle]) * mass[particle]
    end
end
