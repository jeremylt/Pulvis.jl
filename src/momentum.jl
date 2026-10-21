# ------------------------------------------------------------------------------
# momentum computation at a single point
# ------------------------------------------------------------------------------

"""
```julia
momentumequationneighbors!(particle, numneighbors, neighbors, typeids, dynamicviscosity, bulkviscosity, dimension, position, radius, mass, density, pressure, velocity, acceleration, modifiedvelocity, modifiendacceleration)
```

Sum the momentum contributions to the particle from neighbor particles

# Arguments:

  - `particle::Int`:               index number of current particle
  - `numneighbors::Int`:           number of neighboring particles
  - `neighbors::Array`:            input array of neighbor particle index numbers
  - `dimension::Int`:              spatial dimension of problem
  - `position::Array`:             input array of particle positions
  - `radius::Array`:               input array of particle radi
  - `mass::Array`:                 input array of particle masses
  - `density::Array`:              input array of particle densities
  - `pressure::Array`:             input array of particle pressures
  - `velocity::Array`:             input array of particle velocities
  - `acceleration::Array`:         input array of particle acceleration
  - `modifiedvelocity::Array`:     output array of modified particle velocities
  - `modifiedacceleration::Array`: output array of modified particle acceleration

# Returns:

  - none

# Example:

```jldoctest
# setup
particle     = 1
numneighbors = 2
neighbors    = [2, 3]
dimension    = 2
position     = [1.0, 1.0, 1.1, 1.2, 1.3, 1.4]
radius       = [1.0, 1.1, 1,2]
mass         = [42.0, 13.0, 14.0]
density      = [1.0, 0.0, 0.0]

# compute density
sumweightedmassneighbors!(particle, numneighbors, neighbors, dimension, position, radius, mass, density)

# verify
@assert density[particle] ≈ 42.0
@assert density[neighbors[1]] = 0.0
@assert density[neighbors[2]] = 0.0

# output

```
"""
function momentumequationneighbors!(particle, numneighbors, neighbors, typeids, dynamicviscosity, bulkviscosity, dimension, position, radius, mass, density, pressure, velocity, acceleration, modifedvelocity, modifiedacceleration)
    for neighborid = 1:numneighbors
        neighbor = neighbors[neighborid]
        particletype = typeids[particle]
        neighbortype = typeids[neighbor]
        rij = view(position, neighbor*dimension:(neighbor+1)*dimension) - view(position, particle*dimension:(particle+1)*dimension)
        absdistance = distancebetweenparticles(view(position, particle*dimension:(particle+1)*dimension), view(position, neighbor*dimension:(neighbor+1)*dimension))

        # specific coefficient
        dWdrij = dwdrij(absdistance, radius[particle], distance)

        speccoeff = specificcoefficient(density[particle], density[neighbor], mass[particle], mass[neighbor], dwdr)

        # pressure gradient
        eij = (-1.0 / absdistance) * rij

        pressuregradient!(density[particle], density[neighbor], pressure[particle], pressure[neighbor], speccoeff, eij, view(acceleration, particle*dimension:(particle+1)*dimension))

        # shear forces
        shearforces!(density[particle], density[neighbor], view(velocity, particle*dimension:(particle+1)*dimension), view(velocity, neighbor*dimension:(neighbor+1)*dimension), dimension + 2, dynamicviscosity[particletype], dynamicviscosity[neighbortype], bulkviscosity[particletype], bulkviscosity[neighbortype], absdistance, speccoeff, eij, view(acceleration, particle*dimension:(particle+1)*dimension))

        # transport velocity
end
