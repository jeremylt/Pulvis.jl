# ------------------------------------------------------------------------------
# equations of state, ideal gas law
# ------------------------------------------------------------------------------

"""
```julia
densitytopressure!(
    particle,
    density,
    initialdensity,
    speedofsound,
    refdensityfactor,
    exponent,
    pressure,
)
```

Convert density to pressure

# Arguments:

  - `particle::Int`:           index number of current particle
  - `density::Array`:          input array of particle densities
  - `initialdenstity::Float`:  initial density
  - `speedofsound::Float`:     speed of sound
  - `refdensityfactor::Float`: reference density factor
  - `exponent::Int`:           exponent
  - `pressure::Array`:         output array of particle pressures

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
function densitytopressure_idealgas!(
    particle,
    density,
    initialdensity,
    speedofsound,
    refdensityfactor,
    exponent,
    pressure,
)
    pressure[particle] = speedofsound^2 * density[particle]
end
