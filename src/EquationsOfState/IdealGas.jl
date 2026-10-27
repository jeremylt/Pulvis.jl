# ------------------------------------------------------------------------------
# Equations of state - ideal gas law
# ------------------------------------------------------------------------------

"""
```julia
densitytopressure!(particleid, ρ, ρ0, c, refdensityfactor, exponent, p)
```

Convert density to pressure

# Arguments:

  - `particleid::Int`:         index number of current particle
  - `ρ::Array`:                input array of particle densities
  - `ρ0::Float`:               initial density
  - `c::Float`:                speed of sound
  - `refdensityfactor::Float`: reference density factor
  - `exponent::Int`:           exponent
  - `p::Array`:                output array of particle pressures

# Returns:

  - none

# Example:

```jldoctest
# setup
particlei = 1
neighbors = [2, 3]
dimension = 2
r = [1.0, 1.0, 1.1, 1.2, 1.3, 1.4]
radius = [1.0, 1.1, 1, 2]
m = [42.0, 13.0, 14.0]
ρ = [1.0, 0.0, 0.0]

# compute density
sumweightedmassneighbors!(particleid, neighbors, dimension, r, radius, m, ρ)

# verify
@assert ρ[particleid] ≈ 42.0
@assert ρ[neighbors[1]] = 0.0
@assert ρ[neighbors[2]] = 0.0

# output

```
"""
function densitytopressure_IdealGas!(particleid, ρ, ρ0, c, refdensityfactor, exponent, p)
    p[particleid] = c^2 * ρ[particleid]
end
