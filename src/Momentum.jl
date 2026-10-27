# ------------------------------------------------------------------------------
# Momentum computation at a single point
# ------------------------------------------------------------------------------

"""
```julia
momentumequationneighbors!(
    particleid,
    neighborids,
    typeids,
    μ_dynamic,
    μ_bulk,
    dimension,
    r,
    radius,
    m,
    ρ,
    p,
    v,
    a,
    v_modified,
    a_modified,
)
```

Sum the momentum contributions to the particle from neighbor particles

# Arguments:

  - `particle::Int`:       index number of current particle
  - `neighborids::Array`:  input array of neighbor particle indices
  - `typeids::Array`:      input array of particle type ids
  - `μ_dynamic::Array`     input array of dynamic viscosity for each particle type
  - `μ_bulk::Array`        input array of bulk viscosity for each particle type
  - `dimension::Int`:      spatial dimension of problem
  - `r::Array`:            input array of particle positions
  - `radius::Array`:       input array of particle radi
  - `m::Array`:            input array of particle masses
  - `ρ::Array`:            input array of particle densities
  - `p::Array`:            input array of particle pressures
  - `v::Array`:            input array of particle velocities
  - `a::Array`:            input array of particle acceleration
  - `v_modified::Array`:   output array of modified particle velocities
  - `a_modified::Array`:   output array of modified particle acceleration

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
@assert ρ[particleid] ≈ 42.0
@assert ρ[neighborids[3]] ≈ 0.0
@assert ρ[neighborids[4]] ≈ 0.0

# output

```
"""
function momentumequationneighbors!(
    particleid::Int,
    neighborids::Array,
    typeids::Array,
    μ_dynamic::Array,
    μ_bulk::Array,
    dimension::Int,
    r::Array,
    radius::Array,
    m::Array,
    ρ::Array,
    p::Array,
    v::Array,
    a::Array,
    v_modified::Array,
    a_modified::Array,
)
    for i = neighborids[particleid]:(neighborids[particleid+1]-1)
        neighborid = neighborids[i]
        particletype = typeids[particleid]
        neighbortype = typeids[neighborid]
        rij =
            view(r, ((neighborid-1)*dimension+1):(neighborid*dimension)) -
            view(r, ((particleid-1)*dimension+1):(particleid*dimension))
        absdistance = norm(rij)

        # specific coefficient
        dWdrij = dwdrij(absdistance, radius[particleid], absdistance)

        speccoeff = specificcoefficient(
            ρ[particleid],
            ρ[neighborid],
            m[particleid],
            m[neighborid],
            dWdrij,
        )

        # pressure gradient
        eij = (-1.0 / absdistance) * rij

        pressuregradient!(
            ρ[particleid],
            ρ[neighborid],
            p[particleid],
            p[neighborid],
            speccoeff,
            eij,
            view(a, ((particleid-1)*dimension+1):(particleid*dimension)),
        )

        # shear forces
        shearforces!(
            ρ[particleid],
            ρ[neighborid],
            view(v, ((particleid-1)*dimension+1):(particleid*dimension)),
            view(v, ((neighborid-1)*dimension+1):(neighborid*dimension)),
            dimension + 2,
            μ_dynamic[particletype],
            μ_dynamic[neighbortype],
            μ_bulk[particletype],
            μ_bulk[neighbortype],
            absdistance,
            speccoeff,
            eij,
            view(a, ((particleid-1)*dimension+1):(particleid*dimension)),
        )

        # transport velocity
    end
    return nothing
end
