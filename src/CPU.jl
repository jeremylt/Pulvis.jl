using Base.Threads

# ------------------------------------------------------------------------------
# Density
# ------------------------------------------------------------------------------

"""
```julia
computedensity_cpu!(neighborids, dimension, p, r, m, ρ)
```

Compute the density of all particles on the CPU

# Arguments:

  - `neighborids::Array`: input array of neighbor particle indices
  - `dimension::Int`:     spatial dimension of problem
  - `r::Array`:           input array of particle positions
  - `radius::Array`:      input array of particle radi
  - `m::Array`:           input array of particle masses
  - `ρ::Array`:           output array of particle densities

# Returns:

  - none
"""
function computedensity_cpu!(
    neighborids::Array,
    dimension::Int,
    r::Array,
    radius::Array,
    m::Array,
    ρ::Array,
)
    @threads for particleid = 1:length(ρ)
        sumweightedmassself!(particleid, radius, m, ρ)
        sumweightedmassneighbors!(particleid, neighborids, dimension, r, radius, m, ρ)
    end
    return nothing
end

# ------------------------------------------------------------------------------
# Pressure
# ------------------------------------------------------------------------------

"""
```julia
computepressure_cpu!(ρ, typeids, ρ0, c, refdensityfactor, exponent, p)
```

Compute the pressure for all particles on the CPU

# Arguments:

  - `ρ::Array`:                input array of particle densities
  - `typeids::Array`:          input array of particle type ids
  - `ρ0::Array`:               input array of initial densities by particle type
  - `c::Float`:                speed of sound
  - `refdensityfactor::Array`: input array of reference density factors by particle type
  - `exponent::Array`:         input array of exponent by particle type
  - `p::Array`:                output array of particle pressures

# Returns:

  - none
"""
function computepressure_cpu!(
    ρ::Array,
    typeids::Array,
    ρ0::Array,
    c,
    refdensityfactor::Array,
    exponent::Array,
    p::Array,
)
    @threads for particleid = 1:length(p)
        typeid = typeids[particleid]

        densitytopressure!(
            particleid,
            ρ,
            ρ0[typeid],
            c,
            refdensityfactor[typeid],
            exponent[typeid],
            p,
        )
    end
    return nothing
end

# ------------------------------------------------------------------------------
# Momentum and forces
# ------------------------------------------------------------------------------

"""
```julia
computeforces_cpu!(
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

Compute the forces for all particles on the CPU

# Arguments:

  - `neighborids::Array`: input array of neighbor particle indices
  - `typeids::Array`:     input array of particle type ids
  - `μ_dynamic::Array`:   input array of dynamic viscosity for each particle type
  - `μ_bulk::Array`:      input array of bulk viscosity for each particle type
  - `dimension::Int`:     spatial dimension of problem
  - `r::Array`:           input array of particle positions
  - `radius::Array`:      input array of particle radi
  - `m::Array`:           input array of particle masses
  - `ρ::Array`:           input array of particle densities
  - `p::Array`:           input array of particle pressures
  - `v::Array`:           input array of particle velocities
  - `a::Array`:           input array of particle acceleration
  - `v_modified::Array`:  output array of modified particle velocities
  - `a_modified::Array`:  output array of modified particle acceleration

# Returns:

  - none
"""
function computeforces_cpu!(
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
    @threads for particleid = 1:(length(a_modified)/dimension)
        momentumequationneighbors!(
            particleid,
            neighborids,
            typeids,
            μ_dynamic,
            μ_bulk,
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
    end
    return nothing
end
