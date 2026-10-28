using CUDA

# ------------------------------------------------------------------------------
# Density
# ------------------------------------------------------------------------------

"""
```julia
computedensity_particle!(neighborids, dimension, r, radius, m, ρ)
```

Compute the density of a single particle on a single GPU thread

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
function computedensity_particle!(
    neighborids::Array,
    dimension::Int,
    r::Array,
    radius::Array,
    m::Array,
    ρ::Array,
)
    particleid = threadIdx().x

    if particleid <= length(ρ)
        sumweightedmassself!(particleid, radius, m, ρ)
        sumweightedmassneighbors!(particleid, neighborids, dimension, r, radius, m, ρ)
    end
    return nothing
end

"""
```julia
computedensity_cuda!(neighborids, dimension, r, radius, m, ρ)
```

Compute the density of all particles on the GPU

# Arguments:

  - `neighborids::Array`: input array of neighbor particle indices
  - `dimension::Int`:     spatial dimension of problem
  - `r::Array`:           input array of particle positions
  - `radius::Array`:      input array of particle radi
  - `m::Array`:           input array of particle masses
  - `ρ::Array`:            output array of particle densities

# Returns:

  - none
"""
function computedensity_cuda!(
    neighborids::Array,
    dimension::Int,
    r::Array,
    radius::Array,
    m::Array,
    ρ::Array,
)
    @cuda threads = length(ρ) computedensity_particle!(
        neighborids,
        dimension,
        r,
        radius,
        m,
        ρ,
    )
    return nothing
end

# ------------------------------------------------------------------------------
# Pressure
# ------------------------------------------------------------------------------

"""
```julia
computepressure_particle!(ρ, typeids, ρ0, c, refdensityfactor, exponent, p)
```

Compute the pressure for single particle on a single GPU thread

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
function computepressure_particle!(
    ρ::Array,
    typeids::Array,
    ρ0::Array,
    c,
    refdensityfactor::Array,
    exponent::Array,
    p::Array,
)
    particleid = threadIdx().x

    if particleid <= length(p)
        typeid = typeids[particleid]

        p[particleid] = densitytopressure!(
            ρ[particleid],
            ρ0[typeid],
            c,
            refdensityfactor[typeid],
            exponent[typeid],
        )
    end
    return nothing
end

"""
```julia
computepressure_cuda!(ρ, typeids, ρ0, c, refdensityfactor, exponent, p)
```

Compute the pressure for all particles on the GPU

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
function computepressure_cuda!(
    ρ::Array,
    typeids::Array,
    ρ0::Array,
    c,
    refdensityfactor::Array,
    exponent::Array,
    p::Array,
)
    @cuda threads = length(p) computepressure_particle!(
        ρ,
        typeids,
        ρ0,
        c,
        refdensityfactor,
        exponent,
        p,
    )
    return nothing
end

# ------------------------------------------------------------------------------
# Momentum and forces
# ------------------------------------------------------------------------------

"""
```julia
computeforces_particle!(
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

Compute the forces for a single particle on a single GPU thread

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
function computeforces_particle!(
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
    particleid = threadIdx().x

    if particleid <= length(a_modified)/dimension
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

"""
```julia
computeforces_cuda!(
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

Compute the forces for all particles on the GPU

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
function computeforces_cuda!(
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
    @cuda threads = (length(a_modified)/dimension) accelerationaddmomentumcontribution_particle!(
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
    return nothing
end
