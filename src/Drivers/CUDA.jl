# ------------------------------------------------------------------------------
# Density
# ------------------------------------------------------------------------------

"""
```julia
computedensity_particle!(particleids, typeids, neighborindexes, dimension, p, r, m, ρ)
```

Compute the density of a single particle on a single GPU thread

# Arguments:

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
"""
function computedensity_particle!(
    particleids::Array,
    typeids::Array,
    neighboridexes::Array,
    dimension::Int,
    r::StateData,
    radius::StateData,
    m::StateData,
    ρ::StateData,
)
    i = threadIdx().x

    if i <= length(particleids)
        particleid = particleids[i]
        typeid = typeids[i]

        sumweightedmassself!(particleid, radius[typeid], m[typeid], ρ[typeid])
        sumweightedmassneighbors!(
            i,
            particleids,
            typeids,
            neighborindexes,
            dimension,
            r,
            radius,
            m,
            ρ,
        )
    end
    return nothing
end

"""
```julia
computedensity_cuda!(particleids, typeids, neighborindexes, dimension, p, r, m, ρ)
```

Compute the density of all particles on the GPU

# Arguments:

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
"""
function computedensity_cuda!(
    particleids::Array,
    typeids::Array,
    neighboridexes::Array,
    dimension::Int,
    r::StateData,
    radius::StateData,
    m::StateData,
    ρ::StateData,
)
    @cuda threads = length(particleids) computedensity_particle!(
        particleids,
        typeids,
        neighborindexes,
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
computepressure_particle!(particleids, typeids, ρ, ρ0, c, refdensityfactor, exponent, p)
```

Compute the pressure for single particle on a single GPU thread

# Arguments:

  - `particleids::Array`:      input array of particle ids by type
  - `typeids::Array`:          input array of particle type ids
  - `ρ::StateData`:            input arrays of particle densities by particle type and id
  - `ρ0::Array`:               input array of initial densities by particle type
  - `c::Float`:                speed of sound
  - `refdensityfactor::Array`: input array of reference density factors by particle type
  - `exponent::Array`:         input array of exponent by particle type
  - `p::StateData`:            output arrays of particle pressures by particle type and id

# Returns:

  - none
"""
function computepressure_particle!(
    particleids::Array,
    typeids::Array,
    ρ::StateData,
    ρ0::Array,
    c,
    refdensityfactor::Array,
    exponent::Array,
    p::StateData,
)
    i = threadIdx().x

    if i <= length(particleids)
        particleid = particleids[i]
        typeid = typeids[i]

        p[typeid][particleid] = densitytopressure!(
            ρ[typeid][particleid],
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
computepressure_cuda!(particleids, typeids, ρ, ρ0, c, refdensityfactor, exponent, p)
```

Compute the pressure for all particles on the GPU

# Arguments:

  - `particleids::Array`:      input array of particle ids by type
  - `typeids::Array`:          input array of particle type ids
  - `ρ::StateData`:            input arrays of particle densities by particle type and id
  - `ρ0::Array`:               input array of initial densities by particle type
  - `c::Float`:                speed of sound
  - `refdensityfactor::Array`: input array of reference density factors by particle type
  - `exponent::Array`:         input array of exponent by particle type
  - `p::StateData`:            output arrays of particle pressures by particle type and id

# Returns:

  - none
"""
function computepressure_cuda!(
    particleids::Array,
    typeids::Array,
    ρ::StateData,
    ρ0::Array,
    c,
    refdensityfactor::Array,
    exponent::Array,
    p::StateData,
)
    @cuda threads = length(particleids) computepressure_particle!(
        particleids,
        typeids,
        ρ,
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
    particleids,
    typeids,
    neighborids,
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

  - `particleids::Array`:     input array of particle ids by type
  - `typeids::Array`:         input array of particle type ids
  - `neighborindexes::Array`: input array of neighbor particle global indexes
  - `μ_dynamic::Array`:       input array of dynamic viscosity by particle type
  - `μ_bulk::Array`:          input array of bulk viscosity by particle type
  - `dimension::Int`:         spatial dimension of problem
  - `r::StateData`:           input arrays of particle positions by particle type and id
  - `radius::StateData`:      input arrays of particle radi by particle type and id
  - `m::StateData`:           input arrays of particle masses by particle type and id
  - `ρ::StateData`:           input arrays of particle densities by particle type and id
  - `p::StateData`:           input arrays of particle pressures by particle type and id
  - `v::StateData`:           input arrays of particle velocities by particle type and id
  - `a::StateData`:           input arrays of particle acceleration by particle type and id
  - `v_modified::StateData`:  output arrays of modified particle velocities by particle type and id
  - `a_modified::StateData`:  output arrays of modified particle acceleration by particle type and id

# Returns:

  - none
"""
function computeforces_particle!(
    particleids::Array,
    typeids::Array,
    neighborindexes::Array,
    μ_dynamic::Array,
    μ_bulk::Array,
    dimension::Int,
    r::StateData,
    radius::StateData,
    m::StateData,
    ρ::StateData,
    p::StateData,
    v::StateData,
    a::StateData,
    v_modified::StateData,
    a_modified::StateData,
)
    i = threadIdx().x

    if i <= length(particleids)
        momentumequationneighbors!(
            i,
            particleids,
            typeids,
            neighborindexes,
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
    particleids,
    typeids,
    neighborids,
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

  - `particleids::Array`:     input array of particle ids by type
  - `typeids::Array`:         input array of particle type ids
  - `neighborindexes::Array`: input array of neighbor particle global indexes
  - `μ_dynamic::Array`:       input array of dynamic viscosity by particle type
  - `μ_bulk::Array`:          input array of bulk viscosity by particle type
  - `dimension::Int`:         spatial dimension of problem
  - `r::StateData`:           input arrays of particle positions by particle type and id
  - `radius::StateData`:      input arrays of particle radi by particle type and id
  - `m::StateData`:           input arrays of particle masses by particle type and id
  - `ρ::StateData`:           input arrays of particle densities by particle type and id
  - `p::StateData`:           input arrays of particle pressures by particle type and id
  - `v::StateData`:           input arrays of particle velocities by particle type and id
  - `a::StateData`:           input arrays of particle acceleration by particle type and id
  - `v_modified::StateData`:  output arrays of modified particle velocities by particle type and id
  - `a_modified::StateData`:  output arrays of modified particle acceleration by particle type and id

# Returns:

  - none
"""
function computeforces_cuda!(
    particleids::Array,
    typeids::Array,
    neighborindexes::Array,
    μ_dynamic::Array,
    μ_bulk::Array,
    dimension::Int,
    r::StateData,
    radius::StateData,
    m::StateData,
    ρ::StateData,
    p::StateData,
    v::StateData,
    a::StateData,
    v_modified::StateData,
    a_modified::StateData,
)
    @cuda threads = length(particleids) accelerationaddmomentumcontribution_particle!(
        particleids,
        typeids,
        neighborindexes,
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
