# ------------------------------------------------------------------------------
# density
# ------------------------------------------------------------------------------

"""
```julia
computedensity_particle!(
    numneighbors,
    neighbors,
    dimension,
    position,
    radius,
    mass,
    density,
)
```

Compute the density of a single particle on a single GPU thread

# Arguments:

  - `numneighbors::Array`: input array of numbers of neighboring particles
  - `neighbors::Array`:    input array of neighbor particle index numbers
  - `dimension::Int`:      spatial dimension of problem
  - `position::Array`:     input array of particle positions
  - `radius::Array`:       input array of particle radi
  - `mass::Array`:         input array of particle masses
  - `density::Array`:      output array of particle densities

# Returns:

  - none
"""
function computedensity_particle!(
    numneighbors::Array,
    neighbors::Array,
    dimension::Int,
    position::Array,
    radius::Array,
    mass::Array,
    density::Array,
)
    particle = threadIdx().x

    if particle <= numparticles
        sumweightedmassself!(particle, radius, mass, density)
        sumweightedmassneighbors!(
            particle,
            numneighbors[particle],
            neighbors,
            dimension,
            position,
            radius,
            mass,
            density,
        )
    end
end

"""
```julia
computedensity_cuda!(numneighbors, neighbors, dimension, position, radius, mass, density)
```

Compute the density of all particles on the GPU

# Arguments:

  - `numneighbors::Array`: input array of numbers of neighboring particles
  - `neighbors::Array`:    input array of neighbor particle index numbers
  - `dimension::Int`:      spatial dimension of problem
  - `position::Array`:     input array of particle positions
  - `radius::Array`:       input array of particle radi
  - `mass::Array`:         input array of particle masses
  - `density::Array`:      output array of particle densities

# Returns:

  - none
"""
function computedensity_cuda!(
    numneighbors::Array,
    neighbors::Array,
    dimension::Int,
    position::Array,
    radius::Array,
    mass::Array,
    density::Array,
)
    @cuda threads = length(density) computedensity_particle!(
        numneighbors,
        neighbors,
        dimension,
        position,
        radius,
        mass,
        density,
    )
end

# ------------------------------------------------------------------------------
# pressure
# ------------------------------------------------------------------------------

"""
```julia
computepressure_particle!(
    density,
    typeids,
    initialdensity,
    speedofsound,
    refdensityfactor,
    exponent,
    pressure,
)
```

Compute the pressure for single particle on a single GPU thread

# Arguments:

  - `density::Array`:          input array of particle densities
  - `typeids::Array`:          input array of particle type ids
  - `initialdenstity::Array`:  input array of initial densities by particle type
  - `speedofsound::Float`:     speed of sound
  - `refdensityfactor::Array`: input array of reference density factors by particle type
  - `exponent::Array`:         input array of exponent by particle type
  - `pressure::Array`:         output array of particle pressures

# Returns:

  - none
"""
function computepressure_particle!(
    density::Array,
    typeids::Array,
    initialdensity::Array,
    speedofsound,
    refdensityfactor::Array,
    exponent::Array,
    pressure::Array,
)
    particle = threadIdx().x

    if particle <= numparticles
        typeid = typeids[particle]

        densitytopressure!(
            particle,
            density,
            initialdensity[typeid],
            speedofsound,
            refdensityfactor[typeid],
            exponent[typeid],
            pressure,
        )
    end
end

"""
```julia
computepressure_cuda!(
    density,
    typeids,
    initialdensity,
    speedofsound,
    refdensityfactor,
    exponent,
    pressure,
)
```

Compute the pressure for all particles on the GPU

# Arguments:

  - `density::Array`:          input array of particle densities
  - `typeids::Array`:          input array of particle type ids
  - `initialdenstity::Array`:  input array of initial densities by particle type
  - `speedofsound::Float`:     speed of sound
  - `refdensityfactor::Array`: input array of reference density factors by particle type
  - `exponent::Array`:         input array of exponent by particle type
  - `pressure::Array`:         output array of particle pressures

# Returns:

  - none
"""
function computepressure_cuda!(
    density::Array,
    typeids::Array,
    initialdensity::Array,
    speedofsound,
    refdensityfactor::Array,
    exponent::Array,
    pressure::Array,
)
    @cuda threads = length(pressure) computepressure_particle!(
        density,
        initialdensity,
        speedofsound,
        refdensityfactor,
        exponent,
        pressure,
    )
end

# ------------------------------------------------------------------------------
# momentum and forces
# ------------------------------------------------------------------------------

"""
```julia
computeforces_particle!(
    numneighbors,
    neighbors,
    typeids,
    dynamicviscosity,
    bulkviscosity,
    dimension,
    position,
    radius,
    mass,
    density,
    pressure,
    velocity,
    acceleration,
    modifiedvelocity,
    modifiedacceleration,
)
```

Compute the forces for a single particle on a single GPU thread

# Arguments:

  - `numneighbors::Array`:         number of neighboring particles
  - `neighbors::Array`:            input array of neighbor particle index numbers
  - `typeids::Array`:              input array of particle type ids
  - `dynamicviscosity::Array`      input array of dynamic viscosity for each particle type
  - `bulkviscosity::Array`         input array of bulk viscosity for each particle type
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
"""
function computeforces_cuda!(
    numneighbors::Int,
    neighbors::Array,
    typeids::Array,
    dynamicviscosity::Array,
    bulkviscosity::Array,
    dimension::Int,
    position::Array,
    radius::Array,
    mass::Array,
    density::Array,
    pressure::Array,
    velocity::Array,
    acceleration::Array,
    modifedvelocity::Array,
    modifiedacceleration::Array,
)
    particle = threadIdx().x

    if particle <= numparticles
        momentumequationneighbors!(
            particle,
            numneighbors[particle],
            neighbors,
            typeids,
            dynamicviscosity,
            bulkviscosity,
            position,
            radius,
            mass,
            density,
            pressure,
            velocity,
            acceleration,
            modifedvelocity,
            modifiedacceleration,
        )
    end
end

"""
```julia
computeforces_cuda!(
    numneighbors,
    neighbors,
    typeids,
    dynamicviscosity,
    bulkviscosity,
    dimension,
    position,
    radius,
    mass,
    density,
    pressure,
    velocity,
    acceleration,
    modifiedvelocity,
    modifiedacceleration,
)
```

Compute the forces for all particles on the GPU

# Arguments:

  - `numneighbors::Array`:         number of neighboring particles
  - `neighbors::Array`:            input array of neighbor particle index numbers
  - `typeids::Array`:              input array of particle type ids
  - `dynamicviscosity::Array`      input array of dynamic viscosity for each particle type
  - `bulkviscosity::Array`         input array of bulk viscosity for each particle type
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
"""
function computeforces_cuda!(
    numneighbors::Array,
    neighbors::Array,
    typeids::Array,
    dynamicviscosity::Array,
    bulkviscosity::Array,
    dimension::Int,
    position::Array,
    radius::Array,
    mass::Array,
    density::Array,
    pressure::Array,
    velocity::Array,
    acceleration::Array,
    modifedvelocity::Array,
    modifiedacceleration::Array,
)
    @cuda threads = (length(modifiedacceleration)/dimension) accelerationaddmomentumcontribution_particle!(
        numneighbors,
        neighbors,
        typeids,
        dynamicviscosity,
        bulkviscosity,
        position,
        radius,
        mass,
        density,
        pressure,
        velocity,
        acceleration,
        modifedvelocity,
        modifiedacceleration,
    )
end
