using Base.Threads

# ------------------------------------------------------------------------------
# Density
# ------------------------------------------------------------------------------

"""
```julia
computedensity_cpu!(numneighbors, neighbors, dimension, position, radius, mass, density)
```

Compute the density of all particles on the CPU

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
function computedensity_cpu!(
    numneighbors::Array,
    neighbors::Array,
    dimension::Int,
    position::Array,
    radius::Array,
    mass::Array,
    density::Array,
)
    @threads for particle = 1:length(density)
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

# ------------------------------------------------------------------------------
# Pressure
# ------------------------------------------------------------------------------

"""
```julia
computepressure_cpu!(
    density,
    typeids,
    initialdensity,
    speedofsound,
    refdensityfactor,
    exponent,
    pressure,
)
```

Compute the pressure for all particles on the CPU

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
function computepressure_cpu!(
    density::Array,
    typeids::Array,
    initialdensity::Array,
    speedofsound,
    refdensityfactor::Array,
    exponent::Array,
    pressure::Array,
)
    @threads for particle = 1:length(pressure)
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

# ------------------------------------------------------------------------------
# Momentum and forces
# ------------------------------------------------------------------------------

"""
```julia
computeforces_cpu!(
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

Compute the forces for all particles on the CPU

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
function computeforces_cpu!(
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
    @threads for particle = 1:(length(modifiedacceleration)/dimension)
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
